import java.util.Arrays;
import java.util.PriorityQueue;
import java.util.Random;
import java.io.*;
import org.apache.commons.math3.distribution.NormalDistribution;


/**
 * Samples reads from the length bins created by $BuildReadLengthBins$, using a
 * bimodal normal distribution given in input.
 */
public class SampleReadsFromLengthBins {
    /**
     * The distribution from which bins are sampled
     */
    private static double[] cumulativeBimodal;
    
    /**
     * The binned set of reads from which sampling takes place. Not all bins are
     * kept in memory at the same time.
     */
    private static String[][] buffers;
    private static int[] buffers_last;  // Last unsampled read
    private static long[] buffers_stringLength;
    private static boolean[] buffers_isLoaded;
    private static long nBpsInRam, maxBpsInRam;
    private static PriorityQueue<Bin> queue;  // LRU over bins
    private static Bin[] bins;  // Objects to be inserted into $queue$.
    
    /**
     * Global source of randomness
     */
    private static Random random;
    
    
    /**
     * @param args 
     *  4: weight of the left normal, in [0..1];
     * 10: the actual RAM usage is at least twice this much, because of the 
     *     quality track in FASTQs.
     */
	public static void main(String[] args) throws IOException {
        final double MEAN_LEFT = Double.parseDouble(args[0]);
        final double STD_LEFT = Double.parseDouble(args[1]);
        final double MEAN_RIGHT = Double.parseDouble(args[2]);
        final double STD_RIGHT = Double.parseDouble(args[3]);
        final double WEIGHT_LEFT = Double.parseDouble(args[4]);
        final int BIN_LENGTH = Integer.parseInt(args[5]);
        final int MAX_READ_LENGTH = Integer.parseInt(args[6]);  // For a bin
		final String BUCKETS_PREFIX = args[7];
        final long GENOME_LENGTH_HAPLOID = Long.parseLong(args[8]);
        final double TARGET_COVERAGE_HAPLOID = Double.parseDouble(args[9]);
        maxBpsInRam=Long.parseLong(args[10]);
        final String OUTPUT_FILE = args[11];
        
        final int N_BINS = MAX_READ_LENGTH/BIN_LENGTH+1;
        final long TARGET_COVERAGE_BP = (long)(GENOME_LENGTH_HAPLOID*(TARGET_COVERAGE_HAPLOID*2));
        final int MAX_SAMPLING_ATTEMPTS = N_BINS*10;  // Arbitrary
        
        int i, j;
        int bin, length;
        long coverage, nBpsInRam;
        BufferedWriter bw;
        int[] histogram;
        String[] tokens;
        
        random = new Random();
        buildCumulativeBimodal(N_BINS,BIN_LENGTH,MEAN_LEFT,STD_LEFT,MEAN_RIGHT,STD_RIGHT,WEIGHT_LEFT);
        buffers = new String[N_BINS][0];
        buffers_last = new int[N_BINS];
        Arrays.fill(buffers_last,Integer.MAX_VALUE);
        buffers_stringLength = new long[N_BINS];
        Arrays.fill(buffers_stringLength,0);
        buffers_isLoaded = new boolean[N_BINS];
        Arrays.fill(buffers_isLoaded,false);
        histogram = new int[N_BINS];
        Arrays.fill(histogram,0);
        bins = new Bin[N_BINS];
        for (i=0; i<N_BINS; i++) bins[i] = new Bin(i,0);
        queue = new PriorityQueue<Bin>();
        coverage=0; nBpsInRam=0;
        bw = new BufferedWriter(new FileWriter(OUTPUT_FILE));
        while (coverage<TARGET_COVERAGE_BP) {
            j=0; bin=-1;
            while (j<MAX_SAMPLING_ATTEMPTS) {
                bin=Arrays.binarySearch(cumulativeBimodal,random.nextDouble());
                if (bin<0) bin=-1-bin;
                j++;
                if (buffers_last[bin]!=-1) {  // The bin cannot be empty
                    loadBin(bin,BUCKETS_PREFIX);
                    if (buffers_last[bin]!=-1) {  // The bin cannot be empty
                        freeSpace(bin);
                        break;
                    }
                }
            }
            if (j==MAX_SAMPLING_ATTEMPTS) {
                System.err.println("ERROR: Not enough reads to sample "+TARGET_COVERAGE_HAPLOID+"x haploid coverage from the given bins and distribution.");
                System.exit(1);
            }
            histogram[bin]++;
            tokens=buffers[bin][buffers_last[bin]].split(BuildReadLengthBins.BIN_FILE_SEPARATOR);
            bw.write(tokens[0]); bw.newLine();
            bw.write(tokens[1]); bw.newLine();
            length=tokens[1].length();
            coverage+=length;
            buffers_stringLength[bin]-=length;
            bw.write(tokens[2]); bw.newLine();
            bw.write(tokens[3]); bw.newLine();
            buffers_last[bin]--;
        }
        bw.close();
        bw = new BufferedWriter(new FileWriter(OUTPUT_FILE+".histogram"));
        for (i=0; i<N_BINS; i++) bw.write(histogram[i]+"\n");
        bw.close();
	}
    
    
    /**
     * @param weightLeft in [0..1].
     */
    private static final void buildCumulativeBimodal(int nBins, int binLength, double meanLeft, double stdLeft, double meanRight, double stdRight, double weightLeft) {
        int i;
        final double weightRight = 1.0-weightLeft;
        final NormalDistribution normalLeft = new NormalDistribution(meanLeft,stdLeft);
        final NormalDistribution normalRight = new NormalDistribution(meanRight,stdRight);
        cumulativeBimodal = new double[nBins];
        for (i=0; i<nBins; i++) cumulativeBimodal[i] = normalLeft.probability(i*binLength,(i+1)*binLength)*weightLeft + normalRight.probability(i*binLength,(i+1)*binLength)*weightRight;
        for (i=1; i<nBins; i++) cumulativeBimodal[i]+=cumulativeBimodal[i-1];
    }
    
    
    /**
     * Loads from disk a (possibly empty) bin file up until its 
     * $buffers_last$-th record, included.
     */
    private static final void loadBin(int bin, String prefix) throws IOException {
        final int CAPACITY = 1000;  // Arbitrary
        int last;
        String str;
        BufferedReader br;
        
        if (buffers_isLoaded[bin]) {
            queue.remove(bins[bin]);
            bins[bin].priority=buffers_stringLength[bin];  // An approximation
            queue.add(bins[bin]);
            return;
        }
        last=-1; buffers_stringLength[bin]=0;
        br = new BufferedReader(new FileReader(prefix+bin+".bin"));
        str=br.readLine();
        while (str!=null) {
            last++;
            if (last>buffers_last[bin]) break;
            if (buffers[bin]==null) buffers[bin] = new String[CAPACITY];
            else if (last>=buffers[bin].length) {
                String[] newArray = new String[buffers[bin].length+CAPACITY];
                System.arraycopy(buffers[bin],0,newArray,0,buffers[bin].length);
                buffers[bin]=newArray;
            }
            buffers[bin][last]=str;
            buffers_stringLength[bin]+=str.length();
            str=br.readLine();
        }
        br.close();
        buffers_isLoaded[bin]=true;
        if (buffers_last[bin]==Integer.MAX_VALUE) buffers_last[bin]=last;
        buffers_stringLength[bin]>>=1;  // Just an approximation
        nBpsInRam+=buffers_stringLength[bin];
        bins[bin].priority=buffers_stringLength[bin];  // Just an approximation
        queue.add(bins[bin]);
        System.err.println("Allocated bin "+bin+" with "+buffers_stringLength[bin]+" bps");
    }
    
    
    /**
     * Keeps removing a least-recently-used bin until $nBpsInRam$ becomes at
     * most $MAX_BPS_IN_RAM$.
     *
     * @param excludeBin do not deallocate this bin.
     */
    private static final void freeSpace(int excludeBin) {
        boolean someDeallocated;
        int id;
        Bin bin;
        
        someDeallocated=false;
        while (nBpsInRam>maxBpsInRam) {
            bin=queue.poll();
            if (bin==null) break;
            id=bin.id;
            if (id==excludeBin) {
                if (queue.size()==1) break;
                else continue;
            }
            buffers_isLoaded[id]=false;
            nBpsInRam-=buffers_stringLength[id];
            buffers_stringLength[id]=0;
            buffers[id]=null;
            System.err.println("Deallocated bin "+id);
            someDeallocated=true;
        }
        if (someDeallocated) System.gc();
    }
    
    
    private static class Bin implements Comparable {
        public int id;
        public long priority;
        
        public Bin(int id, long priority) {
            this.id=id; this.priority=priority;
        }
        
        public boolean equals(Object other) {
            Bin otherBin = (Bin)other;
            return otherBin.id==id;
        }
        
        public int compareTo(Object other) {
            Bin otherBin = (Bin)other;
            if (priority<otherBin.priority) return -1;
            else if (priority>otherBin.priority) return 1;
            return 0;
        }
    }
    
}