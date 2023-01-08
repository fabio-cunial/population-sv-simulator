import java.util.Arrays;
import java.util.PriorityQueue;
import java.util.Random;
import java.io.*;
import org.apache.commons.math3.distribution.NormalDistribution;


/**
 * Given the length bin files created by $BuildReadLengthBins$, the program 
 * samples reads using a bimodal normal distribution specified in input.
 *
 * If the desired coverage cannot be achieved, the program returns error code 3.
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
    private static Fastq[][] buffers;
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
     *  4: weight of the left normal distribution, in [0..1];
     * 10: max characters in RAM (in billions); the actual RAM usage is at least
     *     twice this much, because of the quality track in a FASTQ file.
     */
	public static void main(String[] args) throws IOException {
        final double MEAN_LEFT = Double.parseDouble(args[0]);
        final double STD_LEFT = Double.parseDouble(args[1]);
        final double MEAN_RIGHT = Double.parseDouble(args[2]);
        final double STD_RIGHT = Double.parseDouble(args[3]);
        final double WEIGHT_LEFT = Double.parseDouble(args[4]);
        final int BIN_LENGTH = Integer.parseInt(args[5]);
        final int MAX_READ_LENGTH = Integer.parseInt(args[6]);  // For a bin
		final String BINS_PREFIX = args[7];
        final long GENOME_LENGTH_ONE_HAPLOTYPE = Long.parseLong(args[8]);
        final double TARGET_COVERAGE_ONE_HAPLOTYPE = Double.parseDouble(args[9]);
        maxBpsInRam=Long.parseLong(args[10])*1000000000;
        final String OUTPUT_FASTQ_FILE = args[11];
        
        final int N_BINS = (MAX_READ_LENGTH+BIN_LENGTH-1)/BIN_LENGTH;
        final long TARGET_COVERAGE_BP = (long)(GENOME_LENGTH_ONE_HAPLOTYPE*(TARGET_COVERAGE_ONE_HAPLOTYPE*2));
        final int MAX_SAMPLING_ATTEMPTS = N_BINS*10;  // Arbitrary
        
        int i, j;
        int bin, length;
        long coverage, nBpsInRam;
        BufferedWriter bw;
        int[] histogram;
        
        random = new Random();
        buildCumulativeBimodal(N_BINS,BIN_LENGTH,MEAN_LEFT,STD_LEFT,MEAN_RIGHT,STD_RIGHT,WEIGHT_LEFT);
        buffers = new Fastq[N_BINS][0];
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
        bw = new BufferedWriter(new FileWriter(OUTPUT_FASTQ_FILE));
        while (coverage<TARGET_COVERAGE_BP) {
 System.err.println("-->0");
            j=0; bin=-1;
            while (j<MAX_SAMPLING_ATTEMPTS) {
                bin=Arrays.binarySearch(cumulativeBimodal,random.nextDouble());
                if (bin<0) bin=-1-bin;
                if (bin==N_BINS) bin=N_BINS-1;
                j++;
                if (buffers_last[bin]!=-1) {  // The bin cannot be empty
                    try { loadBin(bin,BINS_PREFIX); }
                    catch (Error e) {
                        System.err.println("loadBin() threw the following error:");
                        e.printStackTrace();
                        System.exit(1);
                    }
                    if (buffers_last[bin]!=-1) {  // The bin cannot be empty
                        freeSpace(bin);
                        break;
                    }
                }
            }
 System.err.println("-->1  j="+j+" MAX_SAMPLING_ATTEMPTS="+MAX_SAMPLING_ATTEMPTS);
            if (j==MAX_SAMPLING_ATTEMPTS) {
                System.err.println("Not enough reads to sample "+TARGET_COVERAGE_ONE_HAPLOTYPE+"x haploid coverage from the given bins and distribution.");
                System.exit(3);
            }
System.err.println("-->2  bin="+bin);
            histogram[bin]++;
System.err.println("-->3  buffers_last[bin]="+buffers_last[bin]);
            buffers[bin][buffers_last[bin]].writeTo(bw);
System.err.println("-->4  buffers_last[bin]="+buffers_last[bin]);
            length=buffers[bin][buffers_last[bin]].sequence.length();
            coverage+=length;
            buffers[bin][buffers_last[bin]]=null;
System.err.println("-->5 bin="+bin);
            buffers_stringLength[bin]-=length;
            buffers_last[bin]--;
System.err.println("-->6  buffers_last[bin]="+buffers_last[bin]+" coverage="+coverage);
        }
System.err.println("-->6.5");        
        bw.close();
System.err.println("-->7");
        bw = new BufferedWriter(new FileWriter(OUTPUT_FASTQ_FILE+".histogram"));
System.err.println("-->8");
        for (i=0; i<N_BINS; i++) bw.write(histogram[i]+"\n");
System.err.println("-->9");
        bw.close();
System.err.println("-->10");
	}
    
    
    /**
     * @param weightLeft in [0..1].
     */
    private static final void buildCumulativeBimodal(int nBins, int binLength, double meanLeft, double stdLeft, double meanRight, double stdRight, double weightLeft) {
        int i;
        double sum;
        final double weightRight = 1.0-weightLeft;
        final NormalDistribution normalLeft = new NormalDistribution(meanLeft,stdLeft);
        final NormalDistribution normalRight = new NormalDistribution(meanRight,stdRight);
        cumulativeBimodal = new double[nBins];
        for (i=0; i<nBins; i++) cumulativeBimodal[i] = normalLeft.probability(i*binLength,(i+1)*binLength)*weightLeft + normalRight.probability(i*binLength,(i+1)*binLength)*weightRight;
        sum=0.0;
        for (i=0; i<nBins; i++) sum+=cumulativeBimodal[i];
        for (i=0; i<nBins; i++) cumulativeBimodal[i]/=sum;
        System.err.println("Sampling from the following distribution:");
        for (i=0; i<nBins; i++) System.err.println((i*binLength)+","+cumulativeBimodal[i]);
        for (i=1; i<nBins; i++) cumulativeBimodal[i]+=cumulativeBimodal[i-1];
    }
    
    
    /**
     * Loads from disk a (possibly empty) bin file up until its 
     * $buffers_last$-th row, included.
     */
    private static final void loadBin(int bin, String prefix) throws IOException {
        final int CAPACITY = 1000;  // Arbitrary
        int last;
        String header, sequence, separator, quality;
        String str;
        BufferedReader br;
        
 System.err.println("-->loadBin 0");
        if (buffers_isLoaded[bin]) {
 System.err.println("-->loadBin 1");
            queue.remove(bins[bin]);
            bins[bin].priority=buffers_last[bin];  // An approximation
            queue.add(bins[bin]);
            return;
        }
 System.err.println("-->loadBin 2");
        last=-1; buffers_stringLength[bin]=0;
        br = new BufferedReader(new FileReader(prefix+bin+".bin"));
        header=br.readLine();
 System.err.println("-->loadBin 3");        
        while (header!=null) {
            sequence=br.readLine(); separator=br.readLine(); quality=br.readLine();
            last++;
            if (last>buffers_last[bin]) break;
            if (buffers[bin]==null) buffers[bin] = new Fastq[CAPACITY];
            else if (last>=buffers[bin].length) {
 System.err.println("-->loadBin 4");
                Fastq[] newArray = new Fastq[buffers[bin].length+CAPACITY];
                System.arraycopy(buffers[bin],0,newArray,0,buffers[bin].length);
                buffers[bin]=newArray;
 System.err.println("-->loadBin 5");
                freeSpace(bin);
 System.err.println("-->loadBin 6");
                System.gc();
 System.err.println("-->loadBin 7");
            }
 System.err.println("-->loadBin 8");
            buffers[bin][last] = new Fastq(header,sequence,separator,quality);
            buffers_stringLength[bin]+=sequence.length();
            header=br.readLine();
 System.err.println("-->loadBin 9");
        }
 System.err.println("-->loadBin 10");
        br.close();
 System.err.println("-->loadBin 11");
        buffers_isLoaded[bin]=true;
        if (buffers_last[bin]==Integer.MAX_VALUE) buffers_last[bin]=last;
        buffers_stringLength[bin]>>=1;  // Just an approximation
        nBpsInRam+=buffers_stringLength[bin];
        bins[bin].priority=buffers_last[bin];
 System.err.println("-->loadBin 12");
        queue.add(bins[bin]);
        System.err.println("Loaded bin "+bin+" with ~"+buffers_stringLength[bin]+" bps");
    }
    
    
    /**
     * Keeps removing the lowest-priority bin until $nBpsInRam$ becomes at
     * most $maxBpsInRam$.
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
                if (queue.size()==0) break;
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