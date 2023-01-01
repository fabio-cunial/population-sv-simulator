import java.util.Arrays;
import java.util.Random;
import java.io.*;


/**
 * Given a list of FASTQ files (for example the flowcells of a sample), the
 * program partitions the union of all their reads into length bins, permutes
 * every bin randomly, and stores every bin to a separate file. Every bin file
 * contains one line per read, which is the concatenation of the four rows of a
 * FASTQ read.
 */
public class BuildReadLengthBins {

    public static final String BIN_FILE_SEPARATOR = "__________";
    private static int BIN_LENGTH, HALF_BIN_LENGTH;
    private static int N_BINS;
    private static final int MIN_READS_IN_LOCAL_MAXIMUM = 500;  // Arbitrary
    private static long GENOME_LENGTH_HAPLOID;


	public static void main(String[] args) throws IOException {
		final String INPUT_FILES_LIST = args[0];
        BIN_LENGTH=Integer.parseInt(args[1]);
        final int MAX_READ_LENGTH = Integer.parseInt(args[2]);  // For a bin
        GENOME_LENGTH_HAPLOID=Long.parseLong(args[3]);
        final String OUTPUT_PREFIX = args[4];
        
        N_BINS=(MAX_READ_LENGTH+BIN_LENGTH-1)/BIN_LENGTH;
        HALF_BIN_LENGTH=BIN_LENGTH>>1;
        
        int i, j;
        int bin, last, length, maxLength, nReads;
        String str, file, header, sequence, separator, quality;
        Random random = new Random();
        BufferedReader br1, br2;
        BufferedWriter bw;
        int[] buffers_last;
        long[] buffers_stringLength;
        String[] tmpArray, tokens;
        BufferedWriter[] buffers;
        
        // Building bins
        buffers = new BufferedWriter[N_BINS];
        for (i=0; i<N_BINS; i++) buffers[i] = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+i+".bin"));
        buffers_last = new int[N_BINS];
        Arrays.fill(buffers_last,-1);
        buffers_stringLength = new long[N_BINS];
        Arrays.fill(buffers_stringLength,0);
        br1 = new BufferedReader(new FileReader(INPUT_FILES_LIST));
        file=br1.readLine();
        while (file!=null) {
            System.err.println("Loading file <"+file+">...");
            br2 = new BufferedReader(new FileReader(file));
            header=br2.readLine(); nReads=0;
            while (header!=null) {
                sequence=br2.readLine(); separator=br2.readLine(); quality=br2.readLine();
                bin=sequence.length()/BIN_LENGTH;
                if (bin>=N_BINS) bin=N_BINS-1;
                buffers_last[bin]++;
                buffers[bin].write(header); buffers[bin].write(BIN_FILE_SEPARATOR);
                buffers[bin].write(sequence); buffers[bin].write(BIN_FILE_SEPARATOR);
                buffers_stringLength[bin]+=sequence.length();
                buffers[bin].write(separator); buffers[bin].write(BIN_FILE_SEPARATOR);
                buffers[bin].write(quality); buffers[bin].newLine();
                nReads++;
                if (nReads%10000==0) System.err.println(nReads+" reads");
                header=br2.readLine();
            }
            br2.close();
            System.err.println("Done");
            file=br1.readLine();
        }
        br1.close();
        for (i=0; i<N_BINS; i++) buffers[i].close();
        bw = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+".histogram"));
        maxLength=0;
        for (i=0; i<N_BINS; i++) {
            length=buffers_last[i]+1;
            if (length>maxLength) maxLength=length;
            bw.write((i*BIN_LENGTH+HALF_BIN_LENGTH)+","+length+"\n");
        }
        bw.close();
        printStats(buffers_last,buffers_stringLength,OUTPUT_PREFIX+".stats");
        
        // Permuting reads in each bin
        System.err.println("Permuting bins...");
        tmpArray = new String[maxLength];
        for (i=0; i<N_BINS; i++) {
            br1 = new BufferedReader(new FileReader(OUTPUT_PREFIX+i+".bin"));
            str=br1.readLine(); last=-1;
            while (str!=null) {
                tmpArray[++last]=str;
                str=br1.readLine();
            }
            br1.close();
            if (last<=1) continue;
            shuffle(tmpArray,last,random);
            bw = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+i+".bin"));
            for (j=0; j<=last; j++) {
                bw.write(tmpArray[j]); bw.newLine();
            }
            bw.close();
        }
        System.err.println("Done");
	}
    
    
    private static final void shuffle(String[] array, int last, Random random) {
        int i, j;
        String tmp;
        
        for (i=last; i>0; i--) {
            j=random.nextInt(i+1);
            tmp=array[j]; array[j]=array[i]; array[i]=tmp;
        }
    }
    
    
    /**
     * Prints the following values to the rows of $outputFile$, one per row:
     *
     * nMaxima  // If < 2 the following lines are not printed
     * meanLeft  // Last-but-one peak
     * stdLeft
     * meanRight  // Last peak
     * stdRight
     * leftMinBin  // Local minimum between the last-but-one and the previous
     *             // max (if any). 0 otherwise.
     * leftRightMinBin  // Local minimum bin between the last two local maxima
     * massLeft  // Probability mass in [leftMinBin..leftRightMinBin]
     * massRight  // Probability mass in [leftRightMinBin+1..]
     * coverageLeft  // Of one haplotype. In [leftMinBin..leftRightMinBin].
     * coverageRight  // Of one haplotype. In [leftRightMinBin+1..].
     *
     * Remark: mean and std are just crude approximations, and should instead be
     * computed by fitting a mixture of Gaussians.
     */
    private static final void printStats(int[] buffers_last, long[] buffers_stringLength, String outputFile) throws IOException {
        int i;
        int nMaxima, lastMaximum, previousMaximum, leftMinBin, leftRightMinBin;
        double sum, mass;
        BufferedWriter bw;
        
        bw = new BufferedWriter(new FileWriter(outputFile));
        nMaxima=0; lastMaximum=-1; previousMaximum=-1;
        for (i=0; i<N_BINS; i++) {
            if (buffers_last[i]+1>=MIN_READS_IN_LOCAL_MAXIMUM && (i==0 || buffers_last[i]>buffers_last[i-1]) && (i==N_BINS-1 || buffers_last[i]>buffers_last[i+1])) {
                nMaxima++;
                previousMaximum=lastMaximum; lastMaximum=i;
            }
        }
        bw.write(nMaxima+"\n");
        if (nMaxima<2) {
            bw.close();
            return;
        }
        leftRightMinBin=-1;
        for (i=previousMaximum+1; i<lastMaximum; i++) {
            if (buffers_last[i]<buffers_last[i-1] && buffers_last[i]<buffers_last[i+1]) {
                leftRightMinBin=i;
                break;
            }
        }
        leftMinBin=0;
        for (i=1; i<previousMaximum; i++) {
            if (buffers_last[i]<buffers_last[i-1] && buffers_last[i]<buffers_last[i+1]) {
                leftMinBin=i;
                break;
            }
        }
        bw.write((previousMaximum*BIN_LENGTH+HALF_BIN_LENGTH)+"\n");
        bw.write(std(leftMinBin,previousMaximum-1,previousMaximum,buffers_last)+"\n");
        bw.write((lastMaximum*BIN_LENGTH+HALF_BIN_LENGTH)+"\n");
        bw.write(std(lastMaximum+1,N_BINS-1,lastMaximum,buffers_last)+"\n");
        bw.write(leftMinBin+"\n");
        bw.write(leftRightMinBin+"\n");
        sum=0.0;
        for (i=0; i<N_BINS; i++) sum+=buffers_last[i]+1;
        mass=0.0;
        for (i=leftMinBin; i<=leftRightMinBin; i++) mass+=buffers_last[i]+1;
        mass/=sum;
        bw.write(mass+"\n");
        mass=0.0;
        for (i=leftRightMinBin+1; i<N_BINS; i++) mass+=buffers_last[i]+1;
        mass/=sum;
        bw.write(mass+"\n");
        sum=0.0;
        for (i=leftMinBin; i<=leftRightMinBin; i++) sum+=buffers_stringLength[i];
        bw.write((sum/((GENOME_LENGTH_HAPLOID)<<1))+"\n");
        sum=0.0;
        for (i=leftRightMinBin+1; i<N_BINS; i++) sum+=buffers_stringLength[i];
        bw.write((sum/((GENOME_LENGTH_HAPLOID)<<1))+"\n");
        bw.close();
    }
    
    
    /**
     * Computes the standard deviation of the normal distribution centered at 
     * $bin$ using only the values in $[from..to]$, which are assumed to be on 
     * one side of a symmetrical distribution.
     *
     * Remark: this is just a crude approximation. We should instead fit a
     * mixture of Gaussians to the histogram.
     */
    private static final double std(int from, int to, int bin, int[] buffers_last) {
        int i, k, n;
        final int mean = bin*BIN_LENGTH+HALF_BIN_LENGTH;
        double std;
        
        std=0.0; n=0;
        for (i=from; i<=to; i++) {
            n+=buffers_last[i]+1;
            k=i*BIN_LENGTH+HALF_BIN_LENGTH-mean;
            std+=2.0*(buffers_last[i]+1)*k*k;
        }
        n+=buffers_last[bin];
        return Math.sqrt(std/n);
    }

}