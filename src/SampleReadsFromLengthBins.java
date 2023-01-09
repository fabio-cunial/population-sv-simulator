import java.util.Arrays;
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
     * @param args 
     *  4: weight of the left normal distribution, in [0..1].
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
        final String OUTPUT_FASTQ_FILE = args[10];
        
        final int N_BINS = (MAX_READ_LENGTH+BIN_LENGTH-1)/BIN_LENGTH;
        final long TARGET_COVERAGE_BP = (long)(GENOME_LENGTH_ONE_HAPLOTYPE*(TARGET_COVERAGE_ONE_HAPLOTYPE*2));
        final int MAX_SAMPLING_ATTEMPTS = N_BINS*10;  // Arbitrary
        
        int i, j, c;
        int bin, length;
        long coverage;
        String str;
        Random random = new Random();
        BufferedWriter bw;
        int[] histogram;
        double[] cumulativeBimodal;
        BufferedReader[] bins;
        
        cumulativeBimodal=buildCumulativeBimodal(N_BINS,BIN_LENGTH,MEAN_LEFT,STD_LEFT,MEAN_RIGHT,STD_RIGHT,WEIGHT_LEFT);
        bins = new BufferedReader[N_BINS];
        for (i=0; i<N_BINS; i++) bins[i] = new BufferedReader(new FileReader(BINS_PREFIX+i+".bin"));
        histogram = new int[N_BINS];
        Arrays.fill(histogram,0);
        coverage=0;
        bw = new BufferedWriter(new FileWriter(OUTPUT_FASTQ_FILE));
        while (coverage<TARGET_COVERAGE_BP) {
            j=0; bin=-1;
            while (j<MAX_SAMPLING_ATTEMPTS) {
                bin=Arrays.binarySearch(cumulativeBimodal,random.nextDouble());
                if (bin<0) bin=-1-bin;
                else if (bin==N_BINS) bin=N_BINS-1;
                bins[bin].mark(10); c=bins[bin].read(); bins[bin].reset();
                if (c!=-1) break;
                j++;
            }
            if (j==MAX_SAMPLING_ATTEMPTS) {
                System.err.println("Not enough reads to sample "+TARGET_COVERAGE_ONE_HAPLOTYPE+"x haploid coverage from the given bins and distribution.");
                System.exit(3);
            }
            histogram[bin]++;
            bw.write(bins[bin].readLine());
            str=bins[bin].readLine(); coverage+=str.length();
            bw.write(str);
            bw.write(bins[bin].readLine());
            bw.write(bins[bin].readLine());
        }
        for (i=0; i<N_BINS; i++) { bins[i].close(); bins[i]=null; }
        bw.close();
        bw = new BufferedWriter(new FileWriter(OUTPUT_FASTQ_FILE+".histogram"));
        for (i=0; i<N_BINS; i++) bw.write(histogram[i]+"\n");
        bw.close();
	}
    
    
    /**
     * @param weightLeft in [0..1].
     */
    private static final double[] buildCumulativeBimodal(int nBins, int binLength, double meanLeft, double stdLeft, double meanRight, double stdRight, double weightLeft) {
        int i;
        double sum;
        final double weightRight = 1.0-weightLeft;
        final NormalDistribution normalLeft = new NormalDistribution(meanLeft,stdLeft);
        final NormalDistribution normalRight = new NormalDistribution(meanRight,stdRight);
        double[] cumulativeBimodal = new double[nBins];
        for (i=0; i<nBins; i++) cumulativeBimodal[i] = normalLeft.probability(i*binLength,(i+1)*binLength)*weightLeft + normalRight.probability(i*binLength,(i+1)*binLength)*weightRight;
        sum=0.0;
        for (i=0; i<nBins; i++) sum+=cumulativeBimodal[i];
        for (i=0; i<nBins; i++) cumulativeBimodal[i]/=sum;
        System.err.println("Sampling from the following distribution:");
        for (i=0; i<nBins; i++) System.err.println((i*binLength)+","+cumulativeBimodal[i]);
        for (i=1; i<nBins; i++) cumulativeBimodal[i]+=cumulativeBimodal[i-1];
        return cumulativeBimodal;
    }
    
}