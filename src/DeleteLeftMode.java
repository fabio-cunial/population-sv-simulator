import java.util.Arrays;
import java.io.*;


/**
 * Given the length bin files created by $BuildReadLengthBins$, the program
 * tries to enforce the histogram to the left of the rightmost peak to be 
 * symmetrical to the histogram to the right of the rightmost peak.
 */
public class DeleteLeftMode {
    
    /**
     * @param args 
     */
	public static void main(String[] args) throws IOException {
        final double MEAN_RIGHT = Double.parseDouble(args[0]);
        final String BINS_PREFIX = args[1];
        final int BIN_LENGTH = Integer.parseInt(args[2]);
        final int MAX_READ_LENGTH = Integer.parseInt(args[3]);  // For a bin
        final String OUTPUT_FASTQ_FILE = args[4];
        
        final int RIGHT_BIN = (int)(MEAN_RIGHT/BIN_LENGTH);
        final int N_BINS = (MAX_READ_LENGTH+BIN_LENGTH-1)/BIN_LENGTH;
        
        int i, j;
        String str;
        BufferedReader br;
        BufferedWriter bw;
        int[] histogram;
        
        histogram = new int[N_BINS];
        Arrays.fill(histogram,0);
        bw = new BufferedWriter(new FileWriter(OUTPUT_FASTQ_FILE));
        for (i=N_BINS-1; i>=RIGHT_BIN; i--) {
            br = new BufferedReader(new FileReader(BINS_PREFIX+i+".bin"));
            str=br.readLine();
            while (str!=null) {
                histogram[i]++;
                bw.write(str); bw.newLine();
                str=br.readLine(); bw.write(str); bw.newLine();
                str=br.readLine(); bw.write(str); bw.newLine();
                str=br.readLine(); bw.write(str); bw.newLine();
                str=br.readLine();
            }
            br.close();
        }
        for (i=RIGHT_BIN-1; i>=0; i--) {
            br = new BufferedReader(new FileReader(BINS_PREFIX+i+".bin"));
            str=br.readLine();
            while (str!=null) {
                br.readLine(); br.readLine(); br.readLine();
                histogram[i]++;
                str=br.readLine();
            }
            br.close();
            j=RIGHT_BIN+(RIGHT_BIN-i);
            if (histogram[i]>histogram[j]) histogram[i]=histogram[j];
            if (histogram[i]==0) continue;
            br = new BufferedReader(new FileReader(BINS_PREFIX+i+".bin"));
            for (j=0; j<histogram[i]; j++) {
                str=br.readLine(); bw.write(str); bw.newLine();
                str=br.readLine(); bw.write(str); bw.newLine();
                str=br.readLine(); bw.write(str); bw.newLine();
                str=br.readLine(); bw.write(str); bw.newLine();
            }
            br.close();
        }
        bw.close();
	}
    
}