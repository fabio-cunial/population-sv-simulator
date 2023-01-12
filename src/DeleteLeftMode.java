import java.util.Arrays;
import java.io.*;


/**
 * Given the length bin files created by $BuildReadLengthBins$, the program
 * tries to enforce the histogram to the left of the rightmost peak to be 
 * symmetrical to the histogram to the right of the rightmost peak.
 */
public class DeleteLeftMode {

	public static void main(String[] args) throws IOException {
        final double MEAN_RIGHT = Double.parseDouble(args[0]);
        final String BINS_PREFIX = args[1];
        final int BIN_LENGTH = Integer.parseInt(args[2]);
        final int MAX_READ_LENGTH = Integer.parseInt(args[3]);  // For a bin
        final String OUTPUT_FASTQ_FILE = args[4];
        
        final int RIGHT_BIN = (int)(MEAN_RIGHT/BIN_LENGTH);
        final int N_BINS = (MAX_READ_LENGTH+BIN_LENGTH-1)/BIN_LENGTH;
        final String READ_PREFIX = "@read";
        
        int i, j;
        int read;
        String str;
        BufferedReader br;
        BufferedWriter bw;
        int[] histogram;
        
        histogram = new int[N_BINS];
        Arrays.fill(histogram,0);
        read=-1;
        bw = new BufferedWriter(new FileWriter(OUTPUT_FASTQ_FILE));
        for (i=N_BINS-1; i>=RIGHT_BIN; i--) {
            br = new BufferedReader(new FileReader(BINS_PREFIX+i+".bin"));
            str=br.readLine();
            while (str!=null) {
                histogram[i]++;
                bw.write(READ_PREFIX+(++read)); bw.newLine();  // Replacing the original name with a short one, to avoid samtools' "query name too long" error.
                bw.write(br.readLine()); bw.newLine();
                bw.write(br.readLine()); bw.newLine();
                bw.write(br.readLine()); bw.newLine();
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
                br.readLine(); 
                bw.write(READ_PREFIX+(++read)); bw.newLine();  // Replacing the original name with a short one, to avoid samtools' "query name too long" error.
                bw.write(br.readLine()); bw.newLine();
                bw.write(br.readLine()); bw.newLine();
                bw.write(br.readLine()); bw.newLine();
            }
            br.close();
        }
        bw.close();
	}
    
}