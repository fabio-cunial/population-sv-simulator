import java.util.Arrays;
import java.io.*;


/**
 * Builds the histogram of the number of reads that fall in each read length 
 * bin, and prints all local maxima of the histogram.
 */
public class Fastq2LengthHistogram {

	public static void main(String[] args) throws IOException {
		final String FASTQ_FILE = args[0];
		final int BIN_LENGTH = Integer.parseInt(args[1]);
        final int MAX_READ_LENGTH = Integer.parseInt(args[2]);
        final String OUTPUT_FILE_HISTOGRAM = args[3];
        final String OUTPUT_FILE_MAX = args[4];
        
        final int N_BINS = (MAX_READ_LENGTH+BIN_LENGTH-1)/BIN_LENGTH;
        
        int i;
        int bin, length;
        String str;
        BufferedReader br;
        BufferedWriter bw;
        int[] histogram;
        
        // Building histogram
        histogram = new int[N_BINS];
        Arrays.fill(histogram,0);
        br = new BufferedReader(new FileReader(FASTQ_FILE));
        str=br.readLine();
        while (str!=null) {
            str=br.readLine(); length=str.length();
            str=br.readLine(); str=br.readLine();
            bin=length/BIN_LENGTH;
            if (bin>=N_BINS) bin=N_BINS-1;
            histogram[bin]++;
            str=br.readLine();
        }
        br.close();
        
        // Finding local maxima
        bw = new BufferedWriter(new FileWriter(OUTPUT_FILE_MAX));
        for (i=0; i<N_BINS; i++) {
            if ((i==0 || histogram[i]>histogram[i-1]) && (i==N_BINS-1 || histogram[i]>histogram[i+1])) bw.write((i*BIN_LENGTH)+","+((i+1)*BIN_LENGTH-1)+"\n");
        }
        bw.close();
	}

}