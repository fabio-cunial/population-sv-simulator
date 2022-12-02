import java.util.Arrays;
import java.io.*;


/**
 * Given a set of FASTQ files, partitions their reads into bins, permutes every 
 * bin randomly, and stores every bin in a separate file.
 */
public class BuildReadLengthBuckets {

	public static void main(String[] args) throws IOException {
		final String INPUT_FILES_LIST = args[0];
        final int BIN_LENGTH = Integer.parseInt(args[1]);
        final int MAX_READ_LENGTH = Integer.parseInt(args[2]);
        final String OUTPUT_PREFIX = args[3];
        
        final int N_BINS = MAX_READ_LENGTH/BIN_LENGTH+1;
        
        int i;
        int bin;
        String file, header, sequence, separator, quality;
        BufferedReader br1, br2;
        BufferedWriter bw;
        int[] bufferLengths;
        BufferedWriter[] buffers;
        
        buffers = new BufferedWriter[N_BINS];
        for (i=0; i<N_BINS; i++) buffers[i] = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+i+".fastq"));
        bufferLengths = new int[N_BINS];
        Arrays.fill(bufferLengths,0);
        br1 = new BufferedReader(new FileReader(INPUT_FILES_LIST));
        file=br1.readLine();
        while (file!=null) {
            br2 = new BufferedReader(new FileReader(file));
            header=br2.readLine();
            while (str2!=null) {
                sequence=br2.readLine(); separator=br2.readLine(); quality=br2.readLine();
                bin=sequence.length()/BIN_LENGTH;
                if (bin>=N_BINS) bin=N_BINS-1;
                bufferLengths[bin]++;
                buffers[bin].write(header); buffers[bin].newLine();
                buffers[bin].write(sequence); buffers[bin].newLine("\n");
                buffers[bin].write(separator); buffers[bin].newLine("\n");
                buffers[bin].write(quality); buffers[bin].newLine("\n");
                header=br2.readLine();
            }
            br2.close();
            file=br1.readLine();
        }
        br1.close();
        for (i=0; i<N_BINS; i++) buffers[i].close();
        bw = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+"histogram.txt"));
        for (i=0; i<N_BINS; i++) bw.write(bufferLengths[i]);
        bw.close();
	}
    
    ----------->
    private static final void shuffle(int bin) {
        int i, j;
        int tmp;
        
        for (i=buffers_last[bin]; i>0; i--) {
            j=random.nextInt(0,i+1);
            tmp=buffers[bin][j];
            buffers[bin][j]=buffers[bin][i];
            buffers[bin][i]=tmp;
        }
    }

}