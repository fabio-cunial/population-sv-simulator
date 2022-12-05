import java.util.Arrays;
import java.util.Random;
import java.io.*;


/**
 * Given a set of FASTQ files, partitions all their reads into bins, permutes 
 * every bin randomly, and stores every bin to a separate file. Every bin file
 * contains one line per read, which is the concatenation of the four rows of a
 * FASTQ read.
 */
public class BuildReadLengthBins {

    public static final String BIN_FILE_SEPARATOR = "|";


	public static void main(String[] args) throws IOException {
		final String INPUT_FILES_LIST = args[0];
        final int BIN_LENGTH = Integer.parseInt(args[1]);
        final int MAX_READ_LENGTH = Integer.parseInt(args[2]);  // For a bin
        final String OUTPUT_PREFIX = args[3];
        
        final int N_BINS = MAX_READ_LENGTH/BIN_LENGTH+1;
        
        int i, j;
        int bin, last, length, maxLength;
        String str, file, header, sequence, separator, quality;
        Random random = new Random();
        BufferedReader br1, br2;
        BufferedWriter bw;
        int[] bufferLengths;
        String[] tmpArray, tokens;
        BufferedWriter[] buffers;
        
        // Building buckets
        buffers = new BufferedWriter[N_BINS];
        for (i=0; i<N_BINS; i++) buffers[i] = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+i+".bin"));
        bufferLengths = new int[N_BINS];
        Arrays.fill(bufferLengths,0);
        br1 = new BufferedReader(new FileReader(INPUT_FILES_LIST));
        file=br1.readLine();
        while (file!=null) {
            br2 = new BufferedReader(new FileReader(file));
            header=br2.readLine();
            while (header!=null) {
                sequence=br2.readLine(); separator=br2.readLine(); quality=br2.readLine();
                bin=sequence.length()/BIN_LENGTH;
                if (bin>=N_BINS) bin=N_BINS-1;
                bufferLengths[bin]++;
                buffers[bin].write(header); buffers[bin].write(BIN_FILE_SEPARATOR);
                buffers[bin].write(sequence); buffers[bin].write(BIN_FILE_SEPARATOR);
                buffers[bin].write(separator); buffers[bin].write(BIN_FILE_SEPARATOR);
                buffers[bin].write(quality); buffers[bin].newLine();
                header=br2.readLine();
            }
            br2.close();
            file=br1.readLine();
        }
        br1.close();
        for (i=0; i<N_BINS; i++) buffers[i].close();
        bw = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+".histogram"));
        maxLength=0;
        for (i=0; i<N_BINS; i++) {
            length=bufferLengths[i];
            if (length>maxLength) maxLength=length;
            bw.write(length+"\n");
        }
        bw.close();
        
        // Permuting reads in each bucket
        tmpArray = new String[maxLength];
        for (i=0; i<N_BINS; i++) {
            br1 = new BufferedReader(new FileReader(OUTPUT_PREFIX+i+".bin"));
            str=br1.readLine(); last=-1;
            while (str!=null) {
                tmpArray[++last]=str;
                str=br1.readLine();
            }
            br1.close();
            shuffle(tmpArray,last,random);
            bw = new BufferedWriter(new FileWriter(OUTPUT_PREFIX+i+".bin"));
            for (j=0; j<=last; j++) {
                bw.write(tmpArray[j]); bw.newLine();
            }
            bw.close();
        }
	}
    
    
    private static final void shuffle(String[] array, int last, Random random) {
        int i, j;
        String tmp;
        
        for (i=last; i>0; i--) {
            j=random.nextInt(i+1);
            tmp=array[j]; array[j]=array[i]; array[i]=tmp;
        }
    }

}