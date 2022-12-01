import java.io.IOException;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.BufferedWriter;
import java.io.FileWriter;


/**
 * Stores to a file the IDs of all reads that fall in each length bucket.
 */
public class Fastq2LengthBuckets {

	public static void main(String[] args) throws IOException {
		final String FASTQ_FILE = args[0];
		final int BIN_LENGTH = Integer.parseInt(args[1]);   // 500 bp say
        final int MAX_READ_LENGTH = Integer.parseInt(args[2]);  // 30kbp say
        final String OUTPUT_FILE = args[3];
        
        final int N_BINS = MAX_READ_LENGTH/BIN_LENGTH+1;
        final int CAPACITY = 1000;  // Arbitrary
        
        int i, j;
        int bin, length, last;
        String str;
        BufferedReader br;
        BufferedWriter bw;
        int[] ids_last;
        int[][] ids;
        
        // Building histogram
        ids_last = new int[N_BINS];
        for (i=0; i<N_BINS; i++) ids_last[i]=-1;
        ids = new int[N_BINS][CAPACITY];
        br = new BufferedReader(new FileReader(FASTQ_FILE));
        str=br.readLine(); i=0;
        while (str!=null) {
            str=br.readLine(); length=str.length();
            str=br.readLine(); str=br.readLine();
            bin=length/BIN_LENGTH;
            ids_last[bin]++;
            if (ids_last[bin]==ids[bin].length) {
                int[] newArray = new int[ids[bin].length<<1];
                System.arraycopy(ids[bin],0,newArray,0,ids[bin].length);
                ids[bin]=newArray;
            }
            ids[bin][ids_last[bin]]=i;
            str=br.readLine(); i++;
        }
        br.close();
        
        // Serializing
        bw = new BufferedWriter(new FileWriter(OUTPUT_FILE));
        for (i=0; i<N_BINS; i++) {
            last=ids_last[i];
            bw.write((last+1)+"");
            for (j=0; j<=last; j++) bw.write(","+ids[i][j]);
            bw.newLine();
        }
        bw.close();
	}

}