import java.io.BufferedWriter;
import java.io.IOException;


public class Fastq {
    public String header, sequence, separator, quality;
    
    public Fastq(String h, String seq, String sep, String q) {
        header=h; sequence=seq; separator=sep; quality=q;
    }
    
    public void writeTo(BufferedWriter bw) throws IOException {
        bw.write(header); bw.newLine();
        bw.write(sequence); bw.newLine();
        bw.write(separator); bw.newLine();
        bw.write(quality); bw.newLine();
    }
}