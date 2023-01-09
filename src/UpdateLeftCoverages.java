import java.io.*;

/**
 * Given a string of coverages separated by "-", the procedure removes all the
 * values greater than X and makes sure X is in the list.
 */
public class UpdateLeftCoverages {

	public static void main(String[] args) {
		String INPUT_COVERAGES = args[0];
        Double NEW_MAX_COVERAGE = Double.parseDouble(args[1]);
        
        int i;
        double value;
        String out;
        String[] tokens;
        
        out="";
        tokens=INPUT_COVERAGES.split("-");
        for (i=0; i<tokens.length; i++) {
            value=Double.parseDouble(tokens[i]);
            if (value<NEW_MAX_COVERAGE) {
                if (out.length()==0) out=tokens[i];
                else out+="-"+tokens[i];
            }
        }
		out+="-"+NEW_MAX_COVERAGE;
        System.out.println(out);
	}

}