import java.io.IOException;


public class PrintHaplotypes {

	public static void main(String[] args) throws IOException {
		final int FIRST_HAPLOTYPE = Integer.parseInt(args[0]);
		final int LAST_HAPLOTYPE = Integer.parseInt(args[1]);
		final String REFERENCE_FILE = args[2];
		final String HAPLOTYPE2VARIANTS_FILE = args[3];
		final String VARIANTS_FILE = args[4];
		final String OUTPUT_DIR = args[5];
		
		SimulatePopulation.loadReference(REFERENCE_FILE);
		SimulatePopulation.deserialize(HAPLOTYPE2VARIANTS_FILE,VARIANTS_FILE);
		StringBuilder buffer = new StringBuilder(SimulatePopulation.referenceLength);
		SimulatePopulation.buildHaplotypes(FIRST_HAPLOTYPE,LAST_HAPLOTYPE,OUTPUT_DIR,buffer);	
	}

}