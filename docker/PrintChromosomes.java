import java.io.IOException;


public class PrintChromosomes {

	public static void main(String[] args) throws IOException {
		final int FIRST_CHROMOSOME = Integer.parseInt(args[0]);
		final int LAST_CHROMOSOME = Integer.parseInt(args[1]);
		final String REFERENCE_FILE = args[2];
		final String CHROMOSOME2VARIANTS_FILE = args[3];
		final String VARIANTS_FILE = args[4];
		final String OUTPUT_DIR = args[5];
		
		SimulatePopulation.loadReference(REFERENCE_FILE);
		SimulatePopulation.deserialize(CHROMOSOME2VARIANTS_FILE,VARIANTS_FILE);
		StringBuilder buffer = new StringBuilder(SimulatePopulation.referenceLength);
		SimulatePopulation.buildChromosomes(FIRST_CHROMOSOME,LAST_CHROMOSOME,OUTPUT_DIR,buffer);	
	}

}