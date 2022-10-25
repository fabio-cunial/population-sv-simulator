import java.io.IOException;


/**
 * Creates the <repeatIntervals.txt> and <segdups.txt> files of the model, for a
 * genome that is different from the one used by gnomAD-SV.
 */
public class BuildModelRepeats {

	public static void main(String[] args) throws IOException {
		final String REPEATMASKER_FILE = args[0];
		final int REPEATMASKER_FILE_NROWS = Integer.parseInt(args[1]);
		final String SIMPLE_REPEATS_FILE = args[2];
		final int SIMPLE_REPEATS_FILE_NROWS = Integer.parseInt(args[3]);
		final String SEGDUPS_FILE = args[4];
		final int SEGDUPS_FILE_NROWS = Integer.parseInt(args[5]);
        final String OUTPUT_DIR = args[6];
        
        BuildModel.loadRepeatMaskerAnnotations(REPEATMASKER_FILE,REPEATMASKER_FILE_NROWS,BuildModel.MAX_DISTANCE_REPEATMASKER,false);
        BuildModel.loadTrfAnnotations(SIMPLE_REPEATS_FILE,SIMPLE_REPEATS_FILE_NROWS,BuildModel.MIN_OVERLAP_TRF,false);
        BuildModel.buildRepeatIntervals(BuildModel.MIN_OVERLAP_INTERVALS);
        BuildModel.loadSegdups(SEGDUPS_FILE,SEGDUPS_FILE_NROWS,BuildModel.MAX_DISTANCE_SEGDUPS);
        BuildModel.serializeRepeatIntervals(OUTPUT_DIR+"/repeatIntervals.txt");
        BuildModel.serializeSegdups(OUTPUT_DIR+"/segdups.txt");
	}

}