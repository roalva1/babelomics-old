package org.bioinfo.babelomics.tools.genomic.copynumber;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;

import org.bioinfo.babelomics.exception.InvalidParameterException;
import org.bioinfo.babelomics.methods.genomic.copynumber.CopyNumberAnalysisExecutor;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;

public class CopyNumberAnalysis extends BabelomicsTool {

	@Override
	public void initOptions() {
		// second input: the typical ped and map files
		options.addOption(OptionFactory.createOption("normalized-file", "Normalized file", true, true));

		options.addOption(OptionFactory.createOption("segmentation-method", "Segmentation method: dnacopy or glad", true, true));

		options.addOption(OptionFactory.createOption("cgh-mcr", "Minimal common region", false, false));
		options.addOption(OptionFactory.createOption("gap-allowed", "Gap allowed (for CGH-MCR)", false, true));
		options.addOption(OptionFactory.createOption("altered-low", "Gap altered low (for CGH-MCR)", false, true));
		options.addOption(OptionFactory.createOption("altered-high", "Gap altered high (for CGH-MCR)", false, true));
		options.addOption(OptionFactory.createOption("recurrence", "Recurrence (for CGH-MCR)", false, true));
	}

	@Override
	protected void execute() {

		// reading data
		//
		updateJobStatus("20", "reading input data");

		String normalizedFilename = commandLine.getOptionValue("normalized-file", null);

		String segmentation = commandLine.getOptionValue("segmentation-method", null);
		boolean cghMcr = commandLine.hasOption("cgh-mcr");
		int gapAllowed = Integer.parseInt(commandLine.getOptionValue("gap-allowed", "500"));
		double alteredLow = Double.parseDouble(commandLine.getOptionValue("altered-low", "0.2"));
		double alteredHigh = Double.parseDouble(commandLine.getOptionValue("altered-high", "0.8"));
		int recurrence = Integer.parseInt(commandLine.getOptionValue("recurrence", "50"));

		if ( normalizedFilename == null ) {
			abort("normalizedfilemissing_execute_copynumberanalysis", "input normalized file missing", "input normalized file missing", "input normalized file missing");
		}

		if ( ! new File(normalizedFilename).exists() ) {
			abort("normalizedfilenotexist_execute_copynumberanalysis", "input normalized file does not exist", "input normalized file  does not exist", "input normalized file  does not exist");
		}

		CopyNumberAnalysisExecutor cpExecutor = new CopyNumberAnalysisExecutor(new File(normalizedFilename));

		if ( "dnacopy".equalsIgnoreCase(segmentation) ) {
			cpExecutor.setSegmentationBinPath(babelomicsHomePath + "/bin/copynumber/DNAcopy.r");
		} else if ( "glad".equalsIgnoreCase(segmentation) ) {
			cpExecutor.setSegmentationBinPath(babelomicsHomePath + "/bin/copynumber/GLAD.r");
		} else {
			abort("unknownanalysis_execute_copynumberanalysis", "copy number analysis unknown", "copy number analysis unknown, valid values are dnacopy and glad", "copy number analysis unknown, valid values are dnacopy and glad");
		}

		if ( cghMcr ) {
			cpExecutor.setCghMcrBinPath(babelomicsHomePath + "/bin/copynumber/cghMCR.r");
			cpExecutor.setAlteredLow(alteredLow);
			cpExecutor.setAlteredHigh(alteredHigh);
			cpExecutor.setGapAllowed(gapAllowed);
			cpExecutor.setRecurrence(recurrence);
		} else {
			cpExecutor.setCghMcrBinPath(null);			
		}

		cpExecutor.setOutdir(outdir);

		// executing segmentation
		//
		updateJobStatus("40", "executing segmentation");

		try {
			cpExecutor.run();			
		} catch (InvalidParameterException e) {
			printError("invalidparameteexception_execute_copynumberanalysis","error executing segmentation", e.getMessage(), e);
		}

		// saving results
		//
		updateJobStatus("90", "saving results");

		File outFile = new File(outdir + "/" + cpExecutor.getSegmentatedFilename());
		if ( outFile.exists() ) {
			result.addOutputItem(new Item("segmentedfile", outFile.getName(), "Segmented " + segmentation.toUpperCase() + " file", TYPE.FILE, new ArrayList<String>(2), new HashMap<String, String>(2), segmentation.toUpperCase() + " segmentation"));						
		} else {
			printError("segmentedfileerror_execute_copynumberanalyis", "error creating segmented file", "error creating segmented file");
		}
		if ( cghMcr ) {
			outFile = new File(outdir + "/" + cpExecutor.getCghFilename());
			if ( outFile.exists() ) {
				result.addOutputItem(new Item("cghfile", outFile.getName(), "CGH file for " + segmentation.toUpperCase() + " segmentation", TYPE.FILE, new ArrayList<String>(2), new HashMap<String, String>(2), segmentation.toUpperCase() + " CGH"));						
			} else {
				printError("cghfileerror_execute_copynumberanalyis", "error creating cgh file", "error creating cgh file");			
			}
		}


	}

}
