package org.bioinfo.babelomics.tools.expression;

import java.io.File;
import java.io.IOException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.ParseException;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.tool.OptionFactory;

public class OutlierLrs extends BabelomicsTool {


	public OutlierLrs(String[] args) {
		initOptions();
	}

	@Override
	public void initOptions() {
		options.addOption(OptionFactory.createOption("dataset", "the data"));
		options.addOption(OptionFactory.createOption("low-end", "Low-end of permutations", false));
		options.addOption(OptionFactory.createOption("up-end", "Up-end variable", false));		
	}

	@Override
	public void execute() {
		try {
			Dataset expresionDataset = new Dataset(new File(commandLine.getOptionValue("dataset")));
			
			String permutation = commandLine.getOptionValue("low-end");
			String percentile = commandLine.getOptionValue("up-end");
			executeOutLRS(expresionDataset,permutation,percentile);
		} catch (IOException e) {
			logger.error("Error opening the dataset", e.toString());
		} 
	}


	private void executeOutLRS(Dataset dataset, String permutation,String percentile ) {
		logger.info("executing outLRS, not implemented yet");
	}
	
}
