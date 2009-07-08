package org.bioinfo.babelomics.tools.genotype;

import java.io.File;
import java.io.IOException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.ParseException;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.tool.OptionFactory;

public class Association extends BabelomicsTool {
	public Association(String[] args) {
		super(args);
		initOptions();
	}

	@Override
	public void initOptions() {
		options.addOption(OptionFactory.createOption("dataset", "the data"));
		options.addOption(OptionFactory.createOption("fisher", "Fisher exact test"));
		options.addOption(OptionFactory.createOption("duplicates", "Remove duplicates"));
	}

	@Override
	public void execute() {
		try {
			
			
			CommandLine cmd = parse(args, true);
			
			Dataset dataset = new Dataset(new File(cmd.getOptionValue("dataset")));
			String fisher = cmd.getOptionValue("fisher");
			String duplicates = cmd.getOptionValue("duplicates");			
			
			System.out.println(dataset.toString()+"\n");
		
			logger.info("Agilent microarray G2518A converter, not yet implemented");
		
			executeAssociation(dataset,fisher,duplicates);
					

		} catch (ParseException e) {
			logger.error("Error parsing command line", e.toString());
			System.out.println("\n");
			printUsage();
		} catch (IOException e) {
			logger.error("Error opening the dataset", e.toString());
		} 
				
	}
	
	private void executeAssociation(Dataset dataset, String fisher,String duplicates ) {
		logger.info("executing association SNP, not implemented yet");
	}
	
	
}
