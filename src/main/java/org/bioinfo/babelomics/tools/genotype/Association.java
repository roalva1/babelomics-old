package org.bioinfo.babelomics.tools.genotype;

import java.io.File;

import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.tool.OptionFactory;

public class Association extends BabelomicsTool {
	public Association(String[] args) {
		initOptions();
	}

	@Override
	public void initOptions() {
		getOptions().addOption(OptionFactory.createOption("dataset", "the data"));
		getOptions().addOption(OptionFactory.createOption("fisher", "Fisher exact test"));
		getOptions().addOption(OptionFactory.createOption("duplicates", "Remove duplicates"));
	}

	@Override
	public void execute() {
		try {
			Dataset dataset = new Dataset(new File(commandLine.getOptionValue("dataset")));
			String fisher = commandLine.getOptionValue("fisher");
			String duplicates = commandLine.getOptionValue("duplicates");			
			
			System.out.println(dataset.toString()+"\n");
		
			logger.info("Agilent microarray G2518A converter, not yet implemented");
		
			executeAssociation(dataset,fisher,duplicates);
		} catch (Exception e) {
			logger.error("Error opening the dataset", e.toString());
		} 				
	}
	
	private void executeAssociation(Dataset dataset, String fisher,String duplicates ) {
		logger.info("executing association SNP, not implemented yet");
	}
	
	
}
