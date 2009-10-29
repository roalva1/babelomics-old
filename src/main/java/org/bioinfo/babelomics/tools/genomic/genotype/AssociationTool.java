package org.bioinfo.babelomics.tools.genomic.genotype;

import java.io.File;
import java.io.IOException;

import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.tool.OptionFactory;

public class AssociationTool extends GenotypeAnalysisTool {
	
	public AssociationTool() {

	}

	@Override
	public void initOptions() {
		super.initOptions();
		
		options.addOption(OptionFactory.createOption("test", "Valid values: assoc, fisher, linear, logistic"));
		options.addOption(OptionFactory.createOption("maf", "Default value: 0.01", false));
		options.addOption(OptionFactory.createOption("log", "Odd ratio logarithm", false, false));
		
		options.addOption(OptionFactory.createOption("stratification", "Just a flag", false, false));
	}

	@Override
	public void execute() {
		try {
			
			String test = commandLine.getOptionValue("test", "assoc");
			double maf = Double.parseDouble(commandLine.getOptionValue("maf", "0.01"));
			
			
			genotypeAnalysis.association(test, maf);
			
			Dataset dataset = new Dataset(new File(commandLine.getOptionValue("dataset")));
			String fisher = commandLine.getOptionValue("fisher");
			String duplicates = commandLine.getOptionValue("duplicates");			
			
			System.out.println(dataset.toString()+"\n");
		
			logger.info("Agilent microarray G2518A converter, not yet implemented");
		
			executeAssociation(dataset,fisher,duplicates);
		} catch (IOException e) {
			logger.error("Error opening the dataset", e.toString());
		} catch (Exception e) {
			logger.error("Error opening the dataset", e.toString());
		} 				
	}
	
	private void executeAssociation(Dataset dataset, String fisher,String duplicates ) {
		logger.info("executing association SNP, not implemented yet");
	}
	
	
}
