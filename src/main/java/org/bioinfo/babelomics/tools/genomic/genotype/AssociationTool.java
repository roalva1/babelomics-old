package org.bioinfo.babelomics.tools.genomic.genotype;

import java.io.IOException;

import org.bioinfo.babelomics.methods.genomic.genotype.GenotypeAnalysis;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.tool.OptionFactory;

public class AssociationTool extends GenotypeAnalysisTool {

	public AssociationTool() {

	}

	@Override
	public void initOptions() {
		options.addOption(OptionFactory.createOption("test", "Valid values: assoc, fisher, linear, logistic"));
		options.addOption(OptionFactory.createOption("maf", "Default value: 0.01", false));
		options.addOption(OptionFactory.createOption("log", "Odd ratio logarithm", false, false));

		options.addOption(OptionFactory.createOption("stratification", "Just a flag", false, false));
	}

	@Override
	public void execute() {
		parseGenotypeCommonOptions();
		try {
			String test = commandLine.getOptionValue("test", "assoc");
			double maf = Double.parseDouble(commandLine.getOptionValue("maf", "0.01"));

			// prepare the GenotypeAnalysis object for execution
			genotypeAnalysis = new GenotypeAnalysis(pedFilePath, mapFilePath);
			genotypeAnalysis.setOutdir(outdir);
			if(commandLine.hasOption("plink-path")) {
				genotypeAnalysis.setPlinkPath(plinkPath);			
			}else {
				genotypeAnalysis.setPlinkPath(babelomicsHomePath+"/bin/genotype/plink64");
			}

//			System.out.println(""+Math.log(1)/Math.log(2));
//			System.out.println(""+(0.0)*(-1));
//			System.out.println(""+(0.0==-0.0));

			logger.debug("executing the test");
			genotypeAnalysis.association(test, maf);



		} catch (IOException e) {
			e.printStackTrace();
			logger.error("Error opening the dataset", e.toString());
		} catch (Exception e) {
			e.printStackTrace();
			logger.error("Error opening the dataset", e.toString());
		} 				
	}

	private void executeAssociation(Dataset dataset, String fisher,String duplicates ) {
		logger.info("executing association SNP, not implemented yet");
	}


}
