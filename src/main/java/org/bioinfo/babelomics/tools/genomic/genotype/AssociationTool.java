package org.bioinfo.babelomics.tools.genomic.genotype;

import java.io.IOException;

import org.bioinfo.babelomics.methods.genomic.genotype.GenotypeAnalysis;
import org.bioinfo.tool.OptionFactory;

public class AssociationTool extends GenotypeAnalysisTool {

	public AssociationTool() {

	}

	@Override
	public void initOptions() {
		options.addOption(OptionFactory.createOption("test", "Valid values: assoc, fisher, linear, logistic, tdt", true, true));
		options.addOption(OptionFactory.createOption("maf", "Minor allele frequency, default value: 0.02", false));
		options.addOption(OptionFactory.createOption("log", "Wether to apply Odd ratio logarithm to the file result", false, false));
	}

	@Override
	public void execute() {
		logger.debug("executing association test");
		try {
			// update status
			jobStatus.addStatusMessage("10", "Parsing parameters");
			
			parseGenotypeCommonOptions();
			// specific options
			String test = commandLine.getOptionValue("test", "assoc");
			double maf = Double.parseDouble(commandLine.getOptionValue("maf", "0.02"));

			// prepare the GenotypeAnalysis object for execution
			genotypeAnalysis = new GenotypeAnalysis(pedFilePath, mapFilePath);
			genotypeAnalysis.setPlinkPath(plinkPath);
			genotypeAnalysis.setOutdir(outdir);
			
			logger.debug("executing: "+plinkPath+" --ped "+pedFilePath+" --map "+mapFilePath+" --out "+outdir+"/plink --maf "+maf + "--"+test);
			jobStatus.addStatusMessage("30", "Executing association test");
			
			genotypeAnalysis.association(test, maf);
			
		} catch (IOException e) {
			e.printStackTrace();
			logger.error("Error opening the dataset", e.toString());
		} catch (Exception e) {
			e.printStackTrace();
			logger.error("Error opening the dataset", e.toString());
		} 				
	}

}
