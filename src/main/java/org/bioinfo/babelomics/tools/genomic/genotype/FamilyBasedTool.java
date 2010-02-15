package org.bioinfo.babelomics.tools.genomic.genotype;

import java.io.IOException;

import org.bioinfo.babelomics.methods.genomic.genotype.GenotypeAnalysis;
import org.bioinfo.tool.OptionFactory;

public class FamilyBasedTool extends GenotypeAnalysisTool {
	
	public FamilyBasedTool() {

	}

	@Override
	public void initOptions() {
		options.addOption(OptionFactory.createOption("test", "Valid values: tdt"));
	}


	@Override
	public void execute() {		
		logger.debug("executing transmission disequilibrium test (TDT)");
		try {
			parseGenotypeCommonOptions();

			// prepare the GenotypeAnalysis object for execution
			genotypeAnalysis = new GenotypeAnalysis(pedFilePath, mapFilePath);
			genotypeAnalysis.setPlinkPath(plinkPath);
			genotypeAnalysis.setOutdir(outdir);
			logger.debug("executing: "+plinkPath+" --ped "+pedFilePath+" --map "+mapFilePath+" --out "+outdir+"/plink --tdt");

		} catch (IOException e) {
			e.printStackTrace();
			logger.error("Error opening the dataset", e.toString());
		} catch (Exception e) {
			e.printStackTrace();
			logger.error("Error opening the dataset", e.toString());
		} 				
	}
	
	
}
