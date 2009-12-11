package org.bioinfo.babelomics.tools.genomic.genotype;

import java.io.File;
import java.io.IOException;
import java.security.InvalidParameterException;

import org.bioinfo.babelomics.methods.genomic.genotype.GenotypeAnalysis;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.tool.OptionFactory;

public class GenotypeAnalysisTool extends BabelomicsTool {

	protected GenotypeAnalysis genotypeAnalysis;
	
	protected String pedFilePath;
	protected String mapFilePath;
	protected String plinkPath;
	
	public GenotypeAnalysisTool() {
		genotypeAnalysis = new GenotypeAnalysis();
	}

	@Override
	public void initOptions() {
		// the typical ped and map files
		options.addOption(OptionFactory.createOption("ped-file", "Just a flag", false));
		options.addOption(OptionFactory.createOption("map-file", "Just a flag", false));
		options.addOption(OptionFactory.createOption("plink", "Just a flag", false));
		
		options.addOption(OptionFactory.createOption("tdt", "Just a flag", false, false));
	}

	protected void parseGenotypeCommonOptions() {
		pedFilePath = commandLine.getOptionValue("ped-file");
		mapFilePath = commandLine.getOptionValue("map-file");
		plinkPath = commandLine.getOptionValue("plink");
	}
	
	@Override
	protected void execute() {
		
		logger.debug("executing the test");
		File pedFile = null;
		File mapFile = null;
		GenotypeAnalysis genotypeAnalysis;

		if(commandLine.hasOption("pedigree-file") && commandLine.hasOption("chp-dir") && commandLine.hasOption("ped-file") && commandLine.hasOption("map-file")) {
			abort("", "", "", "");
		}

		if(commandLine.hasOption("create-ped-map")) {
			if(commandLine.hasOption("pedigree-file") && commandLine.hasOption("chp-dir")) {
				try {
//					AffyGenotypeUtils.affyCallResultsToPedAndMap(commandLine.getOptionValue("chp-dir"), commandLine.getOptionValue("pedigree-file"), outdir);
					pedFile = new File(outdir+"/plink.ped");
					mapFile = new File(outdir+"/plink.map");
				}catch (InvalidParameterException e) {
					printError("", "", "");
				}
//				catch (IOException e) {
//					printError("", "", "");
//				}
			}else {
				abort("", "", "", "");
			}
		}else {
			if(commandLine.hasOption("ped-file") && commandLine.hasOption("map-file")) {
				pedFile = new File(commandLine.getOptionValue("ped-file"));
				mapFile = new File(commandLine.getOptionValue("map-file"));
			}else {
				abort("", "", "", "");
			}

		}

		if(pedFile == null || !pedFile.exists() || mapFile == null || !mapFile.exists()) {
			abort("", "", "", "");
		}

		genotypeAnalysis = new GenotypeAnalysis(pedFile, mapFile);
		genotypeAnalysis.setPlinkPath(babelomicsHomePath+"/bin/genotype/plink64");
		genotypeAnalysis.setOutdir(outdir);

		if(commandLine.hasOption("stratification")) {
			try {
				genotypeAnalysis.stratification();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

		if(commandLine.hasOption("assoc")) {
			try {
				genotypeAnalysis.association(commandLine.getOptionValue("assoc"));
			} catch (IOException e) {
				e.printStackTrace();
			}
		}


	}


}
