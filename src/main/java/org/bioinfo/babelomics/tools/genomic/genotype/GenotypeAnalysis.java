package org.bioinfo.babelomics.tools.genomic.genotype;

import java.io.File;

import org.bioinfo.babelomics.methods.genomic.genotype.GenotypeUtils;
import org.bioinfo.babelomics.methods.genomic.genotype.PlinkTestsExecutor;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.tool.OptionFactory;

public class GenotypeAnalysis extends BabelomicsTool {

	public GenotypeAnalysis() {
		
	}
	
	@Override
	public void initOptions() {
		// first input: result coming from apt-probeset-genotype and a pedigree like format, this option creates .ped and .map files 
		options.addOption(OptionFactory.createOption("pedigree-file", "Just a flag", false));
		options.addOption(OptionFactory.createOption("chp-dir", "Just a flag", false));
		
		// second input: the typical ped and map files
		options.addOption(OptionFactory.createOption("ped-file", "Just a flag", false));
		options.addOption(OptionFactory.createOption("map-file", "Just a flag", false));

		
		options.addOption(OptionFactory.createOption("create-ped-map", "Just a flag", false, false));
		options.addOption(OptionFactory.createOption("stratification", "Just a flag", false, false));
		options.addOption(OptionFactory.createOption("tdt", "Just a flag", false, false));
		options.addOption(OptionFactory.createOption("assoc", "Either is a assoc or fisher", false));
		options.addOption(OptionFactory.createOption("assoc-odd-ratio-log", "Just a flag", false, false));
	}

	@Override
	protected void execute() {
		File pedFile = null;
		File mapFile = null;
		PlinkTestsExecutor plinkTestsExecutor;
		
		if(commandLine.hasOption("pedigree-file") && commandLine.hasOption("chp-dir") && commandLine.hasOption("ped-file") && commandLine.hasOption("map-file")) {
			abort("", "", "", "");
		}
		
		if(commandLine.hasOption("create-ped-map")) {
			
			pedFile = new File(outdir+"/plink.ped");
			mapFile = new File(outdir+"/plink.map");
			GenotypeUtils.affyToPedAndMap("");
			
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
		
		plinkTestsExecutor = new PlinkTestsExecutor(pedFile, mapFile);
		plinkTestsExecutor.setPlinkPath(babelomicsHomePath+"/bin/genotype/plink64");
		plinkTestsExecutor.setOutdir(outdir);
		
		if(commandLine.hasOption("stratification")) {
			plinkTestsExecutor.stratification();
		}
		
		if(commandLine.hasOption("assoc")) {
			plinkTestsExecutor.association(commandLine.getOptionValue("assoc"));
		}
		
		
	}

	
}
