package org.bioinfo.babelomics.tools.genomic.genotype;

import java.io.File;
import java.io.IOException;
import java.security.InvalidParameterException;

import org.bioinfo.babelomics.methods.genomic.genotype.AffyGenotypeUtils;
import org.bioinfo.babelomics.methods.genomic.genotype.GenotypeAnalysis;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.tool.OptionFactory;

@Deprecated
public class GenotypeAnalysisTool extends BabelomicsTool {

	public GenotypeAnalysisTool() {

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
		GenotypeAnalysis genotypeAnalysis;

		if(commandLine.hasOption("pedigree-file") && commandLine.hasOption("chp-dir") && commandLine.hasOption("ped-file") && commandLine.hasOption("map-file")) {
			abort("", "", "", "");
		}

		if(commandLine.hasOption("create-ped-map")) {
			if(commandLine.hasOption("pedigree-file") && commandLine.hasOption("chp-dir")) {
				try {
					AffyGenotypeUtils.affyToPedAndMap(commandLine.getOptionValue("chp-dir"), commandLine.getOptionValue("pedigree-file"), outdir);
					pedFile = new File(outdir+"/plink.ped");
					mapFile = new File(outdir+"/plink.map");
				}catch (InvalidParameterException e) {
					printError("", "", "");
				} catch (IOException e) {
					printError("", "", "");
				}
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
