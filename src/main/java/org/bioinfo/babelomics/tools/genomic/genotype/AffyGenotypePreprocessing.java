package org.bioinfo.babelomics.tools.genomic.genotype;

import java.io.File;
import java.io.IOException;

import org.bioinfo.babelomics.methods.genomic.genotype.AffyGenotypeUtils;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.commons.exec.Command;
import org.bioinfo.commons.exec.SingleProcess;
import org.bioinfo.commons.io.utils.FileUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.tool.OptionFactory;

public class AffyGenotypePreprocessing extends BabelomicsTool {

	public AffyGenotypePreprocessing() {

	}

	@Override
	public void initOptions() {
		options.addOption(OptionFactory.createOption("data-dir", "the data", true));

		options.addOption(OptionFactory.createOption("calls", "the data", false));
		options.addOption(OptionFactory.createOption("chip-type", "Either is a 500K or 6.0", true));
		options.addOption(OptionFactory.createOption("cdf-file", "the data", false));
		options.addOption(OptionFactory.createOption("chrX-snps-file", "the data", false));

		options.addOption(OptionFactory.createOption("chp-to-ped-map", "the data", false));
		options.addOption(OptionFactory.createOption("pedigree-file", "Just a flag", false));
	}

	@Override
	protected void execute() {
//		File dataDirFile = new File(commandLine.getOptionValue("data-dir"));
//		if(dataDirFile == null || !dataDirFile.exists()) {
//			abort("", "", "", "");
//		}

		if(!commandLine.hasOption("calls") || !commandLine.hasOption("chp-to-ped-map")) {
			printUsage("./babelomics.sh");
			abort("", "", "", "");
		}


		if(commandLine.hasOption("calls")) {
			// check that we have the parameters we need
			if(!commandLine.hasOption("cdf-file") || !commandLine.hasOption("chrX-snps")) {
				printUsage("./babelomics.sh");
				abort("", "", "", "");
			}
			File cdfFile = new File(commandLine.getOptionValue("cdf-file"));
			File chrXSnpsFile = new File(commandLine.getOptionValue("chrX-snps"));
			File dataDirFile = new File(commandLine.getOptionValue("data-dir"));
			try {
				if(!commandLine.getOptionValue("chip-type").equalsIgnoreCase("500k") || !commandLine.getOptionValue("chip-type").equalsIgnoreCase("6.0")) {
					abort("", "", "", "");
					logger.error("No valid chip type: 500k or 6.0");
				}
				
				AffyGenotypeUtils.affyGenotypeNormalization(babelomicsHomePath+"/bin/apt/apt-probeset-genotype", commandLine.getOptionValue("chip-type"), cdfFile.getAbsolutePath(), chrXSnpsFile.getAbsolutePath(), outdir, dataDirFile.getAbsolutePath());
				
			}catch (IOException e) {
				abort("", "", "", "");
				logger.error("Error in calls", StringUtils.getStackTrace(e));
			}
			// update dataDirFile with chp dir
			dataDirFile = new File(outdir+"/chp-files");
		}


		if(commandLine.hasOption("chp-to-ped-map")) {
			//			AffyGenotypeUtils.affyToPedAndMap(chpDirPath, pedigreFilePath, outdir)
		}


	}

}
