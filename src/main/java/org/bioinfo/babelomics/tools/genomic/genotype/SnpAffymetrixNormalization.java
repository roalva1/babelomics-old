package org.bioinfo.babelomics.tools.genomic.genotype;

import java.io.File;

import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.commons.exec.Command;
import org.bioinfo.commons.exec.SingleProcess;
import org.bioinfo.commons.io.utils.FileUtils;
import org.bioinfo.tool.OptionFactory;

public class SnpAffymetrixNormalization extends BabelomicsTool {

	public SnpAffymetrixNormalization() {
	}

	@Override
	public void initOptions() {
		options.addOption(OptionFactory.createOption("chip-type", "Either is a 500K or 6.0", true));
		options.addOption(OptionFactory.createOption("cdf-file", "the data", true));
		options.addOption(OptionFactory.createOption("chrX-snps", "the data", true));
		options.addOption(OptionFactory.createOption("data-dir", "the data", true));
	}

	@Override
	protected void execute() {
		StringBuilder aptCommandLine = new StringBuilder();
		aptCommandLine.append("apt-probeset-genotype ");

		// cdf file
		File cdfFile = new File(commandLine.getOptionValue("cdf-file"));
		if(!cdfFile.exists()) {
			abort("", "", "", "");
		}
		logger.debug("cdf-file valid");
		aptCommandLine.append(" -c " + cdfFile.getAbsolutePath());

		// chrX-snp file
		File chrXSnpsFile = new File(commandLine.getOptionValue("chrX-snps"));
		// to finish
		if(!chrXSnpsFile.exists()) {
			abort("", "", "", "");
		}
		logger.debug("chrX-snps valid");
		aptCommandLine.append(" --chrX-snps " + chrXSnpsFile.getAbsolutePath());

		if(commandLine.getOptionValue("chip-type").equalsIgnoreCase("6.0")) {
			aptCommandLine.append(" -a birdseed ");
		}else {
			if(commandLine.getOptionValue("chip-type").equalsIgnoreCase("500k")) {
				
			}else {
				logger.error("No valid chip type: 500k or 6.0");
			}
		}

		aptCommandLine.append(" --cc-chp-output ");
		aptCommandLine.append(" -o " + new File(outdir).getAbsolutePath());

		File celFiles = new File(commandLine.getOptionValue("data-dir"));
		File[] files = FileUtils.listFiles(celFiles, ".+.cel", true);
		if(files.length < 2) {
			// print error
			logger.error("data-dir does not exist");
		}
		for(int i=0; i<files.length; i++) {
			aptCommandLine.append(" " + files[i].getAbsolutePath());
		}

		Command aptCommand = new Command(aptCommandLine.toString());
		SingleProcess sp = new SingleProcess(aptCommand);
		logger.debug("executing... " + aptCommand.getCommandLine());
		sp.runSync();
		logger.debug("done!");
	}

}
