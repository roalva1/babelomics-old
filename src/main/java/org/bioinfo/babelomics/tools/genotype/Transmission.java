package org.bioinfo.babelomics.tools.genotype;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.ParseException;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.tool.OptionFactory;

public class Transmission extends BabelomicsTool {
	public Transmission(String[] args) {
		initOptions();
	}

	@Override
	public void initOptions() {
		
		getOptions().addOption(OptionFactory.createOption("dataset", "the data"));
		
	}


	@Override
	public void execute() {		
		logger.info("executing transmid¡ssion SNP, not implemented yet");
	}
	
	
}
