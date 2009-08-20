package org.bioinfo.babelomics.tools.genotype;

import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.tool.OptionFactory;

public class Linkage extends BabelomicsTool {
	public Linkage(String[] args) {
		initOptions();
	}

	@Override
	public void initOptions() {
		options.addOption(OptionFactory.createOption("dataset", "the data"));
	}

	@Override
	public void execute() {
		logger.info("executing transmid¡ssion SNP, not implemented yet");		
	}
}
