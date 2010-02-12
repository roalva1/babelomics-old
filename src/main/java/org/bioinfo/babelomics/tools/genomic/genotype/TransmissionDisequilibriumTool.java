package org.bioinfo.babelomics.tools.genomic.genotype;

import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.tool.OptionFactory;

public class TransmissionDisequilibriumTool extends BabelomicsTool {
	
	public TransmissionDisequilibriumTool() {

	}

	@Override
	public void initOptions() {
		options.addOption(OptionFactory.createOption("tdt", "Just a flag", false, false));
		
		getOptions().addOption(OptionFactory.createOption("dataset", "the data"));
		
	}


	@Override
	public void execute() {		
		logger.info("executing transmid¡ssion SNP, not implemented yet");
	}
	
	
}
