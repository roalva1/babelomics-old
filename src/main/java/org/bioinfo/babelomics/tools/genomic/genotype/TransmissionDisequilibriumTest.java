package org.bioinfo.babelomics.tools.genomic.genotype;

import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.tool.OptionFactory;

public class TransmissionDisequilibriumTest extends BabelomicsTool {
	
	public TransmissionDisequilibriumTest() {

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
