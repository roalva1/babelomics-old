package org.bioinfo.babelomics.tools.genotype;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.ParseException;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.tool.OptionFactory;

public class Transmission extends BabelomicsTool {
	public Transmission(String[] args) {
		super(args);
		initOptions();
	}

	@Override
	public void initOptions() {
		
		options.addOption(OptionFactory.createOption("dataset", "the data"));
		
	}


	@Override
	public void execute() {
		try {
			CommandLine cmd = parse(args, true);
		} catch (ParseException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		// TODO Auto-generated method stub
		
		// TODO Auto-generated method stub
		
		logger.info("executing transmid¡ssion SNP, not implemented yet");
		
	}
	
	
}
