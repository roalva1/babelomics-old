package org.bioinfo.babelomics.tools.genotype;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.ParseException;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.tool.OptionFactory;

public class Linkage extends BabelomicsTool {
	public Linkage(String[] args) {
		super(args);
		initOptions();
		
		
		
	}

	@Override
	public void initOptions() {
		options.addOption(OptionFactory.createOption("dataset", "the data"));
	}

	@Override
	public void execute() {
		// TODO Auto-generated method stub
		
		try {
			CommandLine cmd = parse(args, true);
		} catch (ParseException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		// TODO Auto-generated method stub
		logger.info("executing transmid¡ssion SNP, not implemented yet");
		
	}
}
