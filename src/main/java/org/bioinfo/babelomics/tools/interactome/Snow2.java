package org.bioinfo.babelomics.tools.interactome;

import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.tool.OptionFactory;

public class Snow2 extends BabelomicsTool {

	@Override
	public void initOptions() {
		options.addOption(OptionFactory.createOption("sif-file", "A file contoening a SIF interactome", true));
	}

	@Override
	protected void execute() {
		System.out.println("Hello world!!");
	}

}
