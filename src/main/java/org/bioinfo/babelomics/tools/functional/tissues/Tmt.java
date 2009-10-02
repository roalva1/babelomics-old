package org.bioinfo.babelomics.tools.functional.tissues;

import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.tool.OptionFactory;

public class Tmt extends BabelomicsTool {

	@Override
	public void initOptions() {
		options.addOption(OptionFactory.createOption("organism", "Organism, valid values are 'human' and 'mouse'"));
		options.addOption(OptionFactory.createOption("list1", "the feature data containig the list #1 of genes, or the feature data file"));
		options.addOption(OptionFactory.createOption("list2", "the feature data containig the list #2 of genes, or the feature data file", false));
		options.addOption(OptionFactory.createOption("tissues", "the list of tissues separated by commas. Enter 'all tissues' to take into account all available tissues"));
	}

	@Override
	protected void execute() {
	}

}
