package org.bioinfo.babelomics.tools;

import org.bioinfo.tool.OptionFactory;

public class Histogram extends BabelomicsTool {

	@Override
	public void initOptions() {
		options.addOption(OptionFactory.createOption("dataset", "the data"));
		options.addOption(OptionFactory.createOption("logarithm-base", "the logarithm base to apply transformation, possible values: e, 2 10 for log base 2, log base 2 and log base 10", false));
		options.addOption(OptionFactory.createOption("merge-replicates", "method to merge replicates, valid values are: mean or median", false));
		options.addOption(OptionFactory.createOption("filter-missing", "minimum percentage of existing values, from 0 to 100", false));
		options.addOption(OptionFactory.createOption("impute-missing", "method to impute missing values, valid values are: zero, mean, median, knn", false));
		options.addOption(OptionFactory.createOption("kvalue", "kvalue for knn impute method, default 15", false));
		options.addOption(OptionFactory.createOption("gene-file-filter", "This option will remove all the patterns of the genes that are not present in this gene file", false));
		//		options.addOption(OptionFactory.createOption("sample-filter", "class variable", false));
		//		options.addOption(OptionFactory.createOption("feature-filter", "class variable", false));
	}

	@Override
	protected void execute() {
		// TODO Auto-generated method stub

	}

}
