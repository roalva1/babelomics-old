package org.bioinfo.babelomics.tools.functional.textmining;

import java.io.File;
import java.io.IOException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.ParseException;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.data.dataset.FeatureData;
import org.bioinfo.tool.OptionFactory;

public class MarmiteScan extends BabelomicsTool {


	public MarmiteScan(String[] args) {
		super(args);
		initOptions();
	}

	@Override
	public void initOptions() {
		options.addOption(OptionFactory.createOption("dataset", "the data"));
		options.addOption(OptionFactory.createOption("bioentity-name", "possibilities: Disease associated words,Chemical products,Word roots"));
		options.addOption(OptionFactory.createOption("bioentity-score-filter", "Minimum number of genes with a score (0-10000)"));
		options.addOption(OptionFactory.createOption("bioentity-number-filter", "Number of bio-entities in results (0-10000)"));
		options.addOption(OptionFactory.createOption("gene-name-list", "gene name in list", false));
		options.addOption(OptionFactory.createOption("partition-number", "Number of partitions "));
		options.addOption(OptionFactory.createOption("significance", "p-value for statistical significance"));
		options.addOption(OptionFactory.createOption("sort", "Sort list"));

	}

	@Override
	public void execute() {
		try {
			FeatureData fd1 = new FeatureData(new File(commandLine.getOptionValue("dat")));
			String className = commandLine.getOptionValue("class");
			String test = commandLine.getOptionValue("test");

			String timeClass = commandLine.getOptionValue("time-class", null);
			String censoredClass = commandLine.getOptionValue("censor-class", null);

			
			Dataset dataset = new Dataset(new File(commandLine.getOptionValue("dataset")));
			String bioEntityName = commandLine.getOptionValue("bioentity-name");
			String bioentityScoreFilter = commandLine.getOptionValue("bioentity-score-filter");
			String bioentityNumberFilter = commandLine.getOptionValue("bioentity-mumber-filter");			
			String geneNameList = commandLine.getOptionValue("gene-name-list", null);
			String partitionNumber = commandLine.getOptionValue("partition-number");	
			String significance = commandLine.getOptionValue("significance");
			String sort = commandLine.getOptionValue("sort");
			System.out.println(dataset.toString()+"\n");
			
			executeMarmiteScan(dataset, bioEntityName, bioentityScoreFilter, bioentityNumberFilter,geneNameList, partitionNumber, significance,  sort);
			
		} catch (Exception e) {
			logger.error("Error opening the dataset", e.toString());
		}
		
		
		
	}
	
	private void executeMarmiteScan(Dataset dataset, String bioEntityName, String bioentityScoreFilter, String bioentityNumberFilter,String geneNameList, String partitionNumber, String significance, String sort) {
		logger.info("executing svm, not implemented yet");
	}
	
}
