package org.bioinfo.babelomics.tools.functional;

import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.tool.OptionFactory;

public class Affymetrix extends BabelomicsTool {


	public Affymetrix(String[] args) {
		initOptions();
	}

	@Override
	public void initOptions() {
		options.addOption(OptionFactory.createOption("dataset", "the data"));
		options.addOption(OptionFactory.createOption("org-filter", "Organism"));
		options.addOption(OptionFactory.createOption("tissue-name", "Tissue selection name"));
		options.addOption(OptionFactory.createOption("list-filter", "Number of list of genes taken account"));
		options.addOption(OptionFactory.createOption("tissue-normalization-filter", "Normalization method"));
		options.addOption(OptionFactory.createOption("tissue-expression-filter", "Multiple probes expression value "));
	}

	@Override
	public void execute() {
//		try {
//
//			Dataset dataset = new Dataset(new File(cmd.getOptionValue("dataset")));
//			String bioEntityName = cmd.getOptionValue("org-filter");
//			String bientityTissueName = cmd.getOptionValue("tissue-name");
//			String bioentityListFilter = cmd.getOptionValue("list-filter");			
//			String bioentityNormalizationFilter = cmd.getOptionValue("tissue-normalization-filter");
//			String bientityExpressionFilter = cmd.getOptionValue("tissue-expression-filter");
//			System.out.println(dataset.toString()+"\n");
//			
//			executeAffyAnalisys(dataset, bioEntityName, bientityTissueName, bioentityListFilter,bioentityNormalizationFilter, bientityExpressionFilter);
//			
//		} catch (ParseException e) {
//			logger.error("Error parsing command line", e.toString());
//			System.out.println("\n");
//			printUsage();
//		} catch (IOException e) {
//			logger.error("Error opening the dataset", e.toString());
//		}		
	}
	
	private void executeAffyAnalisys(Dataset dataset, String bioEntityName,	String bientityTissueName, String bioentityListFilter, String bioentityNormalizationFilter, String bientityExpressionFilter) {
		// TODO Auto-generated method stub
		
	}

	
	
}
