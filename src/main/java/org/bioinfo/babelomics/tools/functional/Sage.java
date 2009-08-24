package org.bioinfo.babelomics.tools.functional;

import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.tool.OptionFactory;

public class Sage extends BabelomicsTool {


	public Sage(String[] args) {
		initOptions();
	}

	@Override
	public void initOptions() {
		options.addOption(OptionFactory.createOption("dataset", "the data"));
		options.addOption(OptionFactory.createOption("org-filter", "Organism"));
		options.addOption(OptionFactory.createOption("tissue-name", "Tissue selection name"));
		options.addOption(OptionFactory.createOption("histology-name", "Histologies selection name"));
		options.addOption(OptionFactory.createOption("list-filter", "Number of list of genes taken account"));
		options.addOption(OptionFactory.createOption("list-tag-filter", "Type of Tag"));
		options.addOption(OptionFactory.createOption("null-values-genes", "Percentage for null values accepted in genes"));
		options.addOption(OptionFactory.createOption("null-values-libraries", "Percentage for null values accepted in libraries"));
		options.addOption(OptionFactory.createOption("cell-lines", "Percentage for null values accepted in libraries",false,false));
	}

	@Override
	public void execute() {
//		try {
//			CommandLine cmd = parse(args);
//
//			Dataset dataset = new Dataset(new File(cmd.getOptionValue("dataset")));
//			String bioEntityName = cmd.getOptionValue("org-filter");
//			String bientityTissueName = cmd.getOptionValue("tissue-name");
//			String bioentityListFilter = cmd.getOptionValue("list-filter");
//			String bioentityHistologyName = cmd.getOptionValue("histology-name");			
//			String bioentityNullValuesGenes = cmd.getOptionValue("null-values-genes");
//			String bioentityNullValuesLibraries = cmd.getOptionValue("null-values-libraries");
//			String bioentityCellLines = cmd.getOptionValue("cell-lines",null);
//			System.out.println(dataset.toString()+"\n");
//			
//			executeSageAnalisys(dataset, bioEntityName, bientityTissueName, bioentityListFilter,bioentityHistologyName,bioentityNullValuesGenes,bioentityNullValuesLibraries,bioentityCellLines);
//			
//		} catch (ParseException e) {
//			logger.error("Error parsing command line", e.toString());
//			System.out.println("\n");
//			printUsage();
//		} catch (IOException e) {
//			logger.error("Error opening the dataset", e.toString());
//		}		
	}
	
	private void executeSageAnalisys(Dataset dataset, String bioEntityName,	String bientityTissueName, String bioentityListFilter,	String bioentityHistologyName, String bioentityNullValuesGenes,	String bioentityNullValuesLibraries, String bioentityCellLines) {
		// TODO Auto-generated method stub
		
	}

	
}
