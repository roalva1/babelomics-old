package org.bioinfo.babelomics.tools.interactome;

import java.io.File;
import java.io.IOException;

import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.collections.exceptions.InvalidColumnIndexException;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.tool.OptionFactory;

public class Snow extends BabelomicsTool {
	public Snow(String[] args) {
		initOptions();
	}

	@Override
	public void initOptions() {
		
		options.addOption(OptionFactory.createOption("dataset", "the data"));		
		options.addOption(OptionFactory.createOption("interactions-check", " if proteins in interactions and list are in same id",false,false));
		options.addOption(OptionFactory.createOption("bioentity-name", "Select interactome"));
		options.addOption(OptionFactory.createOption("list-genes-prots", "genes/proteins list ", false));
		options.addOption(OptionFactory.createOption("interactions-file", "submit a file with interactome ",false));
		options.addOption(OptionFactory.createOption("genes-prots-file", "ubmit a file with genes/proteins list",false));
		options.addOption(OptionFactory.createOption("label-genes-prots", "Label for list of genes/proteins"));
		options.addOption(OptionFactory.createOption("nature-filter", "Nature of your lists"));
		options.addOption(OptionFactory.createOption("interactions-number", "Nature of your lists"));
		
	}


	@Override
	public void execute() {
		try {

			Dataset dataset = new Dataset(new File(commandLine.getOptionValue("dataset")));
			String bioEntityName = commandLine.getOptionValue("bioentity-name");
			String bioentityInteractionsCheck = commandLine.getOptionValue("interactions-check", null);
			String bioentityListGenesProts = commandLine.getOptionValue("list-genes-prots");			
			String bioentityInteracionsFile = commandLine.getOptionValue("interactions-file", null);
			String bioentityGenesProtsFile = commandLine.getOptionValue("genes-prots-file");	
			String bioentityLabelGenProts = commandLine.getOptionValue("label-genes-prots");
			String bioentityNatureFilter = commandLine.getOptionValue("nature-filter");
			String bioentityInteracionsNumber = commandLine.getOptionValue("interactions-number");			
			
			System.out.println(dataset.toString()+"\n");
			
			executeSnow(dataset, bioEntityName, bioentityInteractionsCheck, bioentityListGenesProts,bioentityInteracionsFile, bioentityGenesProtsFile, bioentityLabelGenProts,  bioentityNatureFilter,bioentityInteracionsNumber);
			
		} catch (IOException e) {
			logger.error("Error opening the dataset", e.toString());
		} catch (InvalidColumnIndexException e) {
			e.printStackTrace();
		}
		
		
		
	}
	
	private void executeSnow(Dataset dataset, String bioEntityName,	String bioentityInteractionsCheck, String bioentityListGenesProts, String bioentityInteracionsFile, String bioentityGenesProtsFile,	String bioentityLabelGenProts, String bioentityNatureFilter, String bioentityInteracionsNumber) {
		// TODO Auto-generated method stub
		
	}

	
	
	
}
