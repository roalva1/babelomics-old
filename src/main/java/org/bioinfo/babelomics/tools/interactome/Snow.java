package org.bioinfo.babelomics.tools.interactome;

import java.io.File;
import java.io.IOException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.ParseException;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.tool.OptionFactory;

public class Snow extends BabelomicsTool {
	public Snow(String[] args) {
		super(args);
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
			CommandLine cmd = parse(args);

			Dataset dataset = new Dataset(new File(cmd.getOptionValue("dataset")));
			String bioEntityName = cmd.getOptionValue("bioentity-name");
			String bioentityInteractionsCheck = cmd.getOptionValue("interactions-check", null);
			String bioentityListGenesProts = cmd.getOptionValue("list-genes-prots");			
			String bioentityInteracionsFile = cmd.getOptionValue("interactions-file", null);
			String bioentityGenesProtsFile = cmd.getOptionValue("genes-prots-file");	
			String bioentityLabelGenProts = cmd.getOptionValue("label-genes-prots");
			String bioentityNatureFilter = cmd.getOptionValue("nature-filter");
			String bioentityInteracionsNumber = cmd.getOptionValue("interactions-number");			
			
			System.out.println(dataset.toString()+"\n");
			
			executeSnow(dataset, bioEntityName, bioentityInteractionsCheck, bioentityListGenesProts,bioentityInteracionsFile, bioentityGenesProtsFile, bioentityLabelGenProts,  bioentityNatureFilter,bioentityInteracionsNumber);
			
		} catch (ParseException e) {
			logger.error("Error parsing command line", e.toString());
			System.out.println("\n");
			printUsage();
		} catch (IOException e) {
			logger.error("Error opening the dataset", e.toString());
		}
		
		
		
	}
	
	private void executeSnow(Dataset dataset, String bioEntityName,	String bioentityInteractionsCheck, String bioentityListGenesProts, String bioentityInteracionsFile, String bioentityGenesProtsFile,	String bioentityLabelGenProts, String bioentityNatureFilter, String bioentityInteracionsNumber) {
		// TODO Auto-generated method stub
		
	}

	
	
	
}
