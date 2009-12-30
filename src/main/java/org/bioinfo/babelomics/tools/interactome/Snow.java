package org.bioinfo.babelomics.tools.interactome;

import java.io.File;
import java.io.IOException;

import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.tool.OptionFactory;

public class Snow extends BabelomicsTool {
	
	public Snow() {
		super();
	}

	@Override
	public void initOptions() {
		options.addOption(OptionFactory.createOption("list1", "the list1"));		
		options.addOption(OptionFactory.createOption("list2", "the list2", false));		
		options.addOption(OptionFactory.createOption("interactome", "Select interactome: join (for human interactome with all ppis), intersect (for human interacionts with only ppis detected by two methods) and own (for your own interactions)", false));
		options.addOption(OptionFactory.createOption("own-interactions", "submit a file with interactome ",false));
		options.addOption(OptionFactory.createOption("check-interactions", "Set this option if proteins in interactions and list are in same id",false,false));
		options.addOption(OptionFactory.createOption("id-nature", "Nature of your lists: genes or proteins", false));
		options.addOption(OptionFactory.createOption("interactions-number", "Maximum number of external proteins introduced: 0, 1, 2 or 3", false));
	}


	@Override
	public void execute() {
		try {
			File f1 = new File(commandLine.getOptionValue("list1"));
			File f2 = commandLine.hasOption("list2") ? new File(commandLine.getOptionValue("list2")) :  null;
			String interactome = commandLine.getOptionValue("interactome", "join");
			File interactionsFile = commandLine.hasOption("own-interactions") ? new File(commandLine.getOptionValue("own-interactions")) :  null;
			boolean checkInteractions = commandLine.hasOption("check-interactions");
			String idNature = commandLine.getOptionValue("id-nature", "proteins");
			int interactionsNumber = Integer.parseInt(commandLine.getOptionValue("interactions-number", "1"));
						
			executeSnow(dataset, bioEntityName, bioentityInteractionsCheck, bioentityListGenesProts,bioentityInteracionsFile, bioentityGenesProtsFile, bioentityLabelGenProts,  bioentityNatureFilter,bioentityInteracionsNumber);
			
		} catch (IOException e) {
			logger.error("Error opening the dataset", e.toString());
		}
	}
	
	private void executeSnow(Dataset dataset, String bioEntityName,	String bioentityInteractionsCheck, String bioentityListGenesProts, String bioentityInteracionsFile, String bioentityGenesProtsFile,	String bioentityLabelGenProts, String bioentityNatureFilter, String bioentityInteracionsNumber) {
		// TODO Auto-generated method stub
c		
	}

	
	
	
}
