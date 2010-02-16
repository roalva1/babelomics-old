package org.bioinfo.babelomics.tools.functional;

import java.io.IOException;

import org.apache.commons.cli.OptionGroup;
import org.apache.commons.cli.ParseException;
import org.bioinfo.data.list.exception.InvalidIndexException;
import org.bioinfo.tool.OptionFactory;

public class GesbapTool extends FunctionalProfilingTool {

	
	public GesbapTool() {
		
	}
	
	
	/* (non-Javadoc)
	 * @see org.bioinfo.babelomics.tools.functional.FunctionalProfilingTool#initOptions()
	 */
	@Override
	public void initOptions() {
		super.initOptions();
		OptionGroup inputaData = new OptionGroup();
		inputaData.setRequired(true);
		inputaData.addOption(OptionFactory.createOption("snp-file", "Two columns file with SNP id and statistics", false, true));
		inputaData.addOption(OptionFactory.createOption("gene-file", "Two columns file with Gene name and statistics", false, true));
		inputaData.addOption(OptionFactory.createOption("ped-file", "PED file path", false));
		inputaData.addOption(OptionFactory.createOption("map-file", "MAP file path", false));
		inputaData.addOption(OptionFactory.createOption("zip-file", "ZIP file containing PED and MAP files", false));
		options.addOptionGroup(inputaData);
		
		options.addOption(OptionFactory.createOption("method", "Gene set analysis method values: fatiscan, logistic. Default value fatiscan", false, true));
		options.addOption(OptionFactory.createOption("partitions", "Set the number of partitions, by default '30'", false, true));
		options.addOption(OptionFactory.createOption("output-format", "Values: short (just most significant partition) or long (term results for all partitions), by deafult 'short'", false, true));
	}



	@Override
	protected void execute() {
		try {
			// update status
			jobStatus.addStatusMessage("10", "Preparing data");
			
			prepare();
			
			
			
		} catch (IOException e) {
			e.printStackTrace();
		} catch (ParseException e) {
			e.printStackTrace();
		} catch (InvalidIndexException e) {
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
