package org.bioinfo.babelomics.tools;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;

import org.apache.commons.cli.ParseException;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.tool.GenericBioTool;
import org.bioinfo.tool.OptionFactory;

public abstract class BabelomicsTool extends GenericBioTool {
	
	protected String species;
	protected String babelomicsHomePath;
	
	
	public BabelomicsTool() {
		babelomicsHomePath = System.getenv("BABELOMICS_HOME");
		initCommonsOptions();
		initOptions();
	}
	
	public abstract void initOptions();
	
	private void initCommonsOptions() {
		getOptions().addOption(OptionFactory.createOption("tool", "to", "tool name", true));
		getOptions().addOption(OptionFactory.createOption("species", "The specie of the ids", false));	

	}

	
	/* (non-Javadoc)
	 * @see org.bioinfo.tool.GenericBioTool#parse(java.lang.String[])
	 */
	@Override
	public void parse(String[] args) throws ParseException, IOException {
		parse(args, false);
	}
	
	
	/* (non-Javadoc)
	 * @see org.bioinfo.tool.GenericBioTool#parse(java.lang.String[], boolean)
	 */
	@Override
	public void parse(String[] args, boolean stopAtNoOption) throws ParseException, IOException {
		super.parse(args, stopAtNoOption);
		
		// must be in commandLine, just in case we initialize...
		this.toolName = commandLine.getOptionValue("tool", "");
		this.species = commandLine.getOptionValue("species", "unknown");
	}
	
	/* (non-Javadoc)
	 * @see org.bioinfo.tool.GenericBioTool#abort(java.lang.String, java.lang.String, java.lang.String, java.lang.String)
	 */
	@Override
	public void abort(String name, String title, String msg, String text) {
		super.abort(name, title, msg, text);
	}
		
	public Dataset initDataset(File file) {
		updateJobStatus("20", "reading dataset");
			
		try {
			// reading data
			//
			return new Dataset(file, true);
		} catch (Exception e) {
			abort("filenotfoundexception_initdataset_babelomicstool", "error reading dataset", e.toString(), StringUtils.getStackTrace(e));
		}
		return null;
	}

	public void updateJobStatus(String progress, String message) {
		logger.debug(message + "...\n");
		try {
			jobStatus.addStatusMessage(progress, message);
		} catch (FileNotFoundException e) {
			abort("filenotfoundexception_updatejobstatus_babelomicstools", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
		}
	}
	
	
	
	/**
	 * @param species the species to set
	 */
	public void setSpecies(String species) {
		this.species = species;
	}

	/**
	 * @return the species
	 */
	public String getSpecies() {
		return species;
	}

	
	
}
