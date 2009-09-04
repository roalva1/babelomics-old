package org.bioinfo.babelomics.tools;

import java.io.IOException;

import org.apache.commons.cli.ParseException;
import org.bioinfo.tool.GenericBioTool;
import org.bioinfo.tool.OptionFactory;

public abstract class BabelomicsTool extends GenericBioTool {
	
	protected String species;
	
	
	public BabelomicsTool() {
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
