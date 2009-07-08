package org.bioinfo.babelomics.tools;

import org.apache.commons.cli.Option;
import org.apache.commons.cli.ParseException;
import org.bioinfo.tool.GenericBioTool;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;

public abstract class BabelomicsTool extends GenericBioTool {

	protected String[] args;
	
	public BabelomicsTool(String[] args) {
		this.args = args;
		initCommonsOptions();
	}
	
	public abstract void initOptions();
	
	
	private void initCommonsOptions() {
		options.addOption(OptionFactory.createOption("outdir", "o",  "outdir to save the results"));
		options.addOption(OptionFactory.createOption("tool", "to", "tool name", false));
		options.addOption(OptionFactory.createOption("log-file", "name of the log file, default: result.log", false));
		options.addOption(OptionFactory.createOption("log-level", "DEBUG -1, INFO -2, WARNING - 3, ERROR - 4, FATAL - 5", false));
		options.addOption(OptionFactory.createOption("species", "The specie of the ids", false));	
	}
	
	public void parse() throws ParseException {
		commandLine = parser.parse( options, args, true);
		
		if(commandLine.hasOption("outdir")) {
			this.outdir = commandLine.getOptionValue("outdir");
		}
		
		if(commandLine.hasOption("log-level")) {
			logger.setLevel(Integer.parseInt(commandLine.getOptionValue("log-level")));
		}
		
//		if(commandLine.hasOption("help")) {
//			printUsage();
//		}
		writeInputData();
	}
	
	protected void writeInputData() {
		for(Option option: commandLine.getOptions()) {
			result.addInputItem(new Item(option.getLongOpt(), option.getValue(), option.getDescription(), TYPE.MESSAGE));
		}
	}
	
	protected void printUsage() {
		printUsage("./babelomics.sh");
	}
	
}
