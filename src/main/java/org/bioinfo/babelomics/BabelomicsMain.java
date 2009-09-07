package org.bioinfo.babelomics;

import java.io.IOException;

import org.apache.commons.cli.ParseException;
import org.bioinfo.babelomics.tools.BabelomicsFactory;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.commons.log.Logger;
import org.bioinfo.commons.utils.StringUtils;

public class BabelomicsMain {
	
	/**
	 * @param args
	 */
	public static void main(String[] args) {
		Logger logger = new Logger();				
				
		if(args.length == 0) {
			logger.println("No options have been provided");
			printUsage();
		}
		
		String toolName = null;
		// looking for the tool name 
		for(int i=0; i < args.length; i++) {
			if(args[i].equals("--tool") || args[i].equals("-tool")) {
				// this fixes the case that --tool is the last option
				// and no name has been provided
				if(i < args.length) {
					toolName = args[i+1];
					break;
				}
			}
		}
		if(toolName == null) {
			System.out.println("No [--tool] option has been provided!");
			printUsage();
			return;
		}
		
		BabelomicsTool tool = BabelomicsFactory.createTool(toolName);
		if(tool == null) {
			System.out.println("tool '" + toolName + "' is not valid!");
			printUsage();
			return;
		}
		tool.setToolName(toolName);
		
		try {
			// parse the command line
			tool.parse(args, false);
			
			// everything is OK with the commandLine
			tool.run();
			
		} catch (ParseException e) {
			logger.println(e.getMessage());
			logger.println("");
			tool.printUsage("./babelomics.sh");
			tool.abort("parseexception_main_babelomics", "ParseException from Main in command line parse", e.toString(), StringUtils.getStackTrace(e));
		} catch (IOException e) {
			tool.printUsage("./babelomics.sh");	
			tool.abort("ioexception_main_babelomics", "IOException from Main in command line parse", e.toString(), StringUtils.getStackTrace(e));
		}
	}

	public static void printUsage() {
		System.out.println();
		System.out.println("Usage: ./babelomics.sh --tool <arg> --outdir <arg> [--species <arg>] [--log-level <arg>] [--log-file <arg>] [--report <arg>] ");
		System.out.println("\t--tool <arg>			The tool to execute: preprocessing, differential-expression, predictor, clustering, fatigo");
		System.out.println("\t--outdir <arg>		");
		System.out.println("\t--species <arg>		The species: hsa, mmu, rno");
		System.out.println("\t--log-level <arg>		Levels: DEBUG: 1, INFO: 2, WARNING: 3, ERROR: 4, FATAL: 5");
		System.out.println("\t--log-file <arg>		The name of the log file, default: result.log");
		System.out.println("\t--report <arg>		The possible values are: pdf, html or txt");
		System.out.println();
	}
}
