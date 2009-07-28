package org.bioinfo.babelomics;

import java.io.IOException;

import org.apache.commons.cli.ParseException;
import org.bioinfo.babelomics.tools.BabelomicsFactory;
import org.bioinfo.babelomics.tools.BabelomicsTool;

public class BabelomicsMain {

	/**
	 * @param args
	 */
	public static void main(String[] args) {
		
		if(args.length == 0) {
			System.out.println("No options provided");
		}
		
		String toolName = null;
		// looking for the tool name and stoping 
		for(int i=0; i < args.length; i++) {
			if(args[i].equals("--tool") || args[i].equals("-tool")) {
				toolName = args[i+1];
				break;
			}
		}
		
		if(toolName == null) {
			System.out.println("no tool provided!!!!");
			return;
		}
		
		BabelomicsTool tool = BabelomicsFactory.createTool(toolName, args);
		if(tool == null) {
			System.out.println("tool no valid!!!! ==> "+toolName);
			return;
		}
		System.out.println("==> "+toolName);
		tool.setToolName(toolName);
		
		
		// if the tool is ok:
		try {
			tool.parse(args, false);
			
			tool.run();
		} catch (ParseException e) {
			tool.printUsage("lkjdaslkj");
			tool.printError("parseexception_main_babelomicsmain", "An error while parsing the command line", e.toString(), e);
		} catch (IOException e) {
			tool.printError("ioexception_main_babelomicsmain", "An error while parsing the command line", e.toString(), e);
		}
		
	}

}
