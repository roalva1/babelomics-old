package org.bioinfo.babelomics;

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
			System.out.println("tool no valid!!!!");
			return;
		}
		
		
		// if the tool is ok:
		try {
			tool.parse();
			
			tool.run();
		} catch (ParseException e) {
			tool.printUsage("lkjdaslkj");
			e.printStackTrace();
		}
		
	}

}
