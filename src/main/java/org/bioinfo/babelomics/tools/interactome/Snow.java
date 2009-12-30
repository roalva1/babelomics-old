package org.bioinfo.babelomics.tools.interactome;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import org.bioinfo.babelomics.exception.InvalidParameterException;
import org.bioinfo.babelomics.methods.interactome.SnowTest;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.commons.io.utils.FileSystemUtils;
import org.bioinfo.commons.io.utils.FileUtils;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;

import sun.management.FileSystem;

public class Snow extends BabelomicsTool {
	
	public Snow() {
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
		
		System.out.println("Snow.java:execute, outdir = " + outdir);
//		try {
			File f1 = new File(commandLine.getOptionValue("list1"));
			File f2 = commandLine.hasOption("list2") ? new File(commandLine.getOptionValue("list2")) :  null;
			String interactome = commandLine.getOptionValue("interactome", "join");
			File interactionsFile = commandLine.hasOption("own-interactions") ? new File(commandLine.getOptionValue("own-interactions")) :  null;
			boolean checkInteractions = commandLine.hasOption("check-interactions");
			String idNature = commandLine.getOptionValue("id-nature", "proteins");
			int interactionsNumber = Integer.parseInt(commandLine.getOptionValue("interactions-number", "1"));
						
			executeSnow(interactome, idNature, interactionsNumber, f1, f2, interactionsFile, checkInteractions);
			
//		} catch (IOException e) {
//			logger.error("Error opening the dataset", e.toString());
//		}
	}
	
	private void executeSnow(String interactome, String idNature, int interactionsNumber, File f1, File f2, File interactionsFile, boolean checkInteractions) {
		SnowTest snowTest = null;
		String snowBinPath = System.getenv("BABELOMICS_HOME") + "/bin/snow/snow.pl";		
		System.out.println("Snow.java:executeSnow, snowBinPath = " + snowBinPath);
		
		// create the snow test object according to the parameters
		//
		if ( f2 == null ) {
			if ( "own".equalsIgnoreCase(interactome) ) {
				// test with one-list and own interactome 
				snowTest = new SnowTest(snowBinPath, f1.getAbsolutePath(), interactionsFile.getAbsolutePath(), checkInteractions, idNature, interactionsNumber, outdir);
			} else {
				// test with one-list
				snowTest = new SnowTest(snowBinPath, f1.getAbsolutePath(), interactome, idNature, interactionsNumber, outdir);
			}
		} else {
			if ( "own".equalsIgnoreCase(interactome) ) {
				// test with one-list and own interactome 
				snowTest = new SnowTest(snowBinPath, f1.getAbsolutePath(), f2.getAbsolutePath(), interactionsFile.getAbsolutePath(), checkInteractions, idNature, interactionsNumber, outdir);
			} else {
				// test with one-list
				snowTest = new SnowTest(snowBinPath, f1.getAbsolutePath(), f2.getAbsolutePath(), interactome, idNature, interactionsNumber, outdir);
			}			
		}
		
		if ( snowTest == null ) {
			abort("snow_executesnow", "invalid parameters", "invalid parameters", "invalid parameters");
		}
		
		try {
			snowTest.run();
			
			String uniqueId = IOUtils.readLines(new File(outdir + "/unique_id.txt")).get(0).trim();
			String baseDir = "/temp/" + uniqueId + "/";
			System.out.println("Snow.java:executeSnow, uniqueId = " + uniqueId);
			
			File file, outFile;
			String fileName, dirName;
			
			// interactome images
			//
			fileName = "interactome_betweenness.png";
			file = new File(outdir + baseDir + fileName);
			if ( file.exists() ) {
				outFile = new File(outdir + "/" + fileName);
				System.out.println("Snow.java:executeSnow, infile = " + file.getAbsolutePath() + ", outfile = " + outFile.getAbsolutePath());
				FileUtils.touch(outFile);
				FileUtils.copy(file, outFile);
				if ( outFile.exists() ) {
					result.addOutputItem(new Item("interactome_betweenness_image", file.getName(), "Interactome betweenness image (png format)", TYPE.IMAGE, new ArrayList<String>(2), new HashMap<String, String>(2), "Interactome images"));
				}
			}

			fileName = "interactome_coefficient.png";
			file = new File(outdir + baseDir + fileName);
			if ( file.exists() ) {
				outFile = new File(outdir + "/" + fileName);
				FileUtils.touch(outFile);
				FileUtils.copy(file, outFile);
				if ( outFile.exists() ) {
					result.addOutputItem(new Item("interactome_coefficient_image", file.getName(), "Interactome coefficient image (png format)", TYPE.IMAGE, new ArrayList<String>(2), new HashMap<String, String>(2), "Interactome images"));
				}
			}

			fileName = "interactome_connections.png";
			file = new File(outdir + baseDir + fileName);
			if ( file.exists() ) {
				outFile = new File(outdir + "/" + fileName);
				FileUtils.touch(outFile);
				FileUtils.copy(file, outFile);
				if ( outFile.exists() ) {
					result.addOutputItem(new Item("interactome_connections_image", file.getName(), "Interactome connections image (png format)", TYPE.IMAGE, new ArrayList<String>(2), new HashMap<String, String>(2), "Interactome images"));
				}
			}

			// network images
			//
			fileName = "network_betweenness.png";
			file = new File(outdir + baseDir + fileName);
			if ( file.exists() ) {
				outFile = new File(outdir + "/" + fileName);
				FileUtils.touch(outFile);
				FileUtils.copy(file, outFile);
				if ( outFile.exists() ) {
					result.addOutputItem(new Item("network_betweenness_image", file.getName(), "Interactome betweenness image (png format)", TYPE.IMAGE, new ArrayList<String>(2), new HashMap<String, String>(2), "Network images"));
				}
			}

			fileName = "network_coefficient.png";
			file = new File(outdir + baseDir + fileName);
			if ( file.exists() ) {
				outFile = new File(outdir + "/" + fileName);
				FileUtils.touch(outFile);
				FileUtils.copy(file, outFile);
				if ( outFile.exists() ) {
					result.addOutputItem(new Item("network_coefficient_image", file.getName(), "Network coefficient image (png format)", TYPE.IMAGE, new ArrayList<String>(2), new HashMap<String, String>(2), "Network images"));
				}
			}

			fileName = "network_connections.png";
			file = new File(outdir + baseDir + fileName);
			if ( file.exists() ) {
				outFile = new File(outdir + "/" + fileName);
				FileUtils.touch(outFile);
				FileUtils.copy(file, outFile);
				if ( outFile.exists() ) {
					result.addOutputItem(new Item("network_connections_image", file.getName(), "Network connections image (png format)", TYPE.IMAGE, new ArrayList<String>(2), new HashMap<String, String>(2), "Network images"));
				}
			}
			
//			result.addOutputItem(new Item("List_file", filename, "MARMITE output file", TYPE.FILE, new ArrayList<String>(), new HashMap<String, String>(), "Results"));

			
		} catch (InvalidParameterException e) {
			printError("snow_executesnow_invalidparameterexception", e.toString(), e.toString(), e);
		} catch (IOException e) {
			printError("snow_executesnow_ioexception", e.toString(), e.toString(), e);		}
	}

	
	
	
}
