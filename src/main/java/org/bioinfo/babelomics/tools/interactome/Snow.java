package org.bioinfo.babelomics.tools.interactome;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;

import org.bioinfo.babelomics.exception.InvalidParameterException;
import org.bioinfo.babelomics.methods.interactome.SnowTest;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.commons.io.utils.FileUtils;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;

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
		File f2 = commandLine.hasOption("list2") ? ("none".equalsIgnoreCase(commandLine.getOptionValue("list2")) ? null : new File(commandLine.getOptionValue("list2"))) :  null;
		String interactome = commandLine.getOptionValue("interactome", "join");
		File interactionsFile = commandLine.hasOption("own-interactions") ? ("none".equalsIgnoreCase(commandLine.getOptionValue("list2")) ? null : new File(commandLine.getOptionValue("own-interactions"))) :  null;
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

			//			File file, outFile;
			//			String fileName, dirName;

			// interactome images
			//
			addResultImage("interactome_betweenness_image", "Interactome betweenness image (png format)", "Interactome images", "interactome_betweenness.png", baseDir);
			addResultImage("interactome_coefficient_image", "Interactome coefficient image (png format)", "Interactome images", "interactome_coefficient.png", baseDir);
			addResultImage("interactome_connections_image", "Interactome connections image (png format)", "Interactome images", "interactome_connections.png", baseDir);

			// network images
			//
			addResultImage("network_betweenness_image", "Network betweenness image (png format)", "Network images", "network_betweenness.png", baseDir);
			addResultImage("network_coefficient_image", "Network coefficient image (png format)", "Network images", "network_coefficient.png", baseDir);
			addResultImage("network_connections_image", "Network connections image (png format)", "Network images", "network_connections.png", baseDir);

			// list #1 files
			//
			addResultFile("list1_bicomp_file", "BiComponents & subnetwork parameters", "List #1 results", "List1_bicomp.txt", baseDir);
			addResultFile("list1_comp_file", "Components & subnetwork parameters", "List #1 results", "List1_comp.txt", baseDir);
			addResultFile("list1_paths_file", "Shortest paths found with a maximun of 1 external proteins introduced", "List #1 results", "List1_paths.txt", baseDir);
			addResultFile("list1_art_file", "BiComponents & articulation points", "List #1 results", "List1_art_points.txt", baseDir);


			if ( f2 != null ) {
				addResultFile("list2_bicomp_file", "BiComponents & subnetwork parameters", "List #2 results", "List2_bicomp.txt", baseDir);
				addResultFile("list2_comp_file", "Components & subnetwork parameters", "List #2 results", "List2_comp.txt", baseDir);
				addResultFile("list2_paths_file", "Shortest paths found with a maximun of 1 external proteins introduced", "List #2 results", "List2_paths.txt", baseDir);
				addResultFile("list2_art_file", "BiComponents & articulation points", "List #2 results", "List2_art_points.txt", baseDir);
			}

			//			result.addOutputItem(new Item("List_file", filename, "MARMITE output file", TYPE.FILE, new ArrayList<String>(), new HashMap<String, String>(), "Results"));
			
			addResultSnowViewer("swnowviewr", "Snow viewer applet (requires Java support)", "Snow viewer", "subnetwork1.xml", baseDir);


		} catch (InvalidParameterException e) {
			printError("snow_executesnow_invalidparameterexception", e.toString(), e.toString(), e);
		} catch (IOException e) {
			printError("snow_executesnow_ioexception", e.toString(), e.toString(), e);		}
	}

	private void addResultImage(String name, String label, String groupName, String fileName, String baseDir) throws IOException {
		File outFile = copyFile(fileName, baseDir);
		if ( outFile != null && outFile.exists() ) {
			result.addOutputItem(new Item(name, outFile.getName(), label, TYPE.IMAGE, new ArrayList<String>(2), new HashMap<String, String>(2), groupName));
		}
	}
	
	private void addResultFile(String name, String label, String groupName, String fileName, String baseDir) throws IOException {
		File outFile = copyFile(fileName, baseDir);
		if ( outFile != null && outFile.exists() ) {
				result.addOutputItem(new Item(name, outFile.getName(), label, TYPE.FILE, new ArrayList<String>(2), new HashMap<String, String>(2), groupName));
		}
	}
	
	private void addResultSnowViewer(String name, String label, String groupName, String fileName, String baseDir) throws IOException {
		File outFile = copyFile(fileName, baseDir);
		if ( outFile != null && outFile.exists() ) {
			String url = "SnowViewer?filename=" + fileName;
			//String url = "http://beta.babelomics.bioinfo.cipf.es/SnowViewer?filename=" + fileName;
			//http://beta.babelomics.bioinfo.cipf.es/SnowViewer?filename=subnetwork1.xml&jobid=262&sessionid=7tySrwJgNRM0tyGOnj9N3eT2BCTsEC8Qw3StPMwwz1WfsFyt1zoIqR6lPa7fnv54
			result.addOutputItem(new Item(name, url, label, TYPE.HTML, StringUtils.toList("SERVER,INCLUDE_REFS", ","), new HashMap<String, String>(2), groupName));
		}
	}

	private File copyFile(String fileName, String baseDir) {
		File outFile = null;
		File file = new File(outdir + baseDir + fileName);
		if ( file.exists() ) {
			outFile = new File(outdir + "/" + fileName);
			try {
				FileUtils.touch(outFile);
				FileUtils.copy(file, outFile);
			} catch (IOException e) {
				outFile = null;
			}
		}
		return outFile;		
	}


}
