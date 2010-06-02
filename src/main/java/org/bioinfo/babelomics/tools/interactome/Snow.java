package org.bioinfo.babelomics.tools.interactome;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;

import org.bioinfo.babelomics.exception.InvalidParameterException;
import org.bioinfo.babelomics.methods.interactome.SnowTest;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.commons.Config;
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
		String snowBinPath = babelomicsHomePath + "/bin/snow/snow.pl";		
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

			String title;
			Config config = new Config(outdir + "/output.properties");
			System.err.println("BET: " + config.getProperty("INTERACTOME.BETWEENNESS.PVALUE"));

						
			// interactome images
			//
			
			if (config==null) {
				abort("exception_execute_snow", "Error", "Internal error accessing to Snow server, please, try later or contact us", "");
			}

			String mode = config.getProperty("MODE");
			String label1 = config.getProperty("LABEL1"), label2 = "Interactome";
			if (!"one-list".equalsIgnoreCase(mode)) {
				label2 = config.getProperty("LABEL2");
			}
			
			title = label1 + ("greater".equalsIgnoreCase(config.getProperty("INTERACTOME.BETWEENNESS.SIDE")) ? " > " : " < ") + label2 + ", p-value = " + config.getProperty("INTERACTOME.BETWEENNESS.PVALUE");
			addResultImage("interactome_betweenness_image", "Betweenness : " + title, "Network parameters evaluation.List's role within interactome of reference", "interactome_betweenness.png", baseDir);

			title = label1 + ("greater".equalsIgnoreCase(config.getProperty("INTERACTOME.CONNECTIONS.SIDE")) ? " > " : " < ") + label2 + ", p-value = " + config.getProperty("INTERACTOME.CONNECTIONS.PVALUE");
			addResultImage("interactome_connections_image", "Connections : " + title, "Network parameters evaluation.List's role within interactome of reference", "interactome_connections.png", baseDir);

			title = label1 + ("greater".equalsIgnoreCase(config.getProperty("INTERACTOME.COEFFICIENT.SIDE")) ? " > " : " < ") + label2 + ", p-value = " + config.getProperty("INTERACTOME.COEFFICIENT.PVALUE");
			addResultImage("interactome_coefficient_image", "Clustering coefficient : " + title, "Network parameters evaluation.List's role within interactome of reference", "interactome_coefficient.png", baseDir);

			// network images
			//
			label2 = "Random";
			if (!"one-list".equalsIgnoreCase(mode)) {
				label2 = config.getProperty("LABEL2");
			}
			
			title = label1 + ("greater".equalsIgnoreCase(config.getProperty("NETWORK.BETWEENNESS.SIDE")) ? " > " : " < ") + label2 + ", p-value = " + config.getProperty("NETWORK.BETWEENNESS.PVALUE");
			addResultImage("network_betweenness_image", "Betweenness : " + title, "Network parameters evaluation.Minimal connected network topological evaluation", "network_betweenness.png", baseDir);

			title = label1 + ("greater".equalsIgnoreCase(config.getProperty("NETWORK.CONNECTIONS.SIDE")) ? " > " : " < ") + label2 + ", p-value = " + config.getProperty("NETWORK.CONNECTIONS.PVALUE");
			addResultImage("network_connections_image", "Network connections : " + title, "Network parameters evaluation.Minimal connected network topological evaluation", "network_connections.png", baseDir);

			title = label1 + ("greater".equalsIgnoreCase(config.getProperty("NETWORK.COEFFICIENT.SIDE")) ? " > " : " < ") + label2 + ", p-value = " + config.getProperty("NETWORK.COEFFICIENT.PVALUE");
			addResultImage("network_coefficient_image", "Clustering coefficient : " + title, "Network parameters evaluation.Minimal connected network topological evaluation", "network_coefficient.png", baseDir);

			// more information
			//
			title = label1 + ": " + config.getProperty("NUMBER_OF_COMPONENTS2") + " [" + config.getProperty("CONF_INTERVAL1") + "]";
			if (f2 != null) {
				title += ", " + label2 + ": " + config.getProperty("NUMBER_OF_COMPONENTS2") + " [" + config.getProperty("CONF_INTERVAL2") + "]";
			}
			result.addOutputItem(new Item("number_of_components", title, "Number of components [95% confidence interval]", Item.TYPE.MESSAGE, new ArrayList<String>(), new HashMap<String,String>(), "Network parameters evaluation.Minimal connected network topological evaluation.More information"));						

			title = label1 + ": " + config.getProperty("NUMBER_OF_COMPONENTS_MORE1");
			if (f2 != null) {
				title += ", " + label2 + ": " + config.getProperty("NUMBER_OF_COMPONENTS_MORE2");				
			}
			result.addOutputItem(new Item("number_of_components_more", title, "Number of components with more than 1 node", Item.TYPE.MESSAGE, new ArrayList<String>(), new HashMap<String,String>(), "Network parameters evaluation.Minimal connected network topological evaluation.More information"));						

			title = label1 + ": " + config.getProperty("NUMBER_OF_BICOMPONENTS1");
			if (f2 != null) {
				title += ", " + label2 + ": " + config.getProperty("NUMBER_OF_BICOMPONENTS2");				
			}
			result.addOutputItem(new Item("number_of_bicomponents", title, "Number of Bicomponents", Item.TYPE.MESSAGE, new ArrayList<String>(), new HashMap<String,String>(), "Network parameters evaluation.Minimal connected network topological evaluation.More information"));						

			title = label1 + ": " + config.getProperty("ARTICULATION_POINTS1");
			if (f2 != null) {
				title += ", " + label2 + ": " + config.getProperty("ARTICULATION_POINTS2");				
			}
			result.addOutputItem(new Item("number_of_articulation_points", title, "Articulation points", Item.TYPE.MESSAGE, new ArrayList<String>(), new HashMap<String,String>(), "Network parameters evaluation.Minimal connected network topological evaluation.More information"));						


			// list #1 files
			//Prots/Genes parameters/info in whole interactome
			
			addResultFile("list1_int_param_file", "Prots/Genes parameters/info in whole interactome", "Topological and functional information.List #1", "List1_int_param.txt", baseDir);			
			addResultFile("list1_input_file", "Prots/Genes from list #1 in network", "Topological and functional information.List #1", "List1_input.txt", baseDir);			
			addResultFile("list1_external_file", "External proteins introduced to the network", "Topological and functional information.List #1", "List1_external.txt", baseDir);			
			addResultFile("list1_paths_file", "Found shortest paths", "Topological and functional information.List #1", "List1_paths.txt", baseDir);
			addResultFile("list1_comp_file", "Components information", "Topological and functional information.List #1", "List1_comp.txt", baseDir);
			addResultFile("list1_bicomp_file", "Bicomponents information", "Topological and functional information.List #1", "List1_bicomp.txt", baseDir);
			addResultFile("list1_art_file", "Articulation points", "Topological and functional information.List #1", "List1_art_points.txt", baseDir);


			if ( f2 != null ) {
				addResultFile("list2_int_param_file", "Prots/Genes parameters/info in whole interactome", "Topological and functional information.List #2", "List2_int_param.txt", baseDir);			
				addResultFile("list2_input_file", "Prots/Genes from list #1 in network", "Topological and functional information.List #2", "List2_input.txt", baseDir);			
				addResultFile("list2_external_file", "External proteins introduced to the network", "Topological and functional information.List #2", "List2_external.txt", baseDir);			
				addResultFile("list2_paths_file", "Found shortest paths", "Topological and functional information.List #2", "List2_paths.txt", baseDir);
				addResultFile("list2_comp_file", "Components information", "Topological and functional information.List #2", "List2_comp.txt", baseDir);
				addResultFile("list2_bicomp_file", "Bicomponents information", "Topological and functional information.List #2", "List2_bicomp.txt", baseDir);
				addResultFile("list2_art_file", "Articulation points", "Topological and functional information.List #2", "List2_art_points.txt", baseDir);
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
