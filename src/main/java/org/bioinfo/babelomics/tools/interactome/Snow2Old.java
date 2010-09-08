package org.bioinfo.babelomics.tools.interactome;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.apache.commons.cli.OptionGroup;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.commons.io.utils.FileUtils;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.graph.SimpleUndirectedGraph;
import org.bioinfo.data.graph.edge.DefaultEdge;
import org.bioinfo.networks.protein.ProteinVertex;
import org.bioinfo.networks.proteinold.ProteinNetworkOld;
import org.bioinfo.tool.OptionFactory;


public class Snow2Old extends BabelomicsTool {

	private File inputFile, inputFileValues;
	private List<String> nodeList1, nodeList2;
	

	@Override
	public void initOptions() {
		options.addOption(OptionFactory.createOption("sif-file", "An input file containing a SIF interactome", false));
		options.addOption(OptionFactory.createOption("file-topo-values", "An input file containing topological values", false, true));

		OptionGroup group = new OptionGroup();
		OptionGroup groupFiles = new OptionGroup();
		
		group.addOption(OptionFactory.createOption("node-list", "A comma separated list of nodes id", false));
		group.addOption(OptionFactory.createOption("node-file1", "An input file containing a node per line. If a list has been provided this option is not read.", false));
		group.addOption(OptionFactory.createOption("node-file2", "An input file containing a node per line. If a list has been provided this option is not read.", false));
		group.addOption(OptionFactory.createOption("randoms", "Number of randoms", false, true));
		options.addOptionGroup(group);
		
		options.addOption(OptionFactory.createOption("random-size", "Size of randoms", false, true));
		options.addOption(OptionFactory.createOption("intermediate", "If there is this argument, it will create the network with 1 intermediate", false, false));
		options.addOption(OptionFactory.createOption("no-relative-betweenness", "If there is this argument, it won't calculate the relative betweenness", false, false));
		options.addOption(OptionFactory.createOption("no-connections", "If there is this argument, it won't calculate the connections", false, false));
		options.addOption(OptionFactory.createOption("no-clustering", "If there is this argument, it won't calculate the clustering", false, false));
		options.addOption(OptionFactory.createOption("no-number-components", "If there is this argument, it won't calculate the number of components", false, false));
		options.addOption(OptionFactory.createOption("no-bicomponents", "If there is this argument, it won't calculate the number of bicomponents", false, false));
		options.addOption(OptionFactory.createOption("dot", "If there is this argument, it will create an output .dot file", false, false));
		options.addOption(OptionFactory.createOption("json", "If there is this argument, it will create an output .json file", false, false));
		options.addOption(OptionFactory.createOption("svg", "If there is this argument, it will create an output .svg and .dot file", false, false));
		options.addOption(OptionFactory.createOption("o-file-components", "If there is this argument, it will create an output .cmp file for the information of each component", false, true));
		groupFiles.addOption(OptionFactory.createOption("o-file-mcn", "Name of the .mcn output file without the extension", false, true));
		groupFiles.addOption(OptionFactory.createOption("o-file-topo", "Name of the .topo output file without the extension", false, true));
		options.addOptionGroup(groupFiles);
		
		options.addOption(OptionFactory.createOption("side", "side for kolmogorov and wilkoxon test. Can be two.sided(by default), less or greater", false, true));
		options.addOption(OptionFactory.createOption("o-stats-file", "Print the statistics file", false, true));
		options.addOption(OptionFactory.createOption("o-images-file", "Print the images for the statistics", false, true));

	}

	@Override
	protected void execute() {
		try {
			boolean intermediate = commandLine.hasOption("intermediate");
			boolean noRelativeBetweenness = commandLine.hasOption("no-relative-betweenness");
			boolean noConnections = commandLine.hasOption("no-connections");
			boolean noClustering = commandLine.hasOption("no-clustering");;
			boolean noNumberOfComponents = commandLine.hasOption("no-number-components");
			boolean noBicomponents = commandLine.hasOption("no-bicomponents");
			boolean json = commandLine.hasOption("json");
			boolean dot = commandLine.hasOption("dot"); 
			boolean svg = commandLine.hasOption("svg");
			String side = "two.sided";

			logger.debug("reading proteins list...");
			if(commandLine.hasOption("node-list")) {
				nodeList1 = StringUtils.toList(commandLine.getOptionValue("node-list"), ",");
				logger.debug("proteins read: " + nodeList1.toString());
			}else {
				if(commandLine.hasOption("node-file1")) {
					FileUtils.checkFile(new File(commandLine.getOptionValue("node-file1")));
					nodeList1 = IOUtils.readLines(commandLine.getOptionValue("node-file1"));
					logger.debug("proteins read: " + nodeList1.toString());
				}
			}
			if(commandLine.hasOption("node-file2")) {
				FileUtils.checkFile(new File(commandLine.getOptionValue("node-file2")));
				nodeList2 = IOUtils.readLines(commandLine.getOptionValue("node-file2"));
				logger.debug("proteins read: " + nodeList2.toString());
			}
			
			if(commandLine.hasOption("sif-file") ){
				inputFile = new File(commandLine.getOptionValue("sif-file"));
				logger.debug("checking if inputFile exist...");
				FileUtils.checkFile(inputFile);
				SimpleUndirectedGraph<ProteinVertex, DefaultEdge> interactomeGraph = ProteinNetworkOld.getInteractomeGraphFilled(inputFile.getAbsolutePath());
				
				// #1:If you have sif-file and you want to create a .topo file from the full sif-file (obviously no randoms)
				if(nodeList1 == null && !commandLine.hasOption("randoms") && commandLine.hasOption("o-file-topo")){
					ProteinNetworkOld.generateParentNetwork(interactomeGraph, commandLine.getOptionValue("o-file-topo"));
					return;
				}
				// #2:If you have sif-file and you want to create a .topo and/or .mcn file random list without a given size
				else if(nodeList1 == null && commandLine.hasOption("randoms") && !commandLine.hasOption("random-size") && ( commandLine.hasOption("o-file-topo") || commandLine.hasOption("o-file-mcn"))){
					ProteinNetworkOld.calculeRandom(interactomeGraph, Integer.parseInt(commandLine.getOptionValue("randoms")), intermediate, noRelativeBetweenness, noConnections, noClustering, 
							noNumberOfComponents, noBicomponents, json, dot, svg, commandLine.getOptionValue("o-file-mcn"), commandLine.getOptionValue("o-file-topo"), commandLine.getOptionValue("o-file-components"));
				}
				// #3:If you have sif-file and you want to create a .topo and/or .mcn file random list with a given size
				else if(nodeList1 == null && commandLine.hasOption("randoms") &&  commandLine.hasOption("random-size") && ( commandLine.hasOption("o-file-topo") || commandLine.hasOption("o-file-mcn"))){
					ProteinNetworkOld.calculeRandom(interactomeGraph, Integer.parseInt(commandLine.getOptionValue("randoms")), Integer.parseInt(commandLine.getOptionValue("random-size")), intermediate, noRelativeBetweenness, noConnections, noClustering, 
							noNumberOfComponents, noBicomponents, json, dot, svg, commandLine.getOptionValue("o-file-mcn"), commandLine.getOptionValue("o-file-topo"), commandLine.getOptionValue("o-file-components"));
				}
				// #4:If you have sif-file and you want to create a .topo and/or .mcn file for the value of your nodes from a given node-list
				else if(nodeList1 != null && (commandLine.hasOption("o-file-topo") || commandLine.hasOption("o-file-mcn") ) && !commandLine.hasOption("o-stats-file")){
					ProteinNetworkOld.calculeRandom(interactomeGraph, nodeList1, intermediate, noRelativeBetweenness, noConnections, noClustering, 
							noNumberOfComponents, noBicomponents, json, dot, svg, commandLine.getOptionValue("o-file-mcn"), commandLine.getOptionValue("o-file-topo"), commandLine.getOptionValue("o-file-components"));
				}
				// #6:If you want to generate statistics a sif-file, topological file and a list of proteins
				else if(commandLine.hasOption("o-stats-file")){
//					if(commandLine.hasOption("o-stats-file"))
//						ProteinNetwork.calculeTopoValues(sif-file, topo-file, list);
					if(commandLine.hasOption("file-topo-values") && nodeList1 != null && commandLine.hasOption("randoms") &&  commandLine.hasOption("random-size") && nodeList2 == null){
						if(commandLine.hasOption("side"))
							side = commandLine.getOptionValue("side");
						FileUtils.checkFile(new File(commandLine.getOptionValue("file-topo-values")));
						inputFileValues = new File(commandLine.getOptionValue("file-topo-values"));
						String outStatsFile = (commandLine.getOptionValue("o-stats-file"));
						String outImagesFile = (commandLine.getOptionValue("o-images-file"));
						ProteinNetworkOld.statisticsTest(inputFileValues, interactomeGraph, nodeList1, Integer.parseInt(commandLine.getOptionValue("randoms")), Integer.parseInt(commandLine.getOptionValue("random-size")),  null, side,  outStatsFile, outImagesFile);
					}
					// #7:If you want to generate statistics a sif-file, topological file and 2 lists of proteins 
					else if(commandLine.hasOption("file-topo-values") && nodeList1 != null && nodeList2 != null && !commandLine.hasOption("randoms") &&  !commandLine.hasOption("random-size")){
						if(commandLine.hasOption("side"))
							side = commandLine.getOptionValue("side");
						FileUtils.checkFile(new File(commandLine.getOptionValue("file-topo-values")));
						inputFileValues = new File(commandLine.getOptionValue("file-topo-values"));
						String outStatsFile = (commandLine.getOptionValue("o-stats-file"));
						String outImagesFile = (commandLine.getOptionValue("o-images-file"));
						ProteinNetworkOld.statisticsTest(inputFileValues, interactomeGraph, nodeList1, Integer.MIN_VALUE, Integer.MIN_VALUE, nodeList2, side, outStatsFile, outImagesFile);
					}
				}
				
			}
			// #5:If you have file-topo-values and you want to create a .mcn file from a given node-list but with the values of the o-file-topo-values
			else if(commandLine.hasOption("file-topo-values") && nodeList1 != null){
				FileUtils.checkFile(new File(commandLine.getOptionValue("file-topo-values")));
				inputFileValues = new File(commandLine.getOptionValue("file-topo-values"));
				ProteinNetworkOld.generateParameters(inputFileValues, nodeList1, commandLine.getOptionValue("o-file-mcn"));
				return;
			}
			
			
		} catch (IOException e) {
			e.printStackTrace();
		}

	}

}
