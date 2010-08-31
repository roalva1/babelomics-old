package org.bioinfo.babelomics.tools.interactome;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.chart.BoxPlotChart;
import org.bioinfo.commons.io.utils.FileUtils;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.data.graph.SimpleUndirectedGraph;
import org.bioinfo.data.graph.Subgraph;
import org.bioinfo.data.graph.alg.Calc;
import org.bioinfo.data.graph.edge.DefaultEdge;
import org.bioinfo.networks.protein.InteractomeParser;
import org.bioinfo.networks.protein.KSTest;
import org.bioinfo.networks.protein.ProteinNetwork;
import org.bioinfo.networks.protein.ProteinNetworkToFile;
import org.bioinfo.networks.protein.ProteinVertex;
import org.bioinfo.networks.protein.WilcoxonTest;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.jfree.chart.plot.PlotOrientation;

public class Snow2  extends BabelomicsTool{

//	private List<String> nodeList1, nodeList2;
	private File inputFile;
	private String outputFileName;
	private ProteinNetwork proteinNetwork;
	private ProteinNetwork subProteinNetwork1, subProteinNetwork2;
	private List<ProteinNetwork> subProteinNetworkRandoms;
//	private String imagesFile;
	private boolean noNumberOfComponents;
	private boolean bicomponents;
	private boolean images;
	private boolean dot;
	private boolean svg;
	private boolean json;

	@Override
	public void initOptions() {
		options.addOption(OptionFactory.createOption("sif-file","s", "An input file containing a SIF interactome", false));
		options.addOption(OptionFactory.createOption("topo-file","t", "An input file containing topological values", false, true));
		options.addOption(OptionFactory.createOption("node-file1","n", "An input file containing a node per line", false));
		options.addOption(OptionFactory.createOption("node-file2", "An input file containing a node per line", false));
		options.addOption(OptionFactory.createOption("randoms", "Number of randoms", false, true));
		options.addOption(OptionFactory.createOption("randoms-size", "Size of randoms", false, true));
		options.addOption(OptionFactory.createOption("intermediate", "If there is this argument, it will create the network with 1 intermediate", false, false));
		options.addOption(OptionFactory.createOption("no-number-components", "If there is this argument, it won't calculate the number of components", false, false));
		options.addOption(OptionFactory.createOption("bicomponents", "If there is this argument, it will calculate the number of bicomponents", false, false));
		options.addOption(OptionFactory.createOption("dot", "It will create an output .dot file", false, true));
		options.addOption(OptionFactory.createOption("json", "It will create an output .json file", false, true));
		options.addOption(OptionFactory.createOption("svg", "It will create an output .svg and .dot file", false, true));
		options.addOption(OptionFactory.createOption("o-sif-topo-file", "Create a full topological file from a sif file", false));
		options.addOption(OptionFactory.createOption("o-file", "If there is this argument, it will create an output .cmp file for the information of each component", false, true));
		options.addOption(OptionFactory.createOption("side", "side for kolmogorov and wilkoxon test. Can be two.sided(by default), less or greater", false, true));
		options.addOption(OptionFactory.createOption("images", "Print the images for the statistics", false, false));
		
		
		
		
		
//		result.addOutputItem(new Item("o-means-file", exp, "Exponential function", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
//		result.addOutputItem(new Item("o-topo-file", mergeMethod, "Merge replicates", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
//		result.addOutputItem(new Item("o-components-file", filterPercentage, "Filter missing values", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
//		result.addOutputItem(new Item("impute_input_param", ("knn".equalsIgnoreCase(imputeMethod) ? ("knn, k-value = " + commandLine.getOptionValue("k-value", "15")) : imputeMethod) , "Impute missing values", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
//		result.addOutputItem(new Item("filenamefilter_input_param", (filterFilename != null && !"none".equalsIgnoreCase(filterFilename) ? "none" : new File(filterFilename).getName()), "ID-file filter", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
//		result.addOutputItem(new Item("extractid_input_param", "" + (extractIds != null && !extractIds.equalsIgnoreCase("none")), "Extract IDs", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
	
	}
	@Override
	protected void execute() {
		
		subProteinNetworkRandoms = new ArrayList<ProteinNetwork>();
		noNumberOfComponents = commandLine.hasOption("no-number-components");
		bicomponents = commandLine.hasOption("bicomponents");
		images = commandLine.hasOption("images");
		dot = commandLine.hasOption("dot");
		svg = commandLine.hasOption("svg");
		json = commandLine.hasOption("json");
		outputFileName = outdir+commandLine.getOptionValue("o-file");
		
		if(json)
			result.addOutputItem(new Item("json", outputFileName+".json", "Json file", Item.TYPE.FILE, Arrays.asList("json"),new HashMap<String,String>(),"Output data"));
		if(dot)
			result.addOutputItem(new Item("dot", outputFileName+".dot", "Dot file", Item.TYPE.FILE, Arrays.asList("dot"),new HashMap<String,String>(),"Output data"));
		if(svg)
			result.addOutputItem(new Item("svg", outputFileName+".svg", "SVG file", Item.TYPE.FILE, Arrays.asList("svg"),new HashMap<String,String>(),"Output data"));
		
		//Analysis one list
		if(commandLine.hasOption("randoms") && !commandLine.hasOption("node-file2")){
			result.addOutputItem(new Item("list_inter_kol", outputFileName+"_list_inter_kol.txt", "List-Inter kolmogorov", Item.TYPE.FILE, Arrays.asList("stats, kolmogorov"),new HashMap<String,String>(),"Output data"));
			result.addOutputItem(new Item("list_inter_wil", outputFileName+"_list_inter_wil.txt", "List-Inter Wilcoxon", Item.TYPE.FILE, Arrays.asList("stats, wilcoxon"),new HashMap<String,String>(),"Output data"));
			result.addOutputItem(new Item("sn_random_kol", outputFileName+"_sn_random_kol.txt", "List-Inter kolmogorov", Item.TYPE.FILE, Arrays.asList("stats, kolmogorov"),new HashMap<String,String>(),"Output data"));
			result.addOutputItem(new Item("sn_random_wil", outputFileName+"_sn_random_wil.txt", "List-Inter Wilcoxon", Item.TYPE.FILE, Arrays.asList("stats, wilcoxon"),new HashMap<String,String>(),"Output data"));
			
			result.addOutputItem(new Item("image_list_inter_relBet", outputFileName+"_list_inter_relBet", "Image list inter relbet", Item.TYPE.FILE, Arrays.asList("image"),new HashMap<String,String>(),"Output data"));
			result.addOutputItem(new Item("image_list_inter_conn", outputFileName+"_list_inter_conn", "Image list inter conn", Item.TYPE.FILE, Arrays.asList("image"),new HashMap<String,String>(),"Output data"));
			result.addOutputItem(new Item("image_list_inter_clust", outputFileName+"_list_inter_clust", "Image list inter clust", Item.TYPE.FILE, Arrays.asList("image"),new HashMap<String,String>(),"Output data"));
			result.addOutputItem(new Item("image_sn_random_relBet", outputFileName+"_sn_random_relBet", "Image subnet random relbet", Item.TYPE.FILE, Arrays.asList("image"),new HashMap<String,String>(),"Output data"));
			result.addOutputItem(new Item("image_sn_random_conn", outputFileName+"_sn_random_conn", "Image subnet random conn", Item.TYPE.FILE, Arrays.asList("image"),new HashMap<String,String>(),"Output data"));
			result.addOutputItem(new Item("image_sn_random_clust", outputFileName+"_sn_random_clust", "Image subnet random clust", Item.TYPE.FILE, Arrays.asList("image"),new HashMap<String,String>(),"Output data"));

		}
		
		//Analysis two lists
		else if(commandLine.hasOption("node-file1") && commandLine.hasOption("node-file2")){
				
		
			result.addOutputItem(new Item("list1_list2_kol", outputFileName+"_list1_list2_kol.txt", "List1-List2 kolmogorov", Item.TYPE.FILE, Arrays.asList("stats, kolmogorov"),new HashMap<String,String>(),"Output data"));
			result.addOutputItem(new Item("list1_list2_wil", outputFileName+"_list1_list2_wil.txt", "List1-List2 Wilcoxon", Item.TYPE.FILE, Arrays.asList("stats, wilcoxon"),new HashMap<String,String>(),"Output data"));
			result.addOutputItem(new Item("sn1_sn2_kol", outputFileName+"_sn1_sn2_kol.txt", "Subnet1-Subnet2 kolmogorov", Item.TYPE.FILE, Arrays.asList("stats, kolmogorov"),new HashMap<String,String>(),"Output data"));
			result.addOutputItem(new Item("sn1_sn2_wil", outputFileName+"_sn1_sn2_wil.txt", "Subnet1-Subnet2 Wilcoxon", Item.TYPE.FILE, Arrays.asList("stats, wilcoxon"),new HashMap<String,String>(),"Output data"));
		
			result.addOutputItem(new Item("image_list1_list2_relBet", outputFileName+"_list1_list2_relBet", "Image list1-list2 inter relbet", Item.TYPE.FILE, Arrays.asList("image"),new HashMap<String,String>(),"Output data"));
			result.addOutputItem(new Item("image_list1_list2_conn", outputFileName+"_list1_list2_conn", "Image list1-list2 conn", Item.TYPE.FILE, Arrays.asList("image"),new HashMap<String,String>(),"Output data"));
			result.addOutputItem(new Item("image_list1_list2_clust", outputFileName+"_list1_list2_clust", "Image list1-list2 clust", Item.TYPE.FILE, Arrays.asList("image"),new HashMap<String,String>(),"Output data"));
			result.addOutputItem(new Item("image_sn1_sn2_relBet", outputFileName+"_sn1_sn2_relBet", "Image subnet1-subnet2 random relbet", Item.TYPE.FILE, Arrays.asList("image"),new HashMap<String,String>(),"Output data"));
			result.addOutputItem(new Item("image_sn1_sn2_conn", outputFileName+"_sn1_sn2_conn", "Image subnet1-subnet2 conn", Item.TYPE.FILE, Arrays.asList("image"),new HashMap<String,String>(),"Output data"));
			result.addOutputItem(new Item("image_sn1_sn2_clust", outputFileName+"_sn1_sn2_clust", "Image subnet1-subnet2 clust", Item.TYPE.FILE, Arrays.asList("image"),new HashMap<String,String>(),"Output data"));
		}
		//
		if(commandLine.hasOption("randoms")){
			int randoms = Integer.parseInt(commandLine.getOptionValue("randoms"));
			result.addOutputItem(new Item("randoms means", outputFileName+"_sn_1-"+(randoms)+"_means.txt", "Random: means", Item.TYPE.FILE, Arrays.asList("randoms","means"),new HashMap<String,String>(),"Output data"));
			result.addOutputItem(new Item("randoms topo", outputFileName+"_sn_1-"+(randoms)+"_topo.txt", "Random: topological values", Item.TYPE.FILE, Arrays.asList("randoms","topological"),new HashMap<String,String>(),"Output data"));
			result.addOutputItem(new Item("randoms components", outputFileName+"_sn_1-"+(randoms)+"_comp.txt", "Randoms: components", Item.TYPE.FILE, Arrays.asList("randoms","components"),new HashMap<String,String>(),"Output data"));
		}
		
		
		String side = "two.sided";
		if(commandLine.hasOption("side"))
			side = commandLine.getOptionValue("side");

		ProteinNetworkToFile file = new ProteinNetworkToFile();

		try {
		
			String folderInteractions = "/opt/babelomics/data/interactions/";
			SimpleUndirectedGraph<ProteinVertex, DefaultEdge> interactomeGraph;
			if(commandLine.hasOption("sif-file") ){
				inputFile = new File(commandLine.getOptionValue("sif-file"));
				logger.debug("checking if inputFile exist...");
				FileUtils.checkFile(inputFile);
				interactomeGraph = InteractomeParser.parseFromSifFile(commandLine.getOptionValue("sif-file") );
				proteinNetwork = new ProteinNetwork(interactomeGraph);
			}
			else{
				
				if(species.equals("sce"))
					interactomeGraph = InteractomeParser.parseFromSifFile(folderInteractions+"sce_alldb_proteins_interactome_nr.sif");
				else if(species.equals("hsa"))
					interactomeGraph = InteractomeParser.parseFromSifFile(folderInteractions+"hsa_alldb_proteins_interactome_nr.sif");
				else{
					System.err.println("Unknown specie");
					return;
				}
				proteinNetwork = new ProteinNetwork(interactomeGraph);
					
			}
			if(commandLine.hasOption("topo-file") ){
				FileUtils.checkFile(new File(commandLine.getOptionValue("topo-file")));
				proteinNetwork.loadTopologicalValues(commandLine.getOptionValue("topo-file"));
			}
			else{
				if(species.equals("sce"))
					proteinNetwork.loadTopologicalValues(folderInteractions+"sce_alldb_proteins_interactome_nr_topo.txt");
				else if(species.equals("hsa"))
					proteinNetwork.loadTopologicalValues(folderInteractions+"hsa_alldb_proteins_interactome_nr_topo.txt");
				else{
					System.err.println("Unknown specie");
					return;
				}
			}
			
			
			
			if(commandLine.hasOption("node-file1")) {
				subProteinNetwork1 = createNodeFile(commandLine.getOptionValue("node-file1"), 1);
			}
			if(commandLine.hasOption("node-file2")) {
				subProteinNetwork2 = createNodeFile(commandLine.getOptionValue("node-file2"), 2);
			}
			if(commandLine.hasOption("o-sif-topo-file")) {
				proteinNetwork.calcTopologicalValues();
				file.toTopologicalFile(commandLine.getOptionValue("o-sif-topo-file")+"_topo.txt", proteinNetwork);
			}
			if(commandLine.hasOption("randoms")){
				int randomSize = 0;
				int randoms = Integer.parseInt(commandLine.getOptionValue("randoms"));
				if(!commandLine.hasOption("randoms-size")){
					int size = this.proteinNetwork.getInteractomeGraph().getVertices().size();
					while(randomSize == 0)
						randomSize = (int) (Math.random() * size);
				}
				else{
					randomSize = Integer.parseInt(commandLine.getOptionValue("randoms-size"));
				}
				createRandoms(randoms, randomSize);
			}
			if(commandLine.hasOption("randoms") && !commandLine.hasOption("node-file2"))
				statsOneListAnalysis(side);
			else if(commandLine.hasOption("node-file1") && commandLine.hasOption("node-file2"))
				statsTwoListsAnalisys(side);
			else
				System.err.println("Not correct arguments for statistic test");
			if(json){
				IOUtils.write(outputFileName+".json", proteinNetwork.getInteractomeGraph().toJson());
			}
			if(dot && !svg){
				IOUtils.write(outputFileName+".dot", proteinNetwork.getInteractomeGraph().toDot());
			}	
			else if(!dot && svg){
				createSVGFile(outputFileName+".dot");
			}
			else if(dot && svg){
				createSVGFile(outputFileName+".dot");
			}
			
				
		} catch (IOException e) {
			e.printStackTrace();
		}
			
	}
	
	private void createSVGFile(String sourceDotFile) throws IOException{
		IOUtils.write(sourceDotFile, proteinNetwork.getInteractomeGraph().toDot());
		IOUtils.write(commandLine.getOptionValue("o-svg-file")+".svg", proteinNetwork.getInteractomeGraph().toSVG(sourceDotFile));
	}
	private void createImages(String fileName, List<Double> list1, String legend1, List<Double> list2, String legend2) throws IOException{
		BoxPlotChart bpc = new BoxPlotChart("", "", "");
		bpc.getLegend().setVisible(true);
		bpc.addSeries(list2, legend2, legend2);
		bpc.addSeries(list1, legend1, legend1);
		bpc.setOrientation(PlotOrientation.HORIZONTAL);
		bpc.save(fileName+".png");
	}
	private void statsOneListAnalysis(String side) throws IOException{
		// 1st Analysis
		
		WilcoxonTest wTest = new WilcoxonTest();
		List<String> list = IOUtils.readLines(commandLine.getOptionValue("node-file1"));
		
		
		List<Double> relBetInter = proteinNetwork.getBetRelList();
		List<Double> relBetList1 = new ArrayList<Double>();
		
		List<Double> connInter = proteinNetwork.getConnList();
		List<Double> connList1 = new ArrayList<Double>();
		
		List<Double> clustInter = proteinNetwork.getClustList();
		List<Double> clustList1 = new ArrayList<Double>();

		createTopoFilterList(list, relBetList1, connList1, clustList1);

		String toWrite = formatStatsFile(relBetList1, relBetInter, connList1, connInter, clustList1, clustInter, side);
		
		IOUtils.write(outputFileName+"_list_inter_kol.txt", toWrite);
		wTest.test(outputFileName+"_list_inter_wil.txt", side, "list1", "inter", relBetList1, relBetInter, connList1, connInter, clustList1, clustInter);
			
		if(images){
			createImages(outputFileName+"_list_inter_relBet", relBetList1, "list1", relBetInter, "inter");
			createImages(outputFileName+"_list_inter_conn", connList1, "list1", connInter, "inter");
			createImages(outputFileName+"_list_inter_clust", clustList1, "list1", clustInter, "inter");
		}
		//2nd Analysis
		List<Double> relBetSubnet1 = subProteinNetwork1.getBetRelList();
		List<Double> relBetRandoms = new ArrayList<Double>();
		
		List<Double> connSubnet1 = subProteinNetwork1.getConnList();
		List<Double> connRandoms = new ArrayList<Double>();
		
		List<Double> clustSubnet1 = subProteinNetwork1.getClustList();
		List<Double> clustRandoms = new ArrayList<Double>();

		for(ProteinNetwork proteinNetwork : subProteinNetworkRandoms){
			relBetRandoms.add(proteinNetwork.getMeanRelBet());
			connRandoms.add(proteinNetwork.getMeanConnections());
			clustRandoms.add(proteinNetwork.getMeanClust());
		}
		
		toWrite = formatStatsFile(relBetSubnet1, relBetRandoms, connSubnet1, connRandoms, clustSubnet1, clustRandoms, side);
		IOUtils.write(outputFileName+"_sn_random_kol.txt", toWrite);		
		wTest.test(outputFileName+"_sn_random_wil.txt", side, "sn1", "rnd", relBetSubnet1, relBetRandoms, connSubnet1, connRandoms, clustSubnet1, clustRandoms);
		if(images){
			createImages(outputFileName+"_sn_random_relBet", relBetSubnet1, "subnet1", relBetRandoms, "randoms");
			createImages(outputFileName+"_sn_random_conn", connSubnet1, "subnet1", connRandoms, "randoms");
			createImages(outputFileName+"_sn_random_clust", clustSubnet1, "subnet1", clustRandoms, "randoms");
		}
	}
	private void statsTwoListsAnalisys(String side) throws IOException{
		// 1st Analysis
		
		WilcoxonTest wTest = new WilcoxonTest();
		List<String> list1 = IOUtils.readLines(commandLine.getOptionValue("node-file1"));
		List<String> list2 = IOUtils.readLines(commandLine.getOptionValue("node-file2"));
		
		List<Double> relBetList1 = new ArrayList<Double>();
		List<Double> connList1 = new ArrayList<Double>();
		List<Double> clustList1 = new ArrayList<Double>();
		List<Double> relBetList2 = new ArrayList<Double>();
		List<Double> connList2 = new ArrayList<Double>();
		List<Double> clustList2 = new ArrayList<Double>();
		
		
		createTopoFilterList(list1, relBetList1, connList1, clustList1);
		createTopoFilterList(list2, relBetList2, connList2, clustList2);
		
		String toWrite = formatStatsFile(relBetList1, relBetList2, connList1, connList2, clustList1, clustList2, side);
		IOUtils.write(outputFileName+"_list1_list2_kol.txt", toWrite);
		wTest.test(outputFileName+"_list1_list2_wil.txt", side, "list1", "list2", relBetList1, relBetList2, connList1, connList2, clustList1, clustList2);
		if(images){
			createImages(outputFileName+"_list1_list2_relBet", relBetList1, "list1", relBetList2, "list2");
			createImages(outputFileName+"_list1_list2_conn", connList1, "list1", connList2, "list2");
			createImages(outputFileName+"_list1_list2_clust", clustList1, "list1", clustList2, "list2");
		}
		
		//2nd Analysis
		List<Double> relBetSubnet1 = subProteinNetwork1.getBetRelList();
		List<Double> relBetSubnet2 = subProteinNetwork2.getBetRelList();
		
		List<Double> connSubnet1 = subProteinNetwork1.getConnList();
		List<Double> connSubnet2 = subProteinNetwork2.getConnList();
		
		List<Double> clustSubnet1 = subProteinNetwork1.getClustList();
		List<Double> clustSubnet2 = subProteinNetwork2.getClustList();
		
		toWrite = formatStatsFile(relBetSubnet1, relBetSubnet2, connSubnet1, connSubnet2, clustSubnet1, clustSubnet2, side);
		IOUtils.write(outputFileName+"_sn1_sn2_kol.txt", toWrite);		
		wTest.test(outputFileName+"_sn1_sn2_wil.txt", side, "sn1", "sn2", relBetSubnet1, relBetSubnet2, connSubnet1, connSubnet2, clustSubnet1, clustSubnet2);
		if(images){
			createImages(outputFileName+"_sn1_sn2_relBet", relBetSubnet1, "subnet1", relBetSubnet2, "subnet2");
			createImages(outputFileName+"_sn1_sn2_conn", connSubnet1, "subnet1", connSubnet2, "subnet2");
			createImages(outputFileName+"_sn1_sn2_clust", clustSubnet1, "subnet1", clustSubnet2, "subnet2");
		}
	
	}
	private String formatStatsFile(List<Double> relBetList1, List<Double> relBetList2, List<Double> connList1, List<Double> connList2, List<Double> clustList1, List<Double> clustList2, String side) throws IOException{
		KSTest kstest = new KSTest();
		StringBuilder sb = new StringBuilder();
		sb.append("#parameter\tpval\tside\n");
		sb.append("betweenness\t"+kstest.resultKolmogorovSmirnov(ListUtils.toArray(relBetList1), ListUtils.toArray(relBetList2), side).getPValue()+"\t"+side+"\n");
		sb.append("connections\t"+kstest.resultKolmogorovSmirnov(ListUtils.toArray(connList1), ListUtils.toArray(connList2), side).getPValue()+"\t"+side+"\n");
		sb.append("coefficient\t"+kstest.resultKolmogorovSmirnov(ListUtils.toArray(clustList1), ListUtils.toArray(clustList2), side).getPValue()+"\t"+side);
		return sb.toString();
	}
	private void createTopoFilterList(List<String> list, List<Double> relBetList, List<Double> connList, List<Double> clustList){
		for(String proteinName : list){
			relBetList.add(proteinNetwork.getRelBetweennessVertex(proteinName));
			connList.add(proteinNetwork.getDegreeVertex(proteinName));
			clustList.add(proteinNetwork.getClusterinVertex(proteinName));
		}
	}
	private void createRandoms(int randoms, int randomSize) throws IOException{
		StringBuilder sbMeans = createMeansHeader();
		StringBuilder sbTopo = createTopoHeader();
		StringBuilder sbComponents =createComponentsHeader();
		
		for(int i=1; i<=randoms; i++){
			SimpleUndirectedGraph<ProteinVertex, DefaultEdge> subgraph = (SimpleUndirectedGraph<ProteinVertex, DefaultEdge>) Subgraph.randomSubgraph(proteinNetwork.getInteractomeGraph(), randomSize);
			if(commandLine.hasOption("intermediate"))
				subgraph = (SimpleUndirectedGraph<ProteinVertex, DefaultEdge>) Subgraph.OneIntermediate(proteinNetwork.getInteractomeGraph(), subgraph);
			
			ProteinNetwork subProteinNetwork = createSubnet(subgraph);
			sbTopo.append(getTopologicalValues(subProteinNetwork, i)).append(System.getProperty("line.separator"));
			sbMeans.append(getTopologicalMeanValues(subProteinNetwork, i)).append(System.getProperty("line.separator"));

			if(!noNumberOfComponents) {
				sbComponents.append(getComponentsValues(subProteinNetwork, i));
			}
			subProteinNetworkRandoms.add(subProteinNetwork);
		}
		randomFilesToString(sbMeans, outputFileName+"_sn_1-"+(randoms)+"_means.txt");
		randomFilesToString(sbTopo,outputFileName+"_sn_1-"+(randoms)+"_topo.txt");
		randomFilesToString(sbComponents,outputFileName+"_sn_1-"+(randoms)+"_comp.txt");
	}
	
	
	private ProteinNetwork createNodeFile(String nodeFile, int node) throws IOException{
		
		FileUtils.checkFile(new File(nodeFile));
		List<String> list = IOUtils.readLines(nodeFile);
		logger.debug("nodes read: " + list.toString());

		StringBuilder sbMeans = createMeansHeader();
		StringBuilder sbTopo = createTopoHeader();
		StringBuilder sbComponents = createComponentsHeader();
		SimpleUndirectedGraph<ProteinVertex, DefaultEdge> subgraph = (SimpleUndirectedGraph<ProteinVertex, DefaultEdge>) Subgraph.randomSubgraph(proteinNetwork.getInteractomeGraph(), toProteinVertex(list));
		
		if(commandLine.hasOption("intermediate"))
			subgraph = (SimpleUndirectedGraph<ProteinVertex, DefaultEdge>) Subgraph.OneIntermediate(proteinNetwork.getInteractomeGraph(), subgraph);
		
		ProteinNetwork subProteinNetwork = createSubnet(subgraph);
		
		IOUtils.write(outputFileName+"_sn_nodeFile"+node+"_topo.txt", sbTopo.toString());
		result.addOutputItem(new Item("sn_nodeFile"+node+"_topo", outputFileName+"_sn_nodeFile"+node+"_topo.txt", "Subnet topo values", Item.TYPE.FILE, Arrays.asList("subnet","topological"),new HashMap<String,String>(),"Output data"));
		
		sbMeans.append(getTopologicalMeanValues(subProteinNetwork, node));
		IOUtils.write(outputFileName+"_sn_nodeFile"+node+"_means.txt", sbMeans.toString());
		result.addOutputItem(new Item("sn_nodeFile"+node+"_means", outputFileName+"_sn_nodeFile"+node+"_means.txt", "Subnet means values", Item.TYPE.FILE, Arrays.asList("subnet","means"),new HashMap<String,String>(),"Output data"));

		
		if(!noNumberOfComponents) {
			sbComponents.append(getComponentsValues(subProteinNetwork, node));
			sbComponents.deleteCharAt(sbComponents.lastIndexOf(System.getProperty("line.separator")));
			IOUtils.write(outputFileName+"_sn_nodeFile"+node+"_comp.txt", sbComponents.toString());
			result.addOutputItem(new Item("sn_nodeFile"+node+"_components", outputFileName+"_sn_nodeFile"+node+"_comp.txt", "Subnet components", Item.TYPE.FILE, Arrays.asList("subnet","components"),new HashMap<String,String>(),"Output data"));

		}
		return subProteinNetwork;
	}
	private ProteinNetwork createSubnet(SimpleUndirectedGraph<ProteinVertex, DefaultEdge> subgraph) {
		ProteinNetwork subProteinNetwork = new ProteinNetwork(subgraph);
		subProteinNetwork.calcTopologicalValues();
		subProteinNetwork.calcTopologicalMeanValues();
		return subProteinNetwork;
	}
	
	private StringBuilder createMeansHeader(){
		StringBuilder sb = new StringBuilder();
		sb.append("#Subnet\tMeanBet\tStdBet\t").append("MeanCon\tStdCon\t").append("MeanCls\tStdCls\t");
		if(!noNumberOfComponents){
			sb.append("Comp\t1Comp");
		}
		if(bicomponents){
			sb.append("\t1BiComp\t");
		}
		sb.append(System.getProperty("line.separator"));
		return sb;
	}
	
	private StringBuilder createTopoHeader(){
		StringBuilder sb = new StringBuilder();
		sb.append("#Subnet\tId\tBet\tClust").append(System.getProperty("line.separator"));
		return sb;
	}
	
	private StringBuilder createComponentsHeader(){
		StringBuilder sb = new StringBuilder();
		sb.append("#Subnet\tComp\tDia\tSize\tNod").append(System.getProperty("line.separator"));
		return sb;
	}
	
	
	private String getTopologicalValues(ProteinNetwork subProteinNetwork, int subnet){
		StringBuilder sb = new StringBuilder();
		for (ProteinVertex proteinVertex : subProteinNetwork.getInteractomeGraph().getVertices()) {
			if(proteinVertex == null)
				continue;
			sb.append("sn"+(subnet)).append("\t");
			sb.append(proteinVertex.getId()).append("\t");
			sb.append(proteinVertex.getRelativeBetweenness()).append("\t").append(subProteinNetwork.getInteractomeGraph().getDegree(proteinVertex)).append("\t").append(proteinVertex.getClusteringCoefficient());
			sb.append(System.getProperty("line.separator"));
		}
		sb.deleteCharAt(sb.length()-1);
		return sb.toString();
	}
	
	private String getTopologicalMeanValues(ProteinNetwork subProteinNetwork, int subnet){
		StringBuilder sb = new StringBuilder();
		SimpleUndirectedGraph<ProteinVertex, DefaultEdge> interactomeGraph = subProteinNetwork.getInteractomeGraph();
		sb.append("sn"+(subnet)).append("\t");
		sb.append(subProteinNetwork.getTopologicalMeanValuesToString());

		if(!noNumberOfComponents){
			int oneComponent = 0;
			int moreOneComponent = 0;
			List<List<ProteinVertex>> listComponents = interactomeGraph.getAllInformationComponents(true);
			for(List<ProteinVertex> listVertex : listComponents){
				if(listVertex.size() == 1)
					oneComponent++;
				else if(listVertex.size() > 1)
					moreOneComponent++;
			}
			sb.append("\t"+moreOneComponent+"\t"+oneComponent);
		}
		if(bicomponents){
			sb.append("\t").append(interactomeGraph.getNumberOfBicomponents());
		}
			
		return sb.toString();
		
	}
	private String getComponentsValues(ProteinNetwork subProteinNetwork, int subnet){
		StringBuilder sb = new StringBuilder();
		SimpleUndirectedGraph<ProteinVertex, DefaultEdge> interactomeGraph = subProteinNetwork.getInteractomeGraph();
		
		if(!noNumberOfComponents){
			List<List<ProteinVertex>> componentsList =  interactomeGraph.getAllInformationComponents(true);
			List<Double> componentsDiameter = Calc.calcDiameter(interactomeGraph, componentsList);
			for(int i=0; i < componentsList.size(); i++){
				sb.append("sn"+(subnet)).append("\t");
				StringBuilder sbNodes = new StringBuilder();
				for (int k = 0; k < componentsList.get(i).size(); k++) {
					sbNodes.append(componentsList.get(i).get(k));
					if(k!=componentsList.get(i).size()-1)
						sbNodes.append(",");
				}
				sb.append((i)+"\t"+componentsDiameter.get(i)+"\t"+componentsList.get(i).size()+"\t"+sbNodes.toString()+ System.getProperty("line.separator"));
			}
		}
		//sb.deleteCharAt(sb.lastIndexOf(System.getProperty("line.separator")));
		return sb.toString();
	}
	
	private void randomFilesToString(StringBuilder sb, String file) throws IOException{
		if(!sb.toString().equals("")){
			sb.deleteCharAt(sb.lastIndexOf(System.getProperty("line.separator")));
			IOUtils.write(file, sb.toString());
		}
	}
	private List<ProteinVertex> toProteinVertex(List<String> vertices){
		List<ProteinVertex> verticesList = new ArrayList<ProteinVertex>(vertices.size());
		for(String proteinName : vertices){
//			if(verticesList.contains(new ProteinVertex(proteinName)))
//				System.out.println("Contiene:" +proteinName);
			if(!proteinName.equals(""))
				verticesList.add(new ProteinVertex(proteinName));
		}
		return verticesList;
	}
}
