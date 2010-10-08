package org.bioinfo.babelomics.tools.interactome;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Set;

import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.chart.BoxPlotChart;
import org.bioinfo.commons.io.TextFileReader;
import org.bioinfo.commons.io.utils.FileUtils;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.commons.utils.StringUtils;
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
import org.bioinfo.networks.protein.files.Dot;
import org.bioinfo.networks.protein.files.Json;
import org.bioinfo.networks.protein.files.Svg;
import org.bioinfo.networks.protein.files.Xml;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;
import org.jfree.chart.plot.PlotOrientation;

public class Snow2  extends BabelomicsTool{

	private File inputFile;
	private String interactome;
	private String outputFileName;
	private String type;
	private ProteinNetwork proteinNetwork;
	private ProteinNetwork subProteinNetwork1, subProteinNetwork2;
	private List<ProteinNetwork> subProteinNetworkRandoms;
	private boolean noNumberOfComponents;
	private boolean bicomponents;
	private boolean images;
	private boolean json;
	private boolean xml;
	private Set<String> intermediatesSub1, intermediatesSub2; 
	private List<List<ProteinVertex>> componentsListSub1, componentsListSub2;
	private int randomSize;

	private String wBinPath;

	@Override
	public void initOptions() {
		options.addOption(OptionFactory.createOption("interactome", "Interactome: hsa, sce or own (for you own interactome)", true));
		options.addOption(OptionFactory.createOption("sif-file", "An input file containing a SIF interactome", false));
		options.addOption(OptionFactory.createOption("topo-file","t", "An input file containing topological values", false, true));
		options.addOption(OptionFactory.createOption("list1", "An input file containing a node per line", false));
		options.addOption(OptionFactory.createOption("list2", "An input file containing a node per line", false));
		options.addOption(OptionFactory.createOption("type", "An argument saying if you want genes or proteins(default)", false, true));
		options.addOption(OptionFactory.createOption("randoms", "Number of randoms", false, true));
//		options.addOption(OptionFactory.createOption("randoms-size", "Size of randoms", false, true));
		options.addOption(OptionFactory.createOption("intermediate", "If there is this argument, it will create the network with 1 intermediate", false, false));
		options.addOption(OptionFactory.createOption("no-number-components", "If there is this argument, it won't calculate the number of components", false, false));
		options.addOption(OptionFactory.createOption("bicomponents", "If there is this argument, it will calculate the number of bicomponents", false, false));
		options.addOption(OptionFactory.createOption("json", "It will create an output .json file", false, false));
		options.addOption(OptionFactory.createOption("o-sif-topo-file", "Create a full topological file from a sif file", false, false));
		options.addOption(OptionFactory.createOption("o-name", "If there is this argument, it will create an output .cmp file for the information of each component", false, true));
		options.addOption(OptionFactory.createOption("side", "side for kolmogorov and wilkoxon test. Can be two.sided(by default), less or greater", false, true));
		options.addOption(OptionFactory.createOption("images", "Print the images for the statistics", false, false));
		options.addOption(OptionFactory.createOption("xml", "Output xml file with the representation of the graph", false, false));

	}
	@Override
	protected void execute() {
		File f = null;

		wBinPath = babelomicsHomePath + "/bin/snow/wilcoxtest.r";
		subProteinNetwork1 = null;
		subProteinNetwork2 = null;

		subProteinNetworkRandoms = new ArrayList<ProteinNetwork>();
		noNumberOfComponents = commandLine.hasOption("no-number-components");
		bicomponents = commandLine.hasOption("bicomponents");
		images = commandLine.hasOption("images");
		json = commandLine.hasOption("json");
		xml = commandLine.hasOption("xml");

		outputFileName = outdir + "/" + commandLine.getOptionValue("o-name");

		String interactomeMsg = "Homo sapiens";
		interactome = commandLine.getOptionValue("interactome");
		
		if ("hsa".equalsIgnoreCase(interactome)) {
			interactomeMsg = "Homo sapiens";
		} else if ("sce".equalsIgnoreCase(interactome)) {
			interactomeMsg = "Saccharomyce cerevisiae";
		} else if ("own".equalsIgnoreCase(interactome)) {
			interactomeMsg = new File(commandLine.getOptionValue("sif-file")).getName();
		} else {
			interactomeMsg = "unknown";
		}
		result.addOutputItem(new Item("interactome_param", interactomeMsg, "Interactome", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));

		String side = "less";
		if(commandLine.hasOption("side")) {
			side = commandLine.getOptionValue("side");
		}
		result.addOutputItem(new Item("side_param", side, "Side for kolmogorov and wilkoxon test", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));

		type = "proteins";		
		if(commandLine.hasOption("type")){
			type = commandLine.getOptionValue("type");
		}
		result.addOutputItem(new Item("type_param", type, "ID type", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));


		ProteinNetworkToFile file = new ProteinNetworkToFile();

		try {

			String folderInteractions = this.babelomicsHomePath + "/conf/interactions/";
			SimpleUndirectedGraph<ProteinVertex, DefaultEdge> interactomeGraph = null;

			
			if (interactome.equals("sce")) {
				if(type.equals("proteins")) {
					interactomeGraph = InteractomeParser.parseFromSifFile(folderInteractions+"sce_alldb_proteins_interactome_nr.sif");
					proteinNetwork = new ProteinNetwork(interactomeGraph);
					proteinNetwork.loadTopologicalValues(folderInteractions+"sce_alldb_proteins_interactome_nr_topo.txt");
				} else if(type.equals("genes")) {
					interactomeGraph = InteractomeParser.parseFromSifFile(folderInteractions+"sce_alldb_genes_interactome_nr.sif");
					proteinNetwork = new ProteinNetwork(interactomeGraph);
					proteinNetwork.loadTopologicalValues(folderInteractions+"sce_alldb_genes_interactome_nr_topo.txt");
				} else{
					System.err.println("Unknown type");
					return;
				}
			} else if (interactome.equals("hsa")) {
				if(type.equals("proteins")) {
					interactomeGraph = InteractomeParser.parseFromSifFile(folderInteractions+"hsa_alldb_proteins_interactome_nr.sif");
					proteinNetwork = new ProteinNetwork(interactomeGraph);
					proteinNetwork.loadTopologicalValues(folderInteractions+"hsa_alldb_proteins_interactome_nr_topo.txt");
				} else if(type.equals("genes")) {
					interactomeGraph = InteractomeParser.parseFromSifFile(folderInteractions+"hsa_alldb_genes_interactome_nr.sif");
					proteinNetwork = new ProteinNetwork(interactomeGraph);
					proteinNetwork.loadTopologicalValues(folderInteractions+"hsa_alldb_genes_interactome_nr_topo.txt");
				} else{
					System.err.println("Unknown type");
					return;
				}
			} else if (interactome.equals("own")){
				if(commandLine.hasOption("sif-file") ){
					inputFile = new File(commandLine.getOptionValue("sif-file"));
					logger.debug("checking if inputFile exist...");
					FileUtils.checkFile(inputFile);
					interactomeGraph = InteractomeParser.parseFromSifFile(commandLine.getOptionValue("sif-file") );
					proteinNetwork = new ProteinNetwork(interactomeGraph);
					
					if(commandLine.hasOption("topo-file"))
						proteinNetwork.loadTopologicalValues(commandLine.getOptionValue("topo-file"));
					else
						proteinNetwork.calcTopologicalValues();
					if(commandLine.hasOption("o-sif-topo-file")) {
						file.toTopologicalFile(outputFileName+"_topo.txt", proteinNetwork);
					}
				} else {
					System.err.println("Missing custom interactome");
					return;					
				}
			} else {
				System.err.println("Unknown interactome");
				return;
			}

			StringBuilder sbMeans = createMeansHeader();
			StringBuilder sbTopo = createTopoHeader();
			StringBuilder sbComponents = createComponentsHeader();
			
			if(commandLine.hasOption("list1")) {
				logger.debug("Starting list1.........");
				String nodeFile = commandLine.getOptionValue("list1");
				int node=1;
				FileUtils.checkFile(new File(nodeFile));
				List<String> list = IOUtils.readLines(nodeFile);
				logger.debug("nodes read: " + list.toString());

				double tInicio = System.currentTimeMillis();
				SimpleUndirectedGraph<ProteinVertex, DefaultEdge> subgraph = (SimpleUndirectedGraph<ProteinVertex, DefaultEdge>) Subgraph.randomSubgraph(proteinNetwork.getInteractomeGraph(), toProteinVertex(list));
				double tFinal = System.currentTimeMillis();
				this.randomSize = subgraph.getVertices().size();
//				System.out.println("Tiempo creando subgrafo["+subgraph.getVertices().size()+","+subgraph.getEdges().size()+"] lista1: "+(tFinal-tInicio)/1000);
				if(commandLine.hasOption("intermediate")){
					tInicio = System.currentTimeMillis();
					intermediatesSub1 = Subgraph.OneIntermediateList(proteinNetwork.getInteractomeGraph(), subgraph);
					tFinal = System.currentTimeMillis();
//					System.out.println("Tiempo creando intermediario subgrafo["+subgraph.getVertices().size()+","+subgraph.getEdges().size()+"] lista1: "+(tFinal-tInicio)/1000);
					if(intermediatesSub1.size()>0)
						result.addOutputItem(new Item("external_nodes_list_"+node, intermediatesSub1.toString().substring(1, intermediatesSub1.toString().length()-1), "External nodes added", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Subnet results.List " + node));
					else
						result.addOutputItem(new Item("external_nodes_list_"+node, "No added", "External nodes added", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Subnet results.List " + node));

				}

				subProteinNetwork1 = createSubnet(subgraph);

				sbTopo.append(getTopologicalValues(subProteinNetwork1, node));
				f = new File(outputFileName+"_sn_nodeFile"+node+"_topo.txt");
				IOUtils.write(f.getAbsoluteFile(), sbTopo.toString());
				result.addOutputItem(new Item("sn_nodeFile"+node+"_topo_param", f.getName(), "Topological values", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Subnet results.List " + node));

				sbMeans.append(getTopologicalMeanValues(subProteinNetwork1, node));
				f = new File(outputFileName+"_sn_nodeFile"+node+"_means.txt");
				IOUtils.write(f.getAbsoluteFile(), sbMeans.toString());
				result.addOutputItem(new Item("sn_nodeFile"+node+"_means_param", f.getName(), "Means values", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Subnet results.List " + node));


				if(!noNumberOfComponents) {
					logger.debug("Starting list1 components.........");
					componentsListSub1 = subProteinNetwork1.getInteractomeGraph().getAllInformationComponents(true);
					
					sbComponents.append(getComponentsValues(subProteinNetwork1, node, componentsListSub1));
					sbComponents.deleteCharAt(sbComponents.lastIndexOf(System.getProperty("line.separator")));
					f = new File(outputFileName+"_sn_nodeFile"+node+"_comp.txt");
					IOUtils.write(f.getAbsoluteFile(), sbComponents.toString());
					result.addOutputItem(new Item("sn_nodeFile"+node+"_components_param", f.getName(), "Component values", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Subnet results.List " + node));

				}
			}

			if(commandLine.hasOption("list2") && !"none".equalsIgnoreCase(commandLine.getOptionValue("list2"))) {
				sbMeans = createMeansHeader();
				sbTopo = createTopoHeader();
				sbComponents = createComponentsHeader();
				logger.debug("Starting list2.........");
				String nodeFile = commandLine.getOptionValue("list2");
				int node=2;
				FileUtils.checkFile(new File(nodeFile));
				List<String> list = IOUtils.readLines(nodeFile);
				logger.debug("nodes read: " + list.toString());

				SimpleUndirectedGraph<ProteinVertex, DefaultEdge> subgraph = (SimpleUndirectedGraph<ProteinVertex, DefaultEdge>) Subgraph.randomSubgraph(proteinNetwork.getInteractomeGraph(), toProteinVertex(list));

				if(commandLine.hasOption("intermediate")){
					intermediatesSub2 = Subgraph.OneIntermediateList(proteinNetwork.getInteractomeGraph(), subgraph);
					if(intermediatesSub1.size()>0)
						result.addOutputItem(new Item("external_nodes_list_"+node, intermediatesSub2.toString().substring(1, intermediatesSub2.toString().length()-1), "External nodes added", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Subnet results.List " + node));
					else
						result.addOutputItem(new Item("external_nodes_list_"+node, "No added", "External nodes added", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Subnet results.List " + node));
				}
				
				/***
				 * if(intermediatesSub1.size()>0)
						result.addOutputItem(new Item("external_nodes_list_"+node, intermediatesSub1.toString().substring(1, intermediatesSub1.toString().length()-1), "External nodes added", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Subnet results.List " + node));
					else
						result.addOutputItem(new Item("external_nodes_list_"+node, "No added", "External nodes added", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Subnet results.List " + node));

				 */
				subProteinNetwork2 = createSubnet(subgraph);

				sbTopo.append(getTopologicalValues(subProteinNetwork2, node));
				f = new File(outputFileName+"_sn_nodeFile"+node+"_topo.txt");
				IOUtils.write(f.getAbsoluteFile(), sbTopo.toString());
				result.addOutputItem(new Item("sn_nodeFile"+node+"_topo_param", f.getName(), "Topografical values", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Subnet results.List " + node));

				sbMeans.append(getTopologicalMeanValues(subProteinNetwork2, node));
				f = new File(outputFileName+"_sn_nodeFile"+node+"_means.txt");
				IOUtils.write(f.getAbsoluteFile(), sbMeans.toString());
				result.addOutputItem(new Item("sn_nodeFile"+node+"_means_param", f.getName(), "Means values", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Subnet results.List " + node));

				if(!noNumberOfComponents) {
					logger.debug("Starting list2 components.........");
					componentsListSub2 = subProteinNetwork2.getInteractomeGraph().getAllInformationComponents(true);
					sbComponents.append(getComponentsValues(subProteinNetwork2, node, componentsListSub2));
					sbComponents.deleteCharAt(sbComponents.lastIndexOf(System.getProperty("line.separator")));
					f = new File(outputFileName+"_sn_nodeFile"+node+"_comp.txt");
					IOUtils.write(f.getAbsoluteFile(), sbComponents.toString());
					result.addOutputItem(new Item("sn_nodeFile"+node+"_components_param", f.getName(), "Component values", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Subnet results.List " + node));
				}
			}
			
			if(commandLine.hasOption("randoms")){
				logger.debug("Starting randoms.........");
				int randoms = Integer.parseInt(commandLine.getOptionValue("randoms"));
//				if(!commandLine.hasOption("randoms-size")){
//					int size = this.proteinNetwork.getInteractomeGraph().getVertices().size();
//					while(randomSize == 0)
//						randomSize = (int) (Math.random() * size);
//				}
//				else{
//					randomSize = Integer.parseInt(commandLine.getOptionValue("randoms-size"));
//				}
				if(randomSize > 0)
					createRandoms(randoms, randomSize);
				result.addOutputItem(new Item("randoms_param", ""+randoms, "Number of randoms", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
				result.addOutputItem(new Item("randoms_size_param", ""+randomSize, "Size of randoms", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
			}

			result.addOutputItem(new Item("no_number_components_param", (commandLine.hasOption("no-number-components") ? "yes" : "no"), "Calculate the number of components", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
			result.addOutputItem(new Item("bicomponents_param", (commandLine.hasOption("bicomponents") ? "yes" : "no"), "Calculate the number of bicomponents", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
			result.addOutputItem(new Item("intermediate_param", (commandLine.hasOption("intermediate") ? "yes" : "no"), "Create the network with 1 intermediate", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));

			if(commandLine.hasOption("randoms") && (!commandLine.hasOption("list2") || "none".equalsIgnoreCase(commandLine.getOptionValue("list2")))) {
				statsOneListAnalysis(side);
			} else if(!commandLine.hasOption("randoms") && commandLine.hasOption("list1") && commandLine.hasOption("list2") && !"none".equalsIgnoreCase(commandLine.getOptionValue("list2")))
				statsTwoListsAnalisys(side);
			else{
				logger.error("Not correct arguments for statistic test");
				return;
			}
			if(json){
				if(subProteinNetwork1 != null){
					createJson(subProteinNetwork1, outputFileName+"_list1.dot", componentsListSub1, intermediatesSub1, 1);
				}
				if(subProteinNetwork2 != null){
					createJson(subProteinNetwork2, outputFileName+"_list2.dot", componentsListSub2, intermediatesSub2, 2);
				}
			}
			if(xml){
				Xml xmlObject = new Xml();
				if(subProteinNetwork1 != null){
					xmlObject.graphToXML(outdir+"/subnetwork1.xml",subProteinNetwork1.getInteractomeGraph(), intermediatesSub1, componentsListSub1);
				}
				if(subProteinNetwork2 != null){
					xmlObject.graphToXML(outdir+"/subnetwork2.xml",subProteinNetwork2.getInteractomeGraph(), intermediatesSub2, componentsListSub2);
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private void createJson(ProteinNetwork proteinNetwork, String sourceDotFile, List<List<ProteinVertex>> componentsListSub, Set<String> intermediatesSub, int node) throws IOException{

		Dot<ProteinVertex, DefaultEdge> dot = new Dot<ProteinVertex, DefaultEdge>();
		Svg svg = new Svg();
		Json<ProteinVertex, DefaultEdge> json = new Json<ProteinVertex, DefaultEdge>();
		List<File> fileList = new ArrayList<File>();
		List<String> layoutsName = new ArrayList<String>();
		
		if(componentsListSub == null)
			logger.error("not components calculated. Remove --no-number-components");

		IOUtils.write(sourceDotFile, dot.toDot(proteinNetwork.getInteractomeGraph()));
		
		File dotFile = new File(outputFileName+"_list"+node+"_dot.svg");
		IOUtils.write(dotFile, svg.toSvg(sourceDotFile, "dot"));
		fileList.add(dotFile);
		layoutsName.add("dot");
		
		File twopiFile = new File(outputFileName+"_list"+node+"_twopi.svg");
		IOUtils.write(twopiFile, svg.toSvg(sourceDotFile, "twopi"));
		fileList.add(twopiFile);
		layoutsName.add("twopi");
		
		IOUtils.write(outputFileName+"_list"+node+".json", json.toJson(fileList, layoutsName, proteinNetwork.getInteractomeGraph(), intermediatesSub, componentsListSub));

		File f = new File(outputFileName+node+".json");
		System.out.println("********************* json file = " + f.getAbsolutePath());
		if (f.exists()) {
			System.out.println("********************* exists JSON !!!");
			String url = "SnowViewer2?filename=" + f.getName();
			result.addOutputItem(new Item("viewer_param", url, "Network viewer", TYPE.HTML, StringUtils.toList("SERVER,INCLUDE_REFS", ","), new HashMap<String, String>(2), "Graph"));			
		}
//		f = new File("result1.json");
//		//String url = "http://beta.babelomics.bioinfo.cipf.es/SnowViewer?filename=" + fileName;
//		//http://beta.babelomics.bioinfo.cipf.es/SnowViewer?filename=subnetwork1.xml&jobid=262&sessionid=7tySrwJgNRM0tyGOnj9N3eT2BCTsEC8Qw3StPMwwz1WfsFyt1zoIqR6lPa7fnv54
//		result.addOutputItem(new Item("graph_param", url, "Graph", TYPE.HTML, StringUtils.toList("SERVER,INCLUDE_REFS", ","), new HashMap<String, String>(2), "Graph"));
//			
//			url = "SnowViewer?filename=" + fileName;
//			result.addOutputItem(new Item(name + "_new_window", url, "Open applet in a new window", TYPE.LINK, StringUtils.toList("SERVER,INCLUDE_REFS", ","), new HashMap<String, String>(2), groupName));
		
	}

	private void createImages(String fileName, List<Double> list1, String legend1, List<Double> list2, String legend2, String itemName, String itemLabel, String itemGroup) throws IOException {
		File f = new File(fileName);
		BoxPlotChart bpc = new BoxPlotChart("", "", "");
		bpc.getLegend().setVisible(true);
		if(!list2.isEmpty())
			bpc.addSeries(list2, legend2, legend2);
		if(!list1.isEmpty())
			bpc.addSeries(list1, legend1, legend1);
		bpc.setOrientation(PlotOrientation.HORIZONTAL);
		bpc.save(f.getAbsolutePath()+".png", 300, 250, "png");

		result.addOutputItem(new Item(itemName, f.getName()+".png", itemLabel, TYPE.IMAGE, new ArrayList<String>(), new HashMap<String, String>(2), itemGroup));		
	}


	private void statsOneListAnalysis(String side) throws IOException{
		// 1st Analysis
		File f = null;
		logger.debug("Starting 1st Analysis..................");
		WilcoxonTest wTest = new WilcoxonTest();
	//	List<String> list = IOUtils.readLines(commandLine.getOptionValue("list1"));
		List<ProteinVertex> list = this.subProteinNetwork1.getInteractomeGraph().getVertices();


		List<Double> relBetInter = proteinNetwork.getBetRelList();
		List<Double> relBetList1 = new ArrayList<Double>();

		List<Double> connInter = proteinNetwork.getConnList();
		List<Double> connList1 = new ArrayList<Double>();

		List<Double> clustInter = proteinNetwork.getClustList();
		List<Double> clustList1 = new ArrayList<Double>();

		createTopoFilterList(list, relBetList1, connList1, clustList1);

		String toWrite = ksTest(relBetList1, relBetInter, connList1, connInter, clustList1, clustInter, side);
		if(!toWrite.equals("")){
			f = new File(outputFileName+"_list_inter_kol.txt");
			IOUtils.write(f.getAbsolutePath(), toWrite);
			result.addOutputItem(new Item("list_inter_kol_param", f.getName(), "List - Interactome", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Statistic results.Kolmogorov-Smirnov test"));
		}
		else
			result.addOutputItem(new Item("list_inter_kol_param", "Empty results", "List - Interactome", Item.TYPE.MESSAGE, new ArrayList<String>(), new HashMap<String,String>(), "Statistic results.Kolmogorov-Smirnov test"));

		
		f = new File(outputFileName+"_list_inter_wil.txt");
		
		if(wTest.test(wBinPath, f.getAbsolutePath(), side, "list1", "inter", relBetList1, relBetInter, connList1, connInter, clustList1, clustInter)){
			result.addOutputItem(new Item("list_inter_wil_param", f.getName(), "List - Interactome", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Statistic results.Wilcoxon test"));
//			List<Double> wilcoxValues = getWilcoxonTestValues(f);
//			result.addOutputItem(new Item("bet_list_inter_wil_param", wilcoxValues.get(0).toString(), "List - Interactome: Betweenness", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Statistic results.Wilcoxon test"));
//			result.addOutputItem(new Item("conn_list_inter_wil_param", wilcoxValues.get(1).toString(), "List - Interactome: Connections", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Statistic results.Wilcoxon test"));
//			result.addOutputItem(new Item("coeff_list_inter_wil_param", wilcoxValues.get(2).toString(), "List - Interactome: Clustering coefficient", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Statistic results.Wilcoxon test"));

		}
		else
			result.addOutputItem(new Item("list_inter_wil_param", "Empty results", "List - Interactome", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Statistic results.Wilcoxon test"));
		
		if(images){
			createImages(outputFileName+"_list_inter_relBet", relBetList1, "list1", relBetInter, "inter", "list_inter_relBet", "relBet", "Images.List vs inter");
			createImages(outputFileName+"_list_inter_conn", connList1, "list1", connInter, "inter", "list_inter_conn", "conn", "Images.List vs inter");
			createImages(outputFileName+"_list_inter_clust", clustList1, "list1", clustInter, "inter", "list_inter_clust", "clust", "Images.List vs inter");
		}
		logger.debug("Finished 1st Analysis..................");
		//2nd Analysis
		logger.debug("Starting 2nd Analysis..................");
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

		toWrite = ksTest(relBetSubnet1, relBetRandoms, connSubnet1, connRandoms, clustSubnet1, clustRandoms, side);
		if(!toWrite.equals("")){
			f = new File(outputFileName+"_sn_random_kol.txt");
			IOUtils.write(f.getAbsolutePath(), toWrite);
			result.addOutputItem(new Item("sn_random_kol_param", f.getName(), "Subnet - Random", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Statistic results.Kolmogorov-Smirnov test"));
		}
		else
			result.addOutputItem(new Item("sn_random_kol_param", "Empty results", "Subnet - Random", Item.TYPE.MESSAGE, new ArrayList<String>(), new HashMap<String,String>(), "Statistic results.Kolmogorov-Smirnov test"));

		
		f = new File(outputFileName+"_sn_random_wil.txt");
		if(wTest.test(wBinPath, f.getAbsolutePath(), side, "sn1", "rnd", relBetSubnet1, relBetRandoms, connSubnet1, connRandoms, clustSubnet1, clustRandoms)){
			result.addOutputItem(new Item("sn_random_wil_param", f.getName(), "Subnet - Random", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Statistic results.Wilcoxon test"));
//			List<Double> wilcoxValues = getWilcoxonTestValues(f);
//			result.addOutputItem(new Item("bet_sn_random_wil_param", wilcoxValues.get(0).toString(), "Subnet - Random: Betweenness", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Statistic results.Wilcoxon test"));
//			result.addOutputItem(new Item("conn_sn_random_wil_param", wilcoxValues.get(1).toString(), "Subnet - Random: Connections", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Statistic results.Wilcoxon test"));
//			result.addOutputItem(new Item("coef_sn_random_wil_param", wilcoxValues.get(2).toString(), "Subnet - Random: Clustering coefficient", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Statistic results.Wilcoxon test"));

		}
		else
			result.addOutputItem(new Item("sn_random_wil_param", "Empty results", "Subnet - Random", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Statistic results.Wilcoxon test"));

		if(images){
			createImages(outputFileName+"_sn_random_relBet", relBetSubnet1, "subnet1", relBetRandoms, "randoms", "sn_random_relBet", "relBet", "Images.Subnet random");
			createImages(outputFileName+"_sn_random_conn", connSubnet1, "subnet1", connRandoms, "randoms", "sn_random_conn", "conn", "Images.Subnet random");
			createImages(outputFileName+"_sn_random_clust", clustSubnet1, "subnet1", clustRandoms, "randoms", "sn_random_clust", "clust", "Images.Subnet random");
		}
		logger.debug("Finished 2nd Analysis..................");
	}
	private void statsTwoListsAnalisys(String side) throws IOException{
		// 1st Analysis
		logger.debug("Starting 1st Analysis..................");
		File f = null;
		WilcoxonTest wTest = new WilcoxonTest();

		List<Double> relBetList1 = new ArrayList<Double>();
		List<Double> connList1 = new ArrayList<Double>();
		List<Double> clustList1 = new ArrayList<Double>();
		List<Double> relBetList2 = new ArrayList<Double>();
		List<Double> connList2 = new ArrayList<Double>();
		List<Double> clustList2 = new ArrayList<Double>();

		List<ProteinVertex> list1 = this.subProteinNetwork1.getInteractomeGraph().getVertices();
		List<ProteinVertex> list2 = this.subProteinNetwork2.getInteractomeGraph().getVertices();


		createTopoFilterList(list1, relBetList1, connList1, clustList1);
		createTopoFilterList(list2, relBetList2, connList2, clustList2);

		String toWrite = ksTest(relBetList1, relBetList2, connList1, connList2, clustList1, clustList2, side);
		if(!toWrite.equals("")){
			f = new File(outputFileName+"_list1_list2_kol.txt");
			IOUtils.write(f.getAbsolutePath(), toWrite);
			result.addOutputItem(new Item("list1_list2_kol", f.getName(), "List1 - List2", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Statistic results.Kolmogorov-Smirnov test"));
		}
		else
			result.addOutputItem(new Item("list1_list2_kol", "Empty results", "List1 - List2", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Statistic results.Kolmogorov-Smirnov test"));

		
		f = new File(outputFileName+"_list1_list2_wil.txt");
		if(wTest.test(wBinPath, f.getAbsolutePath(), side, "list1", "list2", relBetList1, relBetList2, connList1, connList2, clustList1, clustList2)){
			result.addOutputItem(new Item("list1_list2_wil", f.getName(), "List1 - List2", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Statistic results.Wilcoxon test"));
//			List<Double> wilcoxValues = getWilcoxonTestValues(f);
//			result.addOutputItem(new Item("bet_list1_list2_wil", wilcoxValues.get(0).toString(), "Subnet - Random: Betweenness", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Statistic results.Wilcoxon test"));
//			result.addOutputItem(new Item("conn_list1_list2_wil", wilcoxValues.get(1).toString(), "Subnet - Random: Connections", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Statistic results.Wilcoxon test"));
//			result.addOutputItem(new Item("coef_list1_list2_wil", wilcoxValues.get(2).toString(), "Subnet - Random: Clustering coefficient", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Statistic results.Wilcoxon test"));
		}
		else
			result.addOutputItem(new Item("list1_list2_wil", "Empty results", "List1 - List2", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Statistic results.Wilcoxon test"));

		if(images){
			createImages(outputFileName+"_list1_list2_relBet", relBetList1, "list1", relBetList2, "list2", "list1_list2_relBet", "relBet", "Images.List1 vs List2");
			createImages(outputFileName+"_list1_list2_conn", connList1, "list1", connList2, "list2", "list1_list2_conn", "conn", "Images.List1 vs List2");
			createImages(outputFileName+"_list1_list2_clust", clustList1, "list1", clustList2, "list2", "list1_list2_clust", "clust", "Images.List1 vs List2");
		}

		//2nd Analysis
		logger.debug("Starting 2nd Analysis..................");
		List<Double> relBetSubnet1 = subProteinNetwork1.getBetRelList();
		List<Double> relBetSubnet2 = subProteinNetwork2.getBetRelList();

		List<Double> connSubnet1 = subProteinNetwork1.getConnList();
		List<Double> connSubnet2 = subProteinNetwork2.getConnList();

		List<Double> clustSubnet1 = subProteinNetwork1.getClustList();
		List<Double> clustSubnet2 = subProteinNetwork2.getClustList();
		
		toWrite = ksTest(relBetSubnet1, relBetSubnet2, connSubnet1, connSubnet2, clustSubnet1, clustSubnet2, side);
		if(!toWrite.equals("")){
			f = new File(outputFileName+"_sn1_sn2_kol.txt");
			IOUtils.write(f.getAbsolutePath(), toWrite);		
			result.addOutputItem(new Item("sn1_sn2_kol", f.getName(), "Subnet1 - Subnet2", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Statistic results.Kolmogorov-Smirnov test"));
		}
		else
			result.addOutputItem(new Item("sn1_sn2_kol", "Empty results", "Subnet1 - Subnet2", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Statistic results.Kolmogorov-Smirnov test"));

		
		f = new File(outputFileName+"_sn1_sn2_wil.txt");
		if(wTest.test(wBinPath, f.getAbsolutePath(), side, "sn1", "sn2", relBetSubnet1, relBetSubnet2, connSubnet1, connSubnet2, clustSubnet1, clustSubnet2)){
			result.addOutputItem(new Item("sn1_sn2_wil", f.getName(), "Subnet1 - Subnet2", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Statistic results.Wilcoxon test"));
//			List<Double> wilcoxValues = getWilcoxonTestValues(f);
//			result.addOutputItem(new Item("bet_sn1_sn2_wil", wilcoxValues.get(0).toString(), "Subnet1 - Subnet2: Betweenness", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Statistic results.Wilcoxon test"));
//			result.addOutputItem(new Item("conn_sn1_sn2_wil", wilcoxValues.get(1).toString(), "Subnet1 - Subnet2: Connections", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Statistic results.Wilcoxon test"));
//			result.addOutputItem(new Item("coef_sn1_sn2_wil", wilcoxValues.get(2).toString(), "Subnet1 - Subnet2: Clustering coefficient", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Statistic results.Wilcoxon test"));
		
		}
		else
			result.addOutputItem(new Item("sn1_sn2_wil", "Empty results", "Subnet1 - Subnet2", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Statistic results.Wilcoxon test"));
		
		if(images){
			createImages(outputFileName+"_sn1_sn2_relBet", relBetSubnet1, "subnet1", relBetSubnet2, "subnet2", "sn1_sn2_relBet", "relBet", "Images.Subnet1 vs Subnet2");
			createImages(outputFileName+"_sn1_sn2_conn", connSubnet1, "subnet1", connSubnet2, "subnet2", "sn1_sn2_conn", "conn", "Images.Subnet1 vs Subnet2");
			createImages(outputFileName+"_sn1_sn2_clust", clustSubnet1, "subnet1", clustSubnet2, "subnet2", "sn1_sn2_clust", "clust", "Images.Subnet1 vs Subnet2");
		}

	}
	private String ksTest(List<Double> relBetList1, List<Double> relBetList2, List<Double> connList1, List<Double> connList2, List<Double> clustList1, List<Double> clustList2, String side) throws IOException{
		if( relBetList1.size() == 0 || relBetList2.size() == 0 || connList1.size() == 0 || connList2.size() == 0 || clustList1.size() == 0 || clustList2.size() == 0)
			return "";
		KSTest kstest = new KSTest();
		StringBuilder sb = new StringBuilder();
		sb.append("#parameter\tpval\tside\n");
		sb.append("betweenness\t"+kstest.resultKolmogorovSmirnov(ListUtils.toArray(relBetList1), ListUtils.toArray(relBetList2), side).getPValue()+"\t"+side+"\n");
		sb.append("connections\t"+kstest.resultKolmogorovSmirnov(ListUtils.toArray(connList1), ListUtils.toArray(connList2), side).getPValue()+"\t"+side+"\n");
		sb.append("coefficient\t"+kstest.resultKolmogorovSmirnov(ListUtils.toArray(clustList1), ListUtils.toArray(clustList2), side).getPValue()+"\t"+side);
		return sb.toString();
	}
	private void createTopoFilterList(List<ProteinVertex> list, List<Double> relBetList, List<Double> connList, List<Double> clustList){
		for(ProteinVertex protein : list){
			if(proteinNetwork.getInteractomeGraph().getVertex(protein.getId()) == null)
				continue;
			relBetList.add(proteinNetwork.getRelBetweennessVertex(protein.getId()));
			connList.add(proteinNetwork.getDegreeVertex(protein.getId()));
			clustList.add(proteinNetwork.getClusterinVertex(protein.getId()));
		}
	}
	private void createRandoms(int randoms, int randomSize) throws IOException{
		
		StringBuilder sbMeans = createMeansHeader();
		StringBuilder sbTopo = createTopoHeader();
		StringBuilder sbComponents =createComponentsHeader();

		for(int i=1; i<=randoms; i++){
			SimpleUndirectedGraph<ProteinVertex, DefaultEdge> subgraph = (SimpleUndirectedGraph<ProteinVertex, DefaultEdge>) Subgraph.randomSubgraph(proteinNetwork.getInteractomeGraph(), randomSize);
			logger.debug("Randoms["+i+"]: V = "+subgraph.getVertices().size()+" E = "+subgraph.getEdges().size());
//			System.out.println("Randoms["+i+"]: V = "+subgraph.getVertices().size()+" E = "+subgraph.getEdges().size());
			if(commandLine.hasOption("intermediate")) {
				double tInicio = System.currentTimeMillis();
				Subgraph.OneIntermediateList(proteinNetwork.getInteractomeGraph(), subgraph);
				double tFinal = System.currentTimeMillis();
//				System.out.println("Randoms["+i+"]: V = "+subgraph.getVertices().size()+" E = "+subgraph.getEdges().size()+" t= "+(tFinal-tInicio)/1000);
				logger.debug("Randoms intermediate["+i+"]: V = "+subgraph.getVertices().size()+" E = "+subgraph.getEdges().size());
			}
			ProteinNetwork subProteinNetwork = createSubnet(subgraph);
			logger.debug("Subnet created");
			sbTopo.append(getTopologicalValues(subProteinNetwork, i)).append(System.getProperty("line.separator"));
			sbMeans.append(getTopologicalMeanValues(subProteinNetwork, i)).append(System.getProperty("line.separator"));
			if(!noNumberOfComponents) {
				List<List<ProteinVertex>> componentsList =  subProteinNetwork.getInteractomeGraph().getAllInformationComponents(true);
				sbComponents.append(getComponentsValues(subProteinNetwork, i, componentsList));
			}
			logger.debug("Components created");
			subProteinNetworkRandoms.add(subProteinNetwork);
		}
		
		File f = new File(outputFileName+"_sn_1-"+(randoms)+"_topo.txt");
		randomFilesToString(sbTopo, f.getAbsolutePath());
		result.addOutputItem(new Item("randoms_topo_param", f.getName(), "Topografical values", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Randoms results"));

		f = new File(outputFileName+"_sn_1-"+(randoms)+"_means.txt");
		randomFilesToString(sbMeans, f.getAbsolutePath());
		result.addOutputItem(new Item("randoms_means_param", f.getName(), "Means values", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Randoms results"));

		f = new File(outputFileName+"_sn_1-"+(randoms)+"_comp.txt");
		randomFilesToString(sbComponents, f.getAbsolutePath());
		result.addOutputItem(new Item("randoms_means_param", f.getName(), "Components values", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Randoms results"));
	}
	private ProteinNetwork createSubnet(SimpleUndirectedGraph<ProteinVertex, DefaultEdge> subgraph) {
		ProteinNetwork subProteinNetwork = new ProteinNetwork(subgraph);
		double tInicio = System.currentTimeMillis();
		subProteinNetwork.calcTopologicalValues();
		double tFinal = System.currentTimeMillis();
//		System.out.println("\tTiempo calculando valores topológicos: "+(tFinal-tInicio)/1000);
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
		sb.append("#Subnet\tId\tBet\tConn\tClust").append(System.getProperty("line.separator"));
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
		if(!sb.toString().equals(""))
			sb.deleteCharAt(sb.lastIndexOf(System.getProperty("line.separator")));
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
	private String getComponentsValues(ProteinNetwork subProteinNetwork, int subnet, List<List<ProteinVertex>> componentsList){
		StringBuilder sb = new StringBuilder();
		SimpleUndirectedGraph<ProteinVertex, DefaultEdge> interactomeGraph = subProteinNetwork.getInteractomeGraph();

		if(!noNumberOfComponents){
			//			componentsList =  interactomeGraph.getAllInformationComponents(true);
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
			if(!proteinName.equals(""))
				verticesList.add(new ProteinVertex(proteinName));
		}
		return verticesList;
	}
	private List<Double> getWilcoxonTestValues(File file) throws IOException{
		List<Double> list = new ArrayList<Double>();
		
		TextFileReader tfr = new TextFileReader(file.getAbsolutePath());
		String line = null;
		String[] fields;
		while((line = tfr.readLine()) != null) {
			if(line.startsWith("#"))
				continue;
			fields = line.split("\t");
			if(fields.length == 3) {
				list.add(Double.parseDouble(fields[1]));
			}
		}
		tfr.close();
		return list;
	}
}
