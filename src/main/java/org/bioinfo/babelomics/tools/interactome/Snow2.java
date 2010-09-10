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
import org.bioinfo.networks.protein.files.Dot;
import org.bioinfo.networks.protein.files.Json;
import org.bioinfo.networks.protein.files.Svg;
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
	private Set<String> intermediatesSub1, intermediatesSub2; 
	private List<List<ProteinVertex>> componentsListSub1, componentsListSub2;

	private String wBinPath;

	@Override
	public void initOptions() {
		options.addOption(OptionFactory.createOption("interactome", "Interactome: hsa, sce or own (for you own interactome)", true));
		options.addOption(OptionFactory.createOption("sif-file", "An input file containing a SIF interactome", false));
		options.addOption(OptionFactory.createOption("topo-file","t", "An input file containing topological values", false, true));
		options.addOption(OptionFactory.createOption("list1", "An input file containing a node per line", false));
		options.addOption(OptionFactory.createOption("list2", "An input file containing a node per line", false));
		options.addOption(OptionFactory.createOption("type", "An argument saying if you want genes or proteins", false, true));
		options.addOption(OptionFactory.createOption("randoms", "Number of randoms", false, true));
		options.addOption(OptionFactory.createOption("randoms-size", "Size of randoms", false, true));
		options.addOption(OptionFactory.createOption("intermediate", "If there is this argument, it will create the network with 1 intermediate", false, false));
		options.addOption(OptionFactory.createOption("no-number-components", "If there is this argument, it won't calculate the number of components", false, false));
		options.addOption(OptionFactory.createOption("bicomponents", "If there is this argument, it will calculate the number of bicomponents", false, false));
		options.addOption(OptionFactory.createOption("json", "It will create an output .json file", false, false));
		options.addOption(OptionFactory.createOption("o-sif-topo-file", "Create a full topological file from a sif file", false));
		options.addOption(OptionFactory.createOption("o-name", "If there is this argument, it will create an output .cmp file for the information of each component", false, true));
		options.addOption(OptionFactory.createOption("side", "side for kolmogorov and wilkoxon test. Can be two.sided(by default), less or greater", false, true));
		options.addOption(OptionFactory.createOption("images", "Print the images for the statistics", false, false));
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
			//String folderInteractions = "/mnt/commons/babelomics/tests/snow2/data/interactions/";
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
					proteinNetwork.calcTopologicalValues();
				} else {
					System.err.println("Missing custom interactome");
					return;					
				}
			} else {
				System.err.println("Unknown interactome");
				return;
			}


			if(commandLine.hasOption("list1")) {
				System.out.println("Starting list1.........");
				String nodeFile = commandLine.getOptionValue("list1");
				//				subProteinNetwork1 = createNodeFile(commandLine.getOptionValue("node-file1"), 1, intermediatesSub1);
				int node=1;
				FileUtils.checkFile(new File(nodeFile));
				List<String> list = IOUtils.readLines(nodeFile);
				logger.debug("nodes read: " + list.toString());

				StringBuilder sbMeans = createMeansHeader();
				StringBuilder sbTopo = createTopoHeader();

				SimpleUndirectedGraph<ProteinVertex, DefaultEdge> subgraph = (SimpleUndirectedGraph<ProteinVertex, DefaultEdge>) Subgraph.randomSubgraph(proteinNetwork.getInteractomeGraph(), toProteinVertex(list));

				if(commandLine.hasOption("intermediate")){
					System.out.println("Starting intermediate............");
					System.out.println("Original: V = "+subgraph.getVertices().size()+", E = "+subgraph.getEdges().size());
					intermediatesSub1 = Subgraph.OneIntermediateList(proteinNetwork.getInteractomeGraph(), subgraph);
					System.out.println("Intermediate: V = "+subgraph.getVertices().size()+", E = "+subgraph.getEdges().size());
				}

				subProteinNetwork1 = createSubnet(subgraph);

				sbTopo.append(getTopologicalValues(subProteinNetwork1, node));
				f = new File(outputFileName+"_sn_nodeFile"+node+"_topo.txt");
				IOUtils.write(f.getAbsoluteFile(), sbTopo.toString());
				result.addOutputItem(new Item("sn_nodeFile"+node+"_topo_param", f.getName(), "Topografical values", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Subnet results.List " + node));



				sbMeans.append(getTopologicalMeanValues(subProteinNetwork1, node));
				f = new File(outputFileName+"_sn_nodeFile"+node+"_means.txt");
				IOUtils.write(f.getAbsoluteFile(), sbMeans.toString());
				result.addOutputItem(new Item("sn_nodeFile"+node+"_means_param", f.getName(), "Means values", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Subnet results.List " + node));


				if(!noNumberOfComponents) {
					System.out.println("Starting list1 components.........");
					componentsListSub1 = subProteinNetwork1.getInteractomeGraph().getAllInformationComponents(true);
					StringBuilder sbComponents = createComponentsHeader();
					sbComponents.append(getComponentsValues(subProteinNetwork1, node, componentsListSub1));
					sbComponents.deleteCharAt(sbComponents.lastIndexOf(System.getProperty("line.separator")));
					f = new File(outputFileName+"_sn_nodeFile"+node+"_comp.txt");
					IOUtils.write(f.getAbsoluteFile(), sbComponents.toString());
					result.addOutputItem(new Item("sn_nodeFile"+node+"_components_param", f.getName(), "Component values", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Subnet results.List " + node));

				}
			}

			if(commandLine.hasOption("list2") && !"none".equalsIgnoreCase(commandLine.getOptionValue("list2"))) {
				System.out.println("Starting list2.........");
				String nodeFile = commandLine.getOptionValue("list2");
				int node=2;
				FileUtils.checkFile(new File(nodeFile));
				List<String> list = IOUtils.readLines(nodeFile);
				logger.debug("nodes read: " + list.toString());

				StringBuilder sbMeans = createMeansHeader();
				StringBuilder sbTopo = createTopoHeader();

				SimpleUndirectedGraph<ProteinVertex, DefaultEdge> subgraph = (SimpleUndirectedGraph<ProteinVertex, DefaultEdge>) Subgraph.randomSubgraph(proteinNetwork.getInteractomeGraph(), toProteinVertex(list));

				if(commandLine.hasOption("intermediate"))
					intermediatesSub2 = Subgraph.OneIntermediateList(proteinNetwork.getInteractomeGraph(), subgraph);

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
					System.out.println("Starting list2 components.........");
					componentsListSub2 = subProteinNetwork2.getInteractomeGraph().getAllInformationComponents(true);
					StringBuilder sbComponents = createComponentsHeader();
					sbComponents.append(getComponentsValues(subProteinNetwork2, node, componentsListSub2));
					sbComponents.deleteCharAt(sbComponents.lastIndexOf(System.getProperty("line.separator")));
					f = new File(outputFileName+"_sn_nodeFile"+node+"_comp.txt");
					IOUtils.write(f.getAbsoluteFile(), sbComponents.toString());
					result.addOutputItem(new Item("sn_nodeFile"+node+"_components_param", f.getName(), "Component values", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Subnet results.List " + node));
				}
			}
			if(commandLine.hasOption("o-sif-topo-file")) {
				proteinNetwork.calcTopologicalValues();
				file.toTopologicalFile(commandLine.getOptionValue("o-sif-topo-file")+"_topo.txt", proteinNetwork);
			}
			if(commandLine.hasOption("randoms")){
				System.out.println("Starting randoms.........");
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
				System.err.println("Not correct arguments for statistic test");
				return;
			}
				


			if(json){
				if(subProteinNetwork1 != null){
					createJson(subProteinNetwork1, outputFileName+"1.dot", componentsListSub1, intermediatesSub1, 1);
					//					createSVGFile(subProteinNetwork1, outputFileName+"1.dot", 1);
					//					createJson(outputFileName+"1.svg", subProteinNetwork1, componentsListSub1, intermediatesSub1, 1);
				}
				if(subProteinNetwork2 != null){
					createJson(subProteinNetwork2, outputFileName+"2.dot", componentsListSub2, intermediatesSub2, 2);
				}
			}
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private void createJson(ProteinNetwork proteinNetwork, String sourceDotFile, List<List<ProteinVertex>> componentsListSub, Set<String> intermediatesSub, int node) throws IOException{
		//		double tInicio = System.currentTimeMillis();
		//		System.out.println("Starting dot file....");

		Dot<ProteinVertex, DefaultEdge> dot = new Dot<ProteinVertex, DefaultEdge>();
		Svg svg = new Svg();
		Json<ProteinVertex, DefaultEdge> json = new Json<ProteinVertex, DefaultEdge>();

		if(componentsListSub == null)
			System.err.println("not components calculated. Remove --no-number-components");

		IOUtils.write(sourceDotFile, dot.toDot(proteinNetwork.getInteractomeGraph()));
		//		System.out.println("Starting svg files....");
		//		svg.toSvg(sourceDotFile, layout);
		File dotFile = new File(outputFileName+node+"_dot.svg");
		IOUtils.write(dotFile, svg.toSvg(sourceDotFile, "dot"));
		//IOUtils.write(outputFileName+node+"_dot.json", json.toJson(dotFile, proteinNetwork.getInteractomeGraph(), intermediatesSub, componentsListSub));

		File twopiFile = new File(outputFileName+node+"_twopi.svg");
		IOUtils.write(twopiFile, svg.toSvg(sourceDotFile, "twopi"));
		//IOUtils.write(outputFileName+node+"_twopi.json", json.toJson(twopiFile, proteinNetwork.getInteractomeGraph(), intermediatesSub, componentsListSub));
		IOUtils.write(outputFileName+node+"_all.json", json.toJson(dotFile, twopiFile, proteinNetwork.getInteractomeGraph(), intermediatesSub, componentsListSub));

	}

	private void createImages(String fileName, List<Double> list1, String legend1, List<Double> list2, String legend2, String itemName, String itemLabel, String itemGroup) throws IOException {
		File f = new File(fileName);
		BoxPlotChart bpc = new BoxPlotChart("", "", "");
		bpc.getLegend().setVisible(true);
		bpc.addSeries(list2, legend2, legend2);
		bpc.addSeries(list1, legend1, legend1);
		bpc.setOrientation(PlotOrientation.HORIZONTAL);
		bpc.save(f.getAbsolutePath());

		result.addOutputItem(new Item(itemName, f.getName(), itemLabel, TYPE.IMAGE, new ArrayList<String>(), new HashMap<String, String>(2), itemGroup));		
	}


	private void statsOneListAnalysis(String side) throws IOException{
		// 1st Analysis
		File f = null;
		System.out.println("Starting 1st Analysis..................");
		WilcoxonTest wTest = new WilcoxonTest();
		List<String> list = IOUtils.readLines(commandLine.getOptionValue("list1"));


		List<Double> relBetInter = proteinNetwork.getBetRelList();
		List<Double> relBetList1 = new ArrayList<Double>();

		List<Double> connInter = proteinNetwork.getConnList();
		List<Double> connList1 = new ArrayList<Double>();

		List<Double> clustInter = proteinNetwork.getClustList();
		List<Double> clustList1 = new ArrayList<Double>();

		createTopoFilterList(list, relBetList1, connList1, clustList1);

		String toWrite = formatStatsFile(relBetList1, relBetInter, connList1, connInter, clustList1, clustInter, side);

		f = new File(outputFileName+"_list_inter_kol.txt");
		IOUtils.write(f.getAbsoluteFile(), toWrite);
		result.addOutputItem(new Item("list_inter_kol_param", f.getName(), "List-inter", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Statistic results.Kolmogorov-Smirnov test"));

		wTest.test(wBinPath, outputFileName+"_list_inter_wil.txt", side, "list1", "inter", relBetList1, relBetInter, connList1, connInter, clustList1, clustInter);
		f = new File(outputFileName+"_list_inter_wil.txt");
		result.addOutputItem(new Item("list_inter_wil_param", f.getName(), "List-inter", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Statistic results.Wilcoxon test"));

		if(images){
			createImages(outputFileName+"_list_inter_relBet", relBetList1, "list1", relBetInter, "inter", "list_inter_relBet", "relBet", "Images.List vs inter");
			createImages(outputFileName+"_list_inter_conn", connList1, "list1", connInter, "inter", "list_inter_conn", "conn", "Images.List vs inter");
			createImages(outputFileName+"_list_inter_clust", clustList1, "list1", clustInter, "inter", "list_inter_clust", "clust", "Images.List vs inter");
		}
		System.out.println("Finished 1st Analysis..................");
		//2nd Analysis
		System.out.println("Starting 2nd Analysis..................");
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
		f = new File(outputFileName+"_sn_random_kol.txt");
		IOUtils.write(f.getAbsoluteFile(), toWrite);
		result.addOutputItem(new Item("sn_random_kol_param", f.getName(), "Subnet random", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Statistic results.Kolmogorov-Smirnov test"));

		f = new File(outputFileName+"_sn_random_wil.txt");
		wTest.test(wBinPath, f.getAbsolutePath(), side, "sn1", "rnd", relBetSubnet1, relBetRandoms, connSubnet1, connRandoms, clustSubnet1, clustRandoms);
		result.addOutputItem(new Item("sn_random_wil_param", f.getName(), "Subnet random", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Statistic results.Wilcoxon test"));
		if(images){
			createImages(outputFileName+"_sn_random_relBet", relBetSubnet1, "subnet1", relBetRandoms, "randoms", "sn_random_relBet", "relBet", "Images.Subnet random");
			createImages(outputFileName+"_sn_random_conn", connSubnet1, "subnet1", connRandoms, "randoms", "sn_random_conn", "conn", "Images.Subnet random");
			createImages(outputFileName+"_sn_random_clust", clustSubnet1, "subnet1", clustRandoms, "randoms", "sn_random_clust", "clust", "Images.Subnet random");
		}
		System.out.println("Finished 2nd Analysis..................");
	}
	private void statsTwoListsAnalisys(String side) throws IOException{
		// 1st Analysis

		WilcoxonTest wTest = new WilcoxonTest();
		List<String> list1 = IOUtils.readLines(commandLine.getOptionValue("list1"));
		List<String> list2 = IOUtils.readLines(commandLine.getOptionValue("list2"));

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
		wTest.test(wBinPath, outputFileName+"_list1_list2_wil.txt", side, "list1", "list2", relBetList1, relBetList2, connList1, connList2, clustList1, clustList2);
		if(images){
			createImages(outputFileName+"_list1_list2_relBet", relBetList1, "list1", relBetList2, "list2", "list1_list2_relBet", "relBet", "Images.List1 vs List2");
			createImages(outputFileName+"_list1_list2_conn", connList1, "list1", connList2, "list2", "list1_list2_conn", "conn", "Images.List1 vs List2");
			createImages(outputFileName+"_list1_list2_clust", clustList1, "list1", clustList2, "list2", "list1_list2_clust", "clust", "Images.List1 vs List2");
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
		wTest.test(wBinPath, outputFileName+"_sn1_sn2_wil.txt", side, "sn1", "sn2", relBetSubnet1, relBetSubnet2, connSubnet1, connSubnet2, clustSubnet1, clustSubnet2);
		if(images){
			createImages(outputFileName+"_sn1_sn2_relBet", relBetSubnet1, "subnet1", relBetSubnet2, "subnet2", "sn1_sn2_relBet", "relBet", "Images.Subnet1 vs Subnet2");
			createImages(outputFileName+"_sn1_sn2_conn", connSubnet1, "subnet1", connSubnet2, "subnet2", "sn1_sn2_conn", "conn", "Images.Subnet1 vs Subnet2");
			createImages(outputFileName+"_sn1_sn2_clust", clustSubnet1, "subnet1", clustSubnet2, "subnet2", "sn1_sn2_clust", "clust", "Images.Subnet1 vs Subnet2");
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
			System.out.println("Randoms["+i+"]: V = "+subgraph.getVertices().size()+" E = "+subgraph.getEdges().size());
			if(commandLine.hasOption("intermediate")) {
				Subgraph.OneIntermediateList(proteinNetwork.getInteractomeGraph(), subgraph);
				System.out.println("Randoms intermediate["+i+"]: V = "+subgraph.getVertices().size()+" E = "+subgraph.getEdges().size());
			}

			ProteinNetwork subProteinNetwork = createSubnet(subgraph);
			sbTopo.append(getTopologicalValues(subProteinNetwork, i)).append(System.getProperty("line.separator"));
			sbMeans.append(getTopologicalMeanValues(subProteinNetwork, i)).append(System.getProperty("line.separator"));

			if(!noNumberOfComponents) {
				List<List<ProteinVertex>> componentsList =  subProteinNetwork.getInteractomeGraph().getAllInformationComponents(true);
				sbComponents.append(getComponentsValues(subProteinNetwork, i, componentsList));
			}
			subProteinNetworkRandoms.add(subProteinNetwork);
		}

		//			result.addOutputItem(new Item("randoms means", outputFileName+"_sn_1-"+(randoms)+"_means.txt", "Random: means", Item.TYPE.FILE, Arrays.asList("randoms","means"),new HashMap<String,String>(),"Output data"));
		//			result.addOutputItem(new Item("randoms topo", outputFileName+"_sn_1-"+(randoms)+"_topo.txt", "Random: topological values", Item.TYPE.FILE, Arrays.asList("randoms","topological"),new HashMap<String,String>(),"Output data"));
		//			result.addOutputItem(new Item("randoms components", outputFileName+"_sn_1-"+(randoms)+"_comp.txt", "Randoms: components", Item.TYPE.FILE, Arrays.asList("randoms","components"),new HashMap<String,String>(),"Output data"));

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
//	private boolean checkErrors(){
//		if (!interactome.equalsIgnoreCase("sce") || !interactome.equalsIgnoreCase("hsa") || interactome.equalsIgnoreCase("own"))
//			return false;
//		if(!type.equalsIgnoreCase("proteins") || !type.equalsIgnoreCase("genes"))
//			return false;
//		
//		return true;
//	}
}
