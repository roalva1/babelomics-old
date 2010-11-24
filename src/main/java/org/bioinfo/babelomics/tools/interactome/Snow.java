package org.bioinfo.babelomics.tools.interactome;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.chart.BoxPlotChart;
import org.bioinfo.commons.io.utils.FileUtils;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.graph.SimpleUndirectedGraph;
import org.bioinfo.data.graph.Subgraph;
import org.bioinfo.data.graph.alg.Calc;
import org.bioinfo.data.graph.edge.DefaultEdge;
import org.bioinfo.infrared.common.dbsql.DBConnector;
import org.bioinfo.infrared.common.feature.FeatureList;
import org.bioinfo.infrared.core.XRef;
import org.bioinfo.infrared.core.dbsql.XRefDBManager;
import org.bioinfo.networks.protein.InteractomeParser;
import org.bioinfo.networks.protein.KSTest;
import org.bioinfo.networks.protein.ProteinNetwork;
import org.bioinfo.networks.protein.ProteinNetworkToFile;
import org.bioinfo.networks.protein.ProteinVertex;
import org.bioinfo.networks.protein.files.Dot;
import org.bioinfo.networks.protein.files.Json;
import org.bioinfo.networks.protein.files.Sif;
import org.bioinfo.networks.protein.files.Svg;
import org.bioinfo.networks.protein.files.Xml;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;
import org.jfree.chart.plot.PlotOrientation;

public class Snow  extends BabelomicsTool{

	private File inputFile;
	private String interactome;
	private String outputFileName;
	private String type;
	private String group;
	private ProteinNetwork proteinNetwork;
	private ProteinNetwork subProteinNetwork1, subProteinNetwork2;
	private List<ProteinNetwork> subProteinNetworkRandoms;
	private boolean components;
	private boolean bicomponents;
	private boolean intermediate;
	private boolean images;
	private boolean json;
	private boolean xml;
	private boolean sif;
	private Set<String> intermediatesSub1, intermediatesSub2; 
	private List<List<ProteinVertex>> componentsListSub1, componentsListSub2;
	private int randomSize;
	private Map<String, String> mapList1, mapList2;
	
	private DBConnector dbConnector;
	private XRefDBManager xrefDBMan;

	//private String wBinPath;

	@Override
	public void initOptions() {
		options.addOption(OptionFactory.createOption("interactome", "Interactome: hsa, sce or own (for you own interactome)", true));
		options.addOption(OptionFactory.createOption("sif-file", "An input file containing a SIF interactome", false));
		options.addOption(OptionFactory.createOption("topo-file","t", "An input file containing topological values", false, true));
		options.addOption(OptionFactory.createOption("list1", "An input file containing a node per line", false));
		options.addOption(OptionFactory.createOption("list2", "An input file containing a node per line", false));
		options.addOption(OptionFactory.createOption("type", "An argument saying if you want genes, proteins(default) or transcripts", false, true));
		

		options.addOption(OptionFactory.createOption("randoms", "Number of randoms", false, true));
		options.addOption(OptionFactory.createOption("intermediate", "If we want the intermediate we put 1, otherwise we put 0", false, true));
		options.addOption(OptionFactory.createOption("components", "If we want the number of components we put 1, otherwise we put 0", false, false));
		options.addOption(OptionFactory.createOption("bicomponents", "If we want the number of bicomponents we put 1, otherwise we put 0", false, false));
		options.addOption(OptionFactory.createOption("json", "It will create an output .json file", false, false));
		options.addOption(OptionFactory.createOption("o-sif-topo-file", "Create a full topological file from a sif file", false, false));
		options.addOption(OptionFactory.createOption("o-name", "If there is this argument, it will create an output .cmp file for the information of each component", false, true));
		options.addOption(OptionFactory.createOption("side", "side for kolmogorov and wilkoxon test. Can be less or greater", false, true));
		options.addOption(OptionFactory.createOption("images", "Print the images for the statistics", false, false));
		options.addOption(OptionFactory.createOption("xml", "Output xml file with the representation of the graph", false, false));
		options.addOption(OptionFactory.createOption("sif", "Output sif file with the representation of the graph", false, false));
		options.addOption(OptionFactory.createOption("group", "Values all(by default) or curated. It is An argument saying whether you want a curated interactome or the whole interactome", true, true));


	}
	@Override
	protected void execute() {
		File f = null;
		
		List<String> listToVertex1 = new ArrayList<String>();
		List<String> listToVertex2 = new ArrayList<String>();
		mapList1 = new HashMap<String, String>();
		mapList2 = new HashMap<String, String>();
		
//		wBinPath = babelomicsHomePath + "/bin/snow/wilcoxtest.r";
		subProteinNetwork1 = null;
		subProteinNetwork2 = null;

		subProteinNetworkRandoms = new ArrayList<ProteinNetwork>();
		components = commandLine.hasOption("components");
		bicomponents = commandLine.hasOption("bicomponents");
		images = commandLine.hasOption("images");
		json = commandLine.hasOption("json");
		xml = commandLine.hasOption("xml");
		sif = commandLine.hasOption("sif");

		//xml = true;
		//json = false;

		outputFileName = outdir + "/" + commandLine.getOptionValue("o-name", "result");
		
		
		interactome = commandLine.getOptionValue("interactome");
		String interactomeMsg = getInteractomeMsg();
		dbConnector = new DBConnector(interactome, new File(babelomicsHomePath + "/conf/infrared.properties"));
		xrefDBMan = new XRefDBManager(dbConnector);
		
		result.addOutputItem(new Item("interactome_param", interactomeMsg, "Interactome", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));

		String side = !commandLine.hasOption("side") ? "less" : commandLine.getOptionValue("side");
		result.addOutputItem(new Item("side_param", side, "Side for kolmogorov and wilkoxon test", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));

		type = !commandLine.hasOption("type") ? "proteins":commandLine.getOptionValue("type");
		result.addOutputItem(new Item("type_param", type, "ID type", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));

		group = !commandLine.hasOption("group") ? "all": commandLine.getOptionValue("group");
		result.addOutputItem(new Item("group_param", group, "Interactome group (curated or all)", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
		
		intermediate = (commandLine.hasOption("intermediate") && !"0".equalsIgnoreCase(commandLine.getOptionValue("intermediate")));

		ProteinNetworkToFile file = new ProteinNetworkToFile();

		try {

			String folderInteractions = this.babelomicsHomePath + "/conf/interactions/";
			SimpleUndirectedGraph<ProteinVertex, DefaultEdge> interactomeGraph = null;

			
			if (interactome.equals("own")){
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
						System.out.println("File created...");
						return;
					}
				} else {
					System.err.println("Missing custom interactome");
					return;					
				}
			} else {
				//you have to be very strict with the name of the files, otherwise the program won't work fine
				String localType = "genes";
				if(!type.equalsIgnoreCase("genes"))
					localType = "proteins";
				String parseFromSifFile = folderInteractions+interactome+"_alldb_"+group+"physical_"+localType+"_interactome_nr.sif";
				interactomeGraph = InteractomeParser.parseFromSifFile(parseFromSifFile);
				System.out.println("Interactome read: "+interactome+"_alldb_"+group+"physical_"+localType+"_interactome_nr.sif");
				proteinNetwork = new ProteinNetwork(interactomeGraph);
				proteinNetwork.loadTopologicalValues(folderInteractions+interactome+"_alldb_"+group+"physical_"+localType+"_interactome_nr_topo.txt");
				System.out.println("Interactome topo values read: "+interactome+"_alldb_"+group+"physical_"+localType+"_interactome_nr_topo.txt");
//				System.err.println("Unknown interactome");
//				return;
			}
			
			
			StringBuilder sbMeans = createMeansHeader();
			StringBuilder sbTopo = createTopoHeader();
			StringBuilder sbComponents = createComponentsHeader();

			if(commandLine.hasOption("list1")) {
				logger.debug("Starting list1.........");
				String nodeFile = commandLine.getOptionValue("list1");
				int node=1;
				FileUtils.checkFile(new File(nodeFile));
				List<String> list1 = IOUtils.readLines(nodeFile);
				listToVertex1 = new ArrayList<String>();
				if(type.equalsIgnoreCase("transcripts")){
					this.mapList1 = transcriptToUniprot(list1);
					listToVertex1.addAll(mapList1.keySet());
				}
				else if(type.equalsIgnoreCase("genes")){
					this.mapList1 = getGenEnsemble(list1);
					listToVertex1.addAll(mapList1.keySet());
				}
				if(this.mapList1.isEmpty())
					listToVertex1 = list1;
				
				logger.debug("nodes read: " + list1.toString());
				System.out.println("nodes read: " + listToVertex1.toString());

				SimpleUndirectedGraph<ProteinVertex, DefaultEdge> subgraph = (SimpleUndirectedGraph<ProteinVertex, DefaultEdge>) Subgraph.randomSubgraph(proteinNetwork.getInteractomeGraph(), toVertex(listToVertex1));
				System.out.println("Before intermediate: "+subgraph.getVertices().size()+" nodes");
				this.randomSize = subgraph.getVertices().size();
				if(intermediate){
					intermediatesSub1 = Subgraph.OneIntermediateList(proteinNetwork.getInteractomeGraph(), subgraph);
					if(intermediatesSub1.size()>0)
						result.addOutputItem(new Item("external_nodes_list_"+node, intermediatesSub1.toString().substring(1, intermediatesSub1.toString().length()-1), "External nodes added", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Subnet results.List " + node));
					else
						result.addOutputItem(new Item("external_nodes_list_"+node, "No added", "External nodes added", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Subnet results.List " + node));
				}
				System.out.println("With intermediate: "+subgraph.getVertices().size()+" nodes");
				subProteinNetwork1 = createSubnet(subgraph);

				sbTopo.append(getTopologicalValues(subProteinNetwork1, node, false));
				f = new File(outputFileName+"_sn_nodeFile"+node+"_topo.txt");
				IOUtils.write(f.getAbsoluteFile(), sbTopo.toString());
				result.addOutputItem(new Item("sn_nodeFile"+node+"_topo_param", f.getName(), "Topological values", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Subnet results.List " + node));

				sbMeans.append(getTopologicalMeanValues(subProteinNetwork1, node));
				f = new File(outputFileName+"_sn_nodeFile"+node+"_means.txt");
				IOUtils.write(f.getAbsoluteFile(), sbMeans.toString());
				result.addOutputItem(new Item("sn_nodeFile"+node+"_means_param", f.getName(), "Means values", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Subnet results.List " + node));


				if(components) {
					logger.debug("Starting list1 components.........");
					System.out.println("Starting list1 components.........");
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
				List<String> list2 = IOUtils.readLines(nodeFile);
				
				listToVertex2 = new ArrayList<String>();
				if(type.equalsIgnoreCase("transcripts")){
					this.mapList2 = transcriptToUniprot(list2);
					listToVertex2.addAll(mapList2.keySet());
				}
				else if(type.equalsIgnoreCase("genes")){
					this.mapList2 = getGenEnsemble(list2);
					listToVertex2.addAll(mapList2.keySet());
				}
				if(this.mapList2.isEmpty())
					listToVertex2 = list2;
				
				logger.debug("nodes read: " + list2.toString());

				SimpleUndirectedGraph<ProteinVertex, DefaultEdge> subgraph = (SimpleUndirectedGraph<ProteinVertex, DefaultEdge>) Subgraph.randomSubgraph(proteinNetwork.getInteractomeGraph(), toVertex(listToVertex2));

				if(intermediate){
					intermediatesSub2 = Subgraph.OneIntermediateList(proteinNetwork.getInteractomeGraph(), subgraph);
					if(intermediatesSub1.size()>0)
						result.addOutputItem(new Item("external_nodes_list_"+node, intermediatesSub2.toString().substring(1, intermediatesSub2.toString().length()-1), "External nodes added", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Subnet results.List " + node));
					else
						result.addOutputItem(new Item("external_nodes_list_"+node, "No added", "External nodes added", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Subnet results.List " + node));
				}
				subProteinNetwork2 = createSubnet(subgraph);

				sbTopo.append(getTopologicalValues(subProteinNetwork2, node, false));
				f = new File(outputFileName+"_sn_nodeFile"+node+"_topo.txt");
				IOUtils.write(f.getAbsoluteFile(), sbTopo.toString());
				result.addOutputItem(new Item("sn_nodeFile"+node+"_topo_param", f.getName(), "Topografical values", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Subnet results.List " + node));

				sbMeans.append(getTopologicalMeanValues(subProteinNetwork2, node));
				f = new File(outputFileName+"_sn_nodeFile"+node+"_means.txt");
				IOUtils.write(f.getAbsoluteFile(), sbMeans.toString());
				result.addOutputItem(new Item("sn_nodeFile"+node+"_means_param", f.getName(), "Means values", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Subnet results.List " + node));

				if(components) {
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
				if(randomSize > 0)
					createRandoms(randoms, randomSize);
				result.addOutputItem(new Item("randoms_param", ""+randoms, "Number of randoms", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
				result.addOutputItem(new Item("randoms_size_param", ""+randomSize, "Size of randoms", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
			}

			result.addOutputItem(new Item("components_param", (commandLine.hasOption("components") ? "yes" : "no"), "Calculate the number of components", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
			result.addOutputItem(new Item("bicomponents_param", (commandLine.hasOption("bicomponents") ? "yes" : "no"), "Calculate the number of bicomponents", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
			result.addOutputItem(new Item("intermediate_param", (intermediate ? "1" : "0"), "Max. number of external proteins introduced", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));

			if(commandLine.hasOption("randoms") && (!commandLine.hasOption("list2") || "none".equalsIgnoreCase(commandLine.getOptionValue("list2")))) {
				statsOneListAnalysis(side);
			} else if(!commandLine.hasOption("randoms") && commandLine.hasOption("list1") && commandLine.hasOption("list2") && !"none".equalsIgnoreCase(commandLine.getOptionValue("list2")))
				statsTwoListsAnalisys(side);
			else{
				logger.error("Not correct arguments for statistic test");
				return;
			}
			if(json){
				File auxFile = new File(outputFileName);;
				List<String> names = new ArrayList<String>();
				if(subProteinNetwork1 != null){
					names.add(auxFile.getName() + "_list1.json");
					createJson(subProteinNetwork1, outputFileName+"_list1.dot", componentsListSub1, intermediatesSub1, 1);
				}
				if(subProteinNetwork2 != null){
					names.add(auxFile.getName() + "_list2.json");
					createJson(subProteinNetwork2, outputFileName+"_list2.dot", componentsListSub2, intermediatesSub2, 2);
				}

			}
			if(xml){
				File xmlFile = null;
				Xml xmlObject;
				if(interactome.equalsIgnoreCase("own"))
					xmlObject = new Xml();
				else
					xmlObject = new Xml(this.dbConnector);
				//Map<String, String> mapList = new HashMap<String,String>();
				if(subProteinNetwork1 != null){
					xmlFile = new File(outdir+"/subnetwork1.xml");
//					if(type.equals("genes"))
//						mapList = this.getGenEnsemble(listToVertex1);
//					else if(type.equals("proteins"))
//						mapList = this.mapList1;
					xmlObject.graphToXML(xmlFile.getAbsolutePath(),subProteinNetwork1.getInteractomeGraph(), intermediatesSub1, componentsListSub1, type, this.mapList1);
					addOutputAppletItem(xmlFile, 1);
				}
				if(subProteinNetwork2 != null){
					xmlFile = new File(outdir+"/subnetwork2.xml");
					//mapList.clear();
//					if(type.equals("genes"))
//						mapList = this.getGenEnsemble(listToVertex2);
//					else if(type.equals("proteins"))
//						mapList = this.mapList2;
					xmlObject.graphToXML(xmlFile.getAbsolutePath(),subProteinNetwork2.getInteractomeGraph(), intermediatesSub2, componentsListSub2, type, this.mapList2);
					addOutputAppletItem(xmlFile, 2);
				}
			}

			if(sif){
				Sif sifObject = new Sif();
				File fSif;
				if(subProteinNetwork1 != null){
					fSif = new File(outputFileName+"_subnetwork1.sif");
					IOUtils.write(fSif.getAbsoluteFile(), sifObject.graphToSif(subProteinNetwork1.getInteractomeGraph()));
					result.addOutputItem(new Item("subnetwork1.sif", fSif.getName(), "subnetwork1 sif file", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Sif network file"));

				}
				if(subProteinNetwork2 != null){
					fSif = new File(outputFileName+"_subnetwork2.sif");
					IOUtils.write(fSif.getAbsoluteFile(), sifObject.graphToSif(subProteinNetwork2.getInteractomeGraph()));
					result.addOutputItem(new Item("subnetwork2.sif", fSif.getName(), "subnetwork2 sif file", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Sif network file"));

				}
			}
			//dbConnector.disconnect();
		} catch (SQLException e) {
			  //e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
			
	}

	private Map<String, String> transcriptToUniprot(List<String> list) {
		Map<String, String> map = new HashMap<String, String>();
		FeatureList<XRef> xrefsEns = new FeatureList<XRef>();
		for(String id : list){
			try {
				xrefsEns  = xrefDBMan.getByDBName(id, "uniprot_swissprot_accession");
				if(xrefsEns != null && !xrefsEns.isEmpty() && !xrefsEns.get(0).getId().equals(id))
					map.put(xrefsEns.get(0).getId(), id);
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		return map;
	}
	private String getInteractomeMsg() {
		String interactomeMsg = "Homo sapiens";
		if ("hsa".equalsIgnoreCase(interactome) ) {
			interactomeMsg = "Homo sapiens";
		} else if ("sce".equalsIgnoreCase(interactome)) {
			interactomeMsg = "Saccharomyce cerevisiae";
		} else if ("bta".equalsIgnoreCase(interactome)) {
			interactomeMsg = "Bos taurus";
		} else if ("dme".equalsIgnoreCase(interactome)) {
			interactomeMsg = "Drosophila melanogaster";
		} else if ("mmu".equalsIgnoreCase(interactome)) {
			interactomeMsg = "Mus musculus";
		} else if ("ath".equalsIgnoreCase(interactome)) {
			interactomeMsg = "Arabidopsis thaliana";
		} else if ("own".equalsIgnoreCase(interactome)) {
			interactomeMsg = new File(commandLine.getOptionValue("sif-file")).getName();
		} else {
			interactomeMsg = "unknown";
		}
			return interactomeMsg;
	}

	private void createJson(ProteinNetwork proteinNetwork, String sourceDotFile, List<List<ProteinVertex>> componentsListSub, Set<String> intermediatesSub, int node) throws IOException{

		Dot<ProteinVertex, DefaultEdge> dot = new Dot<ProteinVertex, DefaultEdge>();
		Svg svg = new Svg();
		Json<ProteinVertex, DefaultEdge> json = new Json<ProteinVertex, DefaultEdge>();
		List<File> fileList = new ArrayList<File>();
		List<String> layoutsName = new ArrayList<String>();

		if(componentsListSub == null)
			logger.error("not components calculated, please, set the --components option");

		IOUtils.write(sourceDotFile, dot.toDot(proteinNetwork.getInteractomeGraph()));

		File dotFile = new File(outputFileName+"_list"+node+"_dot.svg");
		IOUtils.write(dotFile, svg.toSvg(sourceDotFile, "dot"));
		fileList.add(dotFile);
		layoutsName.add("dot");

		File twopiFile = new File(outputFileName+"_list"+node+"_twopi.svg");
		IOUtils.write(twopiFile, svg.toSvg(sourceDotFile, "twopi"));
		fileList.add(twopiFile);
		layoutsName.add("twopi");

		File f = new File(outputFileName+"_list"+node+".json");
		IOUtils.write(f.getAbsoluteFile(), json.toJson(fileList, layoutsName, proteinNetwork.getInteractomeGraph(), intermediatesSub, componentsListSub));
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
		List<ProteinVertex> list = this.subProteinNetwork1.getInteractomeGraph().getVertices();


		List<Double> relBetInter = proteinNetwork.getBetRelList();
		List<Double> relBetList1 = new ArrayList<Double>();

		List<Double> connInter = proteinNetwork.getConnList();
		List<Double> connList1 = new ArrayList<Double>();

		List<Double> clustInter = proteinNetwork.getClustList();
		List<Double> clustList1 = new ArrayList<Double>();

		//createTopoFilterList(list, relBetList1, connList1, clustList1);
		createTopoFilterListNoIntermediates(list, relBetList1, connList1, clustList1);
		
		String toWrite = ksTest(relBetList1, relBetInter, connList1, connInter, clustList1, clustInter, side);
		if(!toWrite.equals("")){
			f = new File(outputFileName+"_list_inter_kol.txt");
			IOUtils.write(f.getAbsolutePath(), toWrite);
			result.addOutputItem(new Item("list_inter_kol_param", f.getName(), "List - Interactome", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Statistic results.Kolmogorov-Smirnov test"));
		}
		else
			result.addOutputItem(new Item("list_inter_kol_param", "Empty results", "List - Interactome", Item.TYPE.MESSAGE, new ArrayList<String>(), new HashMap<String,String>(), "Statistic results.Kolmogorov-Smirnov test"));

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
			int  randomVertex = (int) (Math.random() * proteinNetwork.getInteractomeGraph().getVertices().size());
			ProteinVertex v = proteinNetwork.getInteractomeGraph().getVertices().get(randomVertex);
			relBetRandoms.add(v.getRelativeBetweenness());
			connRandoms.add(Double.parseDouble(Integer.toString(proteinNetwork.getInteractomeGraph().getDegree(v))));
			clustRandoms.add(v.getClusteringCoefficient());
		}

		toWrite = ksTest(relBetSubnet1, relBetRandoms, connSubnet1, connRandoms, clustSubnet1, clustRandoms, side);
		if(!toWrite.equals("")){
			f = new File(outputFileName+"_sn_random_kol.txt");
			IOUtils.write(f.getAbsolutePath(), toWrite);
			result.addOutputItem(new Item("sn_random_kol_param", f.getName(), "Subnet - Random", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Statistic results.Kolmogorov-Smirnov test"));
		}
		else
			result.addOutputItem(new Item("sn_random_kol_param", "Empty results", "Subnet - Random", Item.TYPE.MESSAGE, new ArrayList<String>(), new HashMap<String,String>(), "Statistic results.Kolmogorov-Smirnov test"));

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
		sb.append("betweenness\t"+kstest.resultKolmogorovSmirnov(ListUtils.toDoubleArray(relBetList1), ListUtils.toDoubleArray(relBetList2), side).getPValue()+"\t"+side+"\n");
		sb.append("connections\t"+kstest.resultKolmogorovSmirnov(ListUtils.toDoubleArray(connList1), ListUtils.toDoubleArray(connList2), side).getPValue()+"\t"+side+"\n");
		sb.append("coefficient\t"+kstest.resultKolmogorovSmirnov(ListUtils.toDoubleArray(clustList1), ListUtils.toDoubleArray(clustList2), side).getPValue()+"\t"+side);
		return sb.toString();
	}
	private void createTopoFilterListNoIntermediates(List<ProteinVertex> list, List<Double> relBetList, List<Double> connList, List<Double> clustList){
		for(ProteinVertex protein : list){
			if( (this.intermediatesSub1 != null && this.intermediatesSub1.contains(protein.getId())) || proteinNetwork.getInteractomeGraph().getVertex(protein.getId()) == null)
				continue;
			relBetList.add(proteinNetwork.getRelBetweennessVertex(protein.getId()));
			connList.add(proteinNetwork.getDegreeVertex(protein.getId()));
			clustList.add(proteinNetwork.getClusterinVertex(protein.getId()));
		}
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

			if(intermediate) {
				Subgraph.OneIntermediateList(proteinNetwork.getInteractomeGraph(), subgraph);
				logger.debug("Randoms intermediate["+i+"]: V = "+subgraph.getVertices().size()+" E = "+subgraph.getEdges().size());
			}
			ProteinNetwork subProteinNetwork = createSubnet(subgraph);
			logger.debug("Subnet created");
			sbTopo.append(getTopologicalValues(subProteinNetwork, i, true)).append(System.getProperty("line.separator"));
			sbMeans.append(getTopologicalMeanValues(subProteinNetwork, i)).append(System.getProperty("line.separator"));

			if(components) {
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
		subProteinNetwork.calcTopologicalValues();
		subProteinNetwork.calcTopologicalMeanValues();
		return subProteinNetwork;
	}

	private StringBuilder createMeansHeader(){
		StringBuilder sb = new StringBuilder();
		sb.append("#Subnet\tMeanBet\tStdBet\t").append("MeanCon\tStdCon\t").append("MeanCls\tStdCls");
		if(components){
			sb.append("\tComp\t1Comp");
		}
		sb.append("\tNodes");
		if(bicomponents){
			sb.append("\tBiComp\t");
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


	private String getTopologicalValues(ProteinNetwork subProteinNetwork, int subnet, boolean random){
		StringBuilder sb = new StringBuilder();
		for (ProteinVertex proteinVertex : subProteinNetwork.getInteractomeGraph().getVertices()) {
			if(proteinVertex == null)
				continue;
			String inputId= "";
			sb.append("sn"+(subnet)).append("\t");
			if( (type.equals("transcripts") || type.equals("genes"))  && random == false){
				inputId = this.getValue(proteinVertex.getId());
			}
			
			sb.append(inputId+proteinVertex.getId()).append("\t");
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

		if(components){
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
		sb.append("\t").append(subProteinNetwork.getInteractomeGraph().getVertices().size());
		if(bicomponents){
			sb.append("\t").append(interactomeGraph.getNumberOfBicomponents());
		}

		return sb.toString();

	}
	private String getComponentsValues(ProteinNetwork subProteinNetwork, int subnet, List<List<ProteinVertex>> componentsList){
		StringBuilder sb = new StringBuilder();
		SimpleUndirectedGraph<ProteinVertex, DefaultEdge> interactomeGraph = subProteinNetwork.getInteractomeGraph();

		if(components){
			List<Double> componentsDiameter = Calc.calcDiameter(interactomeGraph, componentsList);
			for(int i=0; i < componentsList.size(); i++){
				sb.append("sn"+(subnet)).append("\t");
				StringBuilder sbNodes = new StringBuilder();
				for (int k = 0; k < componentsList.get(i).size(); k++) {
					String nodeId = "";
					if(type.equals("transcripts") || type.equals("genes")){
						nodeId = getValue(componentsList.get(i).get(k).getId());
					}
					sbNodes.append(nodeId+componentsList.get(i).get(k).getId());
					if(k!=componentsList.get(i).size()-1)
						sbNodes.append(",");
				}
				sb.append((i)+"\t"+componentsDiameter.get(i)+"\t"+componentsList.get(i).size()+"\t"+sbNodes.toString()+ System.getProperty("line.separator"));
			}
		}
		return sb.toString();
	}
	private String getValue(String idNode){
		String transcriptId="";
		if(this.mapList1.get(idNode) != null)
			transcriptId = this.mapList1.get(idNode)+":";
		else if(this.mapList2.get(idNode) != null)
			transcriptId = this.mapList2.get(idNode)+":";
		return transcriptId;
	}
	private void randomFilesToString(StringBuilder sb, String file) throws IOException{
		if(!sb.toString().equals("")){
			sb.deleteCharAt(sb.lastIndexOf(System.getProperty("line.separator")));
			IOUtils.write(file, sb.toString());
		}
	}
	private List<ProteinVertex> toVertex(List<String> vertices){
		List<ProteinVertex> verticesList = new ArrayList<ProteinVertex>(vertices.size());
//		if(type.equals("genes")){
//			FeatureList<XRef> xrefsEns = new FeatureList<XRef>();
//			for(String proteinName : vertices){
//				try {
//					xrefsEns  = xrefDBMan.getByDBName(proteinName, "ensembl_gene");
//					if(xrefsEns != null && !xrefsEns.isEmpty()){
//						verticesList.add(new ProteinVertex(xrefsEns.get(0).getId()));
//					}
//				} catch (Exception e) {
//				//	e.printStackTrace();
//				}
//			}
//			return verticesList;
//		}
		
		for(String proteinName : vertices){
			if(!proteinName.equals(""))
				verticesList.add(new ProteinVertex(proteinName));
		}
		return verticesList;
	}

	private Map<String,String> getGenEnsemble(List<String> vertices){
		Map<String, String> listGenEnsembl = new HashMap<String, String>();
		FeatureList<XRef> xrefsEns = new FeatureList<XRef>();
		for(String proteinName : vertices){
			try {
				xrefsEns  = xrefDBMan.getByDBName(proteinName, "ensembl_gene");
				if(xrefsEns != null && !xrefsEns.isEmpty() && !xrefsEns.get(0).getId().equals(proteinName))
					listGenEnsembl.put(xrefsEns.get(0).getId(), proteinName);
			} catch (Exception e) {
			//	e.printStackTrace();
			}
		}
		return listGenEnsembl;
	}
	public void addOutputAppletItem(File xmlFile, int index) {
		if (xmlFile.exists()) {
			String url = "SnowViewer?filename=" + xmlFile.getName() + "&width=600&height=600";
			result.addOutputItem(new Item("viewer" + index + "_param", url, "Viewer for network #" + index, TYPE.HTML, StringUtils.toList("SERVER,INCLUDE_REFS", ","), new HashMap<String, String>(2), "Network viewer"));

			url = "SnowViewer?filename=" + xmlFile.getName();
			result.addOutputItem(new Item("viewer" + index + "_param_new_window", url, "Open applet for network #" + index + " in a new window", TYPE.LINK, StringUtils.toList("SERVER,INCLUDE_REFS", ","), new HashMap<String, String>(2), "Network viewer"));
		}
	}

}
