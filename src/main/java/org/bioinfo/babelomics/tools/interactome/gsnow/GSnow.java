package org.bioinfo.babelomics.tools.interactome.gsnow;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.HashSet;
import java.util.List;
import java.util.Map;
import java.util.Set;

import org.apache.commons.math.stat.descriptive.rank.Max;
import org.bioinfo.babelomics.tools.interactome.RandomsSnowManager;
import org.bioinfo.babelomics.tools.interactome.SnowPrinter;
import org.bioinfo.babelomics.tools.interactome.SnowTool;
import org.bioinfo.babelomics.tools.interactome.gsnow.GSnowPreprocessing.ListInfo;
import org.bioinfo.babelomics.tools.interactome.gsnow.GSnowPreprocessing.Node;
import org.bioinfo.commons.io.TextFileReader;
import org.bioinfo.commons.io.utils.FileUtils;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.graph.SimpleUndirectedGraph;
import org.bioinfo.data.graph.Subgraph;
import org.bioinfo.data.graph.edge.DefaultEdge;
import org.bioinfo.infrared.common.DBConnector;
import org.bioinfo.infrared.core.XRefDBManager;
import org.bioinfo.infrared.core.feature.XRef;
import org.bioinfo.infrared.core.feature.XRef.XRefItem;
import org.bioinfo.infrared.core.funcannot.GO;
import org.bioinfo.infrared.funcannot.AnnotationDBManager;
import org.bioinfo.infrared.funcannot.GODBManager;
import org.bioinfo.math.result.KolmogorovSmirnovTestResult;
import org.bioinfo.math.util.MathUtils;
import org.bioinfo.networks.protein.ProteinNetwork;
import org.bioinfo.networks.protein.ProteinVertex;
import org.bioinfo.networks.protein.files.Sif;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;


public class GSnow extends SnowTool{

	private ListInfo listInfo;
	GSnowItem significantItem;
	Map<Integer, GSnowItem> gsnowItems;
	private DBConnector dbConnector;
	private XRefDBManager xrefDBMan;
	private GODBManager gof;
	private int numberOfStartingNodes = 10; // Here we indicate the number of nodes from where we will start the GSnow
	private static int numberOfMaxNodes = 200; // Here we indicate the max number of nodes that we will analyse
	private double significantValue; // Here we indicate which is the minimum significant value that we are interested in.
	private static double defaultSignificantValue = 0.05; // Here we indicate the default significant value that we are interested in.
	private String order;
	private Integer numberItems;
	private double cutOff;
	private String selectOption;
	
	//Important info mapNames<String,String> = <Ensenbl_gen, normal_name>
	private Map<String, String> mapNames; 

	private String decimalFormat;
	@Override
	public void initOptions() {
		
		options.addOption(OptionFactory.createOption("list", "An input file containing a node per line", false, true));
		options.addOption(OptionFactory.createOption("order", "Here we indicate wether we want to order the list ascendant(default) or descendant", false, true));
		options.addOption(OptionFactory.createOption("number-items", "Here we indicate how many nodes we want to process", false, true));
		options.addOption(OptionFactory.createOption("cut-off", "Here we indicate where we cut the list from (considering the order)", false, true));
		options.addOption(OptionFactory.createOption("significant-value", "Here we indicate the significant value", false, true));
		options.addOption(OptionFactory.createOption("select-mcn", "Here we indicate how we can select the chosen network. By the most important value or by the most inportant value considering the 'significant-value'(abs-min, rel-min)", false, true));

		//		options.addOption(OptionFactory.createOption("filter-list", "Here we indicate if we won't filter(default) or filter by number of items or filter by cutting from a value", false, true));

		// This next two options are only used when we want to generate the randoms for the test, only from time to time(maybe twice a year).
		options.addOption(OptionFactory.createOption("size-min", "Minimum size for randoms", false, true));
		options.addOption(OptionFactory.createOption("size-max", "Maximum size for randoms", false, true));
	}

	@Override
	protected void execute() {
		this.decimalFormat = "#.####";
		initExecute();
		if(commandLine.hasOption("list") && !(commandLine.hasOption("size-min") && commandLine.hasOption("size-max")))
			executeGSnow();
		else if(options.hasOption("size-min") && options.hasOption("size-max"))
			executeGenerateRandoms();
		else
			System.err.println("[Error] Some parameters missing.");
	}
	
	
	private void executeGSnow(){
		try {
			
			SnowPrinter snowPrinter = new SnowPrinter();
			String interactomeMsg = getInteractomeMsg();
//			filterList = commandLine.hasOption("filter-list") ? commandLine.getOptionValue("filter-list") : "wholeList";
			order = commandLine.hasOption("order") ? commandLine.getOptionValue("order") : "ascendant";
			numberItems = commandLine.hasOption("number-items") ? Integer.parseInt(commandLine.getOptionValue("number-items")) : numberOfMaxNodes;
			cutOff = commandLine.hasOption("cut-off") ? Double.parseDouble(commandLine.getOptionValue("cut-off")) : Double.NaN;
			significantValue = commandLine.hasOption("significant-value") ? Double.parseDouble(commandLine.getOptionValue("significant-value")) : defaultSignificantValue;

			selectOption = commandLine.hasOption("select-mcn") ? commandLine.getOptionValue("select-mcn") : "abs-min";

			
			dbConnector = new DBConnector(interactome, new File(babelomicsHomePath + "/conf/infrared.properties"));
			xrefDBMan = new XRefDBManager(dbConnector);
			gof = new GODBManager(dbConnector);
			logger.debug("Starting list.........");
			String nodeFile = commandLine.getOptionValue("list");

			
			GSnowPreprocessing preprocessing = new GSnowPreprocessing(proteinNetwork, type, xrefDBMan, /*numberOfStartingNodes,*/ numberOfMaxNodes, order, numberItems, cutOff);
			
			result.addOutputItem(new Item("interactome_param", interactomeMsg, "Interactome", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
			result.addOutputItem(new Item("group_param", group, "Interactome group (curated or all)", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
			result.addOutputItem(new Item("intermediate_param", (intermediate ? "1" : "0"), "Number of external proteins introduced", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
//			result.addOutputItem(new Item("select_option", selectOption, "Option", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));

			String folder = loadFile();
			System.out.println(folder);
			
			FileUtils.checkFile(new File(nodeFile));
			listInfo = preprocessing.preprocess(nodeFile);
			
			mapNames = new HashMap<String, String>();
			
			for(Node n : listInfo.getNodes()){
				mapNames.put(n.getId(), n.getOriginalId());
			}
			updateJobStatus("5","List preprocessed");
			
			File f = new File(outputFileName+"_nodes.txt");
			IOUtils.write(f.getAbsoluteFile(), snowPrinter.printNodesList(listInfo.getNodes()));
			result.addOutputItem(new Item("nodes_file", f.getName(), "Final list", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Summary"));
			
			f = new File(outputFileName+"_not_matched_nodes.txt");
			IOUtils.write(f.getAbsoluteFile(), snowPrinter.printNodes(listInfo.getNotMatchNodes()));
			result.addOutputItem(new Item("list_info_not_matched_nodes", f.getName(), "Not matched nodes", Item.TYPE.FILE, new ArrayList<String>(), new HashMap<String,String>(), "Summary"));
			
			result.addOutputItem(new Item("list_info_repeated_nodes", (listInfo.getRepeatedNodes().size()==0) ? "0" : listInfo.getRepeatedNodes().toString(), "Duplicated nodes", Item.TYPE.MESSAGE, new ArrayList<String>(), new HashMap<String,String>(), "Summary"));
			//result.addOutputItem(new Item("list_info_cutted_nodes", listInfo.getCuttedNodes().toString(), "Not considered nodes because the list is too big", Item.TYPE.MESSAGE, new ArrayList<String>(), new HashMap<String,String>(), "Summary"));
			
			if(listInfo.getNodes().size() < this.numberOfStartingNodes){
				result.addOutputItem(new Item("small_size_list_error", "List to small after preprocessing it. The size should have at least "+numberOfStartingNodes+" after preprocessing ", "An error occurred", Item.TYPE.MESSAGE, Arrays.asList("ERROR"), new HashMap<String, String>(), "Error")); 
				return;
			}
			// 1º Here the map is (size of list, nodes of list)
			gsnowItems = getGSnowItems();
			updateJobStatus("10","List preprocessed");
		
			// 2º Here the map is (size of list, value of this list)
			getDataMatrixList();
			updateJobStatus("15","List preprocessed");
			
			// 3º Here the map is (size of list, value of this list of each random)
			Map<Integer, List<Double>> dataMatrixRandoms = getDataMatrixRandoms(folder);
			updateJobStatus("20","List preprocessed");
			
			// 4º Here the map is (size of list, value of this list compared with randoms)
			compareDataMatrixListRandoms(dataMatrixRandoms);
			updateJobStatus("25","List preprocessed");
			
			// 5º Here we calculate the significant item from gsnowItems the 
			if(selectOption.equalsIgnoreCase("abs-min"))
				significantItem = getSignificatValueAbsMinOption();
			else if(selectOption.equalsIgnoreCase("rel-min"))
				significantItem = getSignificatValueRelMinOption();
			
			
			f = new File(outputFileName+"_all.txt");
			IOUtils.write(f.getAbsoluteFile(), snowPrinter.printGsnowItems(gsnowItems, numberOfStartingNodes));
			result.addOutputItem(new Item("all_results", f.getName(), "All results", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Results"));
			
			
			
			if(significantItem.getComparedValue() > significantValue){
				result.addOutputItem(new Item("significant_value", " There is no significant minimal connected network", "Significant value", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Results"));
				return;
			}
			
			System.out.println("Significant value:"+significantItem.getComparedValue());
			System.out.println("Significant size:"+significantItem.getNodes().size());
			String significantStringItem = StringUtils.decimalFormat(significantItem.getComparedValue(), decimalFormat);
			
			if(significantItem.getComparedValue() == 0)
				significantStringItem = "0.0001";
			result.addOutputItem(new Item("significant_value", "<"+significantStringItem, "pval of MCN chosen", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Results"));
			result.addOutputItem(new Item("significant_size", Integer.toString(significantItem.getNodes().size()), "size of MCN chosen", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Results"));
			
			f = new File(outputFileName+"_size_pvalue.json");
			IOUtils.write(f.getAbsoluteFile(), snowPrinter.getJsonSizePValue(gsnowItems, numberOfStartingNodes));
			List<String> tags = new ArrayList<String>();
			tags.add("NETWORKMINER_SIZE_PVALUE");
			result.addOutputItem(new Item("plot_size_pvalue", f.getName(), "Plot", TYPE.IMAGE, tags, new HashMap<String, String>(2), "Results"));

			// starting SNOW analyses
			this.executeSnow();
			
		} catch (IOException e) {
			e.printStackTrace();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	
	// 1º Here the map is (size of list, nodes of list)
	// This method will have each size[0...listInfo.getNodes().size()] of listInfo with its GSnowItems
	private Map<Integer, GSnowItem> getGSnowItems() {
		
		Map<Integer, GSnowItem> gsnowItems = new HashMap<Integer, GSnowItem>();
		List<Node> nodes = new ArrayList<Node>();
		GSnowItem gsnowItem;
		List<Node> auxNodes;
		for(int i = 0; i < this.listInfo.getNodes().size(); i++){
			nodes.add(this.listInfo.getNodes().get(i));
			auxNodes = new ArrayList<Node>();
			auxNodes.addAll(nodes);
			gsnowItem = new GSnowItem();
			gsnowItem.setNodes(auxNodes);
			gsnowItems.put(i+1, gsnowItem);
//			System.out.println(gsnowItems.get(i+1));
		}
		return gsnowItems;
	}
	
	// 2º Here the map is (size of list, value of this list)
	// This method will fill Map<Integer, GSnowItem> gsnowItems with its raw value
	private void getDataMatrixList(){
		List<ProteinVertex> nodes = new ArrayList<ProteinVertex>();
		List<List<ProteinVertex>> componentsList;
		List<Double> componentsSize;
		double rawValue;
		for(int i=0; i < numberOfStartingNodes; i++){
			nodes.add(new ProteinVertex(listInfo.getNodes().get(i).getId()));
		}
		SimpleUndirectedGraph<ProteinVertex, DefaultEdge> subgraph;
		
		//this.listInfo.getNodes().size() is always <= numberOfMaxNodes, because I cut the list (if it is longer than numberOfMaxNodes) in the preprocessing
		for(int i=numberOfStartingNodes; i <= this.listInfo.getNodes().size(); i++){
			subgraph = (SimpleUndirectedGraph<ProteinVertex, DefaultEdge>) Subgraph.randomSubgraph(proteinNetwork.getInteractomeGraph(), nodes);
			if(intermediate){
//				System.out.println("Nodes size before intermediate: "+subgraph.getVertices().size());
				Subgraph.OneIntermediateList(proteinNetwork.getInteractomeGraph(), subgraph);
//				System.out.println("Nodes size after intermediate: "+subgraph.getVertices().size());
			}
			componentsSize = new ArrayList<Double>();
			componentsList = subgraph.getAllInformationComponents(true);
			for(List<ProteinVertex> component : componentsList)
				componentsSize.add((double)component.size());
			rawValue = (double)subgraph.getVertices().size()/(double)componentsList.size();
			gsnowItems.get(i).setRawValue(rawValue);
			gsnowItems.get(i).setComponentsSize(componentsSize);
			//System.out.println("Size: "+gsnowItems.get(i).getNodes().size()+":"+rawValue);
			if(i < listInfo.getNodes().size())
				nodes.add(new ProteinVertex(listInfo.getNodes().get(i).getId()));
		}
	}
	
	// 3º Here the map is (size of list, value of this list of each random)
	// Here we will load the data of the randoms in Map<Integer, List<Double>>
	public Map<Integer, List<Double>> getDataMatrixRandoms(String folder) throws IOException{
		
		Map<Integer, List<Double>> values1 = new HashMap<Integer, List<Double>>();
		List<Double> values2 = null;
		TextFileReader tfr = new TextFileReader(folder);
		String line = null;
		String[] fields;
		int size=0;;
		while((line = tfr.readLine()) != null) {
			if(!line.startsWith("#")){
				values2 = new ArrayList<Double>();
				fields = line.split("\t");
				size = Integer.parseInt(fields[0]);
				for(int i=1; i <fields.length; i++){
					values2.add(Double.parseDouble(fields[i]));
				}
				values1.put(size, values2);
			}
		}
		tfr.close();
		return values1;
	}
	
	// 4º Here the map is (size of list, value of this list comparing with randoms)
	// Here we will get the compared value against the randoms
	private void compareDataMatrixListRandoms( Map<Integer, List<Double>> dataMatrixRandoms){
		for(int i = numberOfStartingNodes; i <= this.listInfo.getNodes().size(); i++){
			double rawValue = gsnowItems.get(i).getRawValue();
//			if(i==134)
//				System.out.println(gsnowItems.get(i).getRawValue());
			List<Double> valueList = dataMatrixRandoms.get(i);
//			System.out.println("("+i+")List["+gsnowItems.get(i).getNodes().size()+"]:"+rndValue);
			List<Double> intValue = new ArrayList<Double>();
			int rawBiggerValue = 0;
			int valueBiggerRaw = 0;
			for(double value : valueList){
				if(rawValue > value){
					intValue.add(1.0);
					rawBiggerValue++;
				}
					
				else{
					intValue.add(0.0);
					valueBiggerRaw++;
				}
					
			}
//			if(i==134){
//				System.out.println("rawBiggerValue:"+rawBiggerValue);
//				System.out.println("valueBiggerRaw:"+valueBiggerRaw);
//				System.out.println(1-MathUtils.mean(ListUtils.toDoubleArray(intValue)));
//			}
				
			gsnowItems.get(i).setComparedValue(1-MathUtils.mean(ListUtils.toDoubleArray(intValue)));
		}
	}
	
	
	//5º Here we get the significant item in significantItem from gsnowItems
	private GSnowItem getSignificatValueAbsMinOption() {
		GSnowItem significantItemLocal = new GSnowItem();
		List<Node> significantNodes = new ArrayList<Node>();
		double significantValue = Double.MAX_VALUE;
		
		List<Double> componentsSizeGSnowItem = null;
		//This is a pointer to the elements 
		List<Integer> pointerBiggerThanInitSizeElements = new ArrayList<Integer>();
		List<Integer> pointerSmallerThanInitSizeElements = new ArrayList<Integer>();
		int initSizeElements = 5;
		double currentValue;
		Max maxUtil = new Max();
		for(int i = numberOfStartingNodes; i <= this.listInfo.getNodes().size(); i++){
			currentValue = gsnowItems.get(i).getComparedValue();
			componentsSizeGSnowItem = gsnowItems.get(i).getComponentsSize();
			if(currentValue <= this.significantValue)
				// La siguiente linea: cogemos el componente con más elemenentos y vemos si es mayor o menor que initSizeElements
				if(maxUtil.evaluate(ListUtils.toDoubleArray(componentsSizeGSnowItem)) < initSizeElements)
					pointerSmallerThanInitSizeElements.add(i);
				else
					pointerBiggerThanInitSizeElements.add(i);
		}
		
		double min = Double.MAX_VALUE;
		//double maxSizeComponent = Double.MIN_VALUE;
		int pointerToMin = 0;
		for(int i : pointerBiggerThanInitSizeElements){
			currentValue = gsnowItems.get(i).getComparedValue();
			componentsSizeGSnowItem = gsnowItems.get(i).getComponentsSize();
			if(currentValue < min /*&& maxUtil.evaluate(ListUtils.toDoubleArray(componentsSizeGSnowItem)) >  maxSizeComponent*/){
				min = currentValue;
				//maxSizeComponent = maxUtil.evaluate(ListUtils.toDoubleArray(componentsSizeGSnowItem));
				pointerToMin = i;
			}
			
		}
		
		if(pointerBiggerThanInitSizeElements.isEmpty()){
			for(int i : pointerSmallerThanInitSizeElements){
				currentValue = gsnowItems.get(i).getComparedValue();
				if(currentValue < min){
					min = currentValue;
					pointerToMin = i;
				}
			}
		}
		if(gsnowItems.containsKey(pointerToMin)){
			significantValue = gsnowItems.get(pointerToMin).getComparedValue();
			significantNodes = gsnowItems.get(pointerToMin).getNodes();
		}
		significantItemLocal.setNodes(significantNodes);
		significantItemLocal.setComparedValue(significantValue);
		return significantItemLocal;
	}
	private GSnowItem getSignificatValueRelMinOption() {
		GSnowItem significantItemLocal = new GSnowItem();
		List<Node> significantNodes = new ArrayList<Node>();
		double significantValue = Double.MAX_VALUE;
		
		List<Double> componentsSizeGSnowItem = null;
		//This is a pointer to the elements 
		List<Integer> pointerBiggerThanInitSizeElements = new ArrayList<Integer>();
		List<Integer> pointerSmallerThanInitSizeElements = new ArrayList<Integer>();
		int initSizeElements = 5;
		double currentValue;
		double nextCurrentValue;
		Max maxUtil = new Max();
		for(int i = numberOfStartingNodes; i <= this.listInfo.getNodes().size(); i++){
			currentValue = gsnowItems.get(i).getComparedValue();
			if(gsnowItems.get(i+1) != null){
				nextCurrentValue = gsnowItems.get(i+1).getComparedValue();
			}
			else{
				//It is the last item 
				nextCurrentValue = currentValue;
			}
			if(nextCurrentValue < currentValue)
				continue;
			componentsSizeGSnowItem = gsnowItems.get(i).getComponentsSize();
			if(currentValue <= this.significantValue)
				if(maxUtil.evaluate(ListUtils.toDoubleArray(componentsSizeGSnowItem)) < initSizeElements)
					pointerSmallerThanInitSizeElements.add(i);
				else
					pointerBiggerThanInitSizeElements.add(i);
		}
		
		int pointerToMin = 0;
		
		if(!pointerBiggerThanInitSizeElements.isEmpty()){
			//we take the first pick
			pointerToMin = pointerBiggerThanInitSizeElements.get(0);
			currentValue = gsnowItems.get(pointerToMin).getComparedValue();
			componentsSizeGSnowItem = gsnowItems.get(pointerToMin).getComponentsSize();
		}
		else if(pointerBiggerThanInitSizeElements.isEmpty()){
			double min = Double.MAX_VALUE;
			for(int i : pointerSmallerThanInitSizeElements){
				currentValue = gsnowItems.get(i).getComparedValue();
				if(currentValue < min){
					min = currentValue;
					pointerToMin = i;
				}
			}
		}
		if(gsnowItems.containsKey(pointerToMin)){
			significantValue = gsnowItems.get(pointerToMin).getComparedValue();
			significantNodes = gsnowItems.get(pointerToMin).getNodes();
		}
		

		significantItemLocal.setNodes(significantNodes);
		significantItemLocal.setComparedValue(significantValue);
		return significantItemLocal;
	}
	
	
	private void executeSnow() throws IOException, SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException{
		
		List<ProteinVertex> nodes = nodes2ProteinVertices(gsnowItems.get(significantItem.getNodes().size()).getNodes());
		Set<String> intermediates = new HashSet<String>();
		List<List<ProteinVertex>> components = new ArrayList<List<ProteinVertex>>();
		SimpleUndirectedGraph<ProteinVertex, DefaultEdge> subgraph = (SimpleUndirectedGraph<ProteinVertex, DefaultEdge>) Subgraph.randomSubgraph(proteinNetwork.getInteractomeGraph(), nodes);
		if(intermediate){
			System.out.println("Antes intermediario ["+subgraph.getVertices().size()+"]: "+subgraph.getVertices().size());
			intermediates = Subgraph.OneIntermediateList(proteinNetwork.getInteractomeGraph(), subgraph);
			for(ProteinVertex vertex : nodes){
				if(intermediates.contains(vertex.getId())){
					//System.out.println("vertex: "+vertex.getId());
					intermediates.remove(vertex.getId());
				}
			}
			for(String id : intermediates){
				mapNames.put(id, id);
			}
		}
		System.out.println("Despues intermediario ["+subgraph.getVertices().size()+"]: "+subgraph.getVertices().size());
		components = subgraph.getAllInformationComponents(true);
		updateJobStatus("35","List preprocessed");
		//Number of components with more than 1 node
		int compMoreThan1Node = componentsMoreThanOneNode(components);
		result.addOutputItem(new Item("comp_more_than_1_node", Integer.toString(compMoreThan1Node), "Number of components with more than 1 node", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Minimun Connected Network selected. Topology description"));
		int bicomponents = subgraph.getNumberOfBicomponents();
		result.addOutputItem(new Item("bicomponents", Integer.toString(bicomponents), "Number of Bicomponents", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Minimun Connected Network selected. Topology description"));

		File f = new File(outputFileName+"_mcn.sif");
		IOUtils.write(f.getAbsoluteFile(), getSignificantMcnSif(subgraph));
		result.addOutputItem(new Item("mcn.sif", f.getName(), "Minimun Connected Network interactions", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Minimun Connected Network selected"));
		

		f = new File(outputFileName+"_mcn_interactors.txt");
		
		ProteinNetwork subProteinNetwork = new ProteinNetwork(subgraph);
		
		subProteinNetwork.calcTopologicalValues();
		subProteinNetwork.calcTopologicalMeanValues();
		
		IOUtils.write(f.getAbsoluteFile(), this.getMcnInteractors(subProteinNetwork, intermediates));
		result.addOutputItem(new Item("mcn_interactors", f.getName(), "Minimun Connected Network interactors", Item.TYPE.FILE,StringUtils.toList("TABLE,NETWORKMINER_TABLE", ",") ,new HashMap<String,String>(),"Minimun Connected Network selected"));
		if(intermediate && intermediates != null){
			SnowPrinter snowPrinter = new SnowPrinter();
			f = new File(outputFileName+"_external_nodes_added.txt");
			IOUtils.write(f.getAbsoluteFile(), snowPrinter.printNodes(intermediates));
			result.addOutputItem(new Item("external_nodes_list", f.getName(), "External nodes added", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Minimun Connected Network selected"));

		}
		//getGraphViewer(subgraph, intermediates, components);
		File auxFile = new File(outputFileName);
		List<String> names = new ArrayList<String>();
		if(subProteinNetwork != null){
			names.add(auxFile.getName() + "_list1.json");
			createJson(subProteinNetwork, outputFileName+"_list1.dot", components, intermediates, 1, this.mapNames);
			
			auxFile = new File(outputFileName+"_list1.json");
			addOutputSvgViewer(auxFile, 1);
		}
		statsTest(subProteinNetwork, components);
		
	}
	private String getMcnInteractors(ProteinNetwork subProteinNetwork, Set<String> intermediates) throws SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException{
		
		StringBuilder sb = new StringBuilder();
		String tab = "\t";
		sb.append("#input_id").append(tab).append("id").append(tab).append("type").append(tab).append("rank").append(tab).append("bet").append(tab).append("clust").append(tab);
		sb.append("conn").append(tab).append("go").append(tab).append("kegg").append(System.getProperty("line.separator"));
		for(ProteinVertex vertex :subProteinNetwork.getInteractomeGraph().getVertices()){
			Node affectedNode = null;
			for(Node node : this.significantItem.getNodes()){
				if(node.getId().equals(vertex.getId())){
					affectedNode = node;
					break;
				}
			}
			
			if(intermediates.contains(vertex.getId())){
				sb.append("-").append(tab);
				sb.append(vertex.getId()).append(tab);
				sb.append("external").append(tab);
			}
			else{
				sb.append(affectedNode.getOriginalId()).append(tab);
				sb.append(vertex.getId()).append(tab);
				sb.append("list").append(tab);
			}
			if(affectedNode != null){
				sb.append(affectedNode.getValue()).append(tab);
			}
			else{
				sb.append("-").append(tab);
			}
			sb.append(vertex.getRelativeBetweenness()).append(tab);
			sb.append(vertex.getClusteringCoefficient()).append(tab);
			sb.append(subProteinNetwork.getInteractomeGraph().getDegree(vertex)).append(tab);
			
			if(xrefDBMan.getDBConnector().getDbConnection().getDatabase() != null){
				List<String> dbNames = new ArrayList<String>();
				dbNames.add("go");
				dbNames.add("kegg");
				AnnotationDBManager annotationMng = new AnnotationDBManager(dbConnector);
				Map<String, String> keggItems  = annotationMng.getAnnotationTermNames("kegg");
				XRef xrefEns  = xrefDBMan.getByDBName(vertex.getId(), dbNames);
				List<XRefItem> xrefitemListGO = xrefEns.getXrefItems().get("go");
				if(xrefitemListGO != null && !xrefitemListGO.isEmpty()){
					for(XRefItem xrefitem : xrefitemListGO){
						GO go = gof.getByAccesion(xrefitem.getDisplayName());
						sb.append(go.getName()).append(",");
					}
					sb = this.deleteLastCh(sb, ",");
				}
				
				List<XRefItem> xrefitemListKegg = xrefEns.getXrefItems().get("kegg");
				if(xrefitemListKegg != null && !xrefitemListKegg.isEmpty()){
					sb.append(tab);
					for(XRefItem xrefitem : xrefitemListKegg){
						if(keggItems.containsKey(xrefitem.getDisplayName())){
							sb.append(keggItems.get(xrefitem.getDisplayName())).append(",");
						}
					}
					sb = this.deleteLastCh(sb, ",");
				}
				
			}
			sb.append(System.getProperty("line.separator"));
		}
		
		return sb.toString();
	}
	
	
	private void statsTest(ProteinNetwork subProteinNetwork, List<List<ProteinVertex>> components) throws IOException{
		//2nd Analysis
		logger.debug("Starting 2nd Analysis..................");
		List<Double> relBetSubnet1 = subProteinNetwork.getBetRelList();
		List<Double> relBetRandoms = new ArrayList<Double>();

		List<Double> connSubnet1 = subProteinNetwork.getConnList();
		List<Double> connRandoms = new ArrayList<Double>();

		List<Double> clustSubnet1 = subProteinNetwork.getClustList();
		List<Double> clustRandoms = new ArrayList<Double>();
		RandomsSnowManager randomsManager = new RandomsSnowManager(this.proteinNetwork, intermediate, true);

//		RandomsSnowManager randomsManager = new RandomsSnowManager(subProteinNetwork, intermediate, true);
		List<ProteinNetwork> subNetworks = randomsManager.createRandoms(randoms, this.significantItem.getNodes().size());
		double []componentsRandoms = new double[randoms];
		int randonNumber = 0;
		for(ProteinNetwork proteinNetwork : subNetworks){
			int  randomVertex = (int) (Math.random() * proteinNetwork.getInteractomeGraph().getVertices().size());
			ProteinVertex v = proteinNetwork.getInteractomeGraph().getVertices().get(randomVertex);
			relBetRandoms.add(v.getRelativeBetweenness());
			connRandoms.add(Double.parseDouble(Integer.toString(proteinNetwork.getInteractomeGraph().getDegree(v))));
			clustRandoms.add(v.getClusteringCoefficient());
			componentsRandoms[randonNumber] = (double)proteinNetwork.getInteractomeGraph().getAllInformationComponents(true).size();
			randonNumber++;
		}
		
		int[] range = getRange(componentsRandoms);
		String rangeString = "["+Integer.toString(range[0])+", "+Integer.toString(range[1])+"]";
		result.addOutputItem(new Item("components_number", components.size()+" "+rangeString, "Number of components [95% confidence interval]", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Minimun Connected Network selected. Topology description"));

		String toWrite = ksTest(relBetSubnet1, relBetRandoms, connSubnet1, connRandoms, clustSubnet1, clustRandoms/*, side*/);
		if(!toWrite.equals("")){
			String symbol;
				
			KolmogorovSmirnovTestResult resultKs;
			
			resultKs = getPValue(relBetSubnet1, relBetRandoms/*, side*/);
			symbol = getSymbol(resultKs.getSide());
			String stringResult = StringUtils.decimalFormat(resultKs.getPValue(),decimalFormat);
			//if(resultKs.getPValue() == 0)
			//	stringResult = "0.0001";
			result.addOutputItem(new Item("sn_random_kol_param_bet", stringResult , "Relative betweenness: Subnet "+symbol+" Random pval ", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Minimun Connected Network selected. Topology description"));
			createImages(outputFileName+"_sn_random_relBet", relBetSubnet1, "subnet1", relBetRandoms, "randoms", "sn_random_relBet", "Plot", "Minimun Connected Network selected. Topology description","Relative Betweenness");
			
			resultKs = getPValue(connSubnet1, connRandoms);
			symbol = getSymbol(resultKs.getSide());
			stringResult = StringUtils.decimalFormat(resultKs.getPValue(),decimalFormat);
			//if(resultKs.getPValue() == 0)
			//	stringResult = "0.0001";
			result.addOutputItem(new Item("sn_random_kol_param_conn", stringResult, "Connection degree: Subnet "+symbol+" Random pval ", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Minimun Connected Network selected. Topology description"));
			createImages(outputFileName+"_sn_random_conn", connSubnet1, "subnet1", connRandoms, "randoms", "sn_random_conn", "Plot", "Minimun Connected Network selected. Topology description","Connections");
			
			resultKs = getPValue(clustSubnet1, clustRandoms);
			symbol = getSymbol(resultKs.getSide());
			stringResult = StringUtils.decimalFormat(resultKs.getPValue(),decimalFormat);
			//if(resultKs.getPValue() == 0)
			//	stringResult = "0.0001";
			result.addOutputItem(new Item("sn_random_kol_param_clu",  stringResult, "Clustering coefficient: Subnet "+symbol+" Random pval ", Item.TYPE.MESSAGE, new ArrayList<String>(),new HashMap<String,String>(),"Minimun Connected Network selected. Topology description"));
			createImages(outputFileName+"_sn_random_clust", clustSubnet1, "subnet1", clustRandoms, "randoms", "sn_random_clust", "Plot", "Minimun Connected Network selected. Topology description","Clustering coeff");

			File f = new File(outputFileName+"_sn_random_kol.txt");
			IOUtils.write(f.getAbsolutePath(), toWrite);
			result.addOutputItem(new Item("sn_random_kol_param", f.getName(), "File", Item.TYPE.FILE, new ArrayList<String>(),new HashMap<String,String>(),"Minimun Connected Network selected. Topology description"));
		}
		else
			result.addOutputItem(new Item("sn_random_kol_param", "Empty results", "Subnet - Random", Item.TYPE.MESSAGE, new ArrayList<String>(), new HashMap<String,String>(), "Minimun Connected Network selected. Topology description"));
		logger.debug("Finished 2nd Analysis..................");
	}
//	private void getGraphViewer(SimpleUndirectedGraph<ProteinVertex, DefaultEdge> subgraph, Set<String> intermediates, List<List<ProteinVertex>> components) throws IOException, SQLException{
//		
//		Xml xmlObject = new Xml(this.dbConnector);
//		File xmlFile = new File(outdir+"/subnetwork1.xml");
//		xmlObject.graphToXML(xmlFile.getAbsolutePath(), subgraph, intermediates, components, type, this.mapNames);
//		addOutputAppletItem(xmlFile, 1);
//	}
	
	private String getSignificantMcnSif(SimpleUndirectedGraph<ProteinVertex, DefaultEdge> subgraph){
		Sif sif = new Sif();
		StringBuilder sb = new StringBuilder();
		sb.append(sif.graphToSif(subgraph));
		return sb.toString();
	}
	
	public void executeGenerateRandoms(){
		try {
			SnowPrinter snowPrinter = new SnowPrinter();
			int sizeMin = Integer.parseInt(commandLine.getOptionValue("size-min"));	
			int sizeMax = Integer.parseInt(commandLine.getOptionValue("size-max"));
			RandomsSnowManager rndManager = new RandomsSnowManager(proteinNetwork, intermediate, false);
			
			
//			List<List<Double>> allGsnowValues = new ArrayList<List<Double>>();
			List<List<List<ProteinVertex>>> components;
			List<Double> gsnowValues;
			List<ProteinNetwork> list;
			
			double tInicio = 0.0;
			double tFinal = 0.0;
			double tTotal = 0.0;
			IOUtils.append(outputFileName, snowPrinter.printGSnowValuesHeader(randoms));
			for(int i=sizeMin; i <= sizeMax; i++){
				tInicio = System.currentTimeMillis();
				System.out.println("Creating randoms with size: "+i);
				list = rndManager.createRandoms(randoms, i);
				components = rndManager.getComponents(list);
				gsnowValues = rndManager.getGSnowValues(list, components);
				IOUtils.append(outputFileName, i+"\t"+snowPrinter.printGSnowValues(gsnowValues));
				System.gc();
//				allGsnowValues.add(rndManager.getGSnowValues(list, components));
				tFinal = System.currentTimeMillis();
				tTotal = (tFinal-tInicio)/1000;
				System.out.println("\t spend time: "+tTotal);
			}
//			gsnowValues2File(allGsnowValues, sizeMin, randoms);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
//	private void gsnowValues2File(List<List<Double>> allGsnowValues, int sizeMin, int randoms) throws IOException{
//		SnowPrinter snowPrinter = new SnowPrinter();
//		String values = snowPrinter.printAllGSnowValues(allGsnowValues, sizeMin, randoms);
//		IOUtils.write(outputFileName, values);
//	}
	
	private String loadFile(){
		String folder = this.babelomicsHomePath + "/conf/interactions/gsnow/";
		String intermediateString = intermediate ? "intermediate":"nointermediate";
		String localType = type;
		if(type.equals("transcripts")){
			localType = "proteins";
			logger.debug("Transcripts recognised, loading "+localType+" randoms");
		}
		folder += interactome+"_"+localType+"_"+group+"db_"+intermediateString+".txt";
		return folder;
		
	}
	
	public class GSnowItem{

		//List<Node> nodes;
		List<Node> nodes;
		//This is the value(vertices/components) before compared with randoms
		double rawValue;
		//This is the value(mean of significant values) after compared with randoms
		double comparedValue;
		List<Double> componentsSize;
		
		public GSnowItem(){
			nodes = new ArrayList<Node>();
			componentsSize = new ArrayList<Double>();
			
		}
		public List<Node> getNodes() {
			return nodes;
		}
		public void setNodes(List<Node> nodes) {
			this.nodes = nodes;
		}
		public double getRawValue() {
			return rawValue;
		}
		public void setRawValue(double rawValue) {
			this.rawValue = rawValue;
		}
		public double getComparedValue() {
			return comparedValue;
		}
		public void setComparedValue(double comparedValue) {
			this.comparedValue = comparedValue;
		}
		
		public List<Double> getComponentsSize() {
			return componentsSize;
		}
		public void setComponentsSize(List<Double> componentsSize) {
			this.componentsSize = componentsSize;
		}
		@Override
		public String toString(){
			return "{"+nodes.toString()+" - rawValue = "+rawValue+" - comparedValue = "+comparedValue+"}";
		}
		public String getNodesIds(){
			StringBuilder sb = new StringBuilder();
			for(Node n : nodes){
				sb.append(n.getId()).append(",");
			}
			if(sb.toString().equals(""))
				return sb.toString();
			return sb.substring(0,sb.toString().length()-1);
		}
	}
	
	
	
	
	
}

