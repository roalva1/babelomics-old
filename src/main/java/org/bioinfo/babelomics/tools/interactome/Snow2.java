package org.bioinfo.babelomics.tools.interactome;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
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
import org.jfree.chart.plot.PlotOrientation;

public class Snow2  extends BabelomicsTool{

//	private List<String> nodeList1, nodeList2;
	private File inputFile;
	private ProteinNetwork proteinNetwork;
	private ProteinNetwork subProteinNetwork1, subProteinNetwork2;
	private List<ProteinNetwork> subProteinNetworkRandoms;
	private String imagesFile;
	private boolean noNumberOfComponents;
	private boolean bicomponents;
//	private int subnet;

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
		
		options.addOption(OptionFactory.createOption("o-dot-file", "It will create an output .dot file", false, true));
		options.addOption(OptionFactory.createOption("o-json-file", "It will create an output .json file", false, true));
		options.addOption(OptionFactory.createOption("o-svg-file", "It will create an output .svg and .dot file", false, true));
		
		options.addOption(OptionFactory.createOption("o-sif-topo-file", "Create a full topological file from a sif file", false));
		options.addOption(OptionFactory.createOption("o-components-file", "If there is this argument, it will create an output .cmp file for the information of each component", false, true));
		options.addOption(OptionFactory.createOption("o-means-file", "Name of the output file with the mean and standard desviation of the topological values (without the extension)", false, true));
		options.addOption(OptionFactory.createOption("o-topo-file", "Name of the .topo output file (without the extension)", false, true));
		options.addOption(OptionFactory.createOption("side", "side for kolmogorov and wilkoxon test. Can be two.sided(by default), less or greater", false, true));
		options.addOption(OptionFactory.createOption("o-stats-file", "Print the statistics file", false, true));
		options.addOption(OptionFactory.createOption("o-images-file", "Print the images for the statistics", false, true));
		
	}
	@SuppressWarnings("unused")
	@Override
	protected void execute() {
		
		subProteinNetworkRandoms = new ArrayList<ProteinNetwork>();
	//	subnet = 0;
		
		noNumberOfComponents = commandLine.hasOption("no-number-components");
		bicomponents = commandLine.hasOption("bicomponents");
	
		imagesFile = "";
		String side = "two.sided";
		if(commandLine.hasOption("side"))
			side = commandLine.getOptionValue("side");

		ProteinNetworkToFile file = new ProteinNetworkToFile();

		try {
		
			if(commandLine.hasOption("sif-file") ){
				
				inputFile = new File(commandLine.getOptionValue("sif-file"));
				logger.debug("checking if inputFile exist...");
				FileUtils.checkFile(inputFile);
				SimpleUndirectedGraph<ProteinVertex, DefaultEdge> interactomeGraph = InteractomeParser.parseFromSifFile(commandLine.getOptionValue("sif-file") );
				proteinNetwork = new ProteinNetwork(interactomeGraph);
			}
			if(commandLine.hasOption("topo-file") ){
				FileUtils.checkFile(new File(commandLine.getOptionValue("topo-file")));
				proteinNetwork.loadTopologicalValues(commandLine.getOptionValue("topo-file"));
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
			if(commandLine.hasOption("o-stats-file") && commandLine.hasOption("topo-file")){
				if(commandLine.hasOption("o-images-file"))
						imagesFile = commandLine.getOptionValue("o-images-file");
				if(!commandLine.hasOption("node-file2") && commandLine.hasOption("randoms"))
					statsOneListAnalysis(side);
				else if(commandLine.hasOption("node-file1") && commandLine.hasOption("node-file2"))
					statsTwoListsAnalisys(side);
				else
					System.err.println("Not correct arguments for statistic test");
			}
			
			if(commandLine.hasOption("o-json-file")){
				IOUtils.write(commandLine.getOptionValue("o-json-file")+".json", proteinNetwork.getInteractomeGraph().toJson());
			}
			if(commandLine.hasOption("o-dot-file") && !commandLine.hasOption("o-svg-file")){
				String dot = proteinNetwork.getInteractomeGraph().toDot();
				IOUtils.write(commandLine.getOptionValue("o-dot-file")+".dot", proteinNetwork.getInteractomeGraph().toDot());
			}	
			else if(!commandLine.hasOption("o-dot-file") && commandLine.hasOption("o-svg-file")){
				createSVGFile(commandLine.getOptionValue("o-svg-file")+".dot");
			}
			else if(commandLine.hasOption("o-dot-file") && commandLine.hasOption("o-svg-file")){
				createSVGFile(commandLine.getOptionValue("o-dot-file")+".dot");
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
		
		IOUtils.write(commandLine.getOptionValue("o-stats-file")+"_list_inter_kol.txt", toWrite);
		wTest.test(commandLine.getOptionValue("o-stats-file")+"_list_inter_wil.txt", side, "list1", "inter", relBetList1, relBetInter, connList1, connInter, clustList1, clustInter);
			
		if(!imagesFile.equals("")){
			createImages(imagesFile+"_list_inter_relBet", relBetList1, "list1", relBetInter, "inter");
			createImages(imagesFile+"_list_inter_conn", connList1, "list1", connInter, "inter");
			createImages(imagesFile+"_list_inter_clust", clustList1, "list1", clustInter, "inter");
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
		IOUtils.write(commandLine.getOptionValue("o-stats-file")+"_sn_random_kol.txt", toWrite);		
		wTest.test(commandLine.getOptionValue("o-stats-file")+"_sn_random_wil.txt", side, "sn1", "rnd", relBetSubnet1, relBetRandoms, connSubnet1, connRandoms, clustSubnet1, clustRandoms);
		if(!imagesFile.equals("")){
			createImages(imagesFile+"_sn_random_relBet", relBetSubnet1, "subnet1", relBetRandoms, "randoms");
			createImages(imagesFile+"_sn_random_conn", connSubnet1, "subnet1", connRandoms, "randoms");
			createImages(imagesFile+"_sn_random_clust", clustSubnet1, "subnet1", clustRandoms, "randoms");
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
		IOUtils.write(commandLine.getOptionValue("o-stats-file")+"_list1_list2_kol.txt", toWrite);
		wTest.test(commandLine.getOptionValue("o-stats-file")+"_list1_list2_wil.txt", side, "list1", "list2", relBetList1, relBetList2, connList1, connList2, clustList1, clustList2);
		if(!imagesFile.equals("")){
			createImages(imagesFile+"_list1_list2_relBet", relBetList1, "list1", relBetList2, "list2");
			createImages(imagesFile+"_list1_list2_conn", connList1, "list1", connList2, "list2");
			createImages(imagesFile+"_list1_list2_clust", clustList1, "list1", clustList2, "list2");
		}
		
		//2nd Analysis
		List<Double> relBetSubnet1 = subProteinNetwork1.getBetRelList();
		List<Double> relBetSubnet2 = subProteinNetwork2.getBetRelList();
		
		List<Double> connSubnet1 = subProteinNetwork1.getConnList();
		List<Double> connSubnet2 = subProteinNetwork2.getConnList();
		
		List<Double> clustSubnet1 = subProteinNetwork1.getClustList();
		List<Double> clustSubnet2 = subProteinNetwork2.getClustList();
		
		toWrite = formatStatsFile(relBetSubnet1, relBetSubnet2, connSubnet1, connSubnet2, clustSubnet1, clustSubnet2, side);
		IOUtils.write(commandLine.getOptionValue("o-stats-file")+"_sn1_sn2_kol.txt", toWrite);		
		wTest.test(commandLine.getOptionValue("o-stats-file")+"_sn1_sn2_wil.txt", side, "sn1", "sn2", relBetSubnet1, relBetSubnet2, connSubnet1, connSubnet2, clustSubnet1, clustSubnet2);
		if(!imagesFile.equals("")){
			createImages(imagesFile+"_sn1_sn2_relBet", relBetSubnet1, "subnet1", relBetSubnet2, "subnet2");
			createImages(imagesFile+"_sn1_sn2_conn", connSubnet1, "subnet1", connSubnet2, "subnet2");
			createImages(imagesFile+"_sn1_sn2_clust", clustSubnet1, "subnet1", clustSubnet2, "subnet2");
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
			
			ProteinNetwork subProteinNetwork = createSubnet(subgraph/*, subnet*/);
			if(commandLine.hasOption("o-topo-file")) {
				sbTopo.append(getTopologicalValues(subProteinNetwork, i)).append(System.getProperty("line.separator"));
			}
			if(commandLine.hasOption("o-means-file")) {
				sbMeans.append(getTopologicalMeanValues(subProteinNetwork, i)).append(System.getProperty("line.separator"));
			}
			if(commandLine.hasOption("o-components-file")) {
				sbComponents.append(getComponentsValues(subProteinNetwork, i));
			}
			IOUtils.write("random1", subProteinNetwork.getSifToString());
			subProteinNetworkRandoms.add(subProteinNetwork);
		}
		randomFilesToString(sbMeans, commandLine.getOptionValue("o-means-file")+"_sn_1-"+(randoms)+"_means.txt");
		randomFilesToString(sbTopo,commandLine.getOptionValue("o-topo-file")+"_sn_1-"+(randoms)+"_topo.txt");
		randomFilesToString(sbComponents,commandLine.getOptionValue("o-components-file")+"_sn_1-"+(randoms)+"_comp.txt");
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
		ProteinNetwork subProteinNetwork = createSubnet(subgraph/*, subnet++*/);
		
		if(commandLine.hasOption("o-topo-file")) {
			sbTopo.append(getTopologicalValues(subProteinNetwork, node));
			IOUtils.write(commandLine.getOptionValue("o-topo-file")+"_sn_nodeFile"+node+"_topo.txt", sbTopo.toString());
		}
		
		if(commandLine.hasOption("o-means-file")) {
			sbMeans.append(getTopologicalMeanValues(subProteinNetwork, node));
			IOUtils.write(commandLine.getOptionValue("o-means-file")+"_sn_nodeFile"+node+"_means.txt", sbMeans.toString());
		}
		
		if(commandLine.hasOption("o-components-file")) {
			sbComponents.append(getComponentsValues(subProteinNetwork, node));
			sbComponents.deleteCharAt(sbComponents.lastIndexOf(System.getProperty("line.separator")));
			IOUtils.write(commandLine.getOptionValue("o-components-file")+"_sn_nodeFile"+node+"_comp.txt", sbComponents.toString());
		}
		return subProteinNetwork;
	}
	private ProteinNetwork createSubnet(SimpleUndirectedGraph<ProteinVertex, DefaultEdge> subgraph/*, int nSubnet*/) {
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
