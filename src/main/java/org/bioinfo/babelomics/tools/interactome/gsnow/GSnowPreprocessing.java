package org.bioinfo.babelomics.tools.interactome.gsnow;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.bioinfo.infrared.core.XRefDBManager;
import org.bioinfo.infrared.core.feature.XRef;
import org.bioinfo.infrared.core.feature.XRef.XRefItem;
import org.bioinfo.networks.protein.ProteinNetwork;

public class GSnowPreprocessing {

//	private ListInfo listInfo;
	private String type;
	private XRefDBManager xrefDBMan;
//	private int minPreprocessingSize;
//	private int maxPreprocessingSize;
	private static int maxSize;// = 200;
	private String order;
//	private boolean wholeList;
	private int numberItems;
	private double cutOff;
	private ProteinNetwork proteinNetwork;
	
	public GSnowPreprocessing(ProteinNetwork proteinNetwork, String type, XRefDBManager xrefDBMan/*, int minPreprocessingSize*/, int maxPreprocessingSize, String order, /*boolean wholeList,*/ int numberItems, double cutOff){
		this.type = type;
		this.xrefDBMan = xrefDBMan;
		//this.minPreprocessingSize = minPreprocessingSize;
		/*this.maxPreprocessingSize = maxPreprocessingSize;*/
		this.order = order;
		maxSize = maxPreprocessingSize;
		if(numberItems > maxSize)
			this.numberItems = maxSize;
		else
			this.numberItems = numberItems;
		this.cutOff = cutOff;
		this.proteinNetwork = proteinNetwork;
	}
	
	public ListInfo preprocess(String filePath) {
		try{
			ListInfo listInfo = new ListInfo();
			List<Node> nodes = parseFile("#", filePath);
			listInfo.setNodes(nodes);
			orderNodes(listInfo);
//			System.out.println(listInfo.getNodes().toString());
			
//			if(listInfo.getNodes().size() < this.minPreprocessingSize)
//				throw new Exception("[ERROR] List size too small, minimun size permitted "+this.minPreprocessingSize);
//			if(listInfo.getNodes().size() > maxSize){
//				System.err.println("[WARNING] List size too long, maximun size permitted "+maxSize+". GSnow will keep running with the "+maxSize+" firsts elements.");
//			}
			
			System.out.println("Matched nodes after parsing file: "+listInfo.getNodes().size());
			listInfo = deleteDuplicates(listInfo);
			System.out.println("Matched nodes after deleting repeteated nodes: "+listInfo.getNodes().size());
			System.out.println("Repeated nodes: "+listInfo.getRepeatedNodes().size());
			listInfo =  matchingDBAndInteractome(listInfo);
			System.out.println("Matched nodes after matching with db: "+listInfo.getNodes().size());
			System.out.println("Not matched nodes: "+listInfo.getNotMatchNodes().size());
			if(listInfo.getNotMatchNodes().size() < 100)
				System.out.println("Not matched nodes: "+listInfo.getNotMatchNodes());
			listInfo = cutList(listInfo);
			if(listInfo.getCuttedNodes().size() > 0){
				System.out.println("Nodes after cutting the list: "+listInfo.getNodes().size());
				System.out.println("List of cutted nodes: "+listInfo.getCuttedNodes().size());
			}
			System.out.println("Repeated nodes: "+listInfo.getRepeatedNodes().size());
			if(!Double.isNaN(cutOff))
				listInfo = cutOff(listInfo);
			return listInfo;
		}catch(Exception e)
		{
			System.err.println("[ERROR] "+e.getMessage());
			e.printStackTrace();
			//System.exit(0);
		}
		return new ListInfo();
		
	}
	
	

	private List<Node> parseFile(String comment, String nodeFile) throws NumberFormatException, IOException{
		BufferedReader br = new BufferedReader(new FileReader(nodeFile));
		List<Node> nodes = new ArrayList<Node>();
		String line = "";
		String fields[];
		boolean length2 = false, length1 = false, error = false;
		double position = 0;
		while((line = br.readLine()) != null) {
			if(!line.startsWith(comment) && !line.trim().equals("")) {
				fields = line.split("\t");
				if(fields.length == 2){
					length2 = true;
					position = Double.parseDouble(fields[1]);
				}
				else if(fields.length == 1){
					length1 = true;
					position = nodes.size();
				}
				else{
					System.err.println("[ERROR] Format problems detected, two many fields in: "+line);
					error = true;
				}
				if(length1 == true && length2 == true){
					System.err.println("[ERROR] Format problems detected: you have some lines with one column and some others with two columns, for instance: "+line);
					error = true;
					length1 = false;
					length2 = false;
				}
				if(error){
					error = false;
					continue;
				}
				// - Adding node
				Node node = new Node(fields[0], position);
				nodes.add(node);
			}
		}
		return nodes;
	}
	
	private ListInfo deleteDuplicates(ListInfo listInfo){
		List<Node> rawNodes = listInfo.getNodes();
		List<Node> nonDuplicatedNodes = new ArrayList<Node>();
		Set<String> repeatedNodes = new HashSet<String>();
		Set<String> idSetString = new HashSet<String>();
		for(Node node : rawNodes){
			if(!idSetString.contains(node.id)){
				nonDuplicatedNodes.add(node);
				idSetString.add(node.id);
			}
			else
				repeatedNodes.add(node.id);
		}
		listInfo.setNodes(nonDuplicatedNodes);
		listInfo.setRepeatedNodes(repeatedNodes);
		return listInfo;
	}
	
	private ListInfo matchingDBAndInteractome(ListInfo listInfo) {
		String dbName = "";
		if(type.equalsIgnoreCase("proteins") || type.equalsIgnoreCase("transcripts"))
			dbName = "uniprot_swissprot_accession";
		else if(type.equalsIgnoreCase("genes") )
			dbName = "ensembl_gene";
//		System.out.println("DBNAME:"+dbName);
		List<Node> nodes = listInfo.getNodes();
		List<Node> curatedListNodes = new ArrayList<Node>();
		Set<String> notMatchNodes = new HashSet<String>();
		ListInfo listInfoLocal = listInfo;
		Node n = null;
		for(Node node : nodes){
			try {
				XRef xrefEns = null;
				List<XRefItem> xrefitemList = new ArrayList<XRefItem>();//xrefDBMan.getDBConnector().getDbConnection().getDatabase == homo_sapiens
				if(xrefDBMan.getDBConnector().getDbConnection().getDatabase() != null){
					xrefEns  = xrefDBMan.getByDBName(node.getId(), dbName);
					xrefitemList = xrefEns.getXrefItems().get(dbName);
				}
				if(xrefEns != null && !xrefitemList.isEmpty()){
					for(XRefItem xrefitem : xrefitemList){
						String nodeId = xrefitem.getDisplayName();
						n = new Node(nodeId, node.getValue(), node.getId());
						if(proteinNetwork.getInteractomeGraph().getVertex(nodeId) != null)
							curatedListNodes.add(n);
						else
							notMatchNodes.add(node.getId());
					}
				}
				else if(proteinNetwork.getInteractomeGraph().getVertex(node.getId()) != null)
					curatedListNodes.add(node);
				else
					notMatchNodes.add(node.getId());
				
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
		listInfoLocal.setNodes(curatedListNodes);
		listInfoLocal.setNotMatchNodes(notMatchNodes);
		return listInfoLocal;
	}
	private void orderNodes(ListInfo listInfo) {
		Collections.sort(listInfo.getNodes());
	}
	
	private ListInfo cutList(ListInfo listInfo){
		List<Node> nodes = new ArrayList<Node>();
		Set<String> cuttedNodes = new HashSet<String>();
		if(listInfo.getNodes().size() < this.numberItems)
			return listInfo;
		for(int i=0; i < this.numberItems; i++){
			nodes.add(listInfo.getNodes().get(i));
		}
		for(int i=this.numberItems; i < listInfo.getNodes().size(); i++){
//			if(listInfo.getNodes().get(i).getId().equalsIgnoreCase("P03989"))
//				System.out.println(listInfo.getNodes().get(i).getId());
//			if(cuttedNodes.contains(listInfo.getNodes().get(i).getId()))
//				System.out.println(listInfo.getNodes().get(i).getId());
			cuttedNodes.add(listInfo.getNodes().get(i).getId());
		}	
		listInfo.setCuttedNodes(cuttedNodes);
		listInfo.setNodes(nodes);
		
		return listInfo;
	}
	
	private ListInfo cutOff(ListInfo listInfo) {
		List<Node> nodes = listInfo.getNodes();
		List<Node> localNodes = new ArrayList<Node>();
		localNodes.addAll(nodes);
		for(Node node : nodes){
			if(order.equalsIgnoreCase("ascendant")){
				if(node.getValue() > cutOff)
					localNodes.remove(node);
			}
			if(order.equalsIgnoreCase("descendant")){
				if(node.getValue() < cutOff)
					localNodes.remove(node);
			}
		}
		listInfo.setNodes(localNodes);
		return listInfo;
	}
	public class ListInfo{
		
		private List<Node> nodes;
		private Set<String> notMatchNodes;
		private Set<String> repeatedNodes;
		private Set<String> cuttedNodes;

		//this is the original list size after the entire preprocessing, we don't count intermediates
//		private int processedNodesSize;
		
		public ListInfo(){
			nodes = new ArrayList<Node>();
			notMatchNodes = new HashSet<String>();
			repeatedNodes = new HashSet<String>();
			cuttedNodes = new HashSet<String>();

		}
		public void setNodes(List<Node> nodes) {
			this.nodes = nodes;
		}
		public List<Node> getNodes() {
			return nodes;
		}
		public Set<String> getNotMatchNodes() {
			return notMatchNodes;
		}
		public Set<String> getRepeatedNodes() {
			return repeatedNodes;
		}
		public void setNotMatchNodes(Set<String> notMatchNodes) {
			this.notMatchNodes = notMatchNodes;
		}
		public void setRepeatedNodes(Set<String> repeatedNodes) {
			this.repeatedNodes = repeatedNodes;
		}
		public Set<String> getCuttedNodes() {
			return cuttedNodes;
		}
		public void setCuttedNodes(Set<String> cuttedNodes) {
			this.cuttedNodes = cuttedNodes;
		}
//		public int getProcessedNodesSize() {
//			return processedNodesSize;
//		}
//		public void setProcessedNodesSize(int processedNodesSize) {
//			this.processedNodesSize = processedNodesSize;
//		}
		
		
	}
	
	public class Node implements Comparable<Node> {
		
		//this id is retrieved by the database
		private String id;
		//this id is the one coming from the user
		private String originalId;
		//this is the position/p-value of the node in the input list
		private double value;

//		private int connections;
//		private double cluster;
//		private double bet;
//		
//		private List<String> go;
//		private List<String> kegg;
		
		public Node(String id, double value){
			this(id, value, id);
		}
		public Node(String id, double value, String originalId){
			this.id = id;
			this.value = value;
			this.originalId = originalId;
//			go = new ArrayList<String>();
//			kegg = new ArrayList<String>();
		}
		
		public String getId() {
			return id;
		}

		public double getValue() {
			return value;
		}
		
		public String getOriginalId() {
			return originalId;
		}
		public void setOriginalId(String originalId) {
			this.originalId = originalId;
		}
		
		@Override
		public String toString(){
			return id+"\t"+originalId+"\t"+value;
		}
		
		@Override
		public int compareTo(Node o) {
			if(value > o.value){
				if(order.equals("ascendant"))
					return 1;
				else 
					return 0;
			}
			else{
				if(order.equals("ascendant"))
					return 0;
				else
					return 1;
			}
		}
	}

}

