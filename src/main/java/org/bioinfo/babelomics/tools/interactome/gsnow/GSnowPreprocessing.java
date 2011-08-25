package org.bioinfo.babelomics.tools.interactome.gsnow;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashSet;
import java.util.List;
import java.util.Set;

import org.bioinfo.infrared.core.XRefDBManager;
import org.bioinfo.infrared.core.feature.DBName;
import org.bioinfo.infrared.core.feature.XRef;
import org.bioinfo.infrared.core.feature.XRef.XRefItem;
import org.bioinfo.networks.protein.ProteinNetwork;

public class GSnowPreprocessing {

	private String type;
	private XRefDBManager xrefDBMan;
	private static int maxSize;// = 200;
	private String order;
	private int numberItems;
	private double cutOff;
	private ProteinNetwork proteinNetwork;
	
	public GSnowPreprocessing(ProteinNetwork proteinNetwork, String type, XRefDBManager xrefDBMan, int maxPreprocessingSize, String order, int numberItems, double cutOff){
		this.type = type;
		this.xrefDBMan = xrefDBMan;
		this.order = order;
		maxSize = maxPreprocessingSize;
		if(numberItems > maxSize)
			this.numberItems = maxSize;
		else
			this.numberItems = numberItems;
		this.cutOff = cutOff;
		this.proteinNetwork = proteinNetwork;
	}
	/** Priority file pFilePath, normal file: filePath **/
	public ListInfo preprocess(String seedFilePath, String filePath) {
		try{
			ListInfo listInfo = new ListInfo();
			listInfo = parseFile("#", seedFilePath, filePath);
			orderNodes(listInfo);

			System.out.println("Correct nodes after parsing file: "+listInfo.getNodes().size());
			
			listInfo = deleteDuplicates(listInfo);
			System.out.println("Matched nodes after deleting repeteated nodes: "+listInfo.getNodes().size());
			System.out.println("Repeated nodes: "+listInfo.getRepeatedNodes().size());
			
			listInfo =  matchingDBAndInteractome(listInfo);
			
			System.out.println("Matched nodes after matching with db: "+listInfo.getNodes().size());
			System.out.println("Not matched nodes: "+listInfo.getNotMatchNodes().size());
			
			if(listInfo.getNotMatchNodes().size() < 20)
				System.out.println("Not matched nodes: "+listInfo.getNotMatchNodes());
			
			listInfo = cutList(listInfo);
			
			if(listInfo.getCuttedNodes().size() > 0){
				System.out.println("Nodes after cutting the list: "+listInfo.getNodes().size());
				System.out.println("List of cutted nodes: "+listInfo.getCuttedNodes().size());
			}
			
			System.out.println("Repeated nodes: "+listInfo.getRepeatedNodes().size());
			if(!Double.isNaN(cutOff))
				listInfo = cutOff(listInfo);
			
			/** Summary **/
			
			listInfo = fillSeedNodes(listInfo);
			return listInfo;
		}catch(Exception e)
		{
			System.err.println("[ERROR] "+e.getMessage());
			e.printStackTrace();
			//System.exit(0);
		}
		return new ListInfo();
		
	}
	
	
	private ListInfo parseFile(String comment, String seedFilePath, String filePath) throws NumberFormatException, IOException{
		ListInfo listInfo = new ListInfo();
		double position = 0;
		BufferedReader br;
		List<Node> nodes = new ArrayList<Node>();
		String line = "";
		String fields[];
		boolean length2 = false;
		boolean length1 = false;
		boolean error = false;
		double minValue = Double.MAX_VALUE;/** These two values tells the range of the filePath file, with this, we can put the pFilePath in the top of the list, depending on the order  **/
		double maxValue = -Double.MAX_VALUE;/** **/
		br = new BufferedReader(new FileReader(filePath));
		
		while((line = br.readLine()) != null) {
			if(!line.startsWith(comment) && !line.trim().equals("")) {
				fields = line.split("\t");
				if(fields.length == 2){
					length2 = true;
					position = Double.parseDouble(fields[1]);
					listInfo.setRanked(true);
					
				}
				else if(fields.length == 1){
					length1 = true;
					position = nodes.size();
					listInfo.setRanked(false);
				}
				else{
					System.err.println("[ERROR] Format problems detected, two many fields in: "+filePath+" line "+line);
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
				if(position > maxValue)
					maxValue = position;
				if(position < minValue)
					minValue = position;
				/** - Adding node **/
				Node node = new Node(fields[0], position);
				nodes.add(node);
			}
		}
		
		/** There is no seed file **/
		if(seedFilePath.equals("")){
			listInfo.setNodes(nodes);
			return listInfo;
		}
		br = new BufferedReader(new FileReader(seedFilePath));
		if(order.equals("ascending")) {
			position = minValue-1;
		}
		if(order.equals("descending")) {
			position = maxValue+1;
		}
		/** Reading seed file **/
		while((line = br.readLine()) != null) {
			if(!line.startsWith(comment) && !line.trim().equals("")) {
				fields = line.split("\t");
				if(fields.length == 1){
					Node node = new Node(fields[0], position);
					nodes.add(node);
					node.setSeed(true);
				}
				else{
					System.err.println("[ERROR] Format problems detected, two many fields in : "+seedFilePath+" line "+ line);
				}
			}
		}
		listInfo.setNodes(nodes);
		return listInfo;
	}
	
	private ListInfo deleteDuplicates(ListInfo listInfo){
		List<Node> rawNodes = listInfo.getNodes();
		List<Node> nonDuplicatedNodes = new ArrayList<Node>();
		Set<String> repeatedNodes = new HashSet<String>();
		Set<String> idSetString = new HashSet<String>();
		/** Aqui es donde pondré si hay duplicados entre la lista de asociados y la lista de rankeados: Como se ordenan nodos en principio no hay que 
		 ** tocar este código **/
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
	/** Here we discover the type of input id, if it is uniprot, ensembl_gene, hgnc_symbol 
	 **/
	private DBName getOriginalDbName(ListInfo listInfo) throws SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException {
		for(Node node : listInfo.getNodes()){
			List<DBName> dbnames = xrefDBMan.getAllDBNamesById(node.getOriginalId());
			if(dbnames.size() == 1){
				System.out.println("DbName: "+dbnames.get(0));
				return dbnames.get(0);
			}
				
		}
		return null;
	}
	private ListInfo matchingDBAndInteractome(ListInfo listInfo) throws SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException {
		
		String dbName = "";
		if(type.equalsIgnoreCase("proteins") || type.equalsIgnoreCase("transcripts"))
			dbName = "uniprot_swissprot_accession";
		else if(type.equalsIgnoreCase("genes") )
			dbName = "ensembl_gene";
//		System.out.println("DBNAME:"+dbName);
		List<Node> nodes = listInfo.getNodes();
		List<Node> curatedListNodes = new ArrayList<Node>();
		Set<String> notMatchNodes = new HashSet<String>();
		
		/** idsAlreadyMatched: Esta es una variable que nos dice si ya se ha recuperado un id de un originalId, entonces no se inserta en la bbdd
		 * x ejemplo si BRCA2 y BRCA1 dan ensemble_gene: ENSG00001, solo se introduce ENSG00001 una vez, correspondiente a BRCA2, que estaba antes **/
		Set<String> idsAlreadyMatched = new HashSet<String>();
		ListInfo listInfoLocal = listInfo;
		listInfoLocal.setOriginalDbName(getOriginalDbName(listInfo));
		
		Node n = null;
		for(Node node : nodes){
			try {
				//XRef xrefEns = null;
				//List<XRefItem> xrefitemList = new ArrayList<XRefItem>();//xrefDBMan.getDBConnector().getDbConnection().getDatabase == homo_sapiens
				List<XRef> xrefList = null;
				if(xrefDBMan.getDBConnector().getDbConnection().getDatabase() != null){
					List<String> list = new ArrayList<String>();
					list.add(node.getId());
					xrefList = xrefDBMan.getAllIdentifiersByIds(list);
					
					//xrefEns  = xrefDBMan.getByDBName(node.getId(), dbName);
					//xrefitemList = xrefEns.getXrefItems().get(dbName);
				}
				if(xrefList != null && !xrefList.isEmpty()){

					boolean matched = false;
					for(XRef xref : xrefList){
						if(xref.getXrefItems().get(dbName).size()>1)
							System.out.println("ori id:"+node.getOriginalId()+" size: "+xref.getXrefItems().get(dbName).size());
						String nodeId = xref.getXrefItems().get(dbName).get(0).getDisplayName();
						
						if(!idsAlreadyMatched.contains(nodeId)){
							matched = true;
							idsAlreadyMatched.add(nodeId);
						}
						else{
							//System.out.println("Nodes con id ya introducido "+node.getId()+": "+nodeId);
							continue;
						}
						n = new Node(nodeId, node.getValue(), node.getId());
						n.setSeed(node.isSeed());
						if(proteinNetwork.getInteractomeGraph().getVertex(nodeId) != null)
							curatedListNodes.add(n);
						else
							notMatchNodes.add(node.getId());
					}
					if(!matched)
						notMatchNodes.add(node.getId());
				}
				else if(proteinNetwork.getInteractomeGraph().getVertex(node.getId()) != null)
					curatedListNodes.add(node);
				else
					notMatchNodes.add(node.getId());
				
			} catch (Exception e) {
				System.out.println("Problems matching "+node+": "+e.getMessage());
				//e.printStackTrace();
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
			if(order.equalsIgnoreCase("ascending")){
				if(node.getValue() > cutOff)
					localNodes.remove(node);
			}
			if(order.equalsIgnoreCase("descending")){
				if(node.getValue() < cutOff)
					localNodes.remove(node);
			}
		}
		listInfo.setNodes(localNodes);
		return listInfo;
	}
	private ListInfo fillSeedNodes(ListInfo listInfo) {
		Set<String> seedNodes = new HashSet<String>();
		for(Node node : listInfo.getNodes()){
			if(node.isSeed())
				seedNodes.add(node.getId());
		}
		listInfo.setSeedNodes(seedNodes);
		return listInfo;
	}
	public class ListInfo{
		
		private List<Node> nodes;
		private Set<String> notMatchNodes;
		private Set<String> repeatedNodes;
		private Set<String> cuttedNodes;
		private Set<String> seedNodes;
		private boolean ranked;
		
		/** Es la base de datos con la que coincide los inputs id, por ejemplo: si la lista de entrada es BRCA2... entonces: originalDbName = 'hgnc_symbol'**/
		private DBName originalDbName;
		
		public ListInfo(){
			nodes = new ArrayList<Node>();
			notMatchNodes = new HashSet<String>();
			repeatedNodes = new HashSet<String>();
			cuttedNodes = new HashSet<String>();
			seedNodes = new HashSet<String>();
			ranked = false;
			originalDbName = null;

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
		public Set<String> getSeedNodes() {
			return seedNodes;
		}
		public void setSeedNodes(Set<String> seedNodes) {
			this.seedNodes = seedNodes;
		}
		public boolean isRanked() {
			return ranked;
		}
		public void setRanked(boolean ranked) {
			this.ranked = ranked;
		}
		public DBName getOriginalDbName() {
			return originalDbName;
		}
		public void setOriginalDbName(DBName originalDbName) {
			this.originalDbName = originalDbName;
		}
		
		
	}
	
	public class Node implements Comparable<Node> {
		
		/** this id is retrieved by the database **/
		private String id;
		
		/** this id is the one coming from the user **/
		private String originalId;
		
		/** this is the position/p-value of the node in the input list **/
		private double value;
		
		/** This attribute says if it is an seed node or not **/
		boolean seed;
		

		public Node(String id, double value){
			this(id, value, id);
		}
		public Node(String id, double value, String originalId){
			this.id = id;
			this.value = value;
			this.originalId = originalId;
			this.seed = false;
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
		public boolean isSeed() {
			return seed;
		}
		public void setSeed(boolean seed) {
			this.seed = seed;
		}
		@Override
		public String toString(){
			return id+"\t"+originalId+"\t"+value;
		}
		
		@Override
		public int compareTo(Node o) {
			if(value > o.value){
				if(order.equals("ascending"))
					return 1;
				else 
					return 0;
			}
			else{
				if(order.equals("ascending"))
					return 0;
				else
					return 1;
			}
		}
	}

}

