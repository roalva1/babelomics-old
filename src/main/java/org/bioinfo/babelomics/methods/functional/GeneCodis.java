package org.bioinfo.babelomics.methods.functional;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Hashtable;
import java.util.Iterator;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;

import org.bioinfo.commons.exec.Command;
import org.bioinfo.commons.exec.SingleProcess;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.infrared.common.dbsql.DBConnector;
import org.bioinfo.infrared.common.feature.FeatureList;
import org.bioinfo.infrared.funcannot.AnnotationItem;
import org.bioinfo.infrared.funcannot.filter.FunctionalFilter;
import org.bioinfo.math.exception.InvalidParameterException;
import org.bioinfo.math.stats.inference.FisherExactTest;

public class GeneCodis {

	
	public static final int REMOVE_NEVER = 0;
	public static final int REMOVE_EACH = 1;	
	public static final int REMOVE_REF = 2;
	public static final int REMOVE_GENOME = 3;
	public static final int REMOVE_ALL = 4;
	
	
	public enum correctionFactor{fdr, permutation, none};
	public enum testFactor{hypergeometric,chiSquare,both};
	public enum analysisFactor{concurrence,singular};

	
	// input params
	private List<String> list1;
	private List<String> list2;
	private FunctionalFilter filter;
	private DBConnector dbConnector;
	private int testMode;
	private int duplicatesMode;
	private boolean isYourAnnotations;
	private int RFactor; //my list
	private int rFactor;//referenced list
	private String binPath;
	// test
	
	// results	
	List<String> results;	
	private FeatureList<AnnotationItem> annotations;
	private String outdir;
	private String name;
	private int support = 3;
	private analysisFactor analysis;
	private int supportRandom = 3;
	private correctionFactor correction;
	private testFactor test;
	
	

//two list
	
	public GeneCodis(String binPath,String outdir,String name, List<String> list1, List<String> list2, FunctionalFilter filter, DBConnector dbConnector,int duplicatesMode, int support,int supportRandom, correctionFactor correction, testFactor test, analysisFactor analysis) {		
		this.filter = filter;
		this.list1 = list1;
		this.list2 = list2;
		this.dbConnector = dbConnector;
		this.duplicatesMode = duplicatesMode;
		this.isYourAnnotations = false;
		this.binPath = binPath;
		this.outdir= outdir;
		this.name=name;
		this.support= support;
		this.supportRandom = supportRandom;
		this.correction = correction;
		this.test=test;
		this.analysis=analysis;
		this.setRFactor(list1.size());
		this.setrFactor(list2.size());
	}

	//one list agains genome
	public GeneCodis(String binPath,String outdir,String name, List<String> list1, FunctionalFilter filter, DBConnector dbConnector, int duplicatesMode, int support, int supportRandom, correctionFactor correction, testFactor test, analysisFactor analysis) {
		this.list1 = list1;
		this.list2 = InfraredUtils.getGenome(dbConnector);
		this.filter = filter;
		this.dbConnector = dbConnector;
		this.duplicatesMode = REMOVE_GENOME;
		this.isYourAnnotations = false;
		this.binPath = binPath;
		this.outdir= outdir;
		this.name=name;
		this.support= support;
		this.supportRandom = supportRandom;
		this.correction = correction;
		this.test=test;
		this.analysis=analysis;
		this.setRFactor(list1.size());
		this.setrFactor(list2.size());
		System.err.println("list2---------------"+list2.size());
	}
	
	// Your anntoations two list constructor
	public GeneCodis(List<String> list1, FeatureList<AnnotationItem> annotations) {
		this.list1 = list1;
		this.list2 = InfraredUtils.getGenome(dbConnector);
		this.annotations = annotations;
		this.duplicatesMode = REMOVE_GENOME;
		this.isYourAnnotations = true;
	}
	

	public void run() throws SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException, InvalidParameterException {
		
		//commander final list to be filled with genes and annot in compact mode
		System.err.println("00");
		List<String> allItems = new ArrayList<String>(list1.size()+list2.size());
		
		List<String> all = new ArrayList<String>(list1.size()+list2.size());
		all.addAll(this.list1);
		all.addAll(this.list2);
		// duplicates managing
		removeDuplicates();
		//fill my hash list idgenes elements
		Hashtable<String, Boolean> myHash  = new Hashtable<String, Boolean>();		
		for (int i=0; i< list1.size(); i++){
			myHash.put(list1.get(i), true);
		}
		
		// fill annotation expanded
		if(!isYourAnnotations) {
			this.annotations = InfraredUtils.getAnnotations(this.dbConnector, all, this.filter);
//			System.err.println("this.filter--------------"+this.filter.getMaxNumberGenes());
//			System.err.println("all.size()--------------"+all.size());
		}
//		System.err.println("this.dbConnector.toString()-------------"+this.dbConnector.toString());
//		System.err.println("this.annotations.size()---------------"+this.annotations.size());
		//create annot hash
		Map<String, List<String>> myHashAnnotation  = new LinkedHashMap<String, List<String>>();
		
		//fill annotation Hash : <id of gene>, <annotation compact>
		for (int i=0; i< this.annotations.size(); i++){
			if(!myHashAnnotation.containsKey(annotations.get(i).getId())) {
				myHashAnnotation.put(annotations.get(i).getId(), new ArrayList<String>());
			}
			myHashAnnotation.get(annotations.get(i).getId()).add(annotations.get(i).getFunctionalTermId());
		}

		//fill commander final list
		String isRef;
		Iterator<String> e = myHashAnnotation.keySet().iterator();
			String key;
			  while (e.hasNext()) {
			     key = (String)e.next();
			     //if id exist in my first hash list, fill 3er column with id of gene
			     isRef = "0";
			     if (myHash.get(key) != null){
			    	 isRef = (String) key;
			    	 }
			     allItems.add((String) key + "\t" + ListUtils.toString(myHashAnnotation.get(key), ",") + "\t\t" + isRef);
			  }
		setResults(allItems);
//		System.err.println("allItems.size()---------"+allItems.size());
		
		doSingleProcess(this.binPath ,this.outdir+ "/"+name+"_WellFormedInput", outdir + "/"+ name+ ".txt");
		
	}
	
	public void doSingleProcess(String binPath, String inputAbsolutPath, String outputAbsolutPath){		
		
		try {
			IOUtils.write(this.outdir+ "/"+name+"_WellFormedInput", getResults());
		} catch (IOException e1) {
			e1.printStackTrace();
		}
		
		int ANALISIS = 1;
		int CORRECTION = 0;
		int TEST= 1;
		
		//ANALISIS
		switch(this.analysis){
		case concurrence:  ANALISIS= 1;break;
		case singular : ANALISIS = 2;break;
		default:
			System.out.println("ANALISIS format not valid");
		}
		
		//CORRECTION
		switch(this.correction){
		case fdr:CORRECTION= -1;break;
		case permutation: CORRECTION= 1;break;
		case none: CORRECTION =0;break;
		default:
			System.out.println("CORRECTION format not valid");
		}
		
		//TEST
		switch(this.test){
		case chiSquare:TEST = 1;break;
		case hypergeometric: TEST = 0;break;
		case both: TEST = 2;break;
		default:
			System.out.println("TEST format not valid");
		}
	
		String cmdStr = binPath +" "+" "+inputAbsolutPath +" "+ this.support+" "+ " -a" + ANALISIS + " -i" + this.supportRandom +" -r"+this.getrFactor()+ " -R"+this.getRFactor() + " -s"+ CORRECTION + " -t" + TEST+ " -o "+outputAbsolutPath;
		
		Command cmd = new Command(cmdStr); 
		SingleProcess sp = new SingleProcess(cmd);
		sp.runAsync();
		sp.waitFor();

		///end_execution
		System.err.println();
		System.err.println("Output result:");
		System.err.println(sp.getRunnableProcess().getOutput());		
		System.err.println("Error:");
		System.err.println(sp.getRunnableProcess().getError());
		System.err.println("Exception:");
		System.err.println(sp.getRunnableProcess().getException());
		System.err.println("=====================> END OF LOCAL genecodis\n");
	}
	
	public void removeDuplicates(){
		// each list
		if(duplicatesMode!=REMOVE_NEVER){
			list1 = ListUtils.unique(list1);
			list2 = ListUtils.unique(list2);
		}
		//complementary
		if(duplicatesMode==REMOVE_REF){
			for (String id:list1) {
				if(list2.contains(id)){
					list2.remove(id);
				}
			}
		}
		//genome
		if(duplicatesMode==REMOVE_GENOME){
			list1 = ListUtils.unique(list1);
			List<String> ensemblList1 = InfraredUtils.toEnsemblId(dbConnector, list1);
			for (String id:ensemblList1) {				
				if(list2.contains(id)){
					list2.remove(id);
				}
			}
		}
		// all
		if(duplicatesMode==REMOVE_ALL){
			list1 = ListUtils.unique(list1);
			list2 = ListUtils.unique(list2);
			for (String id:list1) {
				if(list2.contains(id)) {
					list1.remove(id);
					list2.remove(id);
				}	
			}
		}
	}
	
	
	/**
	 * @return the annotations
	 */
	public FeatureList<AnnotationItem> getAnnotations() {
		return annotations;
	}


	/**
	 * @param annotations the annotations to set
	 */
	public void setAnnotations(FeatureList<AnnotationItem> annotations) {
		this.annotations = annotations;
	}


	/**
	 * @return the results
	 */
	public List<String> getResults() {
		return results;
	}


	/**
	 * @param results the results to set
	 */
	public void setResults(List<String> results) {
		this.results = results;
	}

	/**
	 * @return the list1
	 */
	public List<String> getList1() {
		return list1;
	}

	/**
	 * @param list1 the list1 to set
	 */
	public void setList1(List<String> list1) {
		this.list1 = list1;
	}

	/**
	 * @return the list2
	 */
	public List<String> getList2() {
		return list2;
	}

	/**
	 * @param list2 the list2 to set
	 */
	public void setList2(List<String> list2) {
		this.list2 = list2;
	}

	public void setRFactor(int rFactor) {
		RFactor = rFactor;
	}

	public int getRFactor() {
		return RFactor;
	}

	public void setrFactor(int rFactor) {
		this.rFactor = rFactor;
	}

	public int getrFactor() {
		return rFactor;
	}
	
}
