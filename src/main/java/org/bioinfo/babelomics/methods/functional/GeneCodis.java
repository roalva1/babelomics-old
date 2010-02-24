package org.bioinfo.babelomics.methods.functional;

import java.lang.reflect.Array;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Enumeration;
import java.util.Hashtable;
import java.util.List;

import org.bioinfo.commons.utils.ArrayUtils;
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
	
	// test
	TwoListFisherTest fisher;
	
	// results	
	List<String> results;	
	FeatureList<AnnotationItem> annotations;

	// Two list constructor
	public GeneCodis(List<String> list1, List<String> list2, FunctionalFilter filter, DBConnector dbConnector, int testMode, int duplicatesMode ) {
		this.list1 = list1;
		this.list2 = list2;
		this.filter = filter;
		this.dbConnector = dbConnector;
		this.testMode = testMode;
		this.duplicatesMode = duplicatesMode;
		this.isYourAnnotations = false;
		System.err.println("filter:::::::::::::::::" + filter.getMaxNumberGenes());
		this.setRFactor(list1.size());
		this.setrFactor(list2.size());
	}

	// One list against Genome constructor
	public GeneCodis(List<String> list1, FunctionalFilter filter, DBConnector dbConnector) {
		this.list1 = list1;
		this.list2 = InfraredUtils.getGenome(dbConnector);
		this.filter = filter;
		this.dbConnector = dbConnector;
		this.testMode = FisherExactTest.GREATER;
		this.duplicatesMode = REMOVE_GENOME;
		this.isYourAnnotations = false;
	}

	// Your annotations two list constructor
	public GeneCodis(List<String> list1, List<String> list2, FeatureList<AnnotationItem> annotations, int testMode, int duplicatesMode ) {
		this.list1 = list1;
		this.list2 = list2;
		this.annotations = annotations;
		this.testMode = testMode;
		this.duplicatesMode = duplicatesMode;
		this.isYourAnnotations = true;
	}
	
	// Your anntoations two list constructor
	public GeneCodis(List<String> list1, FeatureList<AnnotationItem> annotations) {
		this.list1 = list1;
		this.list2 = InfraredUtils.getGenome(dbConnector);
		this.annotations = annotations;
		this.testMode = FisherExactTest.GREATER;
		this.duplicatesMode = REMOVE_GENOME;
		this.isYourAnnotations = true;
	}
	
	
	
	public void run() throws SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException, InvalidParameterException {
		
		//commander final list to be filled with genes and annot in compact mode
		List<String> allItems = new ArrayList<String>(list1.size()+list2.size());
		System.err.println("1111111");
		//
		List<String> all = new ArrayList<String>(list1.size()+list2.size());
		all.addAll(list1);
		all.addAll(list2);
		
		// duplicates managing
		removeDuplicates();
		
		//fill my hash list idgenes elements
		Hashtable<String, Boolean> myHash  = new Hashtable<String, Boolean>();		
		for (int i=0; i< list1.size(); i++){
			myHash.put(list1.get(i), true);
		}
		System.err.println("filter---------"+filter.size());
		// fill annotation expanded
		if(!isYourAnnotations) annotations = InfraredUtils.getAnnotations(dbConnector, all, filter);
		System.err.println("3333333"+filter.size());
		
		//create annot hash
		Hashtable<String, String> myHashAnnotation  = new Hashtable<String, String>();
		
		//fill annotation Hash : <id of gene>, <annotation compact>
		for (int i=0; i< annotations.size(); i++){
			if (myHashAnnotation.get(annotations.get(i).getId()) == null){
				//create key value pair
				System.err.println(annotations.get(i).getId().toString());
				myHashAnnotation.put(annotations.get(i).getId(), annotations.get(i).getFunctionalTermId());
			}else {
				//concat value of created key
				myHashAnnotation.put(annotations.get(i).getId(),myHashAnnotation.get(i).concat("," + annotations.get(i).getFunctionalTermId()));
			}
		}
		System.err.println("2222222");
		//fill commander final list
		String isRef = "0";
		Enumeration e = myHashAnnotation.keys();
		
			Object obj;
			  while (e.hasMoreElements()) {
				  
			     obj = e.nextElement();
			     System.out.println("clave "+ obj +": " + myHashAnnotation.get(obj));
			     //if id exist in my first hash list, fill 3er column with id of gene
			     if (myHash.get(obj) != null){
			    	 isRef = (String) obj;
			    	 }
			     allItems.add((String) obj + "/t" + myHashAnnotation.get(obj) + "/t/t" + isRef);
			     isRef = "0";
			  }
			 
		setResults(allItems);
		
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
