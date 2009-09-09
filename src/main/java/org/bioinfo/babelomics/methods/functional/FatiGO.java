package org.bioinfo.babelomics.methods.functional;

import java.sql.SQLException;
import java.util.ArrayList;
import java.util.List;

import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.infrared.common.dbsql.DBConnector;
import org.bioinfo.infrared.common.feature.FeatureList;
import org.bioinfo.infrared.funcannot.AnnotationItem;
import org.bioinfo.infrared.funcannot.dbsql.AnnotationDBManager;
import org.bioinfo.infrared.funcannot.filter.Filter;
import org.bioinfo.infrared.funcannot.filter.GOFilter;
import org.bioinfo.infrared.funcannot.filter.KeggFilter;
import org.bioinfo.math.exception.InvalidParameterException;
import org.bioinfo.math.result.FisherTestResult;
import org.bioinfo.math.result.TestResultList;
import org.bioinfo.math.stats.inference.FisherExactTest;

public class FatiGO {

	public static final int REMOVE_NEVER = 0;
	public static final int REMOVE_EACH = 1;	
	public static final int REMOVE_REF = 2;
	public static final int REMOVE_GENOME = 2;
	public static final int REMOVE_ALL = 3;
	
	// input params
	private List<String> list1;
	private List<String> list2;
	private Filter filter;
	DBConnector dbConnector;
	private int testMode;
	private int duplicatesMode;
	

	// results	
	List<TwoListFisherTestResult> results;
	FeatureList<AnnotationItem> annotations;

	// Two list constructor
	public FatiGO(List<String> list1, List<String> list2, Filter filter, DBConnector dbConnector, int testMode, int duplicatesMode ) {
		this.list1 = list1;
		this.list2 = list2;
		this.filter = filter;
		this.dbConnector = dbConnector;
		this.testMode = testMode;
		this.duplicatesMode = duplicatesMode;
	}
	
	// One list against Genome constructor
	public FatiGO(List<String> list1, Filter filter, DBConnector dbConnector) {
		this.list1 = list1;
		this.list2 = InfraredUtils.getGenome(dbConnector);
		this.filter = filter;
		this.dbConnector = dbConnector;
		this.testMode = FisherExactTest.GREATER;
		this.duplicatesMode = REMOVE_GENOME;
	}
		
	public void run() throws SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException, InvalidParameterException {

		
		List<String> all = new ArrayList<String>(list1);
		all.addAll(list2);	
		
		// duplicates managing
		removeDuplicates(list1, list2, duplicatesMode);
		
		// annotation
		annotations = InfraredUtils.getAnnotations(dbConnector, all, filter);

		// run test
		TwoListFisherTest fisher = new TwoListFisherTest();
		fisher.test(list1, list2, annotations, testMode);
		results = fisher.getResults();

				
	}

	
	private void removeDuplicates(List<String> list1, List<String> list2, int duplicatesMode){
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
	public List<TwoListFisherTestResult> getResults() {
		return results;
	}


	/**
	 * @param results the results to set
	 */
	public void setResults(List<TwoListFisherTestResult> results) {
		this.results = results;
	}
	
}
