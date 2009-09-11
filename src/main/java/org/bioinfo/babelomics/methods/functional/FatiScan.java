package org.bioinfo.babelomics.methods.functional;

import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.bioinfo.collections.exceptions.InvalidColumnIndexException;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.data.dataset.FeatureData;
import org.bioinfo.infrared.common.dbsql.DBConnector;
import org.bioinfo.infrared.common.feature.FeatureList;
import org.bioinfo.infrared.funcannot.AnnotationItem;
import org.bioinfo.infrared.funcannot.filter.Filter;

public class FatiScan {

	public static final int ASCENDING_SORT = 1;
	public static final int DESCENDING_SORT = 2;
	public static final int SHORT_FORMAT = 1;
	public static final int LONG_FORMAT = 2;
	public final static int DEFAULT_NUMBER_OF_PARTITIONS = 30;
	
	// input params
	private List<String> idList;	
	private List<Double> statistic;
	private FeatureData rankedList;
	private Filter filter;
	private DBConnector dbConnector;
	private int testMode;	
	private int numberOfPartitions;
	private int outputFormat;
	private int order;
	
	
	// test
	TwoListFisherTest fisher;
	
	// results
	List<TwoListFisherTestResult> results;	
	FeatureList<AnnotationItem> annotations;

	// Two list constructor
	public FatiScan(FeatureData rankedList, Filter filter, DBConnector dbConnector, int numberOfPartitions, int testMode, int outputFormat, int order) {
		this.rankedList = rankedList;
		this.filter = filter;
		this.dbConnector = dbConnector;		
		this.numberOfPartitions = numberOfPartitions;
		this.testMode = testMode;
		this.outputFormat = outputFormat;
		this.order = order;
	}
	
	public void prepare() throws InvalidColumnIndexException{
		
		// id list
		idList = rankedList.getDataFrame().getColumn(0);
		// statistic
		statistic = ListUtils.toList(rankedList.getDataFrame().getColumnAsDoubleArray(1));
		
		// order ranked list
		int[] sortIndex = ListUtils.order(statistic,order==DESCENDING_SORT);
		ListUtils.ordered(idList,sortIndex);
		ListUtils.ordered(statistic,sortIndex);
		
	}
	
	public void run() throws InvalidColumnIndexException, SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException {
				
		// prepare list
		prepare();
		
		// annotation		
		annotations = InfraredUtils.getAnnotations(dbConnector, idList, filter);		

		results = new ArrayList<TwoListFisherTestResult>();
		
		double inc = (double)(statistic.get(0)-statistic.get(statistic.size()-1))/(numberOfPartitions+1);
		double acum;
		if(order==DESCENDING_SORT) {
			inc*=-1;
			acum = statistic.get(0) + inc;
		} else {
			acum = statistic.get(statistic.size()-1) + inc;
		}
		
		int thresholdPosition;
				
		for(int i=0; i<numberOfPartitions; i++){
			
			thresholdPosition = getThresholdPosition(acum);
			
			System.err.print(i + ": threshold = " + acum + " (" + thresholdPosition + ")");
			
			List<String> list1 = idList.subList(0, thresholdPosition);
			List<String> list2 = idList.subList(thresholdPosition + 1, idList.size()-1);
			
			System.err.println("l1.size: " + list1.size() + "l2.size: " + list2.size());
			
			// run test
			fisher = new TwoListFisherTest();
			fisher.test(list1, list2, annotations, testMode);
		
			// get result
			results.addAll(fisher.getResults());
						
			acum+=inc;
		}

		if(outputFormat == SHORT_FORMAT) {
			HashMap<String,TwoListFisherTestResult> resultsMap = new HashMap<String,TwoListFisherTestResult>();
			// unique term
			for(TwoListFisherTestResult testResult: results){
				if(resultsMap.containsKey(testResult.getTerm())){
					if(resultsMap.get(testResult.getTerm()).getAdjPValue()>testResult.getAdjPValue()){
						resultsMap.remove(testResult.getTerm());
						resultsMap.put(testResult.getTerm(), testResult);
					}
				} else {
					resultsMap.put(testResult.getTerm(), testResult);
				}
			}
			// update results
			results.clear();
			results.addAll(resultsMap.values());
			
		}
				
		System.err.println("final results.size: " + results.size());	
		
	}
		
	
	
	private int getThresholdPosition(double acum){
		int position = 0;
		for(int i=0; i<statistic.size(); i++){
			if( (order==ASCENDING_SORT && statistic.get(i)>=acum) || (order==DESCENDING_SORT && statistic.get(i)<acum) ) {
				position = i;
				break;
			}			
		}
		return position;
	}

	public List<TwoListFisherTestResult> getSignificant(){
		return getSignificant(TwoListFisherTest.DEFAULT_PVALUE_THRESHOLD);
	}
	
	public List<TwoListFisherTestResult> getSignificant(double threshold){
		if(fisher!=null) return fisher.getSignificantResults(threshold);
		return null;
	}
	
	/**
	 * @return the idList
	 */
	public List<String> getIdList() {
		return idList;
	}

	/**
	 * @param idList the idList to set
	 */
	public void setIdList(List<String> idList) {
		this.idList = idList;
	}

	/**
	 * @return the statistic
	 */
	public List<Double> getStatistic() {
		return statistic;
	}

	/**
	 * @param statistic the statistic to set
	 */
	public void setStatistic(List<Double> statistic) {
		this.statistic = statistic;
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

}
