package org.bioinfo.babelomics.methods.functional;

import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;

import org.bioinfo.babelomics.exception.EmptyAnnotationException;
import org.bioinfo.commons.utils.ArrayUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.data.dataset.FeatureData;
import org.bioinfo.data.list.exception.InvalidIndexException;
import org.bioinfo.infrared.common.dbsql.DBConnector;
import org.bioinfo.infrared.common.feature.FeatureList;
import org.bioinfo.infrared.funcannot.AnnotationItem;
import org.bioinfo.infrared.funcannot.filter.FunctionalFilter;


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
	private FunctionalFilter filter;
	private DBConnector dbConnector;
	private int testMode;	
	private int numberOfPartitions;
	private int outputFormat;
	private int order;
	private boolean isYourAnnotations;
	
	
	// test
	TwoListFisherTest fisher;
	
	// results
	List<GeneSetAnalysisTestResult> results;	
	FeatureList<AnnotationItem> annotations;

	// Two list constructor
	public FatiScan(FeatureData rankedList, FunctionalFilter filter, DBConnector dbConnector, int numberOfPartitions, int testMode, int outputFormat, int order) {
		this.rankedList = rankedList;
		this.filter = filter;
		this.dbConnector = dbConnector;		
		this.numberOfPartitions = numberOfPartitions;
		this.testMode = testMode;
		this.outputFormat = outputFormat;
		this.order = order;
		this.isYourAnnotations = false;
	}
	
	public FatiScan(FeatureData rankedList, FeatureList<AnnotationItem> annotations, int numberOfPartitions, int testMode, int outputFormat, int order) {
		this.rankedList = rankedList;
		this.annotations = annotations;
		this.numberOfPartitions = numberOfPartitions;
		this.testMode = testMode;
		this.outputFormat = outputFormat;
		this.order = order;
		this.isYourAnnotations = true;
	}
	
	
	public void prepare() throws InvalidIndexException{
		
		// id list
		idList = rankedList.getDataFrame().getRowNames();//.getColumn(0);
		// statistic		
		statistic = ArrayUtils.toList(rankedList.getDataFrame().getColumnAsDoubleArray(0));
		// order ranked list		
		int[] sortIndex = ListUtils.order(statistic);
		statistic = ListUtils.ordered(statistic,sortIndex);		
		idList = ListUtils.ordered(idList,sortIndex);
		//if(order==ASCENDING_SORT){
		if(order==DESCENDING_SORT){
			Collections.reverse(idList);
			Collections.reverse(statistic);
		}		
		
	}
	
	public void run() throws InvalidIndexException, SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException, EmptyAnnotationException {
				
		// prepare list
		prepare();
		
		// annotation		
		if(!isYourAnnotations) annotations = InfraredUtils.getAnnotations(dbConnector, idList, filter);
		if(annotations==null || annotations.size()==0) throw new EmptyAnnotationException();
		System.err.println("annotations:" + annotations.size());
		results = new ArrayList<GeneSetAnalysisTestResult>();
		
		double inc = -(double)(statistic.get(0)-statistic.get(statistic.size()-1))/(numberOfPartitions+1);
		System.err.println("inc: " + inc);
		double acum = statistic.get(0) + inc;
		
		int thresholdPosition;
		List<String> list1,list2;
				
		for(int i=0; i<numberOfPartitions; i++){
			
			thresholdPosition = getThresholdPosition(acum);
			
			System.err.println(i + ": threshold = " + acum + " (" + thresholdPosition + ") ");
			
			list1 = idList.subList(0, thresholdPosition);
			if(thresholdPosition<(idList.size()-1)){
				list2 = idList.subList(thresholdPosition + 1, idList.size()-1);
				
				//System.err.println("l1.size: " + list1.size() + "l2.size: " + list2.size());
				
				// run test
				fisher = new TwoListFisherTest();
				fisher.test(list1, list2, annotations, testMode);
			
				// get result
				results.addAll(toGeneSetAnalysisTestResult(fisher.getResults()));
			}
						
			acum+=inc;
			
		}

		if(outputFormat == SHORT_FORMAT) {
			HashMap<String,GeneSetAnalysisTestResult> resultsMap = new HashMap<String,GeneSetAnalysisTestResult>();
			// unique term
			for(GeneSetAnalysisTestResult testResult: results){
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
		
	private List<GeneSetAnalysisTestResult> toGeneSetAnalysisTestResult(List<TwoListFisherTestResult> raw){
		List<GeneSetAnalysisTestResult> result = new ArrayList<GeneSetAnalysisTestResult>(raw.size());
		for(TwoListFisherTestResult test: raw){			
			GeneSetAnalysisTestResult gseaTest = new GeneSetAnalysisTestResult(test);
			result.add(gseaTest);
		}
		return result;
	}
	
	
	private int getThresholdPosition(double acum){
		int position = 0;
		//System.err.println("acum: " + acum + " statistic.size: " + statistic.size());
		for(int i=0; i<statistic.size(); i++){
			//System.err.println("position: " + i + " acum: " + acum + " value: " + statistic.get(i) + " order: " + order);
			if( (order==ASCENDING_SORT && statistic.get(i)>=acum) || (order==DESCENDING_SORT && statistic.get(i)<acum) ) {				
				position = i;
				break;
			}			
		}
		return position;
	}

	public List<GeneSetAnalysisTestResult> getSignificant(){
		return getSignificant(TwoListFisherTest.DEFAULT_PVALUE_THRESHOLD);
	}
	
	public List<GeneSetAnalysisTestResult> getSignificant(double threshold){
		List<GeneSetAnalysisTestResult> significant = new ArrayList<GeneSetAnalysisTestResult>();
		for(GeneSetAnalysisTestResult result: this.results){			
			if(result.getAdjPValue()<threshold) significant.add(result);
		}
		return significant;
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
	public List<GeneSetAnalysisTestResult> getResults() {
		return results;
	}

	/**
	 * @param results the results to set
	 */
	public void setResults(List<GeneSetAnalysisTestResult> results) {
		this.results = results;
	}

}
