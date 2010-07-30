package org.bioinfo.babelomics.methods.functional;

import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.bioinfo.babelomics.exception.EmptyAnnotationException;
import org.bioinfo.commons.utils.ArrayUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.data.dataset.FeatureData;
import org.bioinfo.data.list.exception.InvalidIndexException;
import org.bioinfo.infrared.common.dbsql.DBConnector;
import org.bioinfo.infrared.common.feature.FeatureList;
import org.bioinfo.infrared.funcannot.AnnotationItem;
import org.bioinfo.infrared.funcannot.filter.FunctionalFilter;


public class FatiScan extends GeneSetAnalysis {


	public static final int SHORT_FORMAT = 1;
	public static final int LONG_FORMAT = 2;
	public final static int DEFAULT_NUMBER_OF_PARTITIONS = 30;
	
	// input params	
	private int testMode;	
	private int numberOfPartitions;
	private int outputFormat;
		
	// test
	private TwoListFisherTest fisher;
		
	
	// Two list constructor
	public FatiScan(FeatureData rankedList, FunctionalFilter filter, DBConnector dbConnector, int numberOfPartitions, int testMode, int outputFormat, int order) {
		this.rankedList = rankedList;
		this.filter = filter;
		this.dbConnector = dbConnector;		
		this.order = order;
		this.isYourAnnotations = false;
		// fatiscan specific
		this.numberOfPartitions = numberOfPartitions;
		this.testMode = testMode;
		this.outputFormat = outputFormat;
		
		// set analysis type id
		this.method = FATISCAN;
	}
	
	public FatiScan(FeatureData rankedList, FeatureList<AnnotationItem> annotations, int numberOfPartitions, int testMode, int outputFormat, int order) {
		
		this.rankedList = rankedList;
		this.annotations = annotations;		
		this.order = order;
		this.isYourAnnotations = true;
		
		// fatiscan specific
		this.numberOfPartitions = numberOfPartitions;
		this.testMode = testMode;
		this.outputFormat = outputFormat;
		
		// set analysis type id
		this.method = FATISCAN;
	}
	
	@Override
	public void run() throws InvalidIndexException, SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException, EmptyAnnotationException {
				
		// prepare list
		prepare();
		
		int thresholdPosition;
		List<String> list1,list2;
		
		double inc = -(double)(statistic.get(0)-statistic.get(statistic.size()-1))/(numberOfPartitions+1);		
		double acum = statistic.get(0) + inc;
		
		// test each partition	
		for(int i=0; i<numberOfPartitions; i++){
			
			thresholdPosition = getThresholdPosition(acum);
			
			System.err.println(i + ": threshold = " + acum + " (" + thresholdPosition + ") ");
			
			// top
			list1 = idList.subList(0, thresholdPosition);
			
			if(thresholdPosition<(idList.size()-1)){
				
				// bottom
				list2 = idList.subList(thresholdPosition + 1, idList.size()-1);
		
				// run test
				fisher = new TwoListFisherTest();
				fisher.test(list1,list2,annotations,testMode,termSizes);
			
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
		
	// from TwoListFisherTest result to GeneSetAnalysis result
	private List<GeneSetAnalysisTestResult> toGeneSetAnalysisTestResult(List<TwoListFisherTestResult> twoListFisherTest){
		List<GeneSetAnalysisTestResult> result = new ArrayList<GeneSetAnalysisTestResult>(twoListFisherTest.size());
		for(TwoListFisherTestResult test: twoListFisherTest){			
			GeneSetAnalysisTestResult gseaTest = new GeneSetAnalysisTestResult(test);
			result.add(gseaTest);
		}
		return result;
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


	@Override	
	public List<String> resultListToStringList(List<GeneSetAnalysisTestResult> resultList, boolean header){
		List<String> results = new ArrayList<String>();
		if(header) results.add(GeneSetAnalysisTestResult.fatiScanHeader());
		for(int i=0; i<resultList.size(); i++){			
			results.add(resultList.get(i).toFatiScanString());
		}
		return results;
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

}
