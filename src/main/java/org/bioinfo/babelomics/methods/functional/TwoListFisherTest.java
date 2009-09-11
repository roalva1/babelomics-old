package org.bioinfo.babelomics.methods.functional;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.bioinfo.infrared.common.feature.FeatureList;
import org.bioinfo.infrared.funcannot.AnnotationItem;
import org.bioinfo.math.data.IntegerMatrix;
import org.bioinfo.math.result.FisherTestResult;
import org.bioinfo.math.result.TestResultList;
import org.bioinfo.math.stats.MultipleTestCorrection;
import org.bioinfo.math.stats.inference.FisherExactTest;

public class TwoListFisherTest extends FunctionalTest {
	
	public final static double DEFAULT_PVALUE_THRESHOLD = 0.05;
	
	private List<TwoListFisherTestResult> results;

	@Override
	public void test(List<String> list1, List<String> list2, FeatureList<AnnotationItem> annotations, int testMode) {
		
		// init counters
		Map<String, Integer> map1 = new HashMap<String, Integer>();
		Map<String, Integer> map2 = new HashMap<String, Integer>();
		Map<String, List<String>> list1Positives = new HashMap<String, List<String>>();
		Map<String, List<String>> list2Positives = new HashMap<String, List<String>>();
		List<String> terms = new ArrayList<String>();
		String term,id;
				
		// count annotations for every gene
		for(AnnotationItem annotation: annotations) {
			term = annotation.getFunctionalTermId();
			id = annotation.getId();
			
			// backup term name
			if(!terms.contains(term)) terms.add(term);
			
			// list 1 count
			  // init entry
			if(!map1.containsKey(term)) {
				map1.put(term, 0);
				list1Positives.put(term, new ArrayList<String>());
			}
			  // add term
			if(list1.contains(id)) {
				map1.put(term, map1.get(term) + 1);
				list1Positives.get(term).add(id);
			}
			
			// list 2 count
			  // init entry
			if(!map2.containsKey(term)) {
				map2.put(term, 0);
				list2Positives.put(term, new ArrayList<String>());
			}
			  // add term
			if(list2.contains(id)) {
				map2.put(term, map2.get(term) + 1);
				list2Positives.get(term).add(id);
			}
		}

		// count summary		
		if(terms!= null && terms.size()>0){			
			IntegerMatrix fisherCounts = new IntegerMatrix(terms.size(), 4);
			for(int i=0 ; i<terms.size() ; i++) {				
				fisherCounts.set(i, 0, map1.get(terms.get(i)));
				fisherCounts.set(i, 1, map2.get(terms.get(i)));
				fisherCounts.set(i, 2, list1.size()-map1.get(terms.get(i)));
				fisherCounts.set(i, 3, list2.size()-map2.get(terms.get(i)));				
			}			
			TestResultList<FisherTestResult> testResult = new FisherExactTest().fisherTest(fisherCounts, testMode);
						
			// p-value adjustment
			MultipleTestCorrection.BHCorrection(testResult);
			
			results = new ArrayList<TwoListFisherTestResult>(testResult.size());
			for(int i=0; i<testResult.size(); i++){
				results.add(new TwoListFisherTestResult(terms.get(i),fisherCounts.get(i,0),fisherCounts.get(i,2),fisherCounts.get(i,1),fisherCounts.get(i,3),list1Positives.get(terms.get(i)),list2Positives.get(terms.get(i)),testResult.get(i).getPValue(),testResult.get(i).getAdjPValue()));
			}
			
		} else {
			//FIXME thrown an exception
			System.out.println("\nannotation is nullllll\n");
			results = null;
		}
			
	}

	public List<TwoListFisherTestResult> getSignificantResults(){
		return getSignificantResults(DEFAULT_PVALUE_THRESHOLD);
	}
	
	public List<TwoListFisherTestResult> getSignificantResults(double threshold){
		List<TwoListFisherTestResult> significant = new ArrayList<TwoListFisherTestResult>();
		for(TwoListFisherTestResult result: this.results){
			if(result.getAdjPValue()<threshold) significant.add(result);
		}
		return significant;
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
