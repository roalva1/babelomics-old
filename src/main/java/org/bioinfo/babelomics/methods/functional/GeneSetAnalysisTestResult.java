package org.bioinfo.babelomics.methods.functional;

import java.util.List;

import org.bioinfo.commons.utils.ListUtils;

public class GeneSetAnalysisTestResult {

	// common
	private String term;	
	private List<String> list1Ids;
	private List<String> list2Ids;	
	private double adjPValue;
	// fatiscan
	private int list1Positives;
	private int list1Negatives;
	private int list2Positives;
	private int list2Negatives;
	private double pValue;
	// logistic
	private int size;
	private boolean converged;
	private double logRatio;
	
	/**
	 * @param term
	 * @param list1Positives
	 * @param list1Negatives
	 * @param list2Positives
	 * @param list2Negatives
	 * @param list1Ids
	 * @param list2Ids
	 * @param value
	 * @param adjPValue
	 */
	public GeneSetAnalysisTestResult(String term, int list1Positives, int list1Negatives, int list2Positives, int list2Negatives, List<String> list1Ids, List<String> list2Ids, double pValue, double adjPValue) {
		this.term = term;
		this.list1Positives = list1Positives;
		this.list1Negatives = list1Negatives;
		this.list2Positives = list2Positives;
		this.list2Negatives = list2Negatives;
		this.list1Ids = list1Ids;
		this.list2Ids = list2Ids;
		this.pValue = pValue;
		this.adjPValue = adjPValue;
	}

	public GeneSetAnalysisTestResult(TwoListFisherTestResult test) {		
		this.term = test.getTerm();
		this.list1Positives = test.getList1Positives();
		this.list1Negatives = test.getList1Negatives();
		this.list2Positives = test.getList2Positives();
		this.list2Negatives = test.getList2Negatives();
		this.list1Ids = test.getList1Ids();
		this.list2Ids = test.getList2Ids();
		this.pValue = test.getPValue();
		this.adjPValue = test.getAdjPValue();
	}

	public GeneSetAnalysisTestResult(String term, List<String> list1Ids, int size, boolean converged, double logRatio, double adjPValue) {
		this.term = term;
		this.list1Ids = list1Ids;		
		this.size = size;
		this.converged = converged;
		this.logRatio = logRatio;
		this.adjPValue = adjPValue;
	}
	
		
	public static String fatiScanHeader(){
		StringBuilder out = new StringBuilder();
		out.append("#term").append("\t");
		out.append("list1_positives").append("\t");
		out.append("list1_negatives").append("\t");
		out.append("list2_positives").append("\t");
		out.append("list2_negatives").append("\t");
		out.append("list1_positive_ids").append("\t");
		out.append("list2_positive_ids").append("\t");		
		out.append("pvalue").append("\t");
		out.append("adj_pvalue");
		return out.toString();
	}
	
	public String toFatiScanString(){
		StringBuilder out = new StringBuilder();
		out.append(this.term).append("\t");
		out.append(this.list1Positives).append("\t");
		out.append(this.list1Negatives).append("\t");
		out.append(this.list2Positives).append("\t");
		out.append(this.list2Negatives).append("\t");
		out.append(ListUtils.toString(this.list1Ids,",")).append("\t");
		out.append(ListUtils.toString(this.list2Ids,",")).append("\t");		
		out.append(this.pValue).append("\t");
		out.append(this.adjPValue);
		return out.toString();
	}
	
	public static String LogisticHeader(){
		StringBuilder out = new StringBuilder();
		out.append("#term").append("\t");		
		out.append("list1_positive_ids").append("\t");		
		out.append("size").append("\t");
		out.append("converged").append("\t");
		out.append("lor").append("\t");
		out.append("adj_pvalue");
		return out.toString();
	}

	public String toLogisticString(){
		StringBuilder out = new StringBuilder();
		out.append(this.term).append("\t");		
		out.append(ListUtils.toString(this.list1Ids,",")).append("\t");	
		out.append(this.size).append("\t");
		out.append(this.converged).append("\t");
		out.append(this.logRatio).append("\t");		
		out.append(this.adjPValue);
		return out.toString();
	}
	
	/**
	 * @return the term
	 */
	public String getTerm() {
		return term;
	}


	/**
	 * @param term the term to set
	 */
	public void setTerm(String term) {
		this.term = term;
	}


	/**
	 * @return the list1Positives
	 */
	public int getList1Positives() {
		return list1Positives;
	}


	/**
	 * @param list1Positives the list1Positives to set
	 */
	public void setList1Positives(int list1Positives) {
		this.list1Positives = list1Positives;
	}


	/**
	 * @return the list1Negatives
	 */
	public int getList1Negatives() {
		return list1Negatives;
	}


	/**
	 * @param list1Negatives the list1Negatives to set
	 */
	public void setList1Negatives(int list1Negatives) {
		this.list1Negatives = list1Negatives;
	}


	/**
	 * @return the list2Positives
	 */
	public int getList2Positives() {
		return list2Positives;
	}


	/**
	 * @param list2Positives the list2Positives to set
	 */
	public void setList2Positives(int list2Positives) {
		this.list2Positives = list2Positives;
	}


	/**
	 * @return the list2Negatives
	 */
	public int getList2Negatives() {
		return list2Negatives;
	}


	/**
	 * @param list2Negatives the list2Negatives to set
	 */
	public void setList2Negatives(int list2Negatives) {
		this.list2Negatives = list2Negatives;
	}


	/**
	 * @return the list1Ids
	 */
	public List<String> getList1Ids() {
		return list1Ids;
	}


	/**
	 * @param list1Ids the list1Ids to set
	 */
	public void setList1Ids(List<String> list1Ids) {
		this.list1Ids = list1Ids;
	}


	/**
	 * @return the list2Ids
	 */
	public List<String> getList2Ids() {
		return list2Ids;
	}


	/**
	 * @param list2Ids the list2Ids to set
	 */
	public void setList2Ids(List<String> list2Ids) {
		this.list2Ids = list2Ids;
	}


	/**
	 * @return the pValue
	 */
	public double getPValue() {
		return pValue;
	}


	/**
	 * @param value the pValue to set
	 */
	public void setPValue(double value) {
		pValue = value;
	}


	/**
	 * @return the adjPValue
	 */
	public double getAdjPValue() {
		return adjPValue;
	}


	/**
	 * @param adjPValue the adjPValue to set
	 */
	public void setAdjPValue(double adjPValue) {
		this.adjPValue = adjPValue;
	}
	
}

