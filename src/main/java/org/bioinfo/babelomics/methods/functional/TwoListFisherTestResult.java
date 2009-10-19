package org.bioinfo.babelomics.methods.functional;

import java.util.List;

import org.bioinfo.commons.utils.ListUtils;

public class TwoListFisherTestResult {

	
	private String term;
	private int list1Positives;
	private int list1Negatives;
	private int list2Positives;
	private int list2Negatives;
	private List<String> list1Ids;
	private List<String> list2Ids;
	private double pValue;
	private double adjPValue;
	
	
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
	public TwoListFisherTestResult(String term, int list1Positives, int list1Negatives, int list2Positives, int list2Negatives, List<String> list1Ids, List<String> list2Ids, double pValue, double adjPValue) {
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

	public String toString(){
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

	public static String header(){
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
