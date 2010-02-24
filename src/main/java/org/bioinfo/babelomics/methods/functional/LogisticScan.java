package org.bioinfo.babelomics.methods.functional;

import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

import org.bioinfo.babelomics.utils.RCommand;
import org.bioinfo.commons.io.utils.FileUtils;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ArrayUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.dataset.FeatureData;
import org.bioinfo.data.list.exception.InvalidIndexException;
import org.bioinfo.infrared.common.dbsql.DBConnector;
import org.bioinfo.infrared.common.feature.FeatureList;
import org.bioinfo.infrared.funcannot.AnnotationItem;
import org.bioinfo.infrared.funcannot.filter.Filter;
import org.bioinfo.infrared.funcannot.filter.FunctionalFilter;

public class LogisticScan {

	public static final String ENVIRONMENT_VARIABLE = "BABELOMICS_HOME";
	public static final String DEFAULT_BINARY = "/bin/logistic_univariate.r";
	public static final String TMP_PATH = "/tmp";
	
	private String scriptPath;
	
	public static final int ASCENDING_SORT = 1;
	public static final int DESCENDING_SORT = 2;
	
	// input params
	private List<String> idList;	
	private List<Double> statistic;
	private FeatureData rankedList;
	private FunctionalFilter filter;
	private DBConnector dbConnector;	
	private int order;
	private boolean isYourAnnotations;

	// files
	private String annotationFileName;
	private String rankedListFileName;
	private String outputFileName;
	
	// results
	FeatureList<AnnotationItem> annotations;
	List<GeneSetAnalysisTestResult> results;

	
	// Infrared annotation constructor
	public LogisticScan(FeatureData rankedList, FunctionalFilter filter, DBConnector dbConnector, int order) {
			
		// params
		this.rankedList = rankedList;
		this.filter = filter;
		this.dbConnector = dbConnector;		
		this.order = order;
		this.isYourAnnotations = false;

		init();
	
	}
	
	// Your annotations constructor
	public LogisticScan(FeatureData rankedList, FeatureList<AnnotationItem> annotations, int order) {
			
		// params
		this.rankedList = rankedList;
		this.annotations = annotations;
		this.order = order;
		this.isYourAnnotations = true;
		
		init();
	
	}
	
	
	// init engine params
	private void init(){
		
		// script
		this.scriptPath = System.getenv(ENVIRONMENT_VARIABLE) + DEFAULT_BINARY;
		
		// files
		this.annotationFileName = TMP_PATH + "/" + StringUtils.randomString(30) + "_annotations.txt";
		this.rankedListFileName = TMP_PATH + "/" + StringUtils.randomString(30) + "_ranked_list.txt";
		this.outputFileName = StringUtils.randomString(30) + "_result.txt";
		
	}
	
	public void prepare() throws InvalidIndexException{
		
		// id list
		idList = rankedList.getDataFrame().getRowNames();
		// statistic
		statistic = ArrayUtils.toList(rankedList.getDataFrame().getColumnAsDoubleArray(0));
		// order ranked list
		int[] sortIndex = ListUtils.order(statistic);
		ListUtils.ordered(idList,sortIndex);
		ListUtils.ordered(statistic,sortIndex);
		if(order==ASCENDING_SORT){
			Collections.reverse(idList);
			Collections.reverse(statistic);			
		}		
		
	}
	
	public void run() throws InvalidIndexException, SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException, IOException {
				
		// prepare list
		prepare();
		
		// annotation		
		if(!isYourAnnotations) annotations = InfraredUtils.getAnnotations(dbConnector, idList, filter);		
	
		// prepare files		
		saveRankedList();
		saveAnnotations();
		
		// RComand
		RCommand rCommand = new RCommand(scriptPath,TMP_PATH);
		// input
		rCommand.addParam("rankingfile",rankedListFileName);
		rCommand.addParam("annotationfile",annotationFileName);
		// output
		rCommand.addParam("outfile", outputFileName);

		// RUN
		rCommand.exec();
		
		// Load results
		loadResults();
			
				
	}
		
	
	private void saveRankedList() throws IOException{
		IOUtils.write(rankedListFileName, rankedList.toString());
	}
	
	private void saveAnnotations() throws IOException{
		IOUtils.write(annotationFileName, annotations.toString());
	}
	
	private void loadResults() throws IOException{
		
		// init results
		results = new ArrayList<GeneSetAnalysisTestResult>();
		
		// load result file
		FileUtils.checkFile(TMP_PATH + "/" + outputFileName);
		String resultsContent = IOUtils.toString(TMP_PATH + "/" + outputFileName);
		List<String> terms = StringUtils.toList(resultsContent,"\n");
		
		// init params
		String[] fields;
		String term;
		double adjPValue,logRatio;
		boolean converged;
		int size;
		List<String>list1Ids;
		List<String>list2Ids;
		
		// run rows
		for(String row: terms){
			System.err.println("row:" + row);
			if(!row.startsWith("#") && row.contains("\t")){
				// split row
				fields = row.split("\t");		
				// params
				term = fields[0];
				size = Integer.parseInt(fields[1]);
				if(fields[2].trim().equals("1")) converged = true;
				else converged = false;
				logRatio = Double.parseDouble(fields[3]);
				adjPValue = Double.parseDouble(fields[4]);				
				GeneSetAnalysisTestResult gseaTest = new GeneSetAnalysisTestResult(term,null,size,converged,logRatio,adjPValue);
				results.add(gseaTest);
			}
		}		
	}
	
//	private HashMap<String,String> getTermHash(){
//		for(int i=0; i<annotations.size(); i++){
//			
//		}
//	}
	
	
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
