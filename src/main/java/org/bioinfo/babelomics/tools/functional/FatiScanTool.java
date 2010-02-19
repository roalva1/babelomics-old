package org.bioinfo.babelomics.tools.functional;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import org.apache.commons.cli.ParseException;
import org.bioinfo.babelomics.exception.EmptyAnnotationException;
import org.bioinfo.babelomics.methods.functional.FatiScan;
import org.bioinfo.babelomics.methods.functional.GeneSetAnalysisTestResult;
import org.bioinfo.babelomics.methods.functional.LogisticScan;
import org.bioinfo.babelomics.methods.functional.TwoListFisherTest;
import org.bioinfo.babelomics.methods.functional.graph.FatiScanGraph;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.data.dataset.FeatureData;
import org.bioinfo.data.list.exception.InvalidIndexException;
import org.bioinfo.infrared.common.dbsql.DBConnector;
import org.bioinfo.infrared.common.feature.FeatureList;
import org.bioinfo.infrared.funcannot.AnnotationItem;
import org.bioinfo.infrared.funcannot.filter.Filter;
import org.bioinfo.math.exception.InvalidParameterException;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;

public class FatiScanTool  extends FunctionalProfilingTool{
		
	// list1
	private enum Method {FatiScan,Logistic};	
	private Method method;
	private FeatureData rankedList;
	private int numberOfPartitions;
	private int outputFormat;
	private int order;
	
	public FatiScanTool(){
		initOptions();
	}

	@Override
	public void initOptions() {
		// parent options
		super.initOptions();
		// mode
		options.addOption(OptionFactory.createOption("method", "Gene set analysis method [fatiscan, logistic]", false, true));
		// list 1
		options.addOption(OptionFactory.createOption("ranked-list", "the feature data containig the ranked list"));
		// test
		options.addOption(OptionFactory.createOption("partitions", "the number of partitions",false,true));
		// extras
		options.addOption(OptionFactory.createOption("output-format", "short (just most significant partition) or long (term results for all partitions) [short]",false,true));
		options.addOption(OptionFactory.createOption("order", "ascend or descend [ascend]",false,true));
		options.addOption(OptionFactory.createOption("higher-label", "label for condition with higher statistical values",false,true));
		options.addOption(OptionFactory.createOption("lower-label", "label for condition with lower statistical values",false,true));
		
	}
		
	@Override
	public void prepare() throws IOException, ParseException, InvalidIndexException {
		super.prepare();

		// method 
		method = Method.Logistic;
		
		if(commandLine.hasOption("method") && commandLine.getOptionValue("method").equalsIgnoreCase("fatiscan")) {
			method = Method.FatiScan;
		}
		
		// ranked list		
		rankedList = new FeatureData(new File(commandLine.getOptionValue("ranked-list")), true);
		// test
		numberOfPartitions = Integer.parseInt(commandLine.getOptionValue("partitions","" + FatiScan.DEFAULT_NUMBER_OF_PARTITIONS));
		// output format
		if(commandLine.hasOption("output-format") && commandLine.getOptionValue("output-format").equalsIgnoreCase("long")) {
			outputFormat = FatiScan.LONG_FORMAT;
		} else {
			outputFormat = FatiScan.SHORT_FORMAT;
		}
		// sort order
		if(commandLine.hasOption("order") && commandLine.getOptionValue("order").equalsIgnoreCase("ascend")) {
			order = FatiScan.ASCENDING_SORT;
		} else {
			order = FatiScan.DESCENDING_SORT;
		}

	}
	
	@Override
	protected void execute() {
		try {
			
			// update status
			jobStatus.addStatusMessage("10", "Preparing data");
			
			// infrared connector			
			DBConnector dbConnector = new DBConnector(species, new File(babelomicsHomePath + "/conf/infrared.properties"));		
			
			// prepare params
			prepare();			
	
			// run fatigo's
			if(filterList.size()==0  && !isYourAnnotations){
				throw new ParseException("No biological database selected (eg. --go-bp)");
			} else {
				
				// save id lists
				IOUtils.write(outdir + "/ranked_list.txt", rankedList.toString());
				result.addOutputItem(new Item("ranked_list","ranked_list.txt","Ranked list",Item.TYPE.FILE,Arrays.asList("RANKED_LIST","CLEAN"),new HashMap<String,String>(),"Input data"));
				
				// ordered list
				FatiScan fatiscan = new FatiScan(rankedList,null,dbConnector,FatiScan.DEFAULT_NUMBER_OF_PARTITIONS,testMode,outputFormat,order);
				fatiscan.prepare();
				IOUtils.write(outdir + "/id_list.txt", fatiscan.getIdList());
				result.addOutputItem(new Item("id_list","id_list.txt","Id list (sorted)",Item.TYPE.FILE,Arrays.asList("IDLIST","SORTED"),new HashMap<String,String>(),"Input data"));
				IOUtils.write(outdir + "/statistic.txt", ListUtils.toStringArray(fatiscan.getStatistic()));
				result.addOutputItem(new Item("statistic","statistic.txt","Statistic (sorted)",Item.TYPE.FILE,Arrays.asList("STATISTIC","SORTED"),new HashMap<String,String>(),"Input data"));
				
				// significant terms
				List<GeneSetAnalysisTestResult> significants = new ArrayList<GeneSetAnalysisTestResult>();				
				List<String> significantOutput = new ArrayList<String>();
					
				
				// do fatiscans
				double progress = 20;
				double inc = 60.0/filterList.size();
				for(Filter filter: filterList) {
					doTest(rankedList,filter,dbConnector,significants,method);
					// update status					
					jobStatus.addStatusMessage("" + progress, "Executing test");
					progress+=inc;
				}
				
				// update status					
				jobStatus.addStatusMessage("90", "Executing test");
				
				if(isYourAnnotations){
					doYourAnnotationsTest(rankedList,yourAnnotations,significants,method);
				}
				
				// update status					
				jobStatus.addStatusMessage("95", "Saving results");
				
				System.err.println("significants: " + significants.size());
				if(method==Method.Logistic){
					significantOutput.addAll(LogisticResultToStringList(significants));
				} else {
					significantOutput.addAll(FatiScanResultToStringList(significants));					
				}		
				
				// significant terms				
				if(significants!=null && significants.size()>0){
					
					// Significant results must appear after than complete tables!!
					result.getOutputItems().add(0, new Item("significant","significant_" + TwoListFisherTest.DEFAULT_PVALUE_THRESHOLD + ".txt","Significant terms",Item.TYPE.FILE,Arrays.asList("TABLE","FATISCAN_TABLE"),new HashMap<String,String>(),"Significant Results"));
					IOUtils.write(outdir + "/significant_" + TwoListFisherTest.DEFAULT_PVALUE_THRESHOLD + ".txt", ListUtils.toString(significantOutput,"\n"));
					
					// FatiScan graph
					if(method==Method.FatiScan){
						FatiScanGraph.fatiScanGraph(significants,"Significant terms graph for all databases",outdir + "/significant_graph.png");
						result.addOutputItem(new Item("graph_alldbs","significant_graph.png","Significant term graph",Item.TYPE.IMAGE,Arrays.asList("FATISCAN"),new HashMap<String,String>(),"Significant results (all databases together)"));
					}
					
				} else {
					result.addOutputItem(new Item("graph_alldbs","No significant terms found","Significant term graph",Item.TYPE.MESSAGE,Arrays.asList("WARNING"),new HashMap<String,String>(),"Significant results (all databases together)"));
				}		
			}
					
		}catch(Exception e){
			e.printStackTrace();
		}
		
	}

	
	
	/*
	 * 
	 * FATISCAN
	 * 
	 * 
	 */
	
	
	private void doTest(FeatureData rankedList,Filter filter,DBConnector dbConnector, List<GeneSetAnalysisTestResult> significant, Method method) throws IOException, SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException, InvalidParameterException, InvalidIndexException{
		
		// db attributes
		String name = getDBName(filter);
		String title = getDBTitle(filter);					
				
		logger.info(title + "...\n");

		if(method==Method.Logistic){
			// init test
			LogisticScan logistic = new LogisticScan(rankedList,filter,dbConnector,order);		
			// run
			logistic.run();
			// save results
			saveLogisticScanResults(logistic,name,title);
			// acum significant values
			significant.addAll(logistic.getSignificant());
		} else {
			// init test
			FatiScan fatiscan = new FatiScan(rankedList,filter,dbConnector,numberOfPartitions,testMode,outputFormat,order);		
			// run
			try{
				fatiscan.run();
				// save results
				saveFatiScanResults(fatiscan,name,title);
				// acum significant values
				significant.addAll(fatiscan.getSignificant());
			} catch (EmptyAnnotationException ene){
				result.addOutputItem(new Item("annot_" + name,"No annotation was found for " + name + " ids","Annotations for " + title,Item.TYPE.MESSAGE,Arrays.asList("WARNING"),new HashMap<String,String>(),"Annotation files"));
			}

		}
				

		
		logger.info("...end of " + title);
		
	}
	
	private void doYourAnnotationsTest(FeatureData rankedList, FeatureList<AnnotationItem> annotations, List<GeneSetAnalysisTestResult> significant, Method method) throws IOException, SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException, InvalidParameterException, InvalidIndexException{
		
		// db attributes
		String name = "your_annotations";
		String title = "Your annotations";		
		
		logger.info(title + "...\n");
		if(method==Method.Logistic){
			// init test
			LogisticScan logistic = new LogisticScan(rankedList, annotations,order);		
			// run
			logistic.run();				
			// save results
			saveLogisticScanResults(logistic,name,title);
			// acum significant values
			significant.addAll(logistic.getSignificant());
		} else {
			// init test
			FatiScan fatiscan = new FatiScan(rankedList, annotations,numberOfPartitions,testMode,outputFormat,order);		
			// run
			try{
				fatiscan.run();				
				// save results
				saveFatiScanResults(fatiscan,name,title);		
				// acum significant values
				significant.addAll(fatiscan.getSignificant());
			} catch (EmptyAnnotationException ene){
				result.addOutputItem(new Item("annot_" + name,"No annotation was found for " + name + " ids","Annotations for " + title,Item.TYPE.MESSAGE,Arrays.asList("WARNING"),new HashMap<String,String>(),"Annotation files"));
			}
		}
		
		logger.info("...end of " + title);
		
	}
	
	private void saveFatiScanResults(FatiScan fatiscan, String name, String title) throws IOException{
		
		String fileName = name + ".txt";
		String annotFileName = name + ".annot";
		
		// save statistic results					
		List<String> testResultOutput = FatiScanResultToStringList(fatiscan.getResults());
		
		IOUtils.write(outdir + "/" + fileName, ListUtils.toString(testResultOutput,"\n"));
		result.addOutputItem(new Item(name,fileName,title,Item.TYPE.FILE,Arrays.asList("TABLE","FATISCAN_TABLE"),new HashMap<String,String>(),"Database tests"));
						
		// save annotation
		IOUtils.write(outdir + "/" + annotFileName, fatiscan.getAnnotations().toString());
		result.addOutputItem(new Item("annot_" + name,annotFileName,"Annotations for " + title,Item.TYPE.FILE,Arrays.asList("ANNOTATION"),new HashMap<String,String>(),"Annotation files"));
		
		// save graph
		List<GeneSetAnalysisTestResult> significants = fatiscan.getSignificant();
		if(significants!=null && significants.size()>0){
			String graphFileName = name + "_graph.png"; 
			FatiScanGraph.fatiScanGraph(significants,"Significant terms graph for " + title,outdir + "/" +  graphFileName);			
			result.addOutputItem(new Item("graph_" + name,graphFileName,"Significant terms graph for " + title,Item.TYPE.IMAGE,Arrays.asList("FATISCAN"),new HashMap<String,String>(),"Database tests"));
		} else {
			result.addOutputItem(new Item("graph_" + name,"No significant terms found","Significant terms graph for " + title,Item.TYPE.MESSAGE,Arrays.asList("WARNING"),new HashMap<String,String>(),"Database tests"));
		}
		
	}
	
	private void saveLogisticScanResults(LogisticScan logistic, String name, String title) throws IOException{
		
		String fileName = name + ".txt";
		String annotFileName = name + ".annot";
		
		// save statistic results					
		List<String> testResultOutput = LogisticResultToStringList(logistic.getResults());
		
		IOUtils.write(outdir + "/" + fileName, ListUtils.toString(testResultOutput,"\n"));
		result.addOutputItem(new Item(name,fileName,title,Item.TYPE.FILE,Arrays.asList("TABLE","FATISCAN_TABLE"),new HashMap<String,String>(),"Database tests"));
						
		// save annotation
		IOUtils.write(outdir + "/" + annotFileName, logistic.getAnnotations().toString());
		result.addOutputItem(new Item("annot_" + name,annotFileName,"Annotations for " + title,Item.TYPE.FILE,Arrays.asList("ANNOTATION"),new HashMap<String,String>(),"Annotation files"));
		
//		// save graph
//		List<GeneSetAnalysisTestResult> significants = logistic.getSignificant();
//		if(significants!=null && significants.size()>0){
//			String graphFileName = name + "_graph.png"; 
//			FatiScanGraph.fatiScanGraph(significants,"Significant terms graph for " + title,outdir + "/" +  graphFileName);			
//			result.addOutputItem(new Item("graph_" + name,graphFileName,"Significant terms graph for " + title,Item.TYPE.IMAGE,Arrays.asList("FATISCAN"),new HashMap<String,String>(),"Database tests"));
//		} else {
//			result.addOutputItem(new Item("graph_" + name,"No significant terms found","Significant terms graph for " + title,Item.TYPE.MESSAGE,Arrays.asList("WARNING"),new HashMap<String,String>(),"Database tests"));
//		}
		
	}
	
	
	protected List<String> FatiScanResultToStringList(List<GeneSetAnalysisTestResult> testResult){
		return FatiScanResultToStringList(testResult,true);
	}	
	protected List<String> FatiScanResultToStringList(List<GeneSetAnalysisTestResult> testResult, boolean header){
		List<String> result = new ArrayList<String>();
		if(header) result.add(GeneSetAnalysisTestResult.fatiScanHeader());
		for(int i=0; i<testResult.size(); i++){			
			result.add(testResult.get(i).toFatiScanString());
		}
		return result;
	}
	
	protected List<String> LogisticResultToStringList(List<GeneSetAnalysisTestResult> testResult){
		return LogisticResultToStringList(testResult,true);
	}	
	protected List<String> LogisticResultToStringList(List<GeneSetAnalysisTestResult> testResult, boolean header){
		List<String> result = new ArrayList<String>();
		if(header) result.add(GeneSetAnalysisTestResult.LogisticHeader());
		for(int i=0; i<testResult.size(); i++){			
			result.add(testResult.get(i).toLogisticString());
		}
		return result;
	}
	
	
}
