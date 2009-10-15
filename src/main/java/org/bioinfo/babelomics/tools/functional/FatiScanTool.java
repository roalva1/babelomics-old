package org.bioinfo.babelomics.tools.functional;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.ParseException;
import org.bioinfo.babelomics.methods.functional.FatiScan;
import org.bioinfo.babelomics.methods.functional.TwoListFisherTest;
import org.bioinfo.babelomics.methods.functional.TwoListFisherTestResult;
import org.bioinfo.babelomics.methods.functional.graph.FatiScanGraph;
import org.bioinfo.collections.exceptions.InvalidColumnIndexException;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.data.dataset.FeatureData;
import org.bioinfo.infrared.common.dbsql.DBConnector;
import org.bioinfo.infrared.common.feature.FeatureList;
import org.bioinfo.infrared.funcannot.AnnotationItem;
import org.bioinfo.infrared.funcannot.filter.Filter;
import org.bioinfo.math.exception.InvalidParameterException;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;

public class FatiScanTool  extends FunctionalProfilingTool{
		
	// list1
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
	public void prepare() throws IOException, ParseException, InvalidColumnIndexException {
		super.prepare();

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
			// infrared connector			
			DBConnector dbConnector = new DBConnector(getSpecies(), new File(System.getenv("BABELOMICS_HOME") + "/conf/infrared.conf"));			
			
			// prepare params
			prepare();			
	
			// run fatigo's
			if(filterList.size()==0  && !isYourAnnotations){
				throw new ParseException("No biological database selected (eg. --go-bp)");
			} else {
				
				// significant terms
				List<TwoListFisherTestResult> significants = new ArrayList<TwoListFisherTestResult>();				
				List<String> significantOutput = new ArrayList<String>();
				
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
								
				// Significant results must appear after than complete tables!!
				result.addOutputItem(new Item("significant","significant_" + TwoListFisherTest.DEFAULT_PVALUE_THRESHOLD + ".txt","Significant terms",Item.TYPE.FILE,Arrays.asList("TABLE","FATISCAN_TABLE"),new HashMap<String,String>(),"Significant Results"));
								
				// do fatiscans
				for(Filter filter: filterList) {
					doFatiScan(rankedList,filter,dbConnector,significants);
				}
				if(isYourAnnotations){
					doFatiScanYourAnnotations(rankedList,yourAnnotations,significants);
				}
				
				significantOutput.addAll(testResultToStringList(significants));
				
				// significant terms
				System.err.println("sig size: " + significants.size());
				if(significants!=null && significants.size()>0){						
					IOUtils.write(outdir + "/significant_" + TwoListFisherTest.DEFAULT_PVALUE_THRESHOLD + ".txt", ListUtils.toString(significantOutput,"\n"));
					FatiScanGraph.fatiScanGraph(significants,"Significant terms graph for all databases",outdir + "/significant_graph.png");			
					result.addOutputItem(new Item("graph_alldbs","significant_graph.png","Significant terms graph for all databases",Item.TYPE.IMAGE,Arrays.asList("FATISCAN"),new HashMap<String,String>(),"Significant Results"));
				} else {
					result.addOutputItem(new Item("graph_alldbs","No significant terms found","Significant terms for all databases",Item.TYPE.MESSAGE,Arrays.asList("WARNING"),new HashMap<String,String>(),"Significant Results"));
				}
				
									
			}
					
		}catch(Exception e){
			e.printStackTrace();
		}
		
	}

	private void doFatiScan(FeatureData rankedList,Filter filter,DBConnector dbConnector, List<TwoListFisherTestResult> significant) throws IOException, SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException, InvalidParameterException, InvalidColumnIndexException{
		
		// db attributes
		String name = getDBName(filter);
		String title = getDBTitle(filter);					
		
		
		logger.info(title + "...\n");

		// init test
		FatiScan fatiscan = new FatiScan(rankedList,filter,dbConnector,numberOfPartitions,testMode,outputFormat,order);
		
		// run
		fatiscan.run();

		// save results
		saveFatiScanResults(fatiscan,name,title);
				
		// acum significant values
		significant.addAll(fatiscan.getSignificant());
		
		logger.info("...end of " + title);
		
	}
	
	private void doFatiScanYourAnnotations(FeatureData rankedList, FeatureList<AnnotationItem> annotations, List<TwoListFisherTestResult> significant) throws IOException, SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException, InvalidParameterException, InvalidColumnIndexException{
		
		// db attributes
		String name = "your_annotations";
		String title = "Your annotations";		
		
		logger.info(title + "...\n");

		// init test
		FatiScan fatiscan = new FatiScan(rankedList, annotations,numberOfPartitions,testMode,outputFormat,order);
		
		// run
		fatiscan.run();
				
		// save results
		saveFatiScanResults(fatiscan,name,title);
		
		// acum significant values
		significant.addAll(fatiscan.getSignificant());
		
		logger.info("...end of " + title);
		
	}
	
	private void saveFatiScanResults(FatiScan fatiscan, String name, String title) throws IOException{
		
		String fileName = name + ".txt";
		String annotFileName = name + ".annot";
		
		// save statistic results					
		List<String> testResultOutput = testResultToStringList(fatiscan.getResults());
		
		IOUtils.write(outdir + "/" + fileName, StringUtils.join(testResultOutput,"\n"));
		result.addOutputItem(new Item(name,fileName,title,Item.TYPE.FILE,Arrays.asList("TABLE","FATISCAN_TABLE"),new HashMap<String,String>(),"Database tests"));
						
		// save annotation
		IOUtils.write(outdir + "/" + annotFileName, fatiscan.getAnnotations().toString());
		result.addOutputItem(new Item("annot_" + name,annotFileName,"Annotations for " + title,Item.TYPE.FILE,Arrays.asList("ANNOTATION"),new HashMap<String,String>(),"Annotation files"));
		
		// save graph
		List<TwoListFisherTestResult> significants = fatiscan.getSignificant();
		if(significants!=null && significants.size()>0){
			String graphFileName = name + "_graph.png"; 
			FatiScanGraph.fatiScanGraph(significants,"Significant terms graph for " + title,outdir + "/" +  graphFileName);			
			result.addOutputItem(new Item("graph_" + name,graphFileName,"Significant terms graph for " + title,Item.TYPE.IMAGE,Arrays.asList("FATISCAN"),new HashMap<String,String>(),"Database tests"));
		} else {
			result.addOutputItem(new Item("graph_" + name,"No significant terms found","Significant terms graph for " + title,Item.TYPE.MESSAGE,Arrays.asList("WARNING"),new HashMap<String,String>(),"Database tests"));
		}

		
	}

}
