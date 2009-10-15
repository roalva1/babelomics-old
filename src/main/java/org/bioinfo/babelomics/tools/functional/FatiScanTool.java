package org.bioinfo.babelomics.tools.functional;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.ParseException;
import org.bioinfo.babelomics.methods.functional.FatiScan;
import org.bioinfo.babelomics.methods.functional.TwoListFisherTest;
import org.bioinfo.babelomics.methods.functional.TwoListFisherTestResult;
import org.bioinfo.collections.exceptions.InvalidColumnIndexException;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.data.dataset.FeatureData;
import org.bioinfo.infrared.common.dbsql.DBConnector;
import org.bioinfo.infrared.funcannot.filter.Filter;
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
		options.addOption(OptionFactory.createOption("partitions", "the number of partitions"));
		// extras
		options.addOption(OptionFactory.createOption("output-format", "short (just most significant partition) or long (term results for all partitions) [short]",false,true));
		options.addOption(OptionFactory.createOption("order", "ascend or descend [ascend]",false,true));
		options.addOption(OptionFactory.createOption("greater-label", "label for condition with greater statistical values",false,true));
		options.addOption(OptionFactory.createOption("lower-label", "label for condition with lower statistical values",false,true));
		
	}
		
	@Override
	public void prepare(CommandLine cmdLine) throws IOException, ParseException, InvalidColumnIndexException {
		super.prepare(cmdLine);

		// list 1		
		rankedList = new FeatureData(new File(cmdLine.getOptionValue("ranked-list")), true);
		// test
		numberOfPartitions = Integer.parseInt(cmdLine.getOptionValue("partitions","" + FatiScan.DEFAULT_NUMBER_OF_PARTITIONS));
		// output format
		if(cmdLine.hasOption("output-format") && cmdLine.getOptionValue("output-format").equalsIgnoreCase("long")) {
			outputFormat = FatiScan.LONG_FORMAT;
		} else {
			outputFormat = FatiScan.SHORT_FORMAT;
		}
		// sort order
		if(cmdLine.hasOption("order") && cmdLine.getOptionValue("order").equalsIgnoreCase("ascend")) {
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
			prepare(commandLine);			
	
			// run fatigo's
			if(filterList.size()==0){
				throw new ParseException("No biological database selected (eg. --go-bp)");
			} else {
				List<String> significant = new ArrayList<String>();
				String name,fileName,annotFileName,title;
				
				// save id lists				
				IOUtils.write(outdir + "/ranked_list.txt", rankedList.toString());
				result.addOutputItem(new Item("ranked_list","ranked_list.txt","Ranked list",Item.TYPE.FILE,Arrays.asList("RANKED_LIST","CLEAN"),new HashMap<String,String>(),"Input data"));
				
				// ordered list
				FatiScan fatiscan = new FatiScan(rankedList,filterList.get(0),dbConnector,FatiScan.DEFAULT_NUMBER_OF_PARTITIONS,testMode,outputFormat,order);
				fatiscan.prepare();
				IOUtils.write(outdir + "/id_list.txt", fatiscan.getIdList());
				result.addOutputItem(new Item("id_list","id_list.txt","Id list (sorted)",Item.TYPE.FILE,Arrays.asList("IDLIST","SORTED"),new HashMap<String,String>(),"Input data"));
				IOUtils.write(outdir + "/statistic.txt", ListUtils.toStringArray(fatiscan.getStatistic()));
				result.addOutputItem(new Item("statistic","statistic.txt","Statistic (sorted)",Item.TYPE.FILE,Arrays.asList("STATISTIC","SORTED"),new HashMap<String,String>(),"Input data"));
				
				
				// Significant results must appear after than complete tables!!
				result.addOutputItem(new Item("significant","significant_" + TwoListFisherTest.DEFAULT_PVALUE_THRESHOLD + ".txt","Significant terms",Item.TYPE.FILE,Arrays.asList("TABLE"),new HashMap<String,String>(),"Significant Results"));
				significant.add(TwoListFisherTestResult.header());
				
				for(Filter filter: filterList) {

					// db attributes
					name = getDBName(filter);
					title = getDBTitle(filter);					
					fileName = name + ".txt";
					annotFileName = name + ".annot";
					
					logger.info(title + "...\n");

					// init test
					fatiscan = new FatiScan(rankedList,filter,dbConnector,numberOfPartitions,testMode,outputFormat,order);
					fatiscan.run();
					
					// save statistic results					
					List<String> testResultOutput = testResultToStringList(fatiscan.getResults());
					
					IOUtils.write(outdir + "/" + fileName, ListUtils.toString(testResultOutput, "\n"));
					result.addOutputItem(new Item(name,fileName,title,Item.TYPE.FILE,Arrays.asList("TABLE"),new HashMap<String,String>(),"Database tests"));
									
					// save annotation
					IOUtils.write(outdir + "/" + annotFileName, fatiscan.getAnnotations().toString());
					result.addOutputItem(new Item("annot_" + name,annotFileName,"Annotations for " + title,Item.TYPE.FILE,Arrays.asList("ANNOTATION"),new HashMap<String,String>(),"Annotation files"));
					
					// acum significant values
					significant.addAll(testResultToStringList(fatiscan.getSignificant(),false));
					
					logger.info("...end of " + title);
					
				}
				
				// significant terms
				IOUtils.write(outdir + "/significant_" + TwoListFisherTest.DEFAULT_PVALUE_THRESHOLD + ".txt", ListUtils.toString(significant,"\n"));
				
			}
					
		}catch(Exception e){
			e.printStackTrace();
		}		
	}


}
