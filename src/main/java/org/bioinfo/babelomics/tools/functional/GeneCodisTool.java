package org.bioinfo.babelomics.tools.functional;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import org.apache.commons.cli.ParseException;

import org.bioinfo.babelomics.methods.functional.GeneSetAnalysisTestResult;
import org.bioinfo.babelomics.methods.functional.LogisticScan;
import org.bioinfo.babelomics.methods.functional.TwoListFisherTest;
import org.bioinfo.commons.exec.Command;
import org.bioinfo.commons.exec.SingleProcess;
import org.bioinfo.data.dataset.FeatureData;
import org.bioinfo.data.list.exception.InvalidIndexException;
import org.bioinfo.infrared.common.dbsql.DBConnector;
import org.bioinfo.infrared.funcannot.filter.Filter;
import org.bioinfo.math.exception.InvalidParameterException;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;

public class GeneCodisTool extends FunctionalProfilingTool{
	
	private FeatureData rankedList;
	private String support;   //Minimum number of genes required that have to be implicated in a rule to take it into account (default: 3).
	private String supportRandom;//i . i support for random:  Min support taked into account when the algorithm is correcting p-values (usually the same than Support).
	private int analysis;//  -a [1,2] Analysis. 1: Concurrence analysis, 2: Singular analysis
	private int test;//-t- [1,2] statistical test, 0: Hypergeometric test (default), 1: Chi Square test
	private int correction; // -s 0: none	any number < 0: FDR method	any number > 0: Number of permutations to correct p-values with permutations method

	
	@Override
	public void initOptions() {
		// parent options
		super.initOptions();
		
		
		options.addOption(OptionFactory.createOption("datalist", "the ranked list"));
		options.addOption(OptionFactory.createOption("support", "Minimum number of genes",true));
		options.addOption(OptionFactory.createOption("support-for-random", "Minimum number of genes for correcting p-values",true));
		options.addOption(OptionFactory.createOption("analysis", "singular_analysis or concurrence analysis",true));
		options.addOption(OptionFactory.createOption("hypergeometric", "",false));
		options.addOption(OptionFactory.createOption("chi-square", "",false));
		options.addOption(OptionFactory.createOption("correction", "default correction",true));
	}

	
	@Override
	public void prepare() throws IOException, ParseException, InvalidIndexException {
		super.prepare();		
		// ranked list		
		rankedList = new FeatureData(new File(commandLine.getOptionValue("datalist")), true);
		support = commandLine.getOptionValue("support");  
		supportRandom = commandLine.getOptionValue("support-for-random"); //i
		analysis = (commandLine.getOptionValue("analysis")=="concurrence_analysis")?1:2; //  -a [1,2] Analysis:
		test = ((commandLine.getOptionValue("hypergeometric")) != null? 0 : 1);//-t- [1,2]
		correction = 0; 
		if (commandLine.getOptionValue("correction") == "correction_fdr"){
			correction=-1;
		} else if (commandLine.getOptionValue("correction") == "correction_permutation"){
			correction=1;
		}

	}
	
		//	public List<GeneSetAnalysisTestResult> getSignificant(){
		//		return getSignificant(TwoListFisherTest.DEFAULT_PVALUE_THRESHOLD);
		//	}
		//	
		//	public List<GeneSetAnalysisTestResult> getSignificant(double threshold){
		//		List<GeneSetAnalysisTestResult> significant = new ArrayList<GeneSetAnalysisTestResult>();
		//		for(GeneSetAnalysisTestResult result: this.results){			
		//			if(result.getAdjPValue()<threshold) significant.add(result);
		//		}
		//		return significant;
		//	}
		//	
	
	@Override
	protected void execute() {
		try {
			
			// update status
			jobStatus.addStatusMessage("10", "Preparing data");
			
			// infrared connector			
			DBConnector dbConnector = new DBConnector(species, new File(babelomicsHomePath + "/conf/infrared.properties"));		
			
			// prepare params
			prepare();
			
			jobStatus.addStatusMessage("" + ("40"), "preparing execution");
			
			// significant terms
//			List<GeneSetAnalysisTestResult> significants = new ArrayList<GeneSetAnalysisTestResult>();				
//			List<String> significantOutput = new ArrayList<String>();
//			LogisticScan logistic = new LogisticScan(rankedList, annotations,order);	
//			
//			significants.addAll(logistic.getSignificant());
			
			
			}
		catch (Exception e) {
			
		}
	}
	
	private void doTest(){
		File inputFile =new File(commandLine.getOptionValue("datalist")).getAbsoluteFile();
		File outputFile = new File(outdir + "/geneCodisOut");
		String cmdStr = System.getenv("BABELOMICS_HOME") + "/bin/genecodis/genecodis_bin "+ inputFile+" "+ support +" -a" + analysis + " -i" + supportRandom + " -s"+ correction + "-t" + test+ "-o "+outputFile.getAbsolutePath();
		
		Command cmd = new Command(cmdStr); 
		SingleProcess sp = new SingleProcess(cmd);
		sp.runSync();
		updateJobStatus("90", "saving results");
		if ( outputFile.exists() ) {
			result.addOutputItem(new Item("genecodis_file", outputFile.getName(), "Genecodis file", TYPE.FILE, Arrays.asList("ANNOTATION"), new HashMap<String, String>(2), "geneCodis data"));
		}	
	}
	
}
