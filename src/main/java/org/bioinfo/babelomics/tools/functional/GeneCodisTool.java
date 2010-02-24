package org.bioinfo.babelomics.tools.functional;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;

import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import org.apache.commons.cli.OptionGroup;
import org.apache.commons.cli.ParseException;
import org.apache.commons.math.genetics.Chromosome;

import org.bioinfo.babelomics.methods.functional.FatiGO;
import org.bioinfo.babelomics.methods.functional.GeneCodis;

import org.bioinfo.commons.exec.Command;
import org.bioinfo.commons.exec.SingleProcess;
import org.bioinfo.commons.io.utils.IOUtils;

import org.bioinfo.data.dataset.FeatureData;
import org.bioinfo.data.list.exception.InvalidIndexException;
import org.bioinfo.infrared.common.dbsql.DBConnector;
import org.bioinfo.infrared.funcannot.filter.Filter;
import org.bioinfo.math.exception.InvalidParameterException;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;

public class GeneCodisTool extends FunctionalProfilingTool{
	
	// list1
	private FeatureData list1;
	
	// list2
	private FeatureData list2;		
	
	
	// other capabilities
	protected int duplicatesMode;
	
	
	private List<Chromosome> chromosomes;
	
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
		
		OptionGroup input2 = new OptionGroup();
		input2.setRequired(false);		
			input2.addOption(OptionFactory.createOption("list2", "the file containig the list #2 of genes, or the feature data file",false));
			input2.addOption(OptionFactory.createOption("genome", "compares list #1 with the rest of genome",false,false));
			input2.addOption(OptionFactory.createOption("chromosomes", "the chromosome region to compare with list #1. Use chromosome-start and chromosome-end options tp delimite the region within the chromosome",false));
		options.addOptionGroup(input2);
		
		// extras
		options.addOption(OptionFactory.createOption("duplicates", "to remove duplicated IDs, valid values: never, each, ref. By default, never", false));
		
	}

	
	@Override
	public void prepare() throws IOException, ParseException, InvalidIndexException {
		super.prepare();		

		
		// list 1		
		list1 = new FeatureData(new File(commandLine.getOptionValue("datalist")), true);
		
		// list 2: list to compare
		if(commandLine.hasOption("list2")){
			list2 = new FeatureData(new File(commandLine.getOptionValue("list2")), true);
		} else if(commandLine.hasOption("genome")) {
			this.setRestOfGenome(true);
		} else if(commandLine.hasOption("chromosomes")) {
			throw new ParseException("chromosome reading not yet implemented");
		} else {
			throw new ParseException("No comparison provided, use list2, genome or chromosome options to set your comparison");
		}
		
		//command parameters
		support = commandLine.getOptionValue("support");  
		supportRandom = commandLine.getOptionValue("support-for-random"); //i
		analysis = (commandLine.getOptionValue("analysis")=="concurrence")?1:2; //  -a [1,2] Analysis:
		test = ((commandLine.getOptionValue("hypergeometric")) != null? 0 : 1);//-t- [1,2]
		correction = 0; 
		if (commandLine.getOptionValue("correction") == "fdr"){
			correction=-1;
		} else if (commandLine.getOptionValue("correction") == "permutation"){
			correction=1;
		}
		
		// extras
		String duplicates = commandLine.getOptionValue("remove-duplicates", "never");
		if(duplicates.equalsIgnoreCase("each")) duplicatesMode = FatiGO.REMOVE_EACH;
		if(duplicates.equalsIgnoreCase("ref")) duplicatesMode = FatiGO.REMOVE_REF;
		if(duplicates.equalsIgnoreCase("all")) duplicatesMode = FatiGO.REMOVE_ALL;		

	}
	
	@Override
	protected void execute() {
		try {
			logger.println("Starting GeneCodis...");
			// update status
			jobStatus.addStatusMessage("10", "Preparing data");
			
			// infrared connector			
			DBConnector dbConnector = new DBConnector(species, new File(babelomicsHomePath + "/conf/infrared.properties"));		
			
			// prepare params
			logger.println("Starting 0...");
			prepare();
			logger.println("Starting 1...");
			GeneCodis genecodis = null;
			String list2label = "List 2 (after duplicates managing)";
			
			jobStatus.addStatusMessage("" + ("40"), "preparing Genecodis execution");
			
			
			// list 1			
			System.err.println(list1.getDataFrame());
			List<String> idList1 = list1.getDataFrame().getRowNames();//list1.getDataFrame().getColumn(0); //InfraredUtils.toEnsemblId(dbConnector, list1.getDataFrame().getColumn(0));
			
			// list 2
			List<String> idList2 = null;
//			if(list2!=null) {
//				logger.println("Starting 3...");
//				System.err.println(list2.getDataFrame());
//				idList2 = list2.getDataFrame().getRowNames();//list2.getDataFrame().getColumn(0); //InfraredUtils.toEnsemblId(dbConnector, list2.getDataFrame().getColumn(0));
//				genecodis = new GeneCodis(idList1, idList2, null, dbConnector, testMode, duplicatesMode);
//				
//			} else if(isRestOfGenome()) {
//				logger.println("Starting 4...");
//				duplicatesMode = GeneCodis.REMOVE_GENOME;
//				genecodis = new GeneCodis(idList1, null, dbConnector);
//				list2label = "Genome (after duplicates managing)";
//			} else if(chromosomes!=null) {
//				throw new ParseException("chromosomes comparison not yet implemented");
//			} else {
//				throw new ParseException("No comparison provided, use list2, genome or chromosomes options to set your comparison");				
//			}
//			logger.println("Starting 1...");
//			// run genecodis
			if(filterList.size()==0 && !isYourAnnotations){
				throw new ParseException("No biological database selected (eg. --go-bp)");
			} else {
				logger.println("filterList.size()::::::::"+filterList.size() + "........toString:........"+filterList.toString());
				for(Filter filter: filterList) {
					logger.println("filter::::::::"+filter.toString());
					doTest(idList1,idList2,filter,dbConnector);
					}
				}
			
			}
		catch (Exception e) {
			
		}
		
		
	}
	
	private void doTest(List<String> idList1, List<String> idList2, Filter filter, DBConnector dbConnector) throws SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException, InvalidParameterException{
		File inputFile = new File(outdir + "/geneInputWellFormed");
		logger.println("Starting doing test...");
		GeneCodis genecodis = null;
		if(list2!=null) genecodis = new GeneCodis(idList1, idList2, filter, dbConnector, testMode, duplicatesMode);
		else if(isRestOfGenome()) genecodis = new GeneCodis(idList1, filter, dbConnector);
		
		//run run...
		logger.println("run...........");
		genecodis.run();
		logger.println("test done...");
		
		try {
			IOUtils.write(inputFile, genecodis.getResults());
			
		} catch (IOException e) {
			e.printStackTrace();
		}
		
		File outputFile = new File(outdir + "/geneCodisOut");
		String cmdStr = System.getenv("BABELOMICS_HOME") + "/bin/genecodis/genecodis_bin "+ inputFile.getAbsolutePath()+" "+ support +" -a" + analysis + " -i" + supportRandom +" -r"+genecodis.getrFactor()+ " -R"+genecodis.getRFactor() + " -s"+ correction + "-t" + test+ "-o "+outputFile.getAbsolutePath();
		
		Command cmd = new Command(cmdStr); 
		SingleProcess sp = new SingleProcess(cmd);
		sp.runSync();
		updateJobStatus("90", "saving results");
		if ( outputFile.exists() ) {
			logger.println("exist...");
			result.addOutputItem(new Item("genecodis_file", outputFile.getName(), "Genecodis file", TYPE.FILE, Arrays.asList("ANNOTATION"), new HashMap<String, String>(2), "geneCodis data"));
		}	
	}
	
}
