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
import org.bioinfo.babelomics.methods.functional.GeneCodis;
import org.bioinfo.babelomics.methods.functional.GeneCodis.analysisFactor;
import org.bioinfo.babelomics.methods.functional.GeneCodis.correctionFactor;
import org.bioinfo.babelomics.methods.functional.GeneCodis.testFactor;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.data.dataset.FeatureData;
import org.bioinfo.data.list.exception.InvalidIndexException;
import org.bioinfo.infrared.common.dbsql.DBConnector;
import org.bioinfo.infrared.funcannot.filter.FunctionalFilter;
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
	
	private int support;   //Minimum number of genes required that have to be implicated in a rule to take it into account (default: 3).
	private int supportRandom;//i . i support for random:  Min support taked into account when the algorithm is correcting p-values (usually the same than Support).
	private analysisFactor analysis;//  -a [1,2] Analysis. 1: Concurrence analysis, 2: Singular analysis
	private testFactor test;//-t- [1,2] statistical test, 0: Hypergeometric test (default), 1: Chi Square test
	private correctionFactor correction; // -s 0: none	any number < 0: FDR method	any number > 0: Number of permutations to correct p-values with permutations method

	
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
		support = Integer.parseInt(commandLine.getOptionValue("support"));  
		supportRandom = Integer.parseInt(commandLine.getOptionValue("support-for-random")); //i
		
		analysis = (commandLine.getOptionValue("analysis")=="concurrence")?analysisFactor.concurrence:analysisFactor.singular; //  -a [1,2] Analysis:
		
			if (((commandLine.getOptionValue("hypergeometric")) != null) && ((commandLine.getOptionValue("chi-square")) != null)){
				test = testFactor.both;
			}else if((commandLine.getOptionValue("hypergeometric")) != null){
				test = testFactor.hypergeometric;
			}else if((commandLine.getOptionValue("chi-square")) != null){
				test = testFactor.chiSquare;
			} 
			
		
		correction = correctionFactor.none; 
		if (commandLine.getOptionValue("correction") == "fdr"){
			correction = correctionFactor.fdr;
			
		} else if (commandLine.getOptionValue("correction") == "permutation"){
			correction = correctionFactor.permutation;
		
		}
		
		// extras
		String duplicates = commandLine.getOptionValue("remove-duplicates", "never");
		if(duplicates.equalsIgnoreCase("each")) duplicatesMode = GeneCodis.REMOVE_EACH;
		if(duplicates.equalsIgnoreCase("ref")) duplicatesMode = GeneCodis.REMOVE_REF;
		if(duplicates.equalsIgnoreCase("all")) duplicatesMode = GeneCodis.REMOVE_ALL;		
		if(duplicates.equalsIgnoreCase("never")) duplicatesMode = GeneCodis.REMOVE_NEVER;

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
		
			prepare();
			GeneCodis genecodis = null;
			String list2label = "List 2 (after duplicates managing)";
			
			jobStatus.addStatusMessage("" + ("40"), "preparing Genecodis execution");
			
			
			// list 1			
			List<String> idList1 = list1.getDataFrame().getRowNames();//list1.getDataFrame().getColumn(0); //InfraredUtils.toEnsemblId(dbConnector, list1.getDataFrame().getColumn(0));
			
			// list 2
			List<String> idList2 = null;
			
			if(filterList.size()==0 && !isYourAnnotations){
				throw new ParseException("No biological database selected (eg. --go-bp)");
			} else {
				if(list2!=null)idList2 = list2.getDataFrame().getRowNames();
				for(FunctionalFilter filter: filterList) {
					doTest(idList1,idList2,filter,dbConnector);
					}
				
				}
			
			}
		catch (Exception e) {
		}
		
		
	}
	
	private void doTest(List<String> idList1, List<String> idList2, FunctionalFilter filter, DBConnector dbConnector) throws SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException, InvalidParameterException, IOException{
		
		// db attributes
		String name = getDBName(filter);
		String title = getDBTitle(filter);		
		
		//File inputFile = new File(outdir + "/"+name+"_WellFormedInput");
		logger.println("Starting doing test...");
		GeneCodis genecodis = null;
		if(list2!=null) genecodis = new GeneCodis(babelomicsHomePath + "/bin/genecodis/genecodis_bin ",outdir,name,idList1, idList2, filter, dbConnector, duplicatesMode, support, supportRandom,correction, test,analysis);
		else if(isRestOfGenome()) genecodis = new GeneCodis(babelomicsHomePath + "/bin/genecodis/genecodis_bin ",outdir,name,idList1, filter, dbConnector, duplicatesMode, support, supportRandom, correction, test,analysis);
                                                 
		
		//run run...
		genecodis.run();

		//save results
		saveGenecodisResults(genecodis, name, title);
	}
	
	
	private void saveGenecodisResults(GeneCodis genecodis,String name,String title) throws IOException{
		String fileName = name + ".txt";
		String annotFileName = name + ".annot";
		updateJobStatus("80", "saving results");
		
		//IOUtils.write(outdir + "/" +fileName, genecodis.getResults());
		result.addOutputItem(new Item("genecodis_file", fileName, "Genecodis for "+name, TYPE.FILE, Arrays.asList("TABLE","GENECODIS_TABLE"), new HashMap<String, String>(2), "Significant terms"));
			
		// save annotation
		IOUtils.write(outdir + "/" + annotFileName, genecodis.getAnnotations().toString());
		result.addOutputItem(new Item("annot_" + name,annotFileName,"Annotations for " + title,Item.TYPE.FILE,Arrays.asList("ANNOTATION"),new HashMap<String,String>(),"Annotation files"));
					
	}
	
}
