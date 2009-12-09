package org.bioinfo.babelomics.tools.functional;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import org.apache.commons.cli.OptionGroup;
import org.apache.commons.cli.ParseException;
import org.apache.commons.math.genetics.Chromosome;
import org.bioinfo.babelomics.methods.functional.FatiGO;
import org.bioinfo.babelomics.methods.functional.TwoListFisherTestResult;
import org.bioinfo.collections.exceptions.InvalidColumnIndexException;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.dataset.FeatureData;
import org.bioinfo.infrared.common.dbsql.DBConnector;
import org.bioinfo.infrared.common.feature.FeatureList;
import org.bioinfo.infrared.funcannot.AnnotationItem;
import org.bioinfo.infrared.funcannot.filter.Filter;
import org.bioinfo.math.exception.InvalidParameterException;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;

public class FatiGOTool extends FunctionalProfilingTool{

	public final static double DEFAULT_PVALUE_THRESHOLD = 0.05;
	
	// list1
	private FeatureData list1;
	
	// list2
	private FeatureData list2;		
	private List<Chromosome> chromosomes;
	
	// other capabilities
	protected int duplicatesMode;
	
	
	public FatiGOTool() {
		initOptions();
	}

	@Override
	public void initOptions() {
		// parent options
		super.initOptions();
		// list 1
		options.addOption(OptionFactory.createOption("list1", "the feature data containig the list #1 of genes, or the feature data file"));
		// list 2
		OptionGroup input2 = new OptionGroup();
		input2.setRequired(false);		
			input2.addOption(OptionFactory.createOption("list2", "the file containig the list #2 of genes, or the feature data file",false));
			input2.addOption(OptionFactory.createOption("genome", "compares list #1 with the rest of genome",false,false));
			input2.addOption(OptionFactory.createOption("chromosomes", "the chromosome to compare with list #1. Use chromosome-start and chromosome-end options tp delimite the region within the chromosome",false));
		options.addOptionGroup(input2);
		// filters
		
		// extras
		options.addOption(OptionFactory.createOption("duplicates", "to remove duplicated IDs, valid values: never, each, ref. By default, never", false));
	}

	
	@Override
	public void prepare() throws IOException, ParseException, InvalidColumnIndexException {
		super.prepare();

		// list 1		
		list1 = new FeatureData(new File(commandLine.getOptionValue("list1")), true);
		
		// list 2
		if(commandLine.hasOption("list2")){
			list2 = new FeatureData(new File(commandLine.getOptionValue("list2")), true);
		} else if(commandLine.hasOption("genome")) {
			this.setRestOfGenome(true);
		} else if(commandLine.hasOption("chromosomes")) {
			throw new ParseException("chromosome reading not yet implemented");
		} else {
			throw new ParseException("No comparison provided, use list2, rest-of-genome or chromosome options to set your comparison");
		}
		
		// extras
		String duplicates = commandLine.getOptionValue("remove-duplicates", "never");
		if(duplicates.equalsIgnoreCase("each")) duplicatesMode = FatiGO.REMOVE_EACH;
		if(duplicates.equalsIgnoreCase("ref")) duplicatesMode = FatiGO.REMOVE_REF;
		if(duplicates.equalsIgnoreCase("all")) duplicatesMode = FatiGO.REMOVE_ALL;		
		
		// your annotations
		if(isYourAnnotations){
			
		}
	}

	@Override
	public void execute() {
		try {
			logger.debug("EXECUTING FATIGOOOOO!!!!");
			
			// infrared connector
			DBConnector dbConnector = new DBConnector(getSpecies(), new File(System.getenv("BABELOMICS_HOME") + "/conf/infrared.conf"));			
			
//			System.out.println("FatiGOTool.java, execute: getSpecies() = " + getSpecies());
//			System.out.println("FatiGOTool.java, execute: from file " + new File(System.getenv("BABELOMICS_HOME") + "/conf/infrared.conf").getAbsolutePath() + ", dbConnector = " + dbConnector);

			// prepare params
			prepare();			
	
			// list 1			
			List<String> idList1 = list1.getDataFrame().getColumn(0); //InfraredUtils.toEnsemblId(dbConnector, list1.getDataFrame().getColumn(0));
			
			// list 2
			List<String> idList2 = null;
			if(list2!=null) {
				idList2 = list2.getDataFrame().getColumn(0); //InfraredUtils.toEnsemblId(dbConnector, list2.getDataFrame().getColumn(0));
			} else if(isRestOfGenome()) {				
				duplicatesMode = FatiGO.REMOVE_GENOME;
			} else if(chromosomes!=null) {
				throw new ParseException("chromosome comparison not yet implemented");
			} else {
				throw new ParseException("No comparison provided, use list2, genome or chromosomes options to set your comparison");				
			}
			
			// run fatigo's
			if(filterList.size()==0 && !isYourAnnotations){
				throw new ParseException("No biological database selected (eg. --go-bp)");
			} else {

				// significant terms
				List<String> significant = new ArrayList<String>();

				FatiGO fatigo = new FatiGO(idList1, idList2, null, dbConnector, testMode, duplicatesMode);
				IOUtils.write(outdir + "/clean_list1.txt", ListUtils.toString(fatigo.getList1(),"\n"));
				result.addOutputItem(new Item("clean_list1","clean_list1.txt","List 1 (after duplicates managing)",Item.TYPE.FILE,Arrays.asList("IDLIST","CLEAN"),new HashMap<String,String>(),"Input data"));
				IOUtils.write(outdir + "/clean_list2.txt", ListUtils.toString(fatigo.getList2(),"\n"));
				result.addOutputItem(new Item("clean_list2","clean_list2.txt","List 2 (after duplicates managing)",Item.TYPE.FILE,Arrays.asList("IDLIST","CLEAN"),new HashMap<String,String>(),"Input data"));
				
				// Significant results must appear after than complete tables!!
				result.addOutputItem(new Item("significant","significant_" + DEFAULT_PVALUE_THRESHOLD + ".txt","Significant terms",Item.TYPE.FILE,Arrays.asList("TABLE","FATIGO_TABLE"),new HashMap<String,String>(),"Significant Results"));
				significant.add(TwoListFisherTestResult.header());
				
				// run fatigo's
				for(Filter filter: filterList) {
					doFatigo(idList1,idList2,filter,dbConnector,significant);					
				}				
				if(isYourAnnotations){					
					doFatigoYourAnnotations(idList1, idList2, yourAnnotations, significant);
				}
				
				// significant terms
				IOUtils.write(outdir + "/significant_" + DEFAULT_PVALUE_THRESHOLD + ".txt", StringUtils.join(significant,"\n"));
				
			}
			
		}catch(Exception e){
			e.printStackTrace();
		}
		
	}

	private void doFatigo(List<String> idList1, List<String> idList2,Filter filter,DBConnector dbConnector, List<String> significant) throws IOException, SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException, InvalidParameterException{
		
		// db attributes
		String name = getDBName(filter);
		String title = getDBTitle(filter);					
				
		logger.info(title + "...\n");

		// init test
		FatiGO fatigo = null;
		if(list2!=null) fatigo = new FatiGO(idList1, idList2, filter, dbConnector, testMode, duplicatesMode);
		else if(isRestOfGenome()) fatigo = new FatiGO(idList1, filter, dbConnector);
		
		// run test
		fatigo.run();
		
		// save results
		saveFatigoResults(fatigo,name,title);
		
		// acum significant values
		significant.addAll(testResultToStringList(fatigo.getSignificant(DEFAULT_PVALUE_THRESHOLD),false));

		logger.info("...end of " + title);
	}
	
	private void doFatigoYourAnnotations(List<String> idList1, List<String> idList2,FeatureList<AnnotationItem> yourAnnotations, List<String> significant) throws IOException, SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException, InvalidParameterException{
		
		// db attributes
		String name = "your_annotations";
		String title = "Your annotations";					
				
		logger.info(title + "...\n");

		// init test
		FatiGO fatigo = null;
		if(list2!=null) fatigo = new FatiGO(idList1, idList2, yourAnnotations, testMode, duplicatesMode);
		else if(isRestOfGenome()) fatigo = new FatiGO(idList1,yourAnnotations);
		
		// run test
		fatigo.run();
		
		// save results
		saveFatigoResults(fatigo,name,title);
	
		// acum significant values
		significant.addAll(testResultToStringList(fatigo.getSignificant(DEFAULT_PVALUE_THRESHOLD),false));
		
		logger.info("...end of " + title);
		
	}
	
	private void saveFatigoResults(FatiGO fatigo,String name,String title) throws IOException{
		
		String fileName = name + ".txt";
		String annotFileName = name + ".annot";
		
		// save statistic results					
		List<String> testResultOutput = testResultToStringList(fatigo.getResults());
		IOUtils.write(outdir + "/" + fileName, StringUtils.join(testResultOutput,"\n"));
		result.addOutputItem(new Item(name,fileName,title,Item.TYPE.FILE,Arrays.asList("TABLE","FATIGO_TABLE"),new HashMap<String,String>(),"Database tests"));
						
		// save annotation
		IOUtils.write(outdir + "/" + annotFileName, fatigo.getAnnotations().toString());
		result.addOutputItem(new Item("annot_" + name,annotFileName,"Annotations for " + title,Item.TYPE.FILE,Arrays.asList("ANNOTATION"),new HashMap<String,String>(),"Annotation files"));
				
	}
	
}
