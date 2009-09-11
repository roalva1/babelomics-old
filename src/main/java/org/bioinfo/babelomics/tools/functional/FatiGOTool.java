package org.bioinfo.babelomics.tools.functional;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionGroup;
import org.apache.commons.cli.ParseException;
import org.apache.commons.math.genetics.Chromosome;
import org.bioinfo.babelomics.methods.functional.FatiGO;
import org.bioinfo.babelomics.methods.functional.TwoListFisherTestResult;
import org.bioinfo.collections.exceptions.InvalidColumnIndexException;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.dataset.FeatureData;
import org.bioinfo.infrared.common.dbsql.DBConnector;
import org.bioinfo.infrared.funcannot.filter.BiocartaFilter;
import org.bioinfo.infrared.funcannot.filter.Filter;
import org.bioinfo.infrared.funcannot.filter.GOFilter;
import org.bioinfo.infrared.funcannot.filter.KeggFilter;
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
		options.addOption(OptionFactory.createOption("remove-duplicates", "to remove duplicated IDs, valid values: never, each, ref. By default, never", false));
	}

	
	@Override
	public void prepare(CommandLine cmdLine) throws IOException, ParseException, InvalidColumnIndexException {
		super.prepare(cmdLine);

		// list 1		
		list1 = new FeatureData(new File(cmdLine.getOptionValue("list1")), true);
		
		// list 2
		if(cmdLine.hasOption("list2")){
			list2 = new FeatureData(new File(cmdLine.getOptionValue("list2")), true);
		} else if(cmdLine.hasOption("genome")) {
			this.setRestOfGenome(true);
		} else if(cmdLine.hasOption("chromosomes")) {
			throw new ParseException("chromosome reading not yet implemented");
		} else {
			throw new ParseException("No comparison provided, use list2, rest-of-genome or chromosome options to set your comparison");
		}
		
		// extras
		String duplicates = cmdLine.getOptionValue("remove-duplicates", "never");
		if(duplicates.equalsIgnoreCase("each")) duplicatesMode = FatiGO.REMOVE_EACH;
		if(duplicates.equalsIgnoreCase("ref")) duplicatesMode = FatiGO.REMOVE_REF;
		if(duplicates.equalsIgnoreCase("all")) duplicatesMode = FatiGO.REMOVE_ALL;		
		
	}

	@Override
	public void execute() {
		try {
			// infrared connector			
			DBConnector dbConnector = new DBConnector(getSpecies(), new File(System.getenv("BABELOMICS_HOME") + "/conf/infrared.conf"));			
			
//			System.out.println("FatiGOTool.java, execute: getSpecies() = " + getSpecies());
//			System.out.println("FatiGOTool.java, execute: from file " + new File(System.getenv("BABELOMICS_HOME") + "/conf/infrared.conf").getAbsolutePath() + ", dbConnector = " + dbConnector);

			// prepare params
			prepare(commandLine);			
	
			// list 1
			List<String> idList1;
			
			idList1 = list1.getDataFrame().getColumn(0); //InfraredUtils.toEnsemblId(dbConnector, list1.getDataFrame().getColumn(0));
			
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
			if(filterList.size()==0){
				throw new ParseException("No biological database selected (eg. --go-bp)");
			} else {
				List<String> significant = new ArrayList<String>();
				String name,fileName,annotFileName,title;

				// save id lists
				FatiGO fatigo = new FatiGO(idList1, idList2, filterList.get(0), dbConnector, testMode, duplicatesMode);
				IOUtils.write(outdir + "/clean_list1.txt", StringUtils.join(fatigo.getList1(),"\n"));
				result.addOutputItem(new Item("clean_list1","clean_list1.txt","List 1 (after duplicates managing)",Item.TYPE.FILE,Arrays.asList("IDLIST","CLEAN"),new HashMap<String,String>(),"Input data"));
				IOUtils.write(outdir + "/clean_list2.txt", StringUtils.join(fatigo.getList2(),"\n"));
				result.addOutputItem(new Item("clean_list2","clean_list2.txt","List 2 (after duplicates managing)",Item.TYPE.FILE,Arrays.asList("IDLIST","CLEAN"),new HashMap<String,String>(),"Input data"));

				
				// Significant results must appear after than complete tables!!
				result.addOutputItem(new Item("significant","significant_" + DEFAULT_PVALUE_THRESHOLD + ".txt","Significant terms",Item.TYPE.FILE,Arrays.asList("TABLE"),new HashMap<String,String>(),"Significant Results"));
				significant.add(TwoListFisherTestResult.header());
				
				for(Filter filter: filterList) {

					// db attributes
					name = getDBName(filter);
					title = getDBTitle(filter);					
					fileName = name + ".txt";
					annotFileName = name + ".annot";
					
					logger.info(title + "...\n");

					// init test
					fatigo = null;
					if(list2!=null) fatigo = new FatiGO(idList1, idList2, filter, dbConnector, testMode, duplicatesMode);
					else if(isRestOfGenome()) fatigo = new FatiGO(idList1, filter, dbConnector);
					
					// run test
					fatigo.run();
					
					// save statistic results					
					List<String> testResultOutput = testResultToStringList(fatigo.getResults());
					IOUtils.write(outdir + "/" + fileName, StringUtils.join(testResultOutput,"\n"));
					result.addOutputItem(new Item(name,fileName,title,Item.TYPE.FILE,Arrays.asList("TABLE"),new HashMap<String,String>(),"Database tests"));
									
					// save annotation
					IOUtils.write(outdir + "/" + annotFileName, fatigo.getAnnotations().toString());
					result.addOutputItem(new Item("annot_" + name,annotFileName,"Annotations for " + title,Item.TYPE.FILE,Arrays.asList("ANNOTATION"),new HashMap<String,String>(),"Annotation files"));
					
					// acum significant values
					significant.addAll(testResultToStringList(fatigo.getSignificant(DEFAULT_PVALUE_THRESHOLD),false));
					
					logger.info("...end of " + title);
					
				}
				
				// significant terms
				IOUtils.write(outdir + "/significant_" + DEFAULT_PVALUE_THRESHOLD + ".txt", StringUtils.join(significant,"\n"));
				
			}
			
		}catch(Exception e){
			e.printStackTrace();
		}
		
	}

	
	
}
