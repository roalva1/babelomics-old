package org.bioinfo.babelomics.tools.functional;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import org.apache.commons.cli.OptionGroup;
import org.apache.commons.cli.ParseException;
import org.apache.commons.math.genetics.Chromosome;
import org.bioinfo.babelomics.methods.functional.FatiGO;
import org.bioinfo.babelomics.methods.functional.TwoListFisherTestResult;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.data.dataset.FeatureData;
import org.bioinfo.data.list.exception.InvalidIndexException;
import org.bioinfo.infrared.common.dbsql.DBConnector;
import org.bioinfo.infrared.common.feature.FeatureList;
import org.bioinfo.infrared.funcannot.AnnotationItem;
import org.bioinfo.infrared.funcannot.filter.FunctionalFilter;
import org.bioinfo.math.exception.InvalidParameterException;
import org.bioinfo.math.stats.inference.FisherExactTest;
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
	private StringBuilder duplicatesReport;
	private StringBuilder annotationReport;
	private String list2label;
	
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
			input2.addOption(OptionFactory.createOption("chromosomes", "the chromosome region to compare with list #1. Use chromosome-start and chromosome-end options tp delimite the region within the chromosome",false));
		options.addOptionGroup(input2);
		// filters
		
		// extras
		options.addOption(OptionFactory.createOption("duplicates", "to remove duplicated IDs, valid values: never, each, ref. By default, never", false));
	}

	
	@Override
	public void prepare() throws IOException, ParseException, InvalidIndexException {
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
		}
//		else {
//			throw new ParseException("No comparison provided, use list2, genome or chromosome options to set your comparison");
//		}
		
		// extras
		String duplicates = commandLine.getOptionValue("duplicates", "never");		
		if(duplicates.equalsIgnoreCase("each")) duplicatesMode = FatiGO.REMOVE_EACH;
		if(duplicates.equalsIgnoreCase("ref")) duplicatesMode = FatiGO.REMOVE_REF;
		if(duplicates.equalsIgnoreCase("all")) duplicatesMode = FatiGO.REMOVE_ALL;
		System.err.println("los duplicates son: " + duplicatesMode);
		// your annotations
		if(isYourAnnotations){
			
		}
	}

	@Override
	public void execute() {
		try {
			logger.println("Starting FatiGO...");
			
			// update status
			jobStatus.addStatusMessage("10", "Preparing data");
			
			// infrared connector
			logger.debug("species: " + getSpecies());
			DBConnector dbConnector = new DBConnector(species, new File(babelomicsHomePath + "/conf/infrared.properties"));
			
			// prepare params
			prepare();			
	
			FatiGO fatigo = null;
			list2label = "List 2";
			
			// list 1			
			System.err.println(list1.getDataFrame());
			List<String> idList1 = list1.getDataFrame().getRowNames();//list1.getDataFrame().getColumn(0); //InfraredUtils.toEnsemblId(dbConnector, list1.getDataFrame().getColumn(0));
			
			// list 2
			List<String> idList2 = null;
			if(list2!=null) {
				idList2 = list2.getDataFrame().getRowNames();//list2.getDataFrame().getColumn(0); //InfraredUtils.toEnsemblId(dbConnector, list2.getDataFrame().getColumn(0));
				fatigo = new FatiGO(idList1, idList2, null, dbConnector, testMode, duplicatesMode);
			} else if(isRestOfGenome()) {				
				duplicatesMode = FatiGO.REMOVE_GENOME;
				fatigo = new FatiGO(idList1, null, dbConnector);
				list2label = "Genome";
			} else if(chromosomes!=null) {
				throw new ParseException("chromosomes comparison not yet implemented");
			} else {
				fatigo = new FatiGO(idList1, yourAnnotations);
				list2label = "Rest of ids from your annotations (complementary list)";
			}
//			else {
//				throw new ParseException("No comparison provided, use list2, genome or chromosomes options to set your comparison");				
//			}
			
			addInputParams();
			
			// save data
			IOUtils.write(outdir + "/clean_list1.txt", ListUtils.toString(fatigo.getList1(),"\n"));
			result.addOutputItem(new Item("clean_list1","clean_list1.txt","List 1 (after duplicates managing)",Item.TYPE.FILE,Arrays.asList("IDLIST","CLEAN"),new HashMap<String,String>(),"Input data"));
			IOUtils.write(outdir + "/clean_list2.txt", ListUtils.toString(fatigo.getList2(),"\n"));
			result.addOutputItem(new Item("clean_list2","clean_list2.txt",list2label + " (after duplicates managing)",Item.TYPE.FILE,Arrays.asList("IDLIST","CLEAN"),new HashMap<String,String>(),"Input data"));
			
			
			// run fatigo's
			if(filterList.size()==0 && !isYourAnnotations){
				throw new ParseException("No biological database selected (eg. --go-bp)");
			} else {

				// significant terms
				List<String> significant = new ArrayList<String>();
				
				// Significant results must appear after than complete tables!!
				result.addOutputItem(new Item("significant","significant_" + DEFAULT_PVALUE_THRESHOLD + ".txt","Significant terms",Item.TYPE.FILE,Arrays.asList("TABLE","FATIGO_TABLE"),new HashMap<String,String>(),"Significant Results"));
				significant.add(TwoListFisherTestResult.header());
				
				annotationReport = new StringBuilder();
				annotationReport.append("#DB").append("\t").append("List1").append("\t").append(list2label).append("\n");
				
				// run fatigo's				
				for(FunctionalFilter filter: filterList) {
					fatigo = doFatigo(idList1,idList2,filter,dbConnector,significant);
					addAnnotationReport(fatigo,getDBTitle(filter));
				}				
				if(isYourAnnotations){					
					fatigo = doFatigoYourAnnotations(idList1, idList2, yourAnnotations, significant);
					addAnnotationReport(fatigo,"Your annotations");
				}
				
				saveDuplicatesReport(fatigo);
				saveAnnotationsReport();
				
				// significant terms
				IOUtils.write(outdir + "/significant_" + DEFAULT_PVALUE_THRESHOLD + ".txt", ListUtils.toString(significant,"\n"));
				
			}
			
		}catch(ParseException pe){
			logger.error(pe.getMessage());
		}catch(Exception e){
			e.printStackTrace();
		}
		
	}

	private void addAnnotationReport(FatiGO fatigo, String dbName){
		DecimalFormat formatter = new DecimalFormat("#######.##");
		double list1Percentage = ((double)(fatigo.getList1AnnotatedCounter())/(double)fatigo.getList1SizeBeforeDuplicates())*100.0;
		double list2Percentage = ((double)(fatigo.getList2AnnotatedCounter())/(double)fatigo.getList2SizeBeforeDuplicates())*100.0;
		String list1Message = fatigo.getList1AnnotatedCounter() + " of " + fatigo.getList1SizeAfterDuplicates() + " (" + formatter.format(list1Percentage) + "%, an average of " + formatter.format(fatigo.getList1MeanAnnotationsPerId()) + " per id)";
		String list2Message = fatigo.getList2AnnotatedCounter() + " of " + fatigo.getList2SizeAfterDuplicates() + " (" + formatter.format(list2Percentage) +"%, an average of " + formatter.format(fatigo.getList2MeanAnnotationsPerId()) + " per id)";
		annotationReport.append(dbName).append("\t").append(list1Message).append("\t").append(list2Message).append("\n");
	}
	private void saveAnnotationsReport() throws IOException{
		IOUtils.write(outdir + "/annotations_per_db.txt", annotationReport.toString());
		result.getOutputItems().add(1,new Item("annotations_per_db","annotations_per_db.txt","Id annotations per DB",Item.TYPE.FILE,Arrays.asList("TABLE","SUMMARY_TABLE"),new HashMap<String,String>(),"Summary"));
	}
	
	
	private void addInputParams(){
		// species
		addInputParam("species", "Species", species);
		// duplicates
		HashMap<Integer,String> duplicates = new HashMap<Integer, String>();
		duplicates.put(FatiGO.REMOVE_NEVER, "Never remove");
		duplicates.put(FatiGO.REMOVE_EACH, "Remove separately on each list");
		duplicates.put(FatiGO.REMOVE_REF, "Remove list 1 ids from list 2 (complementary list)");
		duplicates.put(FatiGO.REMOVE_GENOME, "Remove list 1 ids from genome");
		duplicates.put(FatiGO.REMOVE_ALL, "Remove all duplicates (owned and shared duplicates)");
		addInputParam("duplicates", "Duplicates management", duplicates.get(duplicatesMode));
		// fisher laterality
		HashMap<Integer,String> fisher = new HashMap<Integer, String>();		
		fisher.put(FisherExactTest.GREATER, "Over represented terms in list 1");
		fisher.put(FisherExactTest.LESS, "Over represented terms in list 2");
		fisher.put(FisherExactTest.TWO_SIDED, "Two tailed");
		addInputParam("fisher", "Fisher exact test", fisher.get(testMode));
	}	
	private void addInputParam(String id, String label, String value){
		result.addOutputItem(new Item(id,value,label,Item.TYPE.MESSAGE,Arrays.asList("INPUT_PARAM"),new HashMap<String,String>(),"Input params"));
	}
	
	private void saveDuplicatesReport(FatiGO fatigo) throws IOException{
		DecimalFormat formatter = new DecimalFormat("##.##");
		// list 1
		int list1Duplicates = fatigo.getList1SizeBeforeDuplicates()-fatigo.getList1SizeAfterDuplicates();
		double list1DuplicatesPercentage = (double)(list1Duplicates)/(double)fatigo.getList1SizeBeforeDuplicates();
		String list1DuplicatesMessage = list1Duplicates + " of " + fatigo.getList1SizeBeforeDuplicates() + " (" + formatter.format(list1DuplicatesPercentage) + "%)";
		// list 2
		int list2Duplicates = fatigo.getList2SizeBeforeDuplicates()-fatigo.getList2SizeAfterDuplicates();
		double list2DuplicatesPercentage = (double)(list2Duplicates)/(double)fatigo.getList2SizeBeforeDuplicates();
		String list2DuplicatesMessage = list2Duplicates + " of " + fatigo.getList2SizeBeforeDuplicates() + " (" + formatter.format(list2DuplicatesPercentage) + "%)";
		// report
		duplicatesReport = new StringBuilder();
		duplicatesReport.append("#Detail").append("\t").append("List 1").append("\t").append(list2label).append("\n");		
		duplicatesReport.append("Number of duplicates").append("\t").append(list1DuplicatesMessage).append("\t").append(list2DuplicatesMessage).append("\n");
		duplicatesReport.append("Number of finally used ids").append("\t").append(fatigo.getList1SizeAfterDuplicates()).append("\t").append(fatigo.getList2SizeAfterDuplicates()).append("\n");		
		IOUtils.write(outdir + "/duplicates.txt", duplicatesReport.toString());
		result.getOutputItems().add(1,new Item("duplicates","duplicates.txt","Duplicates management",Item.TYPE.FILE,Arrays.asList("TABLE","SUMMARY_TABLE"),new HashMap<String,String>(),"Summary"));		
	}
		
	private FatiGO doFatigo(List<String> idList1, List<String> idList2,FunctionalFilter filter,DBConnector dbConnector, List<String> significant) throws IOException, SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException, InvalidParameterException{
		
		// db attributes
		FunctionalDbDescriptor filterInfo= new FunctionalDbDescriptor(filter);
//		String name = getDBName(filter);
//		String title = getDBTitle(filter);
				
		logger.println(filterInfo.getTitle() + "...");

		// init test
		FatiGO fatigo = null;
		if(list2!=null) fatigo = new FatiGO(idList1, idList2, filter, dbConnector, testMode, duplicatesMode);
		else if(isRestOfGenome()) fatigo = new FatiGO(idList1, filter, dbConnector);
		
		fatigo.setLogger(logger);
		
		// run test
		fatigo.run();
		
		// save results
		saveFatigoResults(fatigo,filterInfo);
		
		// acum significant values
		significant.addAll(testResultToStringList(fatigo.getSignificant(DEFAULT_PVALUE_THRESHOLD),false));

		//logger.println("...end of " + title);
		logger.println("...finished");
		
		return fatigo;
		
	}
	
	private FatiGO doFatigoYourAnnotations(List<String> idList1, List<String> idList2,FeatureList<AnnotationItem> yourAnnotations, List<String> significant) throws IOException, SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException, InvalidParameterException{
		
		// db attributes
		FunctionalDbDescriptor filterInfo = new FunctionalDbDescriptor("your_annotation","Your annotations", "your_annotations","Your annotations");
//		String name = "your_annotations";
//		String title = "Your annotations";					
				
		logger.println(filterInfo.getTitle() + "...");

		// init test
		FatiGO fatigo = null;
		if(list2!=null) fatigo = new FatiGO(idList1, idList2, yourAnnotations, testMode, duplicatesMode);
		//else if(isRestOfGenome()) fatigo = new FatiGO(idList1,yourAnnotations);
		else fatigo = new FatiGO(idList1,yourAnnotations);
		
		// run test
		fatigo.run();
		
		// save results
		saveFatigoResults(fatigo,filterInfo);
	
		// acum significant values
		significant.addAll(testResultToStringList(fatigo.getSignificant(DEFAULT_PVALUE_THRESHOLD),false));
		
		//logger.println("...end of " + title);
		logger.println("...finished");
		
		return fatigo;
	}
		
	private void saveFatigoResults(FatiGO fatigo, FunctionalDbDescriptor filterInfo) throws IOException{
		
		String fileName = filterInfo.getName() + ".txt";
		String annotFileName = filterInfo.getName() + ".annot";
		
		// save statistic results					
		List<String> testResultOutput = testResultToStringList(fatigo.getResults());
		IOUtils.write(outdir + "/" + fileName, ListUtils.toString(testResultOutput,"\n"));
		result.addOutputItem(new Item(filterInfo.getName(),fileName,filterInfo.getTitle(),Item.TYPE.FILE,Arrays.asList("TABLE","FATIGO_TABLE",filterInfo.getPrefix().toUpperCase() + "_TERM"),new HashMap<String,String>(),"Database tests"));
		result.addOutputItem(new Item(filterInfo.getName() + "_description",filterInfo.getDescription(),filterInfo.getTitle(),Item.TYPE.MESSAGE,Arrays.asList("MINI_COMMENT"),new HashMap<String,String>(),"Database tests"));
						
		// save annotation
		IOUtils.write(outdir + "/" + annotFileName, fatigo.getAnnotations().toString());
		result.addOutputItem(new Item("annot_" + filterInfo.getName(),annotFileName,"Annotations for " + filterInfo.getTitle(),Item.TYPE.FILE,Arrays.asList("ANNOTATION"),new HashMap<String,String>(),"Annotation files"));
				
	}
	
}
