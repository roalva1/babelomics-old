package org.bioinfo.babelomics.tools.functional;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.ParseException;
import org.bioinfo.babelomics.methods.functional.TwoListFisherTestResult;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.dataset.FeatureData;
import org.bioinfo.data.list.exception.InvalidIndexException;
import org.bioinfo.infrared.common.feature.FeatureList;
import org.bioinfo.infrared.funcannot.AnnotationItem;
import org.bioinfo.infrared.funcannot.filter.BiocartaFilter;
import org.bioinfo.infrared.funcannot.filter.Filter;
import org.bioinfo.infrared.funcannot.filter.FunctionalFilter;
import org.bioinfo.infrared.funcannot.filter.GOFilter;
import org.bioinfo.infrared.funcannot.filter.InterproFilter;
import org.bioinfo.infrared.funcannot.filter.JasparFilter;
import org.bioinfo.infrared.funcannot.filter.KeggFilter;
import org.bioinfo.infrared.funcannot.filter.OregannoFilter;
import org.bioinfo.infrared.funcannot.filter.ReactomeFilter;
import org.bioinfo.math.stats.inference.FisherExactTest;
import org.bioinfo.tool.OptionFactory;


public abstract class FunctionalProfilingTool extends BabelomicsTool {

//	// JT.2009.09.07
//	private String removeDuplicates;
	
	// input data		
	protected boolean restOfGenome;

	// test
	protected int testMode;
		
	// dbs
	protected List<FunctionalFilter> filterList;
	
	// your annotations
	protected boolean isYourAnnotations;
	protected FeatureList<AnnotationItem> yourAnnotations;

	@Override
	public void initOptions() {

		// commons options
		getOptions().addOption(OptionFactory.createOption("fisher", "the Fisher test mode, valid values: less, greater, two_sided. By default, two_sided", false));
		
		// GO biological process options
		addGOOptions("bp");
		addGOOptions("cc");
		addGOOptions("mf");

		addGenericOptions("interpro");
		addGenericOptions("kegg");
		addGenericOptions("reactome");
		addGenericOptions("biocarta");
		addGenericOptions("jaspar");
		addGenericOptions("oreganno");
		
//		// kegg
//		getOptions().addOption(OptionFactory.createOption("kegg", "Kegg database",false,false));
//		getOptions().addOption(OptionFactory.createOption("kegg-min-num-genes", "Kegg min number of genes",false));
//		getOptions().addOption(OptionFactory.createOption("kegg-max-num-genes", "Kegg max number of genes",false));
//		getOptions().addOption(OptionFactory.createOption("kegg-count-genes-from-genome", "computes the number of annotated genes from all genome, otherwise from you input list",false,false));
//		
//		// biocarta
//		getOptions().addOption(OptionFactory.createOption("biocarta", "Biocarta database",false,false));
//		getOptions().addOption(OptionFactory.createOption("biocarta-min-num-genes", "Biocarta min number of genes",false));
//		getOptions().addOption(OptionFactory.createOption("biocarta-max-num-genes", "Biocarta max number of genes",false));
//		getOptions().addOption(OptionFactory.createOption("biocarta-count-genes-from-genome", "computes the number of annotated genes from all genome, otherwise from you input list",false,false));
//
		// your annotations
		getOptions().addOption(OptionFactory.createOption("annotations", "Your own annotations",false,true));
//		getOptions().addOption(OptionFactory.createOption("annotation-enrichment", "rest of databases will used to enrich your annotations",false,false));
		
	}


	private void addGOOptions(String namespace){
		String namespaceTitle = "biological process";
		if(namespace.equalsIgnoreCase("cc")) namespaceTitle = "cellular component";
		if(namespace.equalsIgnoreCase("mf")) namespaceTitle = "molecular function";
		getOptions().addOption(OptionFactory.createOption("go-" + namespace, "GO " + namespaceTitle + " database",false,false));
		getOptions().addOption(OptionFactory.createOption("go-" + namespace + "-inclusive", "GO " + namespaceTitle + ", inclusive analysis (one per GO level), otherwise joins GO levels",false,false));
		getOptions().addOption(OptionFactory.createOption("go-" + namespace + "-min-level", "GO " + namespaceTitle + ", min go level to take into account, default 3",false));
		getOptions().addOption(OptionFactory.createOption("go-" + namespace + "-max-level", "GO " + namespaceTitle + ", max GO level to take into account, default 15",false));
		getOptions().addOption(OptionFactory.createOption("go-" + namespace + "-min-num-genes", "GO " + namespaceTitle + ", min number of genes filter",false));
		getOptions().addOption(OptionFactory.createOption("go-" + namespace + "-max-num-genes", "GO " + namespaceTitle + ", max number of genes filter",false));
		getOptions().addOption(OptionFactory.createOption("go-" + namespace + "-all-genome", "GO " + namespaceTitle + ", computes the number of annotated genes from all genome, otherwise from you input list",false,false));
		getOptions().addOption(OptionFactory.createOption("go-" + namespace + "-keywords", "GO " + namespaceTitle + ", keywords filter",false));
		getOptions().addOption(OptionFactory.createOption("go-" + namespace + "-keywords-logic", "GO " + namespaceTitle + ", keywords filter logic: all or any",false));
	}
	
	protected void addGenericOptions(String db){
		getOptions().addOption(OptionFactory.createOption(db, db, false, false));
		getOptions().addOption(OptionFactory.createOption(db + "-min-num-genes", db + " min number of genes filter",false));
		getOptions().addOption(OptionFactory.createOption(db + "-max-num-genes", db + " max number of genes filter",false));
		getOptions().addOption(OptionFactory.createOption(db + "-nannot-domain", db + " computes the number of annotated genes from all genome, otherwise from you input list",false,false));
	}
	
	public void prepare() throws IOException, ParseException, InvalidIndexException{

		// species
		setSpecies(commandLine.getOptionValue("species","hsa"));
		
		// fisher test
		String testMode = commandLine.getOptionValue("test-mode", "two-tailed");
		if(testMode.equals("less")) setTestMode(FisherExactTest.LESS);
		if(testMode.equals("greater")) setTestMode(FisherExactTest.GREATER);
		if(testMode.equals("two-tailed")) setTestMode(FisherExactTest.TWO_SIDED);
		
		filterList = new ArrayList<FunctionalFilter>();
		
		// go bp
		if(commandLine.hasOption("go-bp")) {
			parseGODb(commandLine,"bp");
		}
		// go cc
		if(commandLine.hasOption("go-cc")) {
			parseGODb(commandLine,"cc");
		}
		// go mf
		if(commandLine.hasOption("go-mf")) {
			parseGODb(commandLine,"mf");
		}

//		// kegg
//		if(commandLine.hasOption("kegg")) {			
//			KeggFilter keggFilter = new KeggFilter(Integer.parseInt(commandLine.getOptionValue("kegg-min-num-genes","5")),Integer.parseInt(commandLine.getOptionValue("kegg-max-num-genes","500")));
//			filterList.add(keggFilter);
//		}
//		// biocarta
//		if(commandLine.hasOption("biocarta")) {
//			BiocartaFilter biocartaFilter = new BiocartaFilter(Integer.parseInt(commandLine.getOptionValue("biocarta-min-num-genes","5")),Integer.parseInt(commandLine.getOptionValue("biocarta-max-num-genes","500")));
//			filterList.add(biocartaFilter);
//		}
		parseGenericDb(commandLine,"interpro");
		parseGenericDb(commandLine,"kegg");
		parseGenericDb(commandLine,"reactome");		
		parseGenericDb(commandLine,"biocarta");
		parseGenericDb(commandLine,"jaspar");
		parseGenericDb(commandLine,"oreganno");
		

		
		// species must be provided if any db is selected
		if(commandLine.hasOption("go-bp") || commandLine.hasOption("go-mf") || commandLine.hasOption("go-cc") || commandLine.hasOption("kegg") || commandLine.hasOption("biocarta")){
			if(!commandLine.hasOption("species")) throw new ParseException("species is mandatory if any database specified");
		}
		
		// your annotations
		if(commandLine.hasOption("annotations") && !commandLine.getOptionValue("annotations").equalsIgnoreCase("") && !commandLine.getOptionValue("annotations").equalsIgnoreCase("none")) {
			isYourAnnotations = true;			
			FeatureData annotations = new FeatureData(new File(commandLine.getOptionValue("annotations")), true);
			List<String> ids = annotations.getDataFrame().getRowNames();
			List<String> terms = annotations.getDataFrame().getColumn(0);
			yourAnnotations = new FeatureList<AnnotationItem>(ids.size());// FeatureData(new File(cmdLine.getOptionValue("annotations")), true);
			for(int i=0; i<ids.size(); i++){
				yourAnnotations.add(new AnnotationItem(ids.get(i),terms.get(i)));
			}
		}
		
	}

	public void parseGODb(CommandLine cmdLine, String namespace){
		if(cmdLine.hasOption("go-" + namespace + "")) {
			if(cmdLine.hasOption("go-" + namespace + "-inclusive")) {				
				int min = Integer.parseInt(cmdLine.getOptionValue("go-" + namespace + "-min-level","5"));
				int max = Integer.parseInt(cmdLine.getOptionValue("go-" + namespace + "-min-level","12"));
				for(int i=min; i<=max; i++){
					GOFilter goBpFilter = parseGOFilter(cmdLine, namespace);
					goBpFilter.setMinLevel(i);
					goBpFilter.setMaxLevel(i);					
					filterList.add(goBpFilter);	
				}
			} else {
				GOFilter goBpFilter = parseGOFilter(cmdLine, namespace);				
				filterList.add(goBpFilter);
			}			
		}		
	}
	
	public GOFilter parseGOFilter(CommandLine cmdLine, String namespace){
		// ontology
		String infraredNamespace = "biological_process";		
		if(namespace.equalsIgnoreCase("cc")) infraredNamespace = "cellular_component";
		if(namespace.equalsIgnoreCase("mf")) infraredNamespace = "molecular_function";
		GOFilter goFilter = new GOFilter(infraredNamespace);
		// levels
		goFilter.setMinLevel(Integer.parseInt(cmdLine.getOptionValue("go-" + namespace + "-min-level","5")));
		goFilter.setMaxLevel(Integer.parseInt(cmdLine.getOptionValue("go-" + namespace + "-max-level","12")));
		// number of annotations
		goFilter.setMinNumberGenes(Integer.parseInt(cmdLine.getOptionValue("go-" + namespace + "-min-num-genes","5")));	
		goFilter.setMaxNumberGenes(Integer.parseInt(cmdLine.getOptionValue("go-" + namespace + "-max-num-genes","500")));
		
		// filter by keywords
		if(cmdLine.hasOption("go-" + namespace + "-keywords")){
			goFilter.addKeywords(StringUtils.toList(cmdLine.getOptionValue("go-" + namespace + "-keywords", "")));
			goFilter.setLogicalOperator(cmdLine.getOptionValue("go-" + namespace + "-keywords-logic","AND"));			
		}
		return goFilter;
	}
	
	public void parseGenericDb(CommandLine cmdLine, String db){
		if(commandLine.hasOption(db)) {
			FunctionalFilter filter = null;
			
			if(db.equalsIgnoreCase("interpro")){
				filter = new InterproFilter();		
			}
			if(db.equalsIgnoreCase("kegg")){
				filter = new KeggFilter();						
			}
			if(db.equalsIgnoreCase("reactome")){
				filter = new ReactomeFilter();						
			}
			if(db.equalsIgnoreCase("biocarta")){
				filter = new BiocartaFilter();						
			}
			if(db.equalsIgnoreCase("jaspar")){
				filter = new JasparFilter();						
			}
			if(db.equalsIgnoreCase("oreganno")){
				filter = new OregannoFilter();						
			}
			if(filter!=null){
				filter.setMinNumberGenes(Integer.parseInt(cmdLine.getOptionValue(db + "-min-num-genes","2")));	
				filter.setMaxNumberGenes(Integer.parseInt(cmdLine.getOptionValue(db + "-max-num-genes","500")));		
				if(cmdLine.getOptionValue(db + "-nannot-domain","genome").equals("genome")) filter.setGenomicNumberOfGenes(true); 
				else filter.setGenomicNumberOfGenes(false);
				filterList.add(filter);
			}
			
		}
	}

	protected String getDBName(Filter filter){
		String name = StringUtils.randomString();
		if(filter instanceof GOFilter) {						
			GOFilter goFilter = (GOFilter) filter.clone();
			name = "go_" + goFilter.getNamespace() + "_" + goFilter.getMinLevel() + "_" + goFilter.getMaxLevel();
		} else if(filter instanceof KeggFilter) {
			name = "kegg";
		} else if(filter instanceof InterproFilter) {
			name = "interPro";
		} else if(filter instanceof ReactomeFilter) {
			name = "reactome";
		} else if(filter instanceof BiocartaFilter) {
			name = "biocarta";
		} else if(filter instanceof JasparFilter) {
			name = "jaspar";
		} else if(filter instanceof OregannoFilter) {
			name = "oreganno";
		}
		return name;
	}
	
	protected String getDBTitle(Filter filter){
		String title = "Untitled",levels;
		if(filter instanceof GOFilter) {						
			GOFilter goFilter = (GOFilter) filter.clone();
			if(goFilter.getMinLevel()==goFilter.getMaxLevel()) {
				levels = "(level " + goFilter.getMinLevel() + ")";
			} else{
				levels = "(levels from " + goFilter.getMinLevel() + " to " + goFilter.getMaxLevel() + ")"; 
			}						
			title = "GO " + goFilter.getNamespace().replace("_", " ") + " " + levels;
		} else if(filter instanceof KeggFilter) {
			title = "Kegg";
		} else if(filter instanceof InterproFilter) {
			title = "InterPro";
		}  else if(filter instanceof ReactomeFilter) {
			title = "Reactome";
		}	else if(filter instanceof BiocartaFilter) {
			title = "Biocarta";
		}	else if(filter instanceof JasparFilter) {
			title = "Jaspar";
		}	else if(filter instanceof OregannoFilter) {
			title = "Oreganno";
		}	
		return title;
	}

	
	protected List<String> testResultToStringList(List<TwoListFisherTestResult> testResult){
		return testResultToStringList(testResult,true);
	}
	
	protected List<String> testResultToStringList(List<TwoListFisherTestResult> testResult, boolean header){
		List<String> result = new ArrayList<String>();
		if(header) result.add(TwoListFisherTestResult.header());
		for(int i=0; i<testResult.size(); i++){			
			result.add(testResult.get(i).toString());
		}
		return result;
	}
		
	
	public void setRestOfGenome(boolean restOfGenome) {
		this.restOfGenome = restOfGenome;
	}

	public boolean isRestOfGenome() {
		return restOfGenome;
	}

	public void setTestMode(int testMode) {
		this.testMode = testMode;
	}

	public int getTestMode() {
		return testMode;
	}


}
