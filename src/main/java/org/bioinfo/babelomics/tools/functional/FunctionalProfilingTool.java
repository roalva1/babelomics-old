package org.bioinfo.babelomics.tools.functional;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.ParseException;
import org.bioinfo.babelomics.methods.functional.TwoListFisherTestResult;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.collections.exceptions.InvalidColumnIndexException;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.dataset.FeatureData;
import org.bioinfo.infrared.common.dbsql.DBConnector;
import org.bioinfo.infrared.common.feature.FeatureList;
import org.bioinfo.infrared.core.Gene;
import org.bioinfo.infrared.core.dbsql.GeneDBManager;
import org.bioinfo.infrared.core.dbsql.XRefDBManager;
import org.bioinfo.infrared.funcannot.AnnotationItem;
import org.bioinfo.infrared.funcannot.dbsql.AnnotationDBManager;
import org.bioinfo.infrared.funcannot.filter.BiocartaFilter;
import org.bioinfo.infrared.funcannot.filter.Filter;
import org.bioinfo.infrared.funcannot.filter.GOFilter;
import org.bioinfo.infrared.funcannot.filter.KeggFilter;
import org.bioinfo.infrared.funcannot.filter.Keywords;
import org.bioinfo.math.stats.inference.FisherExactTest;
import org.bioinfo.tool.OptionFactory;


public abstract class FunctionalProfilingTool extends BabelomicsTool {

//	// JT.2009.09.07
//	private String removeDuplicates;
	
	// input data		
	protected boolean restOfGenome;

	// test
	protected int testMode;
		
// GO biological process
//	protected List<GOFilter> goBpFilterList;
//	protected List<GOFilter> goMfFilterList;
//	protected List<GOFilter> goCcFilterList;
//	protected KeggFilter keggFilter;
//	protected BiocartaFilter biocartaFilter;
	
	protected List<Filter> filterList;
	
//	private int goBpMinLevel, goBpMaxLevel, goBpMinNumberGenes, goBpMaxNumberGenes;
//	private String goBpDBKeywords, goBpDBKeywordsLogic;
//	private boolean goBpAllGenome, goBpInclusive;

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

		// kegg
		getOptions().addOption(OptionFactory.createOption("kegg", "Kegg database",false,false));
		getOptions().addOption(OptionFactory.createOption("kegg-min-num-genes", "Kegg min number of genes",false));
		getOptions().addOption(OptionFactory.createOption("kegg-max-num-genes", "Kegg max number of genes",false));
		getOptions().addOption(OptionFactory.createOption("kegg-count-genes-from-genome", "computes the number of annotated genes from all genome, otherwise from you input list",false,false));
		
		// biocarta
		getOptions().addOption(OptionFactory.createOption("biocarta", "Biocarta database",false,false));
		getOptions().addOption(OptionFactory.createOption("biocarta-min-num-genes", "Biocarta min number of genes",false));
		getOptions().addOption(OptionFactory.createOption("biocarta-max-num-genes", "Biocarta max number of genes",false));
		getOptions().addOption(OptionFactory.createOption("biocarta-count-genes-from-genome", "computes the number of annotated genes from all genome, otherwise from you input list",false,false));

		// your annotations
		getOptions().addOption(OptionFactory.createOption("annotations", "Your own annotations",false,true));
		getOptions().addOption(OptionFactory.createOption("annotation-enrichment", "rest of databases will used to enrich your annotations",false,false));
		
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
	
	public void prepare() throws IOException, ParseException, InvalidColumnIndexException{

		// species
		setSpecies(commandLine.getOptionValue("species","hsa"));
		
		// fisher test
		String testMode = commandLine.getOptionValue("test-mode", "two-tailed");
		if(testMode.equals("less")) setTestMode(FisherExactTest.LESS);
		if(testMode.equals("greater")) setTestMode(FisherExactTest.GREATER);
		if(testMode.equals("two-tailed")) setTestMode(FisherExactTest.TWO_SIDED);
		
		filterList = new ArrayList<Filter>();
		
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
		// kegg
		if(commandLine.hasOption("kegg")) {			
			KeggFilter keggFilter = new KeggFilter(Integer.parseInt(commandLine.getOptionValue("kegg-min-num-genes","5")),Integer.parseInt(cmdLine.getOptionValue("kegg-max-num-genes","500")));
			filterList.add(keggFilter);
		}
		// biocarta
		if(commandLine.hasOption("biocarta")) {
			BiocartaFilter biocartaFilter = new BiocartaFilter(Integer.parseInt(commandLine.getOptionValue("biocarta-min-num-genes","5")),Integer.parseInt(cmdLine.getOptionValue("biocarta-max-num-genes","500")));
			filterList.add(biocartaFilter);
		}
		// your annotations
		if(cmdLine.hasOption("annotations") && !cmdLine.getOptionValue("annotations").equalsIgnoreCase("") && !cmdLine.getOptionValue("annotations").equalsIgnoreCase("none")) {
			isYourAnnotations = true;
			System.err.println("annotation file: " + cmdLine.getOptionValue("annotations"));
			FeatureData annotations = new FeatureData(new File(cmdLine.getOptionValue("annotations")), true);
			List<String> ids = annotations.getDataFrame().getColumn(0);
			List<String> terms = annotations.getDataFrame().getColumn(1);			
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
		String infraredNamespace = "biological_process";		
		if(namespace.equalsIgnoreCase("cc")) infraredNamespace = "cellular_component";
		if(namespace.equalsIgnoreCase("mf")) infraredNamespace = "molecular_function";
		GOFilter goFilter = new GOFilter(infraredNamespace);
		goFilter.setMinLevel(Integer.parseInt(cmdLine.getOptionValue("go-" + namespace + "-min-level","5")));
		goFilter.setMaxLevel(Integer.parseInt(cmdLine.getOptionValue("go-" + namespace + "-max-level","12")));
		goFilter.setMinNumberGenes(Integer.parseInt(cmdLine.getOptionValue("go-" + namespace + "-min-num-genes","5")));	
		goFilter.setMaxNumberGenes(Integer.parseInt(cmdLine.getOptionValue("go-" + namespace + "-max-num-genes","500")));
		//goFilter.addKeywords(arg0);
		goFilter.setLogicalOperator(cmdLine.getOptionValue("go-" + namespace + "-keywords-logic","AND"));
		return goFilter;
	}
	

	protected String getDBName(Filter filter){
		String name = StringUtils.randomString();
		if(filter instanceof GOFilter) {						
			GOFilter goFilter = (GOFilter) filter.clone();
			name = "go_" + goFilter.getNamespace() + "_" + goFilter.getMinLevel() + "_" + goFilter.getMaxLevel();
		} else if(filter instanceof KeggFilter) {
			name = "kegg";
		}
		else if(filter instanceof BiocartaFilter) {
			name = "biocarta";
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
		} else if(filter instanceof BiocartaFilter) {
			title = "Biocarta";
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
