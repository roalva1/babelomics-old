package org.bioinfo.babelomics.tools.functional;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.ParseException;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.collections.exceptions.InvalidColumnIndexException;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.data.dataset.FeatureData;
import org.bioinfo.infrared.common.dbsql.DBConnector;
import org.bioinfo.infrared.common.feature.FeatureList;
import org.bioinfo.infrared.core.Gene;
import org.bioinfo.infrared.core.dbsql.GeneDBManager;
import org.bioinfo.infrared.core.dbsql.XRefDBManager;
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


	@Override
	public void initOptions() {

		// commons options
		getOptions().addOption(OptionFactory.createOption("test-mode", "the Fisher test mode, valid values: less, greater, two_sided. By default, two_sided", false));
		
		// GO biological process options
		getOptions().addOption(OptionFactory.createOption("go-bp", "GO biological process database",false,false));
		getOptions().addOption(OptionFactory.createOption("go-bp-inclusive", "GO biological process, inclusive analysis (one per GO level), otherwise joins GO levels",false,false));
		getOptions().addOption(OptionFactory.createOption("go-bp-min-level", "GO biological process, min go level to take into account, default 3",false));
		getOptions().addOption(OptionFactory.createOption("go-bp-max-level", "GO biological process, max GO level to take into account, default 15",false));
		getOptions().addOption(OptionFactory.createOption("go-bp-min-num-genes", "GO biological process, min number of genes filter",false));
		getOptions().addOption(OptionFactory.createOption("go-bp-max-num-genes", "GO biological process, max number of genes filter",false));
		getOptions().addOption(OptionFactory.createOption("go-bp-all-genome", "GO biological process, computes the number of annotated genes from all genome, otherwise from you input list",false,false));
		getOptions().addOption(OptionFactory.createOption("go-bp-keywords", "GO biological process, keywords filter",false));
		getOptions().addOption(OptionFactory.createOption("go-bp-keywords-logic", "GO biological process, keywords filter logic: all or any",false));

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



	}

	public void prepare(CommandLine cmdLine) throws IOException, ParseException, InvalidColumnIndexException{

		// species
		setSpecies(cmdLine.getOptionValue("species","hsa"));
		
		// fisher test
		String testMode = cmdLine.getOptionValue("test-mode", "two-tailed");
		if(testMode.equals("less")) setTestMode(FisherExactTest.LESS);
		if(testMode.equals("greater")) setTestMode(FisherExactTest.GREATER);
		if(testMode.equals("two-tailed")) setTestMode(FisherExactTest.TWO_SIDED);
		
		filterList = new ArrayList<Filter>();
		
		// go bp
		if(cmdLine.hasOption("go-bp")) {
			//goBpFilterList = new ArrayList<GOFilter>();
			if(cmdLine.hasOption("go-bp-inclusive")) {
				throw new ParseException("inclusive analysis not yet implemented");
			} else {
				GOFilter goBpFilter = new GOFilter("biological_process");
				goBpFilter.setMinLevel(Integer.parseInt(cmdLine.getOptionValue("go-bp-min-level","5")));
				goBpFilter.setMaxLevel(Integer.parseInt(cmdLine.getOptionValue("go-bp-max-level","12")));
				goBpFilter.setMinNumberGenes(Integer.parseInt(cmdLine.getOptionValue("go-bp-min-num-genes","5")));	
				goBpFilter.setMaxNumberGenes(Integer.parseInt(cmdLine.getOptionValue("go-bp-max-num-genes","500")));				
				goBpFilter.setLogicalOperator(cmdLine.getOptionValue("go-bp-keywords-logic","AND"));
				filterList.add(goBpFilter);
			}			
		}
		if(cmdLine.hasOption("kegg")) {
			System.err.println("entrando en kegg");
			KeggFilter keggFilter = new KeggFilter(Integer.parseInt(cmdLine.getOptionValue("kegg-min-num-genes","5")),Integer.parseInt(cmdLine.getOptionValue("kegg-max-num-genes","500")));
			filterList.add(keggFilter);
		}
		if(cmdLine.hasOption("biocarta")) {
			BiocartaFilter biocartaFilter = new BiocartaFilter(Integer.parseInt(cmdLine.getOptionValue("biocarta-min-num-genes","5")),Integer.parseInt(cmdLine.getOptionValue("biocarta-max-num-genes","500")));
			filterList.add(biocartaFilter);
		}
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
