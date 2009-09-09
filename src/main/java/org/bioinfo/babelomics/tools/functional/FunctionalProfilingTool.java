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

	// JT.2009.09.07
	private String removeDuplicates;
	
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
		getOptions().addOption(OptionFactory.createOption("go-bp-db", "GO biological process database",false,false));
		getOptions().addOption(OptionFactory.createOption("go-bp-inclusive", "GO biological process, inclusive analysis (one per GO level), otherwise joins GO levels",false,false));
		getOptions().addOption(OptionFactory.createOption("go-bp-min-level", "GO biological process, min go level to take into account, default 3",false));
		getOptions().addOption(OptionFactory.createOption("go-bp-max-level", "GO biological process, max GO level to take into account, default 15",false));
		getOptions().addOption(OptionFactory.createOption("go-bp-min-num-genes", "GO biological process, min number of genes filter",false));
		getOptions().addOption(OptionFactory.createOption("go-bp-max-num-genes", "GO biological process, max number of genes filter",false));
		getOptions().addOption(OptionFactory.createOption("go-bp-all-genome", "GO biological process, computes the number of annotated genes from all genome, otherwise from you input list",false,false));
		getOptions().addOption(OptionFactory.createOption("go-bp-keywords", "GO biological process, keywords filter",false));
		getOptions().addOption(OptionFactory.createOption("go-bp-keywords-logic", "GO biological process, keywords filter logic: all or any",false));

		// kegg
		getOptions().addOption(OptionFactory.createOption("kegg-db", "Kegg database",false,false));
		getOptions().addOption(OptionFactory.createOption("kegg-min-num-genes", "Kegg min number of genes",false));
		getOptions().addOption(OptionFactory.createOption("kegg-max-num-genes", "Kegg max number of genes",false));
		getOptions().addOption(OptionFactory.createOption("kegg-count-genes-from-genome", "computes the number of annotated genes from all genome, otherwise from you input list",false,false));
		
		// biocarta
		getOptions().addOption(OptionFactory.createOption("biocarta-db", "Biocarta database",false,false));
		getOptions().addOption(OptionFactory.createOption("biocarta-min-num-genes", "Biocarta min number of genes",false));
		getOptions().addOption(OptionFactory.createOption("biocarta-max-num-genes", "Biocarta max number of genes",false));
		getOptions().addOption(OptionFactory.createOption("biocarta-count-genes-from-genome", "computes the number of annotated genes from all genome, otherwise from you input list",false,false));



	}

	public void prepare(CommandLine cmdLine) throws IOException, ParseException, InvalidColumnIndexException{

		// species
		setSpecies(cmdLine.getOptionValue("species","hsa"));
		
		// fisher test
		String testMode = cmdLine.getOptionValue("test-mode", "two-sided");
		if(testMode.equals("less")) setTestMode(FisherExactTest.LESS);
		if(testMode.equals("greater")) setTestMode(FisherExactTest.GREATER);
		if(testMode.equals("two-sided")) setTestMode(FisherExactTest.TWO_SIDED);
		
		filterList = new ArrayList<Filter>();
		
		// go bp
		if(cmdLine.hasOption("go-bp-db")) {
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
		if(cmdLine.hasOption("kegg-db")) {
			KeggFilter keggFilter = new KeggFilter(Integer.parseInt(cmdLine.getOptionValue("kegg-min-num-genes","5")),Integer.parseInt(cmdLine.getOptionValue("kegg-max-num-genes","500")));
			filterList.add(keggFilter);
		}
		if(cmdLine.hasOption("biocarta-db")) {
			BiocartaFilter biocartaFilter = new BiocartaFilter(Integer.parseInt(cmdLine.getOptionValue("biocarta-min-num-genes","5")),Integer.parseInt(cmdLine.getOptionValue("biocarta-max-num-genes","500")));
			filterList.add(biocartaFilter);
		}
	}


	public static List<String> getGenes(DBConnector dbConnector) {
		try {
			return new GeneDBManager(dbConnector).getAllEnsemblIds();
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}

	public static List<String> getGenes(DBConnector dbConnector, String chromosome, int start, int end) {
		try {
			FeatureList<Gene> genes = new GeneDBManager(dbConnector).getAllByLocation(chromosome, start, end);
			List<String> list = new ArrayList<String> (genes.size());
			for (Gene gene: genes) {
				list.add(gene.getId());
			}
			return list;

		} catch (Exception e) {
			e.printStackTrace();
			return null;
		}
	}

	public static List<String> toEnsembId(DBConnector dbConnector, List<String> ids) {
		try {
			
			//List<FeatureList<Gene>> list = new GeneDBManager(dbConnector).getAllByExternalIds(ids);
			System.err.println("dbconnector" + dbConnector);
			List<List<String>> list = new XRefDBManager(dbConnector).getIdsByDBName(ids, "ensembl_gene");
			
			List<String> ensemblIds = new ArrayList<String>();
			if ( list == null ) {
				System.out.println(">>>>> getAllByExternalIds returns null");
				ensemblIds = null;
			} else {
				
				for(List<String> stringList: list) {					
					if(stringList!=null){
						for(String ensemblId: stringList) {						
							ensemblIds.add(ensemblId);
						}
					}
//				for(FeatureList<Gene> featureList: list) {
//					ensemblIds.addAll(featureList.getFeaturesIds());
//					for(Gene gene: featureList) {
//						ensemblIds.add(gene.getStableId());
//					}
				}
				ensemblIds = ListUtils.unique(ensemblIds);
			}
			return ensemblIds;
		} catch (Exception e) {
			e.printStackTrace();
			return null;
		}
	}
	
	public static Map<String, List<String>> getEnsemblMap(DBConnector dbConnector, List<String> ids) {
		try {
			Map<String, List<String>> map = new HashMap<String, List<String>>();
			
			List<List<String>> lists = new XRefDBManager(dbConnector).getIdsByDBName(ids, "ensembl_gene");

//			System.out.println("ids size = " + ids.size() + ", list size " + lists.size());
			List<String> ensemblIds = new ArrayList<String>();
			if ( lists == null ) {
				System.out.println(">>>>> getAllByExternalIds returns null");
				map = null;
			} else {
				for(int i=0 ; i<lists.size() ; i++) {
					if ( lists.get(i) == null || lists.get(i).size() <= 0) {
						map.put(ids.get(i), null);
					} else {
//						if ( lists.get(i).size() > 1 ) {
//							System.out.println("-------------->>>> gene " + ids.get(i) + ", more than one ensemblId : " + ListUtils.toString(lists.get(i)));
//						}
						if ( !map.containsKey(ids.get(i)) ) {
							map.put(ids.get(i), new ArrayList<String>());
						}
						map.get(ids.get(i)).addAll(lists.get(i));
					}
				}
			}
			return map;
		} catch (Exception e) {
			e.printStackTrace();
			return null;
		}
	}

//	public Filter getFilter() {
//		if ( isGoBpDB() ) {
//			logger.info("go biological process");
//
//			// creating GO biological process filter
//			//
//			GOFilter gof = new GOFilter("biological_process", getGoBpMinLevel(), getGoBpMaxLevel(), getGoBpMinNumberGenes(), getGoBpMaxNumberGenes());
//			gof.setLogicalOperator(getGoBpDBKeywordsLogic());
//			if ( getGoBpDBKeywords() != null && !getGoBpDBKeywords().equalsIgnoreCase("")) {
//				logger.info("adding the keywords: " + getGoBpDBKeywords());
//				gof.addKeyword(getGoBpDBKeywords());
//			}
//			return gof;
//		}
//		return null;
//	}
	
//	public String getDbTarget() {
//		if ( isGoBpDB() ) {
//			return FuncTest.GO;
//		}
//		return null;
//	}

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

	public void setRemoveDuplicates(String removeDuplicates) {
		this.removeDuplicates = removeDuplicates;
	}

	public String getRemoveDuplicates() {
		return removeDuplicates;
	}


}
