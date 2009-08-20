package org.bioinfo.babelomics.tools.functional;

import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.ParseException;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.data.dataset.FeatureData;
import org.bioinfo.infrared.common.dbsql.DBConnector;
import org.bioinfo.infrared.common.feature.FeatureList;
import org.bioinfo.infrared.core.Gene;
import org.bioinfo.infrared.core.dbsql.GeneDBManager;
import org.bioinfo.infrared.core.dbsql.XRefDBManager;
import org.bioinfo.infrared.funcannot.filter.Filter;
import org.bioinfo.infrared.funcannot.filter.GOFilter;
import org.bioinfo.math.stats.inference.FisherExactTest;
import org.bioinfo.tool.OptionFactory;


public abstract class FunctionalProfilingTool extends BabelomicsTool {

	// commons
	//
	private FeatureData featureData;
	private FeatureData featureData2;

	private String species;

	private boolean restOfGenome;
	private String chromosome;
	private int chromosomeStart;
	private int chromosomeEnd;

	private int testMode;
	private String removeDuplicates;
	
	private int NumberOfPartitions;
	
	// GO biological process
	//
	private boolean goBpDB;
	private int goBpMinLevel, goBpMaxLevel, goBpMinNumberGenes, goBpMaxNumberGenes;
	private String goBpDBKeywords, goBpDBKeywordsLogic;
	private boolean goBpAllGenome, goBpInclusive;

	public static String REMOVE_NEVER = "never";
	public static String REMOVE_EACH = "each";
	public static String REMOVE_REF = "ref";

	public static String FISHER_LESS = "less";
	public static String GREATER = "greater";
	public static String TWO_SIDED = "two_sided";

//	public FunctionalProfilingTool(String[] args) {
//	}

	@Override
	public void initOptions() {

		// commons options
		//
		options.addOption(OptionFactory.createOption("test-mode", "the Fisher test mode, valid values: less, greater, two_sided. By default, two_sided", false));
		
		// GO biological process options
		//
		options.addOption(OptionFactory.createOption("go-bp-db", "GO biological process database",false,false));
		options.addOption(OptionFactory.createOption("go-bp-inclusive", "GO biological process, inclusive analysis (one per GO level), otherwise joins GO levels",false,false));
		options.addOption(OptionFactory.createOption("go-bp-min-level", "GO biological process, min go level to take into account, default 3",false));
		options.addOption(OptionFactory.createOption("go-bp-max-level", "GO biological process, max GO level to take into account, default 15",false));
		options.addOption(OptionFactory.createOption("go-bp-min-num-genes", "GO biological process, min number of genes filter",false));
		options.addOption(OptionFactory.createOption("go-bp-max-num-genes", "GO biological process, max number of genes filter",false));
		options.addOption(OptionFactory.createOption("go-bp-all-genome", "GO biological process, computes the number of annotated genes from all genome, otherwise from you input list",false,false));
		options.addOption(OptionFactory.createOption("go-bp-keywords", "GO biological process, keywords filter",false));
		options.addOption(OptionFactory.createOption("go-bp-keywords-logic", "GO biological process, keywords filter logic: all or any",false));

		// GO molecular function options
		//


		// GO cellular component options
		//


	}

	public void prepare(CommandLine cmdLine) throws IOException, ParseException {

		// commons
		//
		setSpecies(cmdLine.getOptionValue("species","hsa"));			

		if  ( FISHER_LESS.equals(cmdLine.getOptionValue("test-mode")) ) {
			setTestMode(FisherExactTest.LESS);
		} else if ( FISHER_LESS.equals(cmdLine.getOptionValue("test-mode")) ) {
			setTestMode(FisherExactTest.GREATER);			
		} else {
			setTestMode(FisherExactTest.TWO_SIDED);			
		}

		// GO biological process database
		//
		if( cmdLine.hasOption("go-bp-db") ) {
			setGoBpDB(cmdLine.hasOption("go-bp-db"));
			setGoBpInclusive(cmdLine.hasOption("go-bp-inclusive"));
			setGoBpMinLevel(Integer.parseInt(cmdLine.getOptionValue("go-bp-min-level","5")));
			setGoBpMaxLevel(Integer.parseInt(cmdLine.getOptionValue("go-bp-max-level","12")));
			setGoBpMinNumberGenes(Integer.parseInt(cmdLine.getOptionValue("go-bp-min-num-genes","5")));	
			setGoBpMaxNumberGenes(Integer.parseInt(cmdLine.getOptionValue("go-bp-max-num-genes","500")));
			setGoBpAllGenome(cmdLine.hasOption("go-bp-all-genome"));
			setGoBpDBKeywords(cmdLine.getOptionValue("go-bp-keywords",""));
			setGoBpDBKeywordsLogic(cmdLine.getOptionValue("go-bp-keywords-logic","AND"));
			setGoBpDBKeywords(cmdLine.getOptionValue("go-bp-keywords",""));
			setGoBpDBKeywordsLogic(cmdLine.getOptionValue("go-bp-keywords-logic","AND"));
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
			
			List<List<String>> list = new XRefDBManager(dbConnector).getIdsByDBName(ids, "ensembl_gene");
			
			List<String> ensemblIds = new ArrayList<String>();
			if ( list == null ) {
				System.out.println(">>>>> getAllByExternalIds returns null");
				ensemblIds = null;
			} else {
				for(List<String> stringList: list) {
					for(String ensemblId: stringList) {
						ensemblIds.add(ensemblId);
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

	public Filter getFilter() {
		if ( isGoBpDB() ) {
			logger.info("go biological process");

			// creating GO biological process filter
			//
			GOFilter gof = new GOFilter("biological_process", getGoBpMinLevel(), getGoBpMaxLevel(), getGoBpMinNumberGenes(), getGoBpMaxNumberGenes());
			gof.setLogicalOperator(getGoBpDBKeywordsLogic());
			if ( getGoBpDBKeywords() != null && !getGoBpDBKeywords().equalsIgnoreCase("")) {
				logger.info("adding the keywords: " + getGoBpDBKeywords());
				gof.addKeyword(getGoBpDBKeywords());
			}
			return gof;
		}
		return null;
	}
	
	public String getDbTarget() {
		if ( isGoBpDB() ) {
			return FuncTest.GO;
		}
		return null;
	}

	public FeatureData getFeatureData() {
		return featureData;
	}

	public void setFeatureData(FeatureData featureData) {
		this.featureData = featureData;
	}

	public String getSpecies() {
		return species;
	}

	public void setSpecies(String species) {
		this.species = species;
	}

	public boolean isGoBpDB() {
		return goBpDB;
	}

	public void setGoBpDB(boolean goBpDB) {
		this.goBpDB = goBpDB;
	}

	public boolean isGoBpAllGenome() {
		return goBpAllGenome;
	}

	public void setGoBpAllGenome(boolean goBpAllGenome) {
		this.goBpAllGenome = goBpAllGenome;
	}

	public int getGoBpMinLevel() {
		return goBpMinLevel;
	}

	public void setGoBpMinLevel(int goBpMinLevel) {
		this.goBpMinLevel = goBpMinLevel;
	}

	public int getGoBpMaxLevel() {
		return goBpMaxLevel;
	}

	public void setGoBpMaxLevel(int goBpMaxLevel) {
		this.goBpMaxLevel = goBpMaxLevel;
	}

	public int getGoBpMinNumberGenes() {
		return goBpMinNumberGenes;
	}

	public void setGoBpMinNumberGenes(int goBpMinNumberGenes) {
		this.goBpMinNumberGenes = goBpMinNumberGenes;
	}

	public int getGoBpMaxNumberGenes() {
		return goBpMaxNumberGenes;
	}

	public void setGoBpMaxNumberGenes(int goBpMaxNumberGenes) {
		this.goBpMaxNumberGenes = goBpMaxNumberGenes;
	}

	public String getGoBpDBKeywords() {
		return goBpDBKeywords;
	}

	public void setGoBpDBKeywords(String goBpDBKeywords) {
		this.goBpDBKeywords = goBpDBKeywords;
	}

	public String getGoBpDBKeywordsLogic() {
		return goBpDBKeywordsLogic;
	}

	public void setGoBpDBKeywordsLogic(String goBpDBKeywordsLogic) {
		this.goBpDBKeywordsLogic = goBpDBKeywordsLogic;
	}

	public boolean isGoBpInclusive() {
		return goBpInclusive;
	}

	public void setGoBpInclusive(boolean goBpInclusive) {
		this.goBpInclusive = goBpInclusive;
	}

	public void setFeatureData2(FeatureData featureData2) {
		this.featureData2 = featureData2;
	}

	public FeatureData getFeatureData2() {
		return featureData2;
	}

	public void setRestOfGenome(boolean restOfGenome) {
		this.restOfGenome = restOfGenome;
	}

	public boolean isRestOfGenome() {
		return restOfGenome;
	}

	public void setChromosome(String chromosome) {
		this.chromosome = chromosome;
	}

	public String getChromosome() {
		return chromosome;
	}

	public void setChromosomeStart(int chromosomeStart) {
		this.chromosomeStart = chromosomeStart;
	}

	public int getChromosomeStart() {
		return chromosomeStart;
	}

	public void setChromosomeEnd(int chromosomeEnd) {
		this.chromosomeEnd = chromosomeEnd;
	}

	public int getChromosomeEnd() {
		return chromosomeEnd;
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

	public void setNumberOfPartitions(int numberOfPartitions) {
		NumberOfPartitions = numberOfPartitions;
	}

	public int getNumberOfPartitions() {
		return NumberOfPartitions;
	}


}
