package org.bioinfo.babelomics.tools.functional;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionGroup;
import org.apache.commons.cli.ParseException;
import org.apache.commons.math.genetics.Chromosome;
import org.bioinfo.collections.exceptions.InvalidColumnIndexException;
import org.bioinfo.commons.Config;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.data.dataset.FeatureData;
import org.bioinfo.infrared.common.dbsql.DBConnector;
import org.bioinfo.infrared.funcannot.filter.Filter;
import org.bioinfo.infrared.funcannot.filter.GOFilter;
import org.bioinfo.math.result.FisherTestResult;
import org.bioinfo.math.result.TestResultList;
import org.bioinfo.math.stats.MultipleTestCorrection;
import org.bioinfo.tool.OptionFactory;

public class FatiGOTool extends FunctionalProfilingTool{

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
			input2.addOption(OptionFactory.createOption("rest-of-genome", "compares list #1 with the rest of genome",false,false));
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
		} else if ( cmdLine.hasOption("rest-of-genome")) {
			this.setRestOfGenome(true);
		} else if(cmdLine.hasOption("chromosomes")) {
			throw new ParseException("chromosome reading not yet implemented");
		} else {
			throw new ParseException("No comparison provided, use list2, rest-of-genome or chromosome options to set your comparison");
		}
		
		// extras
		setRemoveDuplicates(cmdLine.getOptionValue("remove-duplicates", "never"));
	}

	@Override
	public void execute() {
		try {
			// infrared connector
			Config infraredConfig = new Config(System.getenv("BABELOMICS_HOME") + "/conf/infrared.conf");
			DBConnector dbConnector = new DBConnector(getSpecies(), new File(System.getenv("BABELOMICS_HOME") + "/conf/infrared.conf"));
			
			// prepare params
			prepare(commandLine);			
	
			// list 1
			List<String> idList1;
			
			idList1 = FunctionalProfilingTool.toEnsembId(dbConnector, list1.getDataFrame().getColumn(0));
			
			// list 2
			List<String> idList2;
			if(list2!=null) {
				idList2 = FunctionalProfilingTool.toEnsembId(dbConnector, list2.getDataFrame().getColumn(0));
			} else if(isRestOfGenome()) {
				idList2 = FunctionalProfilingTool.getGenes(dbConnector);
			} else if(chromosomes!=null) {
				throw new ParseException("chromosome comparison not yet implemented");
			} else {
				throw new ParseException("No comparison provided, use list2, rest-of-genome or chromosome options to set your comparison");				
			}
			
			
			if(filterList.size()==0){
				throw new ParseException("No biological database selected (eg. --go-db-bp)");
			} else {
				for(Filter filter: filterList) {
					// GO
					if(filter instanceof GOFilter) {						
						GOFilter goFilter = (GOFilter) filter.clone();
						logger.info(goFilter.getNamespace() + " (" + goFilter.getMinLevel() + "," + goFilter.getMaxLevel() + ")...\n");
						
						FatigoTest fatigo = new FatigoTest(idList1, idList2, goFilter, dbConnector, testMode, duplicatesMode);						
						fatigo.run();
						
						TestResultList<FisherTestResult> testResult = fatigo.getResult();
						
						logger.info(goFilter.getNamespace() + " (" + goFilter.getMinLevel() + "," + goFilter.getMaxLevel() + ")\n" + testResult.toString());
						logger.info("...end of GO " + goFilter.getNamespace());				
					}
				}
			}
			
		}catch(Exception e){
			e.printStackTrace();
		}
		
	}
}
