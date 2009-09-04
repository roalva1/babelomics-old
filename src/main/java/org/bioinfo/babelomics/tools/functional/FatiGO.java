package org.bioinfo.babelomics.tools.functional;

import java.io.File;
import java.io.IOException;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionGroup;
import org.apache.commons.cli.ParseException;
import org.apache.commons.math.genetics.Chromosome;
import org.bioinfo.collections.exceptions.InvalidColumnIndexException;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.log.Logger;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.data.dataset.FeatureData;
import org.bioinfo.infrared.common.dbsql.DBConnector;
import org.bioinfo.infrared.funcannot.filter.GOFilter;
import org.bioinfo.tool.OptionFactory;

public class FatiGO extends FunctionalProfilingTool{

	private FeatureData list1;
	private FeatureData list2;
	
	private List<Chromosome> chromosomes;
	
	public FatiGO() {
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
		input2.setRequired(true);
			input2.addOption(OptionFactory.createOption("list2", "the file containig the list #2 of genes, or the feature data file"));
			input2.addOption(OptionFactory.createOption("rest-of-genome", "compares list #1 with the rest of genome"));
			input2.addOption(OptionFactory.createOption("chromosomes", "the chromosome to compare with list #1. Use chromosome-start and chromosome-end options tp delimite the region within the chromosome"));
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
			DBConnector dbConnector = new DBConnector(getSpecies());
			
			// prepare params
			prepare(commandLine);			
	
			// list 1
			List<String> idList1;
			System.err.println("list1.getDataFrame().getColumn(0): " + dbConnector);
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

			// extra
			if(REMOVE_EACH.equalsIgnoreCase(this.getRemoveDuplicates())) {
				idList1 = ListUtils.unique(idList1);
				idList2 = ListUtils.unique(idList2);
			} else if(REMOVE_REF.equalsIgnoreCase(this.getRemoveDuplicates())) {
				for (String s: idList1) {
					if(idList2.contains(s)) idList2.remove(s);
				}
			}
			
			// go bp
			if(goBpFilterList!=null) {
				logger.info("GO biological process...");
				FuncTest test = new FuncTest(idList1, idList2, FuncTest.GO, dbConnector);
				for(GOFilter goFilter: goBpFilterList) {					
					test.setFilter(goFilter);
					System.err.println(goFilter.getNamespace() + " (" + goFilter.getMinLevel() + "," + goFilter.getMaxLevel() + ")\n" + test.run().toString());					
				}
				logger.info("...end of GO bp execution");				
			}
			
		}catch(Exception e){
			e.printStackTrace();
		}
		
	}
}
