package org.bioinfo.babelomics.tools.functional;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.OptionGroup;
import org.apache.commons.cli.ParseException;
import org.bioinfo.data.dataset.FeatureData;
import org.bioinfo.infrared.common.dbsql.DBConnector;
import org.bioinfo.infrared.funcannot.filter.GOFilter;
import org.bioinfo.math.exception.InvalidParameterException;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.utils.ListUtils;

public class FatiGO extends FunctionalProfilingTool{

	public FatiGO(String[] args) {
		initOptions();
	}
	
	@Override
	public void initOptions() {
		super.initOptions();
		
		getOptions().addOption(OptionFactory.createOption("list1", "the feature data containig the list #1 of genes, or the feature data file"));

		OptionGroup input2 = new OptionGroup();
		input2.setRequired(true);
		input2.addOption(OptionFactory.createOption("list2", "the file containig the list #2 of genes, or the feature data file"));
		input2.addOption(OptionFactory.createOption("rest-of-genome", "compares list #1 with the rest of genome"));
		input2.addOption(OptionFactory.createOption("chromosome", "the chromosome to compare with list #1. Use chromosome-start and chromosome-end options tp delimite the region within the chromosome"));
		getOptions().addOptionGroup(input2);

		getOptions().addOption(OptionFactory.createOption("chromosome-start", "the region start position within the chromosome", false));
		getOptions().addOption(OptionFactory.createOption("chromosome-end", "the region end position within the chromosome, -1 means the last position of the chromosome", false));

		getOptions().addOption(OptionFactory.createOption("filter-name1", "the feature filter name, requires a feature data file in list #1", false));
		getOptions().addOption(OptionFactory.createOption("filter-value1", "the feature filter value, requires a feature data file in list #1", false));

		getOptions().addOption(OptionFactory.createOption("filter-name2", "the feature filter name, requires a feature data file in list #2", false));
		getOptions().addOption(OptionFactory.createOption("filter-value2", "the feature filter value, requires a feature data file in list #2", false));

		getOptions().addOption(OptionFactory.createOption("remove-duplicates", "to remove duplicated IDs, valid values: never, each, ref. By default, never", false));
	}

	@Override
	public void prepare(CommandLine cmdLine) throws IOException, ParseException {
		super.prepare(cmdLine);

		setFeatureData(new FeatureData(new File(cmdLine.getOptionValue("list1")), true));
		if ( cmdLine.hasOption("filter-name1") && cmdLine.hasOption("filter-value1") ) {
			setFeatureData(getFeatureData().getSubFeatureData(cmdLine.getOptionValue("filter-name1"), cmdLine.getOptionValue("filter-value1"))); 
		}

		if ( cmdLine.hasOption("list2") ) {
			setFeatureData2(new FeatureData(new File(cmdLine.getOptionValue("list2")), true));
			if ( cmdLine.hasOption("filter-name2") && cmdLine.hasOption("filter-value2") ) {
				setFeatureData2(getFeatureData2().getSubFeatureData(cmdLine.getOptionValue("filter-name2"), cmdLine.getOptionValue("filter-value2"))); 
			}
		} else if ( cmdLine.hasOption("rest-of-genome") ) {
			this.setRestOfGenome(true);
		} else if ( cmdLine.hasOption("chromosome") ) {
			setChromosome(cmdLine.getOptionValue("chromosome"));
			setChromosomeStart(Integer.parseInt(cmdLine.getOptionValue("chromosome-start", "0")));
			setChromosomeEnd(Integer.parseInt(cmdLine.getOptionValue("chromosome-end", "-1")));
		} else {
			throw new ParseException("No comparison provided, use list2, rest-of-genome or chromosome options to set your comparison");
		}

		setRemoveDuplicates(cmdLine.getOptionValue("remove-duplicates", "never"));
	}
	
	@Override
	public void execute() {
		try {
//			CommandLine cmdLine = parse(args);
			prepare(commandLine);			

			DBConnector dbConnector = new DBConnector(getSpecies());
			logger.info("db connector (" + dbConnector.toString() + ")");

			// setting the lists to compare
			//
			List<String> list1, list2;
			list1 = FunctionalProfilingTool.toEnsembId(dbConnector, getFeatureData().getDataFrame().getColumn(0));
			if ( getFeatureData2() != null ) {
				list2 = FunctionalProfilingTool.toEnsembId(dbConnector, getFeatureData2().getDataFrame().getColumn(0));
			} else if ( isRestOfGenome() ) {
				list2 = FunctionalProfilingTool.getGenes(dbConnector);
			} else if ( getChromosome() != null ) {
				list2 = FunctionalProfilingTool.getGenes(dbConnector, getChromosome(), getChromosomeStart(), getChromosomeEnd());
			} else {
				throw new ParseException("No comparison provided, use list2, rest-of-genome or chromosome options to set your comparison");				
			}
			
			//System.out.println("before removing duplictes ********* list1 (" + list1.size() + "):\n" + StringUtils.arrayToString(list1, "\n"));		
			//System.out.println("before removing duplicates ********* list2 (" + list2.size() + "):\n" + StringUtils.arrayToString(list2, "\n"));		

			System.out.println("before removing duplicates ********* list1 size (" + list1.size() + ")\n");		
			System.out.println("before removing duplicates ********* list2 size (" + list2.size() + ")\n");		
			System.out.println("removing duplicates ********* (" + getRemoveDuplicates() + ")\n");		

			// and now removing duplicates if necessary
			//
			if ( REMOVE_EACH.equalsIgnoreCase(this.getRemoveDuplicates()) ) {
				list1 = ListUtils.unique(list1);
				list2 = ListUtils.unique(list2);
			} else if ( REMOVE_REF.equalsIgnoreCase(this.getRemoveDuplicates()) ) {
				for (String s: list1) {
					if ( list2.contains(s) ) {
						list2.remove(s);
					}
				}				
			}
			
			System.out.println("after removing duplicates ********* list1 size (" + list1.size() + ")\n");		
			System.out.println("after removing duplicates ********* list2 size (" + list2.size() + ")\n");		
			//System.out.println("after removing duplictes ********* list1 (" + list1.size() + ") :\n" + StringUtils.arrayToString(list1, "\n"));		
			//System.out.println("after removing duplicates ********* list2 (" + list2.size() + "):\n" + StringUtils.arrayToString(list2, "\n"));		

			// GO biological process database
			//
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
				
				FuncTest test = new FuncTest(list1, list2, FuncTest.GO, dbConnector);
				test.setFilter(gof);

				if ( isGoBpInclusive() ) {
					logger.info("inclusive, creating the annotation file for go biological process");
										
					for(int i = getGoBpMinLevel() ; i<=getGoBpMaxLevel() ; i++) {
						gof.setMinLevel(i);
						gof.setMaxLevel(i);
						
						test.setFilter(gof);
						logger.info("\ninclusive, level " + i + ", fisher, result :\n" + test.run().toString());
						
					}
				} else {
					logger.info("\nno inclusive, fisher, result :\n" + test.run().toString());					
				}
				logger.info("go-bp-db executed");
			}

			
			
		} catch (ParseException e) {
			logger.error("Error parsing command line", e.toString());
			System.out.println("\n");
			this.printUsage("script !!!");
		} catch (IOException e) {
			logger.error("Error opening the feature data", e.toString());
		} catch (IndexOutOfBoundsException e) {
			e.printStackTrace();
		} catch (SQLException e) {
			e.printStackTrace();
		} catch (IllegalAccessException e) {
			e.printStackTrace();
		} catch (ClassNotFoundException e) {
			e.printStackTrace();
		} catch (InstantiationException e) {
			e.printStackTrace();
		} catch (InvalidParameterException e) {
			e.printStackTrace();
		}

	}
}
