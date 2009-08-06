package org.bioinfo.babelomics.tools.functional;

import java.io.File;
import java.io.IOException;
import java.sql.SQLException;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.ParseException;
import org.bioinfo.data.dataset.FeatureData;
import org.bioinfo.infrared.common.dbsql.DBConnector;
import org.bioinfo.infrared.funcannot.filter.Filter;
import org.bioinfo.math.exception.InvalidParameterException;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.utils.StringUtils;

public class FatiScan extends FunctionalProfilingTool{

	public FatiScan(String[] args) {
		initOptions();
	}

	@Override
	public void initOptions() {
		super.initOptions();

		getOptions().addOption(OptionFactory.createOption("list", "the feature data containig the list of genes and its corresponding statistic"));
		getOptions().addOption(OptionFactory.createOption("nb-partitions", "Number of partitions, 30",false));
	}

	@Override
	public void prepare(CommandLine cmdLine) throws IOException, ParseException {
		super.prepare(cmdLine);
		
		setFeatureData(new FeatureData(new File(cmdLine.getOptionValue("list")), true));
		setNumberOfPartitions(Integer.parseInt(cmdLine.getOptionValue("nb-partitions", "30")));
	}

	@Override
	public void execute() {
		try {
			logger.info("begin of fatiscan");
			prepare(commandLine);			

			DBConnector dbConnector = new DBConnector(getSpecies());
			logger.info("db connector (" + dbConnector.toString() + ")");

			List<String> list1, list2;
			int partitionSize = getFeatureData().getDataFrame().getRowNumber() / getNumberOfPartitions();

			Filter filter = getFilter();
			String dbTarget = getDbTarget();

			for (int p=0 ; p<getNumberOfPartitions()-1 ; p++) {
				logger.info("************** partition " + p + ":");
				list1 = getFeatureData().getDataFrame().getColumn(0).subList(0, (p*partitionSize) + partitionSize - 1); 
				list2 = getFeatureData().getDataFrame().getColumn(0).subList((p*partitionSize) + partitionSize, getFeatureData().getDataFrame().getRowNumber() - 1);

				logger.info("list1: " + StringUtils.arrayToString(list1, "---"));
				logger.info("list2: " + StringUtils.arrayToString(list2, "---"));
				
				FuncTest test = new FuncTest(list1, list2, dbTarget, dbConnector);
				test.setFilter(filter);
				if ( test.run() != null ) {
					logger.info("\nno inclusive, fisher, result :\n" + test.getResult().toString());
				} else {
					logger.info("\nno inclusive, fisher, result is NULL\n");					
				}
			}
			logger.info("end of fatiscan");
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
		} catch (ParseException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
}
