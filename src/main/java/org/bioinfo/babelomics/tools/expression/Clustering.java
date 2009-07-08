package org.bioinfo.babelomics.tools.expression;

import java.io.File;
import java.io.IOException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.ParseException;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.tool.OptionFactory;

public class Clustering extends BabelomicsTool {


	public Clustering(String[] args) {
		super(args);
		initOptions();
	}

	@Override
	public void initOptions() {
		options.addOption(OptionFactory.createOption("dataset", "the data"));
		options.addOption(OptionFactory.createOption("method", "the method, possible values: upgma, sota, som, kmeans"));
		options.addOption(OptionFactory.createOption("distance", "the distance, possible values: euclidean, spearman, pearson. Default value: euclidean, false"));
		options.addOption(OptionFactory.createOption("kvalue", "k-value for kmeans clustering. Default value: 15", false));
		options.addOption(OptionFactory.createOption("sample-filter", "class variable", false));
		options.addOption(OptionFactory.createOption("feature-filter", "class variable", false));

	}

	@Override
	public void execute() {
		try {
			CommandLine cmd = parse(args);

			Dataset dataset = new Dataset(new File(cmd.getOptionValue("dataset")));
			String method = cmd.getOptionValue("method");

			String distance = cmd.getOptionValue("distance", "euclidean");
			String kvalue = cmd.getOptionValue("time-class", "15");
				
			if(cmd.hasOption("sample-filter") || cmd.hasOption("feature-filter")) {
				dataset = dataset.getSubDataset(cmd.getOptionValue("sample-filter"), "4", cmd.getOptionValue("feature-filter"), ""); 
			}
			System.out.println(dataset.toString()+"\n");
			
			if ( "upgma".equalsIgnoreCase(method) ) {
				executeUpgma(dataset, distance);
				return;
			}
			
			if ( "sota".equalsIgnoreCase(method) ) {
				executeSota(dataset, distance);
				return;
			}

			if ( "som".equalsIgnoreCase(method) ) {
				executeSom(dataset, distance);
				return;
			}

			if ( "kmeans".equalsIgnoreCase(method) ) {
				int k;
				try {
					k = Integer.parseInt(kvalue);
					executeKmeans(dataset, distance, k);
				} catch (NumberFormatException e ) {
					logger.error("Invalid k-value: " + kvalue);
				}
				return;
			}

			logger.warn("que raroo....");
		} catch (ParseException e) {
			logger.error("Error parsing command line", e.toString());
			System.out.println("\n");
			printUsage();
		} catch (IOException e) {
			logger.error("Error opening the dataset", e.toString());
		} 
	}


	private void executeUpgma(Dataset dataset, String distance) {
		logger.info("executing upgma, not implemented yet");
	}

	private void executeSota(Dataset dataset, String distance) {
		logger.info("executing data adaptive, not implemented yet");
	}
	
	private void executeSom(Dataset dataset, String distance) {
		logger.info("executing som, not implemented yet");
	}

	private void executeKmeans(Dataset dataset, String distance, int kvalue) {
		logger.info("executing kmeans, not implemented yet");
	}
}
