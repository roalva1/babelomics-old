package org.bioinfo.babelomics.tools.expression;

import java.io.File;
import java.io.IOException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.ParseException;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.tool.OptionFactory;

public class TimeDosageSeries extends BabelomicsTool {
	//Time/dosage series
	
	public TimeDosageSeries(String[] args) {
		super(args);
		initOptions();
	}

	@Override
	public void initOptions() {
		// TODO Auto-generated method stub
		options.addOption(OptionFactory.createOption("dataset", "the data"));
		options.addOption(OptionFactory.createOption("time", "time variable"));
		options.addOption(OptionFactory.createOption("series", "series variable"));
		options.addOption(OptionFactory.createOption("names", "names variable"));
		options.addOption(OptionFactory.createOption("degree", "class variable"));
		options.addOption(OptionFactory.createOption("comparison", "Multiple testing adjustment", false));
		options.addOption(OptionFactory.createOption("significance", "Significance Level for model variable(0-0.9)", false));
		options.addOption(OptionFactory.createOption("cluster_comparison", "Cluster method", false));
		options.addOption(OptionFactory.createOption("kmeans_options", "Number of clusters", false));
	}

	@Override
	public void execute() {
		// TODO Auto-generated method stub
		
		try {
			CommandLine cmd = parse(this.args, true);
			Dataset dataset = new Dataset(new File(cmd.getOptionValue("dataset")));
			String time = cmd.getOptionValue("time");
			String series = cmd.getOptionValue("series");
			String names = cmd.getOptionValue("names");
			String degree = cmd.getOptionValue("degree");
			String comparison = cmd.getOptionValue("comparison", null);
			String significance = cmd.getOptionValue("significance", null);
			String clusterComparison = cmd.getOptionValue("cluster_comparison", null);
			String kmeansOptions = cmd.getOptionValue("kmeans_options", null);
			
			executeDetds(dataset, time,series ,names,degree,comparison,significance,clusterComparison,kmeansOptions);
			
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (ParseException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
		
		
		
	}
	
	private void executeDetds(Dataset dataset, String time,String series ,String names,String degree,String comparison,String significance,String clusterComparison,String kmeansOptions) {
		logger.info("executing TimeDosageSeries, not implemented yet");
	}
	

}
