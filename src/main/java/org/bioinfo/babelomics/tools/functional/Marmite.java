package org.bioinfo.babelomics.tools.functional;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.ParseException;
import org.apache.commons.math.MathException;
import org.apache.commons.math.linear.MatrixIndexException;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.collections.array.NamedArrayList;
import org.bioinfo.collections.matrix.DataFrame;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.data.dataset.FeatureData;
import org.bioinfo.graphics.canvas.Canvas;
import org.bioinfo.graphics.canvas.feature.ScoreFeature;
import org.bioinfo.graphics.canvas.panel.GridPanel;
import org.bioinfo.graphics.canvas.track.GridTrack;
import org.bioinfo.math.result.CorrelationTestResult;
import org.bioinfo.math.result.TestResultList;
import org.bioinfo.math.stats.correlation.CorrelationTest;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;
import org.bioinfo.utils.ListUtils;

public class Marmite extends BabelomicsTool {


	public Marmite(String[] args) {
		initOptions();
	}

	@Override
	public void initOptions() {
		options.addOption(OptionFactory.createOption("list1", "the feature data containig the list #1 of genes, or the feature data file"));
		options.addOption(OptionFactory.createOption("list2", "the feature data containig the list #2 of genes, or the feature data file"));
		options.addOption(OptionFactory.createOption("bioentity-name", "Valid values are: disease (for disease associated words), products (for chemical products) or roots (for word roots)"));
		options.addOption(OptionFactory.createOption("bioentity-score-filter", "Minimum number of genes with a score (0-10000)"));
		options.addOption(OptionFactory.createOption("bioentity-number-filter", "Number of bio-entities in results (0-10000)"));
	}

	@Override
	public void execute() {
		try {
			FeatureData fd1 = new FeatureData(new File(commandLine.getOptionValue("list1")));
			FeatureData fd2 = new FeatureData(new File(commandLine.getOptionValue("list2")));

			String bioentity = commandLine.getOptionValue("bioentity-name");
			int scoreFilter = Integer.parseInt(commandLine.getOptionValue("bioentity-score-filter", "5"));
			int numberFilter = Integer.parseInt(commandLine.getOptionValue("bioentity-number-filter", "50"));

			executeMarmite(fd1, fd2, bioentity, scoreFilter, numberFilter);
			
		} catch (IOException e) {
			logger.error("Error opening the dataset", e.toString());
		}
	}
	
	private void executeMarmite(FeatureData fd1, FeatureData fd2, String bioentity, int scoreFilter, int numberFilter) {
		
		try {
			
			// reading data
			//
			jobStatus.addStatusMessage("20", "reading data");
			logger.debug("reading data...\n");
			
			
			// marmite test
			//
			jobStatus.addStatusMessage("40", "computing marmite functional enrichment");
			logger.debug("computing marmite functional enrichment...\n");
									
			
			
			// generating boxplots
			//
			jobStatus.addStatusMessage("60", "generating boxplots");
			logger.debug("generating boxplots...\n");
					

			// saving data
			//
			jobStatus.addStatusMessage("80", "saving results");
			logger.debug("saving results...");		
//			result.addOutputItem(new Item("pearson_correlation_file", getOutdir() + "/pearson_correlation.txt", "The pearson correlation file is: ", TYPE.FILE));
//			result.addOutputItem(new Item("pearson_correlation_heatmap", getOutdir() + "/pearson_heatmap.png", "The pearson correlation heatmap is: ", TYPE.IMAGE));

			
			// done
			//
			jobStatus.addStatusMessage("100", "done");
			logger.debug("marmite funcitonal enrichment done\n");
						
		} catch (java.security.InvalidParameterException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (MatrixIndexException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (MathException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
}
