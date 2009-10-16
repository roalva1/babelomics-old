package org.bioinfo.babelomics.tools.expression.differential;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.collections.matrix.DataFrame;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.data.dataset.FeatureData;
import org.bioinfo.graphics.canvas.Canvas;
import org.bioinfo.math.regression.SimpleRegressionTest;
import org.bioinfo.math.result.CorrelationTestResult;
import org.bioinfo.math.result.SimpleRegressionTestResult;
import org.bioinfo.math.result.TestResultList;
import org.bioinfo.math.stats.correlation.CorrelationTest;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;

public class Correlation extends BabelomicsTool {

	Dataset dataset;
	String test;
	String className;
	List<String> vars;
	List<Double> doubleVars;
	
	public Correlation() {
		initOptions();
	}

	public void initOptions() {
		options.addOption(OptionFactory.createOption("dataset", "the data"));
		options.addOption(OptionFactory.createOption("test", "the test, possible values: pearson, spearman, regression"));
		options.addOption(OptionFactory.createOption("class", "class variable", false));
	}

	@Override
	public void execute() {

		// reading dataset
		//
		dataset = initDataset(new File(commandLine.getOptionValue("dataset")));

		test = commandLine.getOptionValue("test", null);
		className = commandLine.getOptionValue("class", null);

		if ( ! "pearson".equalsIgnoreCase(test) && ! "spearman".equalsIgnoreCase(test) && ! "regression".equalsIgnoreCase(test) ) {
			abort("unknowntest_execute_correlation", "unknown test (" + test + ")", "unknown test (" + test + ")", "unknown test (" + test + ")");
		}

		vars = dataset.getVariables().getByName(className).getValues();
		doubleVars = new ArrayList<Double>(vars.size());
		for(String str: vars) {
			doubleVars.add(Double.parseDouble(str));
		}
		
		
		if ( "regression".equalsIgnoreCase(test) ) {
			executeRegression();
		} else {
			executeCorrelation();
		}
	}
		
	
	private void executeCorrelation() {
		
		// running test (pearson or spearman)
		//
		updateJobStatus("40", "computing " + test + " test");		
		TestResultList<CorrelationTestResult> res = null;
		try {
			CorrelationTest corrTest = new CorrelationTest(dataset.getDoubleMatrix(), doubleVars, test);
			res = corrTest.compute();
		} catch (Exception e) {
			abort("exception_executecorrelation_correlation", "error running " + test + " test", e.toString(), StringUtils.getStackTrace(e));
		}
		
		int[] columnOrder = ListUtils.order(vars);
		int[] rowOrder = ListUtils.order(ListUtils.toList(res.getCorrelations()), true);

		// generating heatmap
		//
		updateJobStatus("60", "generating heatmap");
		Canvas heatmap = DiffExpressionUtils.generateHeatmap(dataset, className, columnOrder, rowOrder, "correlation", res.getCorrelations(), "adj. p-value", res.getAdjPValues());
		String heatmapFilename = getOutdir() + "/" + test + "_heatmap";
		try {
			heatmap.save(heatmapFilename);
		} catch (IOException e) {
			printError("ioexception_executecorrelation_correlation", "error generating heatmap", e.toString(), e);
		}

		// saving data
		//
		updateJobStatus("80", "saving results");
		DataFrame dataFrame = new DataFrame(dataset.getFeatureNames().size(), 0);
		dataFrame.setRowNames(ListUtils.ordered(dataset.getFeatureNames(), rowOrder));

		try {
			dataFrame.addColumn("statistic", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getStatistics()), rowOrder)));
			dataFrame.addColumn("correlation", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getCorrelations()), rowOrder)));
			dataFrame.addColumn("p-value", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getPValues()), rowOrder)));
			dataFrame.addColumn("adj. p-value", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getAdjPValues()), rowOrder)));

			File outFile = new File(getOutdir() + "/" + test + ".txt");
			FeatureData featureData = new FeatureData(dataFrame);
			featureData.write(outFile);
			result.addOutputItem(new Item(test + "_file", outFile.getName(), test.toUpperCase() + " output file", TYPE.FILE, new ArrayList<String>(2), new HashMap<String, String>(2), "Output files"));

			outFile = new File(getOutdir() + "/" + test + "_table.txt");
			IOUtils.write(outFile, dataFrame.toString(true, true));			
			List<String> tags = new ArrayList<String>();
			tags.add("TABLE");
			result.addOutputItem(new Item(test + "_table", outFile.getName(), test.toUpperCase() + " output table", TYPE.FILE, tags, new HashMap<String, String>(2), "Output files"));
		} catch (Exception e) {
			printError("ioexception_executecorrelation_correlation", "error saving " + test + " results", e.toString(), e);
		}

		if ( new File(heatmapFilename + ".png").exists() ) {
			result.addOutputItem(new Item(test + "_heatmap", test + "_heatmap.png", test.toUpperCase() + " heatmap", TYPE.IMAGE, new ArrayList<String>(2), new HashMap<String, String>(2), "Heatmap image"));
		}
	}




	private void executeRegression() {
		
		// running test (regression)
		//
		updateJobStatus("40", "computing " + test + " test");		
		TestResultList<SimpleRegressionTestResult> res = null;
		try {
			SimpleRegressionTest regression = new SimpleRegressionTest();
			res = regression.compute(dataset.getDoubleMatrix(), doubleVars);
		} catch (Exception e) {
			abort("exception_executecorrelation_correlation", "error running " + test + " test", e.toString(), StringUtils.getStackTrace(e));
		}
		
		int[] columnOrder = ListUtils.order(vars);
		int[] rowOrder = ListUtils.order(ListUtils.toList(res.getStatistics()), true);

		// generating heatmap
		//
		updateJobStatus("60", "generating heatmap");
		Canvas heatmap = DiffExpressionUtils.generateHeatmap(dataset, className, columnOrder, rowOrder, "slope", res.getSlopes(), "adj. p-value", res.getAdjPValues());
		String heatmapFilename = getOutdir() + "/" + test + "_heatmap";
		try {
			heatmap.save(heatmapFilename);
		} catch (IOException e) {
			printError("ioexception_executecorrelation_correlation", "error generating heatmap", e.toString(), e);
		}

		// saving data
		//
		updateJobStatus("80", "saving results");
		DataFrame dataFrame = new DataFrame(dataset.getFeatureNames().size(), 0);
		dataFrame.setRowNames(ListUtils.ordered(dataset.getFeatureNames(), rowOrder));

		try {
			dataFrame.addColumn("statistic", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getStatistics()), rowOrder)));
			dataFrame.addColumn("slope", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getSlopes()), rowOrder)));
			dataFrame.addColumn("intercept", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getIntercepts()), rowOrder)));
			dataFrame.addColumn("p-value", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getPValues()), rowOrder)));
			dataFrame.addColumn("adj. p-value", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getAdjPValues()), rowOrder)));

			File outFile = new File(getOutdir() + "/" + test + ".txt");
			FeatureData featureData = new FeatureData(dataFrame);
			featureData.write(outFile);
			result.addOutputItem(new Item(test + "_file", outFile.getName(), test.toUpperCase() + " output file", TYPE.FILE, new ArrayList<String>(2), new HashMap<String, String>(2), "Output files"));

			outFile = new File(getOutdir() + "/" + test + "_table.txt");
			IOUtils.write(outFile, dataFrame.toString(true, true));			
			List<String> tags = new ArrayList<String>();
			tags.add("TABLE");
			result.addOutputItem(new Item(test + "_table", outFile.getName(), test.toUpperCase() + " output table", TYPE.FILE, tags, new HashMap<String, String>(2), "Output files"));
		} catch (Exception e) {
			printError("ioexception_executecorrelation_correlation", "error saving " + test + " results", e.toString(), e);
		}

		if ( new File(heatmapFilename + ".png").exists() ) {
			result.addOutputItem(new Item(test + "_heatmap", test + "_heatmap.png", test.toUpperCase() + " heatmap", TYPE.IMAGE, new ArrayList<String>(2), new HashMap<String, String>(2), "Heatmap image"));
		}		
	}
}
