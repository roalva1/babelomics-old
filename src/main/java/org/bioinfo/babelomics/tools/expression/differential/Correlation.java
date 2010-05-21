package org.bioinfo.babelomics.tools.expression.differential;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ArrayUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.data.dataset.FeatureData;
import org.bioinfo.data.list.DataFrame;
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
	String correction;
	private double pValue = 0.05;

	public Correlation() {
	}

	public void initOptions() {
		options.addOption(OptionFactory.createOption("dataset", "the data"));
		options.addOption(OptionFactory.createOption("test", "the test, possible values: pearson, spearman, regression"));
		options.addOption(OptionFactory.createOption("classname", "class variable", false));
		options.addOption(OptionFactory.createOption("correction", "Multiple-test correction: fdr, bh, by, bonferroni, hochberg, hold"));
	}

	@Override
	public void execute() {

		// reading dataset
		//
		dataset = initDataset(new File(commandLine.getOptionValue("dataset")));

		test = commandLine.getOptionValue("test", null);
		className = commandLine.getOptionValue("classname", null);
		correction = commandLine.getOptionValue("correction", null);

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

			// apply multiple test correction according to input correction
			DiffExpressionUtils.multipleTestCorrection(res, correction);			
		} catch (Exception e) {
			abort("exception_executecorrelation_correlation", "error running " + test + " test", e.toString(), StringUtils.getStackTrace(e));
		}

		int[] columnOrder = ListUtils.order(vars);
		int[] rowOrder = ListUtils.order(ArrayUtils.toList(res.getCorrelations()), true);

		// generating heatmap
		//
		Canvas heatmap = null, sigHeatmap = null;
		try {
			updateJobStatus("60", "generating heatmap");
			heatmap = DiffExpressionUtils.generateHeatmap(dataset, className, columnOrder, rowOrder, "correlation", res.getCorrelations(), "adj. p-value", res.getAdjPValues());
			sigHeatmap = DiffExpressionUtils.generateSigHeatmap(dataset, className, columnOrder, "correlation", res.getCorrelations(), "adj. p-value", res.getAdjPValues(), pValue);
		} catch (Exception e) {
			printError("ioexception_executecorrelation_correlation", "error generating heatmaps for " + test + " test", e.toString(), e);
		}

		// saving data
		//
		updateJobStatus("80", "saving results");

		DataFrame dataFrame = new DataFrame(dataset.getFeatureNames().size(), 0);

		try {
			dataFrame.addColumn("statistic", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getStatistics()), rowOrder)));
			dataFrame.addColumn("correlation", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getCorrelations()), rowOrder)));
			dataFrame.addColumn("p-value", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getPValues()), rowOrder)));
			dataFrame.addColumn("adj. p-value", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getAdjPValues()), rowOrder)));

			dataFrame.setRowNames(ListUtils.ordered(dataset.getFeatureNames(), rowOrder));

			File outFile = new File(getOutdir() + "/" + test + ".txt");
			FeatureData featureData = new FeatureData(dataFrame);
			featureData.save(outFile);
//			result.addOutputItem(new Item(test + "_file", outFile.getName(), test.toUpperCase() + " output file", TYPE.FILE, new ArrayList<String>(2), new HashMap<String, String>(2), "Output files"));
			result.addOutputItem(new Item(test + "_file", outFile.getName(), test.toUpperCase() + " output file", TYPE.FILE, StringUtils.toList("TABLE,CORRELATION_TABLE", ","), new HashMap<String, String>(2), "Output files"));

//			outFile = new File(getOutdir() + "/" + test + "_table.txt");
//			IOUtils.write(outFile, dataFrame.toString(true, true));			
//			result.addOutputItem(new Item(test + "_table", outFile.getName(), test.toUpperCase() + " output table", TYPE.FILE, StringUtils.toList("TABLE,CORRELATION_TABLE", ","), new HashMap<String, String>(2), "Output files"));

			DiffExpressionUtils.addOutputLists(dataFrame, test, "correlation", result, outdir);
		} catch (Exception e) {
			printError("ioexception_executecorrelation_correlation", "error saving " + test + " results", e.toString(), e);
		}

		String sigHeatmapFilename = getOutdir() + "/" + test + "_heatmap_significative.png";
		if ( sigHeatmap == null ) {
			result.addOutputItem(new Item(test + "_heatmap_significative", "None significative terms", test.toUpperCase() + " heatmap with significative terms (p-value = " + pValue + ")", TYPE.TEXT, new ArrayList<String>(2), new HashMap<String, String>(2), "Heatmap images"));
		} else {
			try {
				sigHeatmap.save(sigHeatmapFilename);
				if ( new File(sigHeatmapFilename).exists() ) {
					result.addOutputItem(new Item(test + "_heatmap_significative", test + "_heatmap_significative.png", test.toUpperCase() + " heatmap with significative terms (p-value = " + pValue + ")", TYPE.IMAGE, new ArrayList<String>(2), new HashMap<String, String>(2), "Heatmap images"));
				}
			} catch (IOException e) {
				printError("ioexception_executet_classcomparison", "error generating heatmap", e.toString(), e);
			}
		}

//		String heatmapFilename = getOutdir() + "/" + test + "_heatmap.png";
//		try {
//			heatmap.save(heatmapFilename);
//			if ( new File(heatmapFilename).exists() ) {
//				result.addOutputItem(new Item(test + "_heatmap", test + "_heatmap.png", test.toUpperCase() + " heatmap", TYPE.FILE, new ArrayList<String>(2), new HashMap<String, String>(2), "Heatmap images"));
//			}
//		} catch (IOException e) {
//			printError("ioexception_executecorrelation_correlation", "error generating heatmap", e.toString(), e);
//		}
	}




	private void executeRegression() {

		// running test (regression)
		//
		updateJobStatus("40", "computing " + test + " test");		
		TestResultList<SimpleRegressionTestResult> res = null;
		try {
			SimpleRegressionTest regression = new SimpleRegressionTest();
			res = regression.compute(dataset.getDoubleMatrix(), doubleVars);

			// apply multiple test correction according to input correction
			DiffExpressionUtils.multipleTestCorrection(res, correction);			
		} catch (Exception e) {
			abort("exception_executecorrelation_correlation", "error running " + test + " test", e.toString(), StringUtils.getStackTrace(e));
		}

		int[] columnOrder = ListUtils.order(vars);
		int[] rowOrder = ListUtils.order(ArrayUtils.toList(res.getStatistics()), true);

		// generating heatmap
		//
		updateJobStatus("60", "generating heatmap");
		Canvas heatmap = null, sigHeatmap = null;
		try {
			heatmap = DiffExpressionUtils.generateHeatmap(dataset, className, columnOrder, rowOrder, "slope", res.getSlopes(), "adj. p-value", res.getAdjPValues());
			sigHeatmap = DiffExpressionUtils.generateSigHeatmap(dataset, className, columnOrder, "slope", res.getSlopes(), "adj. p-value", res.getAdjPValues(), pValue);
		} catch (Exception e) {
			printError("ioexception_executecorrelation_correlation", "error generating heatmaps for " + test + " test", e.toString(), e);
		}

		// saving data
		//
		updateJobStatus("80", "saving results");
		DataFrame dataFrame = new DataFrame(dataset.getFeatureNames().size(), 0);

		try {
			dataFrame.addColumn("statistic", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getStatistics()), rowOrder)));
			dataFrame.addColumn("slope", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getSlopes()), rowOrder)));
			dataFrame.addColumn("intercept", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getIntercepts()), rowOrder)));
			dataFrame.addColumn("p-value", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getPValues()), rowOrder)));
			dataFrame.addColumn("adj. p-value", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getAdjPValues()), rowOrder)));

			dataFrame.setRowNames(ListUtils.ordered(dataset.getFeatureNames(), rowOrder));

			File outFile = new File(getOutdir() + "/" + test + ".txt");
			FeatureData featureData = new FeatureData(dataFrame);
			featureData.save(outFile);
//			result.addOutputItem(new Item(test + "_file", outFile.getName(), test.toUpperCase() + " output file", TYPE.FILE, new ArrayList<String>(2), new HashMap<String, String>(2), "Output files"));
			result.addOutputItem(new Item(test + "_file", outFile.getName(), test.toUpperCase() + " output file", TYPE.FILE, StringUtils.toList("TABLE,REGRESSION_TABLE", ","), new HashMap<String, String>(2), "Output files"));

//			outFile = new File(getOutdir() + "/" + test + "_table.txt");
//			IOUtils.write(outFile, dataFrame.toString(true, true));			
//			result.addOutputItem(new Item(test + "_table", outFile.getName(), test.toUpperCase() + " output table", TYPE.FILE, StringUtils.toList("TABLE,REGRESSION_TABLE", ","), new HashMap<String, String>(2), "Output files"));

			DiffExpressionUtils.addOutputLists(dataFrame, test, "statistic", result, outdir);

		} catch (Exception e) {
			printError("ioexception_executecorrelation_correlation", "error saving " + test + " results", e.toString(), e);
		}

		String sigHeatmapFilename = getOutdir() + "/" + test + "_heatmap_significative.png";
		if ( sigHeatmap == null ) {
			result.addOutputItem(new Item(test + "_heatmap_significative", "None significative terms", test.toUpperCase() + " heatmap with significative terms (p-value = " + pValue + ")", TYPE.TEXT, new ArrayList<String>(2), new HashMap<String, String>(2), "Heatmap images"));
		} else {
			try {
				sigHeatmap.save(sigHeatmapFilename);
				if ( new File(sigHeatmapFilename).exists() ) {
					result.addOutputItem(new Item(test + "_heatmap_significative", test + "_heatmap_significative.png", test.toUpperCase() + " heatmap with significative terms (p-value = " + pValue + ")", TYPE.IMAGE, new ArrayList<String>(2), new HashMap<String, String>(2), "Heatmap images"));
				}
			} catch (IOException e) {
				printError("ioexception_executet_classcomparison", "error generating heatmap", e.toString(), e);
			}
		}

//		String heatmapFilename = getOutdir() + "/" + test + "_heatmap.png";
//		try {
//			heatmap.save(heatmapFilename);
//			if ( new File(heatmapFilename).exists() ) {
//				result.addOutputItem(new Item(test + "_heatmap", test + "_heatmap.png", test.toUpperCase() + " heatmap", TYPE.IMAGE, new ArrayList<String>(2), new HashMap<String, String>(2), "Heatmap images"));
//			}		
//		} catch (IOException e) {
//			printError("ioexception_executecorrelation_correlation", "error generating heatmap", e.toString(), e);
//		}
	}
}
