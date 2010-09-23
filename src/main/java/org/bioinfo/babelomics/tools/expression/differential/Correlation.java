package org.bioinfo.babelomics.tools.expression.differential;

import java.io.File;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.commons.utils.ArrayUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.data.dataset.FeatureData;
import org.bioinfo.data.list.DataFrame;
import org.bioinfo.math.regression.SimpleRegressionTest;
import org.bioinfo.math.result.CorrelationTestResult;
import org.bioinfo.math.result.SimpleRegressionTestResult;
import org.bioinfo.math.result.TestResultList;
import org.bioinfo.math.stats.correlation.CorrelationTest2;
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
	
	private int maxDisplay = 500;	

	public Correlation() {
		initOptions();
	}

	public void initOptions() {
		options.addOption(OptionFactory.createOption("dataset", "the data"));
		options.addOption(OptionFactory.createOption("test", "the test, possible values: pearson, spearman, regression"));
		options.addOption(OptionFactory.createOption("class-name", "class variable", false));
		options.addOption(OptionFactory.createOption("correction", "Multiple-test correction: fdr, bh, by, bonferroni, hochberg, holm", false));
		options.addOption(OptionFactory.createOption("p-value", "p-value for significative genes", false));
	}

	@Override
	public void execute() {

		// reading dataset
		//
		dataset = initDataset(new File(commandLine.getOptionValue("dataset")));

		test = commandLine.getOptionValue("test", null);
		className = commandLine.getOptionValue("class-name", null);
		correction = commandLine.getOptionValue("correction", null);
		String pValueParam = commandLine.getOptionValue("p-value", "0.05");

		try {
			vars = dataset.getVariables().getByName(className).getValues();
		} catch (Exception e) {
			vars = new ArrayList<String>();
		}

		try {
			pValue = Double.parseDouble(pValueParam);
			if (pValue > 1 || pValue < 0) {
				pValue = 0.05;
			}
		} catch (NumberFormatException e) {
			pValue = 0.05;
		}			

		// input parameters
		//
		result.addOutputItem(new Item("test_input_param", test, "Test", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
		result.addOutputItem(new Item("correction_input_param", correction, "Multiple-test correction", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
		result.addOutputItem(new Item("class_input_param", className + " [" + ListUtils.toString(vars, ",") + "]", "Class", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
		result.addOutputItem(new Item("pvalue_input_param", pValueParam, "Adjusted p-value", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));

		if ( ! "pearson".equalsIgnoreCase(test) && ! "spearman".equalsIgnoreCase(test) && ! "regression".equalsIgnoreCase(test) ) {
			abort("unknowntest_execute_correlation", "unknown test (" + test + ")", "unknown test (" + test + ")", "unknown test (" + test + ")");
		}

		String aux = "";
		if (vars.size()>0) {
			try {
				doubleVars = new ArrayList<Double>(vars.size());
				for(int i=0 ; i<vars.size() ; i++) {
					aux = vars.get(i);
					doubleVars.add(Double.parseDouble(aux));
				}
			} catch (NumberFormatException e) {
				abort("exception_execute_correlation", "Error", "Invalid class value: " + aux + ". Class values must be numeric for correlation tests", "");
			}
		} else {
			abort("exception_execute_correlation", "Error", "Missing values for class " + className, "");			
		}


		if ( "regression".equalsIgnoreCase(test) ) {
			executeRegression();
		} else {
			executeCorrelation();
		}
	}


	private void executeCorrelation() {
		try {
			updateJobStatus("40", "computing " + test + "correlation");
			
			// apply regression and multiple test correction according
			//
			
//			DoubleMatrix dm = new DoubleMatrix(dataset.getDoubleMatrix().getRowDimension(), 1);
//			dm.setColumn(0, dataset.getDoubleMatrix().getColumn(0));
//			dataset.setDoubleMatrix(dm);
			//CorrelationTest corrTest = new CorrelationTest(dataset.getDoubleMatrix(), doubleVars, test);
			
			CorrelationTest2 corrTest = new CorrelationTest2(dataset.getDoubleMatrix(), doubleVars, test);
			TestResultList<CorrelationTestResult> res = corrTest.compute();
			DiffExpressionUtils.multipleTestCorrection(res, correction);			

			// create output file
			//
			//int[] columnOrder = ListUtils.order(vars);
			int[] columnOrder = ListUtils.order(dataset.getVariables().getByName(className).getValues());
			int[] rowOrder = ListUtils.order(ArrayUtils.toList(res.getStatistics()), true);

			DataFrame dataFrame = new DataFrame(dataset.getFeatureNames().size(), 0);
			dataFrame.addColumn("statistic", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getStatistics()), rowOrder)));
			dataFrame.addColumn("correlation", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getCorrelations()), rowOrder)));
			dataFrame.addColumn("p-value", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getPValues()), rowOrder)));
			dataFrame.addColumn("adj. p-value", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getAdjPValues()), rowOrder)));
			dataFrame.setRowNames(ListUtils.ordered(dataset.getFeatureNames(), rowOrder));

			FeatureData featureData = new FeatureData(dataFrame);
			File file = new File(outdir + "/" + test + ".txt");
			featureData.save(file);
			if ( file.exists() ) {
				result.addOutputItem(new Item(test + "_correlation_file", file.getName(), "Correlation output file", TYPE.FILE, new ArrayList<String>(), new HashMap<String, String>(), "Correlation output files"));							
			}

			// getting significative genes
			//
			DiffExpressionUtils.addSignificativeResults(dataset, test, "statistic", res.getStatistics(), "adj. p-value", res.getAdjPValues(), "p-value", res.getPValues(), "correlation", res.getCorrelations(), null, null, className, columnOrder, pValue, maxDisplay, this);
			DiffExpressionUtils.createFatiScanRedirection(dataFrame, test, "statistic", result, outdir);			
		} catch (Exception e) {
			e.printStackTrace();
			abort("exception_execute_" + test + "_correlation", "ERROR", "Error running " + test + " correlation", "");
		}		

//////////////////////////////////////////////		
//		// running test (pearson or spearman)
//		//
//		updateJobStatus("40", "computing " + test + " test");		
//		TestResultList<CorrelationTestResult> res = null;
//		try {
//			CorrelationTest corrTest = new CorrelationTest(dataset.getDoubleMatrix(), doubleVars, test);
//			res = corrTest.compute();
//
//			// apply multiple test correction according to input correction
//			DiffExpressionUtils.multipleTestCorrection(res, correction);			
//		} catch (Exception e) {
//			abort("exception_executecorrelation_correlation", "error running " + test + " test", e.toString(), StringUtils.getStackTrace(e));
//		}
//
//		int[] columnOrder = ListUtils.order(vars);
//		int[] rowOrder = ListUtils.order(ArrayUtils.toList(res.getCorrelations()), true);
//
//		// generating heatmap
//		//
//		Canvas heatmap = null, sigHeatmap = null;
//		try {
//			updateJobStatus("60", "generating heatmap");
//			//heatmap = DiffExpressionUtils.generateHeatmap(dataset, className, columnOrder, rowOrder, "correlation", res.getCorrelations(), "adj. p-value", res.getAdjPValues());
//			sigHeatmap = DiffExpressionUtils.generateSigHeatmap(dataset, className, columnOrder, "correlation", res.getCorrelations(), "adj. p-value", res.getAdjPValues(), pValue);
//		} catch (Exception e) {
//			printError("ioexception_executecorrelation_correlation", "error generating heatmaps for " + test + " test", e.toString(), e);
//		}
//
//		// saving data
//		//
//		updateJobStatus("80", "saving results");
//
//		DataFrame dataFrame = new DataFrame(dataset.getFeatureNames().size(), 0);
//
//		try {
//			dataFrame.addColumn("statistic", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getStatistics()), rowOrder)));
//			dataFrame.addColumn("correlation", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getCorrelations()), rowOrder)));
//			dataFrame.addColumn("p-value", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getPValues()), rowOrder)));
//			dataFrame.addColumn("adj. p-value", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getAdjPValues()), rowOrder)));
//
//			dataFrame.setRowNames(ListUtils.ordered(dataset.getFeatureNames(), rowOrder));
//
//			File outFile = new File(getOutdir() + "/" + test + ".txt");
//			FeatureData featureData = new FeatureData(dataFrame);
//			featureData.save(outFile);
//			//			result.addOutputItem(new Item(test + "_file", outFile.getName(), test.toUpperCase() + " output file", TYPE.FILE, new ArrayList<String>(2), new HashMap<String, String>(2), "Output files"));
//			result.addOutputItem(new Item(test + "_file", outFile.getName(), test.toUpperCase() + " output file", TYPE.FILE, StringUtils.toList("TABLE,CORRELATION_TABLE", ","), new HashMap<String, String>(2), "Output files"));
//
//			//			outFile = new File(getOutdir() + "/" + test + "_table.txt");
//			//			IOUtils.write(outFile, dataFrame.toString(true, true));			
//			//			result.addOutputItem(new Item(test + "_table", outFile.getName(), test.toUpperCase() + " output table", TYPE.FILE, StringUtils.toList("TABLE,CORRELATION_TABLE", ","), new HashMap<String, String>(2), "Output files"));
//
//			DiffExpressionUtils.addOutputLists(dataFrame, test, "correlation", result, outdir);
//		} catch (Exception e) {
//			printError("ioexception_executecorrelation_correlation", "error saving " + test + " results", e.toString(), e);
//		}
//
//		String sigHeatmapFilename = getOutdir() + "/" + test + "_heatmap_significative.png";
//		if ( sigHeatmap == null ) {
//			result.addOutputItem(new Item(test + "_heatmap_significative", "None significative terms", test.toUpperCase() + " heatmap with significative terms (p-value = " + pValue + ")", TYPE.TEXT, new ArrayList<String>(2), new HashMap<String, String>(2), "Heatmap images"));
//		} else {
//			try {
//				sigHeatmap.save(sigHeatmapFilename);
//				if ( new File(sigHeatmapFilename).exists() ) {
//					result.addOutputItem(new Item(test + "_heatmap_significative", test + "_heatmap_significative.png", test.toUpperCase() + " heatmap with significative terms (p-value = " + pValue + ")", TYPE.IMAGE, new ArrayList<String>(2), new HashMap<String, String>(2), "Heatmap images"));
//				}
//			} catch (IOException e) {
//				printError("ioexception_executet_classcomparison", "error generating heatmap", e.toString(), e);
//			}
//		}
//
//		//		String heatmapFilename = getOutdir() + "/" + test + "_heatmap.png";
//		//		try {
//		//			heatmap.save(heatmapFilename);
//		//			if ( new File(heatmapFilename).exists() ) {
//		//				result.addOutputItem(new Item(test + "_heatmap", test + "_heatmap.png", test.toUpperCase() + " heatmap", TYPE.FILE, new ArrayList<String>(2), new HashMap<String, String>(2), "Heatmap images"));
//		//			}
//		//		} catch (IOException e) {
//		//			printError("ioexception_executecorrelation_correlation", "error generating heatmap", e.toString(), e);
//		//		}
//////////////////////////////////////////////
	}


	private void executeRegression() {
		try {
			updateJobStatus("40", "computing regression");
			
			// apply regression and multiple test correction according
			//
			SimpleRegressionTest regression = new SimpleRegressionTest();
			TestResultList<SimpleRegressionTestResult> res = regression.compute(dataset.getDoubleMatrix(), doubleVars);
			DiffExpressionUtils.multipleTestCorrection(res, correction);			

			// create output file
			//
			//int[] columnOrder = ListUtils.order(vars);
			int[] columnOrder = ListUtils.order(dataset.getVariables().getByName(className).getValues());
			int[] rowOrder = ListUtils.order(ArrayUtils.toList(res.getStatistics()), true);

			DataFrame dataFrame = new DataFrame(dataset.getFeatureNames().size(), 0);
			dataFrame.addColumn("statistic", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getStatistics()), rowOrder)));
			dataFrame.addColumn("slope", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getSlopes()), rowOrder)));
			dataFrame.addColumn("intercept", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getIntercepts()), rowOrder)));
			dataFrame.addColumn("p-value", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getPValues()), rowOrder)));
			dataFrame.addColumn("adj. p-value", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getAdjPValues()), rowOrder)));
			dataFrame.setRowNames(ListUtils.ordered(dataset.getFeatureNames(), rowOrder));

			FeatureData featureData = new FeatureData(dataFrame);
			File file = new File(outdir + "/" + test + ".txt");
			featureData.save(file);
			if ( file.exists() ) {
				result.addOutputItem(new Item(test + "_file", file.getName(), "Regression output file", TYPE.FILE, new ArrayList<String>(), new HashMap<String, String>(), "Regression output files"));							
			}

			// getting significative genes
			//
			DiffExpressionUtils.addSignificativeResults(dataset, test, "statistic", res.getStatistics(), "adj. p-value", res.getAdjPValues(), "p-value", res.getPValues(), "slope", res.getSlopes(), "intercept", res.getIntercepts(), className, columnOrder, pValue, maxDisplay, this);
			DiffExpressionUtils.createFatiScanRedirection(dataFrame, test, "statistic", result, outdir);			
		} catch (Exception e) {
			e.printStackTrace();
			abort("exception_execute_regression", "ERROR", "Error running simple regression", "");
		}		
	}
	
	


//	private void executeRegression() {
//
//		// running test (regression)
//		//
//		updateJobStatus("40", "computing " + test + " test");		
//		TestResultList<SimpleRegressionTestResult> res = null;
//		try {
//			SimpleRegressionTest regression = new SimpleRegressionTest();
//			res = regression.compute(dataset.getDoubleMatrix(), doubleVars);
//
//			// apply multiple test correction according to input correction
//			DiffExpressionUtils.multipleTestCorrection(res, correction);			
//		} catch (Exception e) {
//			abort("exception_executecorrelation_correlation", "error running " + test + " test", e.toString(), StringUtils.getStackTrace(e));
//		}
//
//		int[] columnOrder = ListUtils.order(vars);
//		int[] rowOrder = ListUtils.order(ArrayUtils.toList(res.getStatistics()), true);
//
//		// generating heatmap
//		//
//		updateJobStatus("60", "generating heatmap");
//		Canvas heatmap = null, sigHeatmap = null;
//		try {
//			//heatmap = DiffExpressionUtils.generateHeatmap(dataset, className, columnOrder, rowOrder, "slope", res.getSlopes(), "adj. p-value", res.getAdjPValues());
//			sigHeatmap = DiffExpressionUtils.generateSigHeatmap(dataset, className, columnOrder, "slope", res.getSlopes(), "adj. p-value", res.getAdjPValues(), pValue);
//		} catch (Exception e) {
//			printError("ioexception_executecorrelation_correlation", "error generating heatmaps for " + test + " test", e.toString(), e);
//		}
//
//		// saving data
//		//
//		updateJobStatus("80", "saving results");
//		DataFrame dataFrame = new DataFrame(dataset.getFeatureNames().size(), 0);
//
//		try {
//			dataFrame.addColumn("statistic", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getStatistics()), rowOrder)));
//			dataFrame.addColumn("slope", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getSlopes()), rowOrder)));
//			dataFrame.addColumn("intercept", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getIntercepts()), rowOrder)));
//			dataFrame.addColumn("p-value", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getPValues()), rowOrder)));
//			dataFrame.addColumn("adj. p-value", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getAdjPValues()), rowOrder)));
//
//			dataFrame.setRowNames(ListUtils.ordered(dataset.getFeatureNames(), rowOrder));
//
//			File outFile = new File(getOutdir() + "/" + test + ".txt");
//			FeatureData featureData = new FeatureData(dataFrame);
//			featureData.save(outFile);
//			//			result.addOutputItem(new Item(test + "_file", outFile.getName(), test.toUpperCase() + " output file", TYPE.FILE, new ArrayList<String>(2), new HashMap<String, String>(2), "Output files"));
//			result.addOutputItem(new Item(test + "_file", outFile.getName(), test.toUpperCase() + " output file", TYPE.FILE, StringUtils.toList("TABLE,REGRESSION_TABLE", ","), new HashMap<String, String>(2), "Output files"));
//
//			//			outFile = new File(getOutdir() + "/" + test + "_table.txt");
//			//			IOUtils.write(outFile, dataFrame.toString(true, true));			
//			//			result.addOutputItem(new Item(test + "_table", outFile.getName(), test.toUpperCase() + " output table", TYPE.FILE, StringUtils.toList("TABLE,REGRESSION_TABLE", ","), new HashMap<String, String>(2), "Output files"));
//
//			DiffExpressionUtils.addOutputLists(dataFrame, test, "statistic", result, outdir);
//
//		} catch (Exception e) {
//			printError("ioexception_executecorrelation_correlation", "error saving " + test + " results", e.toString(), e);
//		}
//
//		String sigHeatmapFilename = getOutdir() + "/" + test + "_heatmap_significative.png";
//		if ( sigHeatmap == null ) {
//			result.addOutputItem(new Item(test + "_heatmap_significative", "None significative terms", test.toUpperCase() + " heatmap with significative terms (p-value = " + pValue + ")", TYPE.TEXT, new ArrayList<String>(2), new HashMap<String, String>(2), "Heatmap images"));
//		} else {
//			try {
//				sigHeatmap.save(sigHeatmapFilename);
//				if ( new File(sigHeatmapFilename).exists() ) {
//					result.addOutputItem(new Item(test + "_heatmap_significative", test + "_heatmap_significative.png", test.toUpperCase() + " heatmap with significative terms (p-value = " + pValue + ")", TYPE.IMAGE, new ArrayList<String>(2), new HashMap<String, String>(2), "Heatmap images"));
//				}
//			} catch (IOException e) {
//				printError("ioexception_executet_classcomparison", "error generating heatmap", e.toString(), e);
//			}
//		}
//
//		//		String heatmapFilename = getOutdir() + "/" + test + "_heatmap.png";
//		//		try {
//		//			heatmap.save(heatmapFilename);
//		//			if ( new File(heatmapFilename).exists() ) {
//		//				result.addOutputItem(new Item(test + "_heatmap", test + "_heatmap.png", test.toUpperCase() + " heatmap", TYPE.IMAGE, new ArrayList<String>(2), new HashMap<String, String>(2), "Heatmap images"));
//		//			}		
//		//		} catch (IOException e) {
//		//			printError("ioexception_executecorrelation_correlation", "error generating heatmap", e.toString(), e);
//		//		}
//	}
}
