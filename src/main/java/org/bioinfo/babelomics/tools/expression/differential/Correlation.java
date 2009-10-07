package org.bioinfo.babelomics.tools.expression.differential;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.bioinfo.babelomics.tools.expression.DifferentialAnalysis;
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

public class Correlation extends DifferentialAnalysis {

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
			result.addOutputItem(new Item(test + "_file", outFile.getName(), test.toUpperCase() + " output file", TYPE.FILE));

			outFile = new File(getOutdir() + "/" + test + "_table.txt");
			IOUtils.write(outFile, dataFrame.toString(true, true));			
			List<String> tags = new ArrayList<String>();
			tags.add("TABLE");
			result.addOutputItem(new Item(test + "_table", outFile.getName(), test.toUpperCase() + " output table", TYPE.FILE, tags, new HashMap<String, String>(2)));
		} catch (Exception e) {
			printError("ioexception_executecorrelation_correlation", "error saving " + test + " results", e.toString(), e);
		}

		if ( new File(heatmapFilename + ".png").exists() ) {
			result.addOutputItem(new Item(test + "_heatmap", test + "_heatmap.png", test.toUpperCase() + " heatmap", TYPE.IMAGE));
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
			result.addOutputItem(new Item(test + "_file", outFile.getName(), test.toUpperCase() + " output file", TYPE.FILE));

			outFile = new File(getOutdir() + "/" + test + "_table.txt");
			IOUtils.write(outFile, dataFrame.toString(true, true));			
			List<String> tags = new ArrayList<String>();
			tags.add("TABLE");
			result.addOutputItem(new Item(test + "_table", outFile.getName(), test.toUpperCase() + " output table", TYPE.FILE, tags, new HashMap<String, String>(2)));
		} catch (Exception e) {
			printError("ioexception_executecorrelation_correlation", "error saving " + test + " results", e.toString(), e);
		}

		if ( new File(heatmapFilename + ".png").exists() ) {
			result.addOutputItem(new Item(test + "_heatmap", test + "_heatmap.png", test.toUpperCase() + " heatmap", TYPE.IMAGE));
		}		
	}





















	//	
	//	@Override
	//	public void execute() {
	//		initDataset();
	//				
	//		String test = commandLine.getOptionValue("test", null);
	//		String className = commandLine.getOptionValue("class", null);
	//
	//		if ( ! "pearson".equalsIgnoreCase(test) && ! "spearman".equalsIgnoreCase(test) && ! "regression".equalsIgnoreCase(test) ) {
	//			abort("unknowntest_execute_correlation", "unknown test (" + test + ")", "unknown test (" + test + ")", "unknown test (" + test + ")");
	//		}
	//		
	//		if ( "none".equalsIgnoreCase(className) ) {
	//			abort("execute_correlation", "Missing independent variable to test", "Missing independent variable to test", "Missing independent variable to test");
	//		}
	//		
	//		List<String> vars = dataset.getVariables().getByName(className).getValues();
	//		List<Double> doubleVars = new ArrayList<Double>(vars.size());
	//		for(String str: vars) {
	//			doubleVars.add(Double.parseDouble(str));
	//		}
	//
	//		// running test
	//		//
	//		logger.debug("computing " + test + "...\n");
	//		try {
	//			jobStatus.addStatusMessage("40", "computing " + test);
	//		} catch (FileNotFoundException e) {
	//			abort("filenotfoundexception_execute_correlation", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
	//		}
	//
	//		if ( "regression".equalsIgnoreCase(test) ) {
	//			executeRegression();
	//			SimpleRegressionTest regression = new SimpleRegressionTest();
	//			TestResultList<SimpleRegressionTestResult> res = regression.compute(dataset.getDoubleMatrix(), doubleVars);			
	//		} else {
	//			CorrelationTest corrTest = new CorrelationTest(dataset.getDoubleMatrix(), doubleVars, test);
	//			TestResultList<CorrelationTestResult> res = corrTest.compute();
	//		}
	//		
	//		
	//		CoxTest coxtest = new CoxTest();
	//		//System.out.println("input matrix = \n" + dataset.getDoubleMatrix().toString());
	//		TestResultList<CoxTestResult> res = null;
	//		try {
	//			res = coxtest.compute(dataset.getDoubleMatrix(), timeVars, censoredVars);
	//		} catch (Exception e) {
	//			abort("exception_run_cox", "error running cox test", e.toString(), StringUtils.getStackTrace(e));
	//		}
	//
	//		int[] columnOrder = ListUtils.order(vars);
	//		int[] rowOrder = ListUtils.order(ListUtils.toList(res.getStatistics()), true);
	//
	//		// generating heatmap
	//		//
	//		logger.debug("generating heatmap...\n");
	//		try {
	//			jobStatus.addStatusMessage("60", "generating heatmap");
	//		} catch (FileNotFoundException e) {
	//			abort("filenotfoundexception_run_cox", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
	//		}
	//
	//		Canvas heatmap = generateHeatmap(dataset, timeClass, columnOrder, rowOrder, "coeff.", res.getCoefs(), "adj. p-value", res.getAdjPValues());
	//		String heatmapFilename = getOutdir() + "/cox_heatmap";
	//		try {
	//			heatmap.save(heatmapFilename);
	//		} catch (IOException e) {
	//			printError("ioexception_cox_cox", "error generating heatmap", e.toString(), e);
	//		}
	//
	//		// saving data
	//		//
	//		logger.debug("saving results...");
	//		try {
	//			jobStatus.addStatusMessage("80", "saving results");
	//		} catch (FileNotFoundException e) {
	//			abort("filenotfoundexception_run_cox", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
	//		}
	//
	//		DataFrame dataFrame = new DataFrame(dataset.getFeatureNames().size(), 0);
	//		dataFrame.setRowNames(ListUtils.ordered(dataset.getFeatureNames(), rowOrder));
	//
	//		//dataFrame.addColumn("id", ListUtils.ordered(dataset.getFeatureNames(), rowOrder));
	//		try {
	//			dataFrame.addColumn("statistic", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getStatistics()), rowOrder)));
	//			dataFrame.addColumn("coeff.", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getCoefs()), rowOrder)));
	//			dataFrame.addColumn("p-value", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getPValues()), rowOrder)));
	//			dataFrame.addColumn("adj. p-value", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getAdjPValues()), rowOrder)));
	//
	//			FeatureData featureData = new FeatureData(dataFrame);
	//			featureData.write(new File(getOutdir() + "/cox.txt"));
	//			result.addOutputItem(new Item("cox_file", "cox.txt", "Cox output file", TYPE.FILE));
	//			
	//			IOUtils.write(new File(getOutdir() + "/cox_table.txt"), dataFrame.toString(true, true));			
	//			List<String> tags = new ArrayList<String>();
	//			tags.add("TABLE");
	//			result.addOutputItem(new Item("cox_table", "cox_table.txt", "Cox table", TYPE.FILE, tags, new HashMap<String, String>(2)));
	//		} catch (Exception e) {
	//			printError("ioexception_cox_cox", "error saving results", e.toString(), e);
	//		}
	//
	//		if ( new File(heatmapFilename + ".png").exists() ) {
	//			result.addOutputItem(new Item("cox_heatmap", "cox_heatmap.png", "Cox heatmap", TYPE.IMAGE));
	//		}
	//	}
	//
	//	
	//	
	//	
	//	
	//	
	//	
	//	
	//	
	//	private void executeCorrelation(Dataset dataset, String className, String test) {
	//	
	//
	//		CorrelationTest corrTest = new CorrelationTest(dataset.getDoubleMatrix(), doubleVars, test);
	//		TestResultList<CorrelationTestResult> res = corrTest.compute();
	//
	//		int[] columnOrder = ListUtils.order(vars);
	//		int[] rowOrder = ListUtils.order(ListUtils.toList(res.getCorrelations()), true);
	//
	//		// generating heatmap
	//		//
	//		jobStatus.addStatusMessage("60", "generating heatmap");
	//		logger.debug("generating heatmap...\n");
	//
	//		Canvas heatmap = generateHeatmap(dataset, className, columnOrder, rowOrder, "correlation", res.getCorrelations(), "adj. p-value", res.getAdjPValues());
	//		heatmap.save(getOutdir() + "/" + test + "_heatmap");
	//		
	//		// saving data
	//		//
	//		jobStatus.addStatusMessage("80", "saving " + test + " results");
	//		logger.debug("saving results...");
	//
	//		DataFrame dataFrame = new DataFrame(dataset.getFeatureNames().size(), 0);
	//		dataFrame.setRowNames(ListUtils.ordered(dataset.getFeatureNames(), rowOrder));
	//
	//		//dataFrame.addColumn("id", ListUtils.ordered(dataset.getFeatureNames(), rowOrder));
	//		dataFrame.addColumn("statistic", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getStatistics()), rowOrder)));
	//		dataFrame.addColumn("correlation", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getCorrelations()), rowOrder)));
	//		dataFrame.addColumn("p-value", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getPValues()), rowOrder)));
	//		dataFrame.addColumn("adj. p-value", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getAdjPValues()), rowOrder)));
	//
	//		FeatureData featureData = new FeatureData(dataFrame);
	//		featureData.write(new File(getOutdir() + "/" + test + "_correlation.txt"));
	//
	//		IOUtils.write(new File(getOutdir() + "/" + test + "_table.txt"), dataFrame.toString(true, true));
	//
	//		Item item;
	//		
	//		item = new Item(test + "_correlation_file", test + "_correlation.txt", "Correlation file", TYPE.FILE);
	//		item.setGroup(test.toUpperCase() + " results");
	//		result.addOutputItem(item);
	//		
	//		item = new Item(test + "_table", test + "_table.txt", "Correlation table", TYPE.FILE);
	//		item.addTag("TABLE");
	//		item.setGroup(test.toUpperCase() +" results");
	//		result.addOutputItem(item);
	//
	//		item = new Item(test + "_heatmap", test + "_heatmap.png", test.toUpperCase() + " Heatmap", TYPE.IMAGE);
	//		item.setGroup(test.toUpperCase() + " results");
	//		result.addOutputItem(item);			
	//	} catch (java.security.InvalidParameterException e) {
	//		logger.error("not valid parameter: execute " + test + " correlation");
	//		e.printStackTrace();
	//	} catch (MatrixIndexException e) {
	//		logger.error("MatrixIndexException:  execute " + test + " correlation");
	//		// TODO Auto-generated catch block
	//		e.printStackTrace();
	//	} catch (MathException e) {
	//		// TODO Auto-generated catch block
	//		logger.error("math exception:  execute " + test + " correlation");
	//		e.printStackTrace();
	//	} catch (IOException e) {
	//		// TODO Auto-generated catch block
	//		logger.error("IOException:  execute " + test + " correlation");
	//		e.printStackTrace();
	//	} catch (InvalidColumnIndexException e) {
	//		// TODO Auto-generated catch block
	//		e.printStackTrace();
	//	}
	//}
	//
	//
	//private void executeRegression(Dataset dataset, List<Double> doubleVars) {
	//
	//		// regression test
	//		//
	//		jobStatus.addStatusMessage("40", "computing regression");
	//		logger.debug("computing regression...\n");
	//
	//		SimpleRegressionTest regression = new SimpleRegressionTest();
	//		TestResultList<SimpleRegressionTestResult> res = regression.compute(dataset.getDoubleMatrix(), doubleVars);
	//
	//		int[] columnOrder = ListUtils.order(doubleVars);
	//		int[] rowOrder = ListUtils.order(ListUtils.toList(res.getStatistics()), true);			
	//		
	//		// generating heatmap
	//		//
	//		jobStatus.addStatusMessage("60", "generating heatmap");
	//		logger.debug("generating heatmap...\n");
	//
	//		Canvas heatmap = generateHeatmap(dataset, className, columnOrder, rowOrder, "slope", res.getSlopes(), "adj. p-value", res.getAdjPValues());
	//		heatmap.save(getOutdir() + "/regression_heatmap");
	//		
	//		// saving data
	//		//
	//		jobStatus.addStatusMessage("80", "saving results");
	//		logger.debug("saving results...");
	//
	//		DataFrame dataFrame = new DataFrame(dataset.getFeatureNames().size(), 0);
	//		dataFrame.setRowNames(ListUtils.ordered(dataset.getFeatureNames(), rowOrder));
	//
	//		//dataFrame.addColumn("id", ListUtils.ordered(dataset.getFeatureNames(), rowOrder));
	//		dataFrame.addColumn("statistic", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getStatistics()), rowOrder)));
	//		dataFrame.addColumn("slope", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getSlopes()), rowOrder)));
	//		dataFrame.addColumn("intercept", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getIntercepts()), rowOrder)));
	//		dataFrame.addColumn("p-value", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getPValues()), rowOrder)));
	//		dataFrame.addColumn("adj. p-value", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getAdjPValues()), rowOrder)));
	//
	//		FeatureData featureData = new FeatureData(dataFrame);
	//		featureData.write(new File(getOutdir() + "/regression.txt"));
	//
	//		IOUtils.write(new File(getOutdir() + "/regression_table.txt"), dataFrame.toString(true, true));
	//
	//		Item item;
	//		
	//		item = new Item("regression_file", "regression.txt", "Regression file", TYPE.FILE);
	//		item.setGroup("Regression results");
	//		result.addOutputItem(item);
	//		
	//		item = new Item("regression_table", "regression_table.txt", "Regression table", TYPE.FILE);
	//		item.addTag("TABLE");
	//		item.setGroup("Regression results");
	//		result.addOutputItem(item);
	//
	//		item = new Item("regression_heatmap", "regression_heatmap.png", "The regression heatmap is: ", TYPE.IMAGE);
	//		item.setGroup("Regression results");
	//		result.addOutputItem(item);
	//	} catch (java.security.InvalidParameterException e) {
	//		logger.error("not valid parameter: execute regression");
	//		e.printStackTrace();
	//	} catch (MatrixIndexException e) {
	//		logger.error("MatrixIndexException: execute regression");
	//		// TODO Auto-generated catch block
	//		e.printStackTrace();
	//	} catch (MathException e) {
	//		// TODO Auto-generated catch block
	//		logger.error("math exception: execute regression");
	//		e.printStackTrace();
	//	} catch (IOException e) {
	//		// TODO Auto-generated catch block
	//		logger.error("IOException: execute regression");
	//		e.printStackTrace();
	//	} catch (InvalidColumnIndexException e) {
	//		// TODO Auto-generated catch block
	//		e.printStackTrace();
	//	}
	//}

}
