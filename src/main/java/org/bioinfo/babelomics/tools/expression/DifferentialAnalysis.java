package org.bioinfo.babelomics.tools.expression;

import java.awt.Color;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.security.InvalidParameterException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math.MathException;
import org.apache.commons.math.linear.MatrixIndexException;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.collections.exceptions.InvalidColumnIndexException;
import org.bioinfo.collections.list.NamedArrayList;
import org.bioinfo.collections.matrix.DataFrame;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.data.dataset.FeatureData;
import org.bioinfo.graphics.canvas.Canvas;
import org.bioinfo.graphics.canvas.feature.ScoreFeature;
import org.bioinfo.graphics.canvas.panel.GridPanel;
import org.bioinfo.graphics.canvas.track.GridTrack;
import org.bioinfo.math.data.DoubleMatrix;
import org.bioinfo.math.regression.SimpleRegressionTest;
import org.bioinfo.math.result.AnovaTestResult;
import org.bioinfo.math.result.CorrelationTestResult;
import org.bioinfo.math.result.CoxTestResult;
import org.bioinfo.math.result.SimpleRegressionTestResult;
import org.bioinfo.math.result.TTestResult;
import org.bioinfo.math.result.TestResultList;
import org.bioinfo.math.stats.MultipleTestCorrection;
import org.bioinfo.math.stats.correlation.CorrelationTest;
import org.bioinfo.math.stats.inference.AnovaTest;
import org.bioinfo.math.stats.inference.TTest;
import org.bioinfo.math.stats.survival.CoxTest;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;

public class DifferentialAnalysis extends BabelomicsTool {


	public DifferentialAnalysis() {
		initOptions();
	}

	@Override
	public void initOptions() {
		options.addOption(OptionFactory.createOption("dataset", "the data"));
		options.addOption(OptionFactory.createOption("test", "the test, possible values: t-test, bayes, sam, fold-change, anova, pearson, spearman, regression, cox"));
		options.addOption(OptionFactory.createOption("class", "class variable", false));
		options.addOption(OptionFactory.createOption("time-class", "class variable", false));
		options.addOption(OptionFactory.createOption("censored-class", "class variable", false));
		options.addOption(OptionFactory.createOption("sample-filter", "class variable", false));
		options.addOption(OptionFactory.createOption("feature-filter", "class variable", false));
	}

	//	
	//	DEBUG_LEVEL = 1;
	//	INFO_LEVEL = 2;;
	//	WARNING_LEVEL = 3;
	//	ERROR_LEVEL = 4;
	//	FATAL_LEVEL = 5;
	@Override
	public void execute() {

		Dataset dataset = null;
		try {
			dataset = new Dataset(new File(commandLine.getOptionValue("dataset")));
		} catch (Exception e) {
			logger.error("Error opening the dataset", e.toString());
			return;
		}
		String test = commandLine.getOptionValue("test");

		String className = commandLine.getOptionValue("class", null);
		String timeClass = commandLine.getOptionValue("time-class", null);
		String censoredClass = commandLine.getOptionValue("censored-class", null);

		if(commandLine.hasOption("sample-filter") || commandLine.hasOption("feature-filter")) {
			try {
				dataset = dataset.getSubDataset(commandLine.getOptionValue("sample-filter"), "4", commandLine.getOptionValue("feature-filter"), "");
			} catch (Exception e) {
				logger.error("Error filtering the dataset", e.toString());
				return;
			}
		}


		System.out.println(dataset.toString()+"\n");
		if(test.equals("t-test")) {
			executeTTest(dataset, className);
			return;
		}
		if(test.equals("bayes")) {
			executeBayes(dataset, className);
			return;
		}
		if(test.equals("data-adaptive")) {
			executeDataAdaptive(dataset, className);
			return;
		}
		if(test.equals("sam")) {
			executeSam(dataset, className);
			return;
		}
		if(test.equals("anova")) {
			executeAnova(dataset, className);
			return;
		}
		if(test.equals("pearson")) {
			executeCorrelation(dataset, className, "pearson");
			return;
		}
		if(test.equals("spearman")) {
			executeCorrelation(dataset, className, "spearman");
			return;
		}
		if(test.equals("regression")) {
			executeRegression(dataset, className);
			return;
		}
		if(test.equals("cox")) {
			executeCox(dataset, timeClass, censoredClass);
			return;
		}

		logger.warn("que raroo....");
	}


	private void executeTTest(Dataset dataset, String className) {
		logger.info("executing t-test");
		
		try {
			// reading data
			//
			jobStatus.addStatusMessage("20", "reading data");
			logger.debug("reading data...\n");

			if(dataset.getVariables().getByName(className).getLabels().size() != 2) { 
				throw new InvalidParameterException("Number of labels distinct of 2");
			}

			String label1 = dataset.getVariables().getByName(className).getLabels().get(0);
			String label2 = dataset.getVariables().getByName(className).getLabels().get(1);

			if ( dataset.getDoubleMatrix() == null ) { 
				try {
					dataset.load();
				} catch (Exception e) {
					logger.error("Error loading the dataset", e.toString());
					return;
				}
				dataset.validate();
			}			
			
//			int[] colOrder = new int[dataset.getColumnDimension()];
			int[] cols = dataset.getColumnIndexesByVariableValue(className, label1);
			DoubleMatrix sample1 = dataset.getSubMatrixByColumns(cols);

			cols = dataset.getColumnIndexesByVariableValue(className, label2);
			DoubleMatrix sample2 = dataset.getSubMatrixByColumns(cols);

//			FileUtils.writeStringToFile(new File("/tmp/ttest/matrix1.txt"), sample1.toString());
//			FileUtils.writeStringToFile(new File("/tmp/ttest/matrix2.txt"), sample2.toString());
//			FileUtils.writeStringToFile(new File("/tmp/ttest/rownames.txt"), ListUtils.toString(dataset.getFeatureNames(), "\n"));

			// t-test
			//
			jobStatus.addStatusMessage("40", "computing t-test");
			logger.debug("computing t-test...\n");

			TTest tTest = new TTest();
			TestResultList<TTestResult> res = tTest.tTest(sample1, sample2);				
			MultipleTestCorrection.BHCorrection(res);

			int[] columnOrder = ListUtils.order(dataset.getVariables().getByName(className).getValues());
			int[] rowOrder = ListUtils.order(ListUtils.toList(res.getStatistics()), true);

			System.out.println("result size = " + res.size());
			
			// generating heatmap
			//
			jobStatus.addStatusMessage("60", "generating heatmap");
			logger.debug("generating heatmap...\n");

			Canvas heatmap = generateHeatmap(dataset, className, columnOrder, rowOrder, "statistic", res.getStatistics(), "adj. p-value", res.getAdjPValues());
			heatmap.save(getOutdir() + "/ttest_heatmap");
			
			// saving data
			//
			jobStatus.addStatusMessage("80", "saving results");
			logger.debug("saving results...");

			DataFrame dataFrame = new DataFrame(dataset.getFeatureNames().size(), 0);
			dataFrame.setRowNames(ListUtils.ordered(dataset.getFeatureNames(), rowOrder));

			//dataFrame.addColumn("id", ListUtils.ordered(dataset.getFeatureNames(), rowOrder));
			dataFrame.addColumn("statistic", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getStatistics()), rowOrder)));
			dataFrame.addColumn("p-value", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getPValues()), rowOrder)));
			dataFrame.addColumn("adj. p-value", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getAdjPValues()), rowOrder)));

			FeatureData featureData = new FeatureData(dataFrame);
			featureData.write(new File(getOutdir() + "/ttest.txt"));

			BufferedWriter bw = new BufferedWriter(new FileWriter(getOutdir() + "/ttest_table.txt"));
			bw.write(dataFrame.toString(true, true));
			bw.close();

			result.addOutputItem(new Item("ttest_file", "ttest.txt", "The t-test file is: ", TYPE.FILE));
			result.addOutputItem(new Item("ttest_heatmap", "ttest_heatmap.png", "The t-test heatmap is: ", TYPE.IMAGE));

			Item item = new Item("ttest_table", "ttest_table.txt", "The t-test table is: ", TYPE.FILE);
			item.addTag("TABLE");
			result.addOutputItem(item);
			
			// done
			//
			jobStatus.addStatusMessage("100", "done");
			logger.debug("t-test done\n");
		} catch (java.security.InvalidParameterException e) {
			logger.error("not valid parameter: execute t-test");
			e.printStackTrace();
		} catch (MatrixIndexException e) {
			logger.error("MatrixIndexException: execute t-test");
			e.printStackTrace();
		} catch (MathException e) {
			logger.error("math exception: execute t-test");
			e.printStackTrace();
		} catch (IOException e) {
			logger.error("IOException: execute t-test");
			e.printStackTrace();
		} catch (InvalidColumnIndexException e) {
			logger.error("IOException: execute t-test");
			e.printStackTrace();
		} catch (org.bioinfo.math.exception.InvalidParameterException e) {
			logger.error("IOException: execute t-test");
			e.printStackTrace();
		}
	}

	private void executeBayes(Dataset dataset, String className) {
		logger.info("executing bayes, not implemented yet");
	}

	private void executeDataAdaptive(Dataset dataset, String className) {
		logger.info("executing data adaptive, not implemented yet");
	}

	private void executeSam(Dataset dataset, String className) {
		logger.info("executing sam, not implemented yet");
	}

	// anova test
    //
	private void executeAnova(Dataset dataset, String className) {
		logger.info("executing anova test");
		List<String> vars = dataset.getVariables().getByName(className).getValues();

		try {
			// reading data
			//
			jobStatus.addStatusMessage("20", "reading data");
			logger.debug("reading data...\n");

			if ( dataset.getDoubleMatrix() == null ) { 
				try {
					dataset.load();
				} catch (Exception e) {
					logger.error("Error loading the dataset", e.toString());
					return;
				}
				dataset.validate();
			}

			// anova test
			//
			jobStatus.addStatusMessage("40", "computing anova test");
			logger.debug("computing anova test...\n");

			AnovaTest anova = new AnovaTest(dataset.getDoubleMatrix(), vars);			
			TestResultList<AnovaTestResult> res = anova.compute();			

			int[] columnOrder = ListUtils.order(vars);
			int[] rowOrder = ListUtils.order(ListUtils.toList(res.getStatistics()), true);

			
			// generating heatmap
			//
			jobStatus.addStatusMessage("60", "generating heatmap");
			logger.debug("generating heatmap...\n");

			Canvas heatmap = generateHeatmap(dataset, className, columnOrder, rowOrder, "statistic", res.getStatistics(), "adj. p-value", res.getAdjPValues());
			heatmap.save(getOutdir() + "/anova_heatmap");
			
			// saving data
			//
			jobStatus.addStatusMessage("80", "saving results");
			logger.debug("saving results...");

			DataFrame dataFrame = new DataFrame(dataset.getFeatureNames().size(), 0);
			dataFrame.setRowNames(ListUtils.ordered(dataset.getFeatureNames(), rowOrder));

			dataFrame.addColumn("statistic", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getStatistics()), rowOrder)));
			dataFrame.addColumn("p-value", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getPValues()), rowOrder)));
			dataFrame.addColumn("adj. p-value", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getAdjPValues()), rowOrder)));

			FeatureData featureData = new FeatureData(dataFrame);
			featureData.write(new File(getOutdir() + "/anova.txt"));

			BufferedWriter bw = new BufferedWriter(new FileWriter(getOutdir() + "/anova_table.txt"));
			bw.write(dataFrame.toString(true, true));
			bw.close();
			
			result.addOutputItem(new Item("anova_file", "anova.txt", "The anova file is: ", TYPE.FILE));
			result.addOutputItem(new Item("anova_heatmap", "anova_heatmap.png", "The anova heatmap is: ", TYPE.IMAGE));

			Item item = new Item("anova_table", "anova_table.txt", "The anova table is: ", TYPE.FILE);
			item.addTag("TABLE");
			result.addOutputItem(item);
			
			// done
			//
			jobStatus.addStatusMessage("100", "done");
			logger.debug("anova test done\n");
		} catch (java.security.InvalidParameterException e) {
			logger.error("not valid parameter: execute anova");
			e.printStackTrace();
		} catch (MatrixIndexException e) {
			logger.error("MatrixIndexException: execute anova");
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (MathException e) {
			// TODO Auto-generated catch block
			logger.error("math exception: execute anova");
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			logger.error("IOException: execute anova");
			e.printStackTrace();
		} catch (InvalidColumnIndexException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private void executeCorrelation(Dataset dataset, String className, String test) {
		logger.info("executing " + test);
		List<String> vars = dataset.getVariables().getByName(className).getValues();
		List<Double> doubleVars = new ArrayList<Double>(vars.size());
		for(String str: vars) {
			doubleVars.add(Double.parseDouble(str));
		}

		try {
			// reading data
			//
			jobStatus.addStatusMessage("20", "reading data");
			logger.debug("reading data...\n");

			if ( dataset.getDoubleMatrix() == null ) { 
				try {
					dataset.load();
				} catch (Exception e) {
					logger.error("Error loading the dataset", e.toString());
					return;
				}
				dataset.validate();
			}

			// spearman test
			//
			jobStatus.addStatusMessage("40", "computing " + test + " correlation");
			logger.debug("computing " + test + " correlation...\n");

			CorrelationTest corrTest = new CorrelationTest(dataset.getDoubleMatrix(), doubleVars, test);
			TestResultList<CorrelationTestResult> res = corrTest.compute();

			int[] columnOrder = ListUtils.order(vars);
			int[] rowOrder = ListUtils.order(ListUtils.toList(res.getCorrelations()), true);

			// generating heatmap
			//
			jobStatus.addStatusMessage("60", "generating heatmap");
			logger.debug("generating heatmap...\n");

			Canvas heatmap = generateHeatmap(dataset, className, columnOrder, rowOrder, "correlation", res.getCorrelations(), "adj. p-value", res.getAdjPValues());
			heatmap.save(getOutdir() + "/" + test + "_heatmap");
			
			// saving data
			//
			jobStatus.addStatusMessage("80", "saving " + test + " results");
			logger.debug("saving results...");

			DataFrame dataFrame = new DataFrame(dataset.getFeatureNames().size(), 0);
			dataFrame.setRowNames(ListUtils.ordered(dataset.getFeatureNames(), rowOrder));

			//dataFrame.addColumn("id", ListUtils.ordered(dataset.getFeatureNames(), rowOrder));
			dataFrame.addColumn("statistic", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getStatistics()), rowOrder)));
			dataFrame.addColumn("correlation", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getCorrelations()), rowOrder)));
			dataFrame.addColumn("p-value", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getPValues()), rowOrder)));
			dataFrame.addColumn("adj. p-value", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getAdjPValues()), rowOrder)));

			FeatureData featureData = new FeatureData(dataFrame);
			featureData.write(new File(getOutdir() + "/" + test + "_correlation.txt"));

			BufferedWriter bw = new BufferedWriter(new FileWriter(getOutdir() + "/" + test + "_table.txt"));
			bw.write(dataFrame.toString(true, true));
			bw.close();
			
			result.addOutputItem(new Item(test + "_correlation_file", test + "_correlation.txt", "The " + test + " correlation file is: ", TYPE.FILE));
			result.addOutputItem(new Item(test + "_correlation_heatmap", test + "_heatmap.png", "The " + test + " correlation heatmap is: ", TYPE.IMAGE));

			Item item = new Item(test + "_table", test + "_table.txt", "The " + test + " table is: ", TYPE.FILE);
			item.addTag("TABLE");
			result.addOutputItem(item);

			// done
			//
			jobStatus.addStatusMessage("100", "done");
			logger.debug(test + " correlation done\n");
		} catch (java.security.InvalidParameterException e) {
			logger.error("not valid parameter: execute " + test + " correlation");
			e.printStackTrace();
		} catch (MatrixIndexException e) {
			logger.error("MatrixIndexException:  execute " + test + " correlation");
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (MathException e) {
			// TODO Auto-generated catch block
			logger.error("math exception:  execute " + test + " correlation");
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			logger.error("IOException:  execute " + test + " correlation");
			e.printStackTrace();
		} catch (InvalidColumnIndexException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}


	private void executeRegression(Dataset dataset, String className) {
		logger.info("executing regression");
		List<String> vars = dataset.getVariables().getByName(className).getValues();
		List<Double> doubleVars = new ArrayList<Double>(vars.size());
		for(String str: vars) {
			doubleVars.add(Double.parseDouble(str));
		}

		try {
			// reading data
			//
			jobStatus.addStatusMessage("20", "reading data");
			logger.debug("reading data...\n");

			if ( dataset.getDoubleMatrix() == null ) { 
				try {
					dataset.load();
				} catch (Exception e) {
					logger.error("Error loading the dataset", e.toString());
					return;
				}
				dataset.validate();
			}

			// spearman test
			//
			jobStatus.addStatusMessage("40", "computing regression");
			logger.debug("computing regression...\n");

			SimpleRegressionTest regression = new SimpleRegressionTest();
			TestResultList<SimpleRegressionTestResult> res = regression.compute(dataset.getDoubleMatrix(), doubleVars);

			int[] columnOrder = ListUtils.order(vars);
			int[] rowOrder = ListUtils.order(ListUtils.toList(res.getStatistics()), true);

			
			
			// generating heatmap
			//
			jobStatus.addStatusMessage("60", "generating heatmap");
			logger.debug("generating heatmap...\n");

			Canvas heatmap = generateHeatmap(dataset, className, columnOrder, rowOrder, "slope", res.getSlopes(), "adj. p-value", res.getAdjPValues());
			heatmap.save(getOutdir() + "/regression_heatmap");
			
			// saving data
			//
			jobStatus.addStatusMessage("80", "saving results");
			logger.debug("saving results...");

			DataFrame dataFrame = new DataFrame(dataset.getFeatureNames().size(), 0);
			dataFrame.setRowNames(ListUtils.ordered(dataset.getFeatureNames(), rowOrder));

			//dataFrame.addColumn("id", ListUtils.ordered(dataset.getFeatureNames(), rowOrder));
			dataFrame.addColumn("statistic", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getStatistics()), rowOrder)));
			dataFrame.addColumn("slope", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getSlopes()), rowOrder)));
			dataFrame.addColumn("intercept", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getIntercepts()), rowOrder)));
			dataFrame.addColumn("p-value", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getPValues()), rowOrder)));
			dataFrame.addColumn("adj. p-value", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getAdjPValues()), rowOrder)));

			FeatureData featureData = new FeatureData(dataFrame);
			featureData.write(new File(getOutdir() + "/regression.txt"));

			BufferedWriter bw = new BufferedWriter(new FileWriter(getOutdir() + "/regression_table.txt"));
			bw.write(dataFrame.toString(true, true));
			bw.close();

			result.addOutputItem(new Item("regression_file", "regression.txt", "The regression file is: ", TYPE.FILE));
			result.addOutputItem(new Item("regression_heatmap", "regression_heatmap.png", "The regression heatmap is: ", TYPE.IMAGE));

			Item item = new Item("regression_table", "regression_table.txt", "The regression table is: ", TYPE.FILE);
			item.addTag("TABLE");
			result.addOutputItem(item);

			// done
			//
			jobStatus.addStatusMessage("100", "done");
			logger.debug("regression done\n");
		} catch (java.security.InvalidParameterException e) {
			logger.error("not valid parameter: execute regression");
			e.printStackTrace();
		} catch (MatrixIndexException e) {
			logger.error("MatrixIndexException: execute regression");
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (MathException e) {
			// TODO Auto-generated catch block
			logger.error("math exception: execute regression");
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			logger.error("IOException: execute regression");
			e.printStackTrace();
		} catch (InvalidColumnIndexException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private void executeCox(Dataset dataset, String timeClass, String censoredClass) {
		logger.info("executing cox");
		
		List<String> vars = dataset.getVariables().getByName(timeClass).getValues();
		List<Double> timeVars = new ArrayList<Double>(vars.size());
		for(String str: vars) {
			timeVars.add(Double.parseDouble(str));
		}
		vars = dataset.getVariables().getByName(censoredClass).getValues();
		List<Double> censoredVars = new ArrayList<Double>(vars.size());
		for(String str: vars) {
			censoredVars.add(Double.parseDouble(str));
		}

		try {
			// reading data
			//
			jobStatus.addStatusMessage("20", "reading data");
			logger.debug("reading data...\n");

			if ( dataset.getDoubleMatrix() == null ) { 
				try {
					dataset.load();
				} catch (Exception e) {
					logger.error("Error loading the dataset", e.toString());
					return;
				}
				dataset.validate();
			}

			// cox test
			//
			jobStatus.addStatusMessage("40", "computing cox");
			logger.debug("computing cox...\n");

			CoxTest coxtest = new CoxTest();
			//System.out.println("input matrix = \n" + dataset.getDoubleMatrix().toString());
			TestResultList<CoxTestResult> res = coxtest.compute(dataset.getDoubleMatrix(), timeVars, censoredVars);

			int[] columnOrder = ListUtils.order(vars);
			int[] rowOrder = ListUtils.order(ListUtils.toList(res.getStatistics()), true);
			
			// generating heatmap
			//
			jobStatus.addStatusMessage("60", "generating heatmap");
			logger.debug("generating heatmap...\n");

			Canvas heatmap = generateHeatmap(dataset, timeClass, columnOrder, rowOrder, "coeff.", res.getCoefs(), "adj. p-value", res.getAdjPValues());
			heatmap.save(getOutdir() + "/cox_heatmap");
			
			// saving data
			//
			jobStatus.addStatusMessage("80", "saving results");
			logger.debug("saving results...");

			DataFrame dataFrame = new DataFrame(dataset.getFeatureNames().size(), 0);
			dataFrame.setRowNames(ListUtils.ordered(dataset.getFeatureNames(), rowOrder));

			//dataFrame.addColumn("id", ListUtils.ordered(dataset.getFeatureNames(), rowOrder));
			dataFrame.addColumn("statistic", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getStatistics()), rowOrder)));
			dataFrame.addColumn("coeff.", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getCoefs()), rowOrder)));
			dataFrame.addColumn("p-value", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getPValues()), rowOrder)));
			dataFrame.addColumn("adj. p-value", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getAdjPValues()), rowOrder)));

			FeatureData featureData = new FeatureData(dataFrame);
			featureData.write(new File(getOutdir() + "/cox.txt"));

			BufferedWriter bw = new BufferedWriter(new FileWriter(getOutdir() + "/cox_table.txt"));
			bw.write(dataFrame.toString(true, true));
			bw.close();
			
			result.addOutputItem(new Item("cox_file", "cox.txt", "The cox file is: ", TYPE.FILE));
			result.addOutputItem(new Item("cox_heatmap", "cox_heatmap.png", "The cox heatmap is: ", TYPE.IMAGE));
			
			Item item = new Item("cox_table", "cox_table.txt", "The cox table is: ", TYPE.FILE);
			item.addTag("TABLE");
			result.addOutputItem(item);

			// done
			//
			jobStatus.addStatusMessage("100", "done");
			logger.debug("regression done\n");
		} catch (java.security.InvalidParameterException e) {
			logger.error("not valid parameter: execute regression");
			e.printStackTrace();
		} catch (MatrixIndexException e) {
			logger.error("MatrixIndexException: execute regression");
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (MathException e) {
			// TODO Auto-generated catch block
			logger.error("math exception: execute regression");
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			logger.error("IOException: execute regression");
			e.printStackTrace();
		} catch (InvalidColumnIndexException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}


	public Canvas generateHeatmap(Dataset dataset, String className, int[] columnOrder, int[] rowOrder, String infoName1, double[] infoList1, String infoName2, double[] infoList2) {
		int xHeatMap = 2;				
		int yHeatMap = 2;
		int rowDimension = dataset.getRowDimension();
		int columnDimension = dataset.getColumnDimension();
		int cellSide = 16;
		int rowLabelsWidth = 70;
		int colLabelsWidth = 70;
		int infoWidth = 140;
		//double min = dataset.getDoubleMatrix().getMinValue();
		//double max = dataset.getDoubleMatrix().getMaxValue();
		double offset, min, max, mean, deviation, standard;
		//System.out.println("heatmap dimensions: (rowDimension, columnDimension) = (" + rowDimension + ", " + columnDimension + ")(min, max) = (" + min + ", " + max + ")");

		Canvas canvas = new Canvas("");
		canvas.setBorderWidth(0);
		canvas.setBorderColor(Color.BLACK);
		canvas.setBackGroundColor(Color.WHITE);

		GridPanel gridPanel = new GridPanel("", (columnDimension * cellSide) + rowLabelsWidth + infoWidth + canvas.getBorderPadding() + canvas.getSpaceSeparator(), ((rowDimension+1) * cellSide) + colLabelsWidth + canvas.getBorderPadding() + canvas.getSpaceSeparator(), xHeatMap, yHeatMap);
		GridTrack gridTrack = new GridTrack(rowDimension, columnDimension, cellSide, cellSide);
		gridTrack.setRowLabels(ListUtils.ordered(dataset.getFeatureNames(), rowOrder));
		List<String> columnLabels = new ArrayList<String>(dataset.getSampleNames().size());
		for(int i=0 ; i <columnOrder.length ; i++) {
			columnLabels.add(dataset.getSampleNames().get(columnOrder[i]) + " (" + dataset.getVariables().getByName(className).getValues().get(columnOrder[i]) + ")");
		}
		gridTrack.setColumnLabels(columnLabels);
		//gridTrack.setName();
		gridTrack.setTopRegion(colLabelsWidth);
		gridTrack.setLeftRegion(rowLabelsWidth);
		gridTrack.setRightRegion(infoWidth);

		double[] values = new double[columnOrder.length]; 
		int row, column;
		
		ScoreFeature feature;
		for(int i=0 ; i<rowOrder.length ; i++) {
			row = rowOrder[i];
			mean = dataset.getDoubleMatrix().getRowMean(row);
			deviation = dataset.getDoubleMatrix().getRowStdDeviation(row);
			min = Double.MAX_VALUE;
			max = Double.MIN_NORMAL;
			for(int j=0 ; j<columnOrder.length ; j++) {
				column = columnOrder[j];
				values[column] = (deviation == 0) ? Double.NaN : (dataset.getDoubleMatrix().get(row, column)-mean)/(deviation);
				if ( min > values[column] ) min = values[column];
				if ( max < values[column] ) max = values[column];
			}
			offset = ( min <= 0 ) ? Math.abs(min) : (-1 * min);
			for(int j=0 ; j<columnOrder.length ; j++) {
				column = columnOrder[j];
				//System.out.print("row, column = " + row + ", " + column + ": value = "); System.out.println(dataset.getDoubleMatrix().get(row, column));
//				feature = new ScoreFeature("name (" + column + ", " + row + ")", "", 0, 0, (dataset.getDoubleMatrix().get(row, column)-min)/(max-min));
				//standard = (deviation == 0) ? Double.NaN : (dataset.getDoubleMatrix().get(row, column)-mean)/(deviation);
				standard = (values[column] + offset) / ( max + offset);
//				System.out.println("(value, standard) = (" + dataset.getDoubleMatrix().get(row, column) + ", " + standard + ")");
				feature = new ScoreFeature("name, " + dataset.getDoubleMatrix().get(row, column), "", 0, 0, standard);
				//feature = new ScoreFeature("name (" + row + ", " + column + ")", "", 0, 0, dataset.getDoubleMatrix().get(row, column));
				//feature.setJsFunction("http://www.cipf.es");
				//				gridTrack.setFeature(row, column, feature);
				gridTrack.setFeature(i, j, feature);
			}
		}

		DecimalFormat df = new DecimalFormat("##0.000000000");
		List<String> aux = new ArrayList<String>(infoList1.length);
		for(int i=0 ; i<infoList1.length ; i++) {
			aux.add(df.format(infoList1[i]));
		}
		gridTrack.addInfo(new NamedArrayList(infoName1, ListUtils.ordered(aux, rowOrder)));
		aux.clear();
		for(int i=0 ; i<infoList2.length ; i++) {
			aux.add(df.format(infoList2[i]));
		}
		gridTrack.addInfo(new NamedArrayList(infoName2, ListUtils.ordered(aux, rowOrder)));

		gridPanel.add(gridTrack);
		canvas.addPanel(gridPanel);		
		canvas.render();
		
		return canvas;
	}

}































//			//generating boxplot
//			XYPlotChart xyplot = new XYPlotChart();
//			xyplot.setMainTitle("correlation Box plot");
//			CommandLine commandLine = parse(args);
//			//XYSeries xyseries = new XYSeries(dataset.getDoubleMatrix());
//			Dataset x = dataset.getSubDataset(commandLine.getOptionValue("class"), "1", "", "");
//			Dataset y = dataset.getSubDataset(commandLine.getOptionValue("class"), "0", "", "");		
//			
//						
//			double [] meanRowX = new double[dataset.getDoubleMatrix().getRowDimension()];
//			List<Double> lx  = new ArrayList<Double>(dataset.getDoubleMatrix().getRowDimension());
//			List<Double> ly  = new ArrayList<Double>(dataset.getDoubleMatrix().getRowDimension());
//			
//			for(int j=0; j<dataset.getDoubleMatrix().getRowDimension(); j++) {
//				lx.add(x.getDoubleMatrix().getRowMean(j));				
//				meanRowX[j] = x.getDoubleMatrix().getRowMean(j);
//			}
//			
//			double [] meanRowY = new double[dataset.getDoubleMatrix().getRowDimension()];
//			for(int j=0; j<dataset.getDoubleMatrix().getRowDimension(); j++) {
//				meanRowY[j] = 0;
//				ly.add(y.getDoubleMatrix().getRowMean(j));				
//				meanRowY[j] = y.getDoubleMatrix().getRowMean(j);
//			}
//			xyplot.addSeries(meanRowX,meanRowY, "male-female:");
//			xyplot.render();
//			File out = new File(outdir+"/correlationimage.png" );
//			logger.debug("saving pearson correlation Box Plot");			
//			jobStatus.addStatusMessage("90", "saving pearson Box PLot");			
//			ImageIO.write(xyplot.getBufferedImage(), "png", out);
//			//int[] order = ArrayUtils.order(lx);
////			List<Double> slx = ListUtils.ordered(lx, order);
////			List<Double> sly = ListUtils.ordered(ly, order);
//			
////			xyplot.addSeries(ListUtils.toArray(slx),ListUtils.toArray(sly), "male-female");
//		
//			
//			
////			for(int row=0 ; row < x.getDoubleMatrix().getRowDimension(); row++)
////			{		
////				rowListMeanX.add(x.getDoubleMatrix().getRowMean(row));
////				rowListMeanY.add(y.getDoubleMatrix().getRowMean(row));
////								
////			}
//			//			
////			xyplot.addSeries(dataset.getSubDataset("only_males", "1", "", ""), dataset.getSubDataset("only_males", "1", "", ""), "") ;




//private void executePearson(Dataset dataset, String className) {
//	logger.info("executing pearson");
//	List<String> vars = dataset.getVariables().getByName(className).getValues();
//	List<Double> doubleVars = new ArrayList<Double>(vars.size());
//	for(String str: vars) {
//		doubleVars.add(Double.parseDouble(str));
//	}
//
//	try {
//
//		// reading data
//		//
//		jobStatus.addStatusMessage("20", "reading data");
//		logger.debug("reading data...\n");
//
//		if ( dataset.getDoubleMatrix() == null ) { 
//			try {
//				dataset.load();
//			} catch (Exception e) {
//				logger.error("Error loading the dataset", e.toString());
//				return;
//			}
//			dataset.validate();
//		}
//
//		// pearson test
//		//
//		jobStatus.addStatusMessage("40", "computing pearson correlation");
//		logger.debug("computing pearson correlation...\n");
//
//		CorrelationTest pearson = new CorrelationTest(dataset.getDoubleMatrix(), doubleVars, "pearson");
//		TestResultList<CorrelationTestResult> res = pearson.compute();
//
//		int[] columnOrder = ListUtils.order(vars);
//		int[] rowOrder = ListUtils.order(ListUtils.toList(res.getCorrelations()), true);
//
//		//			System.out.println("corr = " + ListUtils.toString(corr));
//		//			System.out.println("sorted corr = " + ListUtils.toString(ListUtils.ordered(corr, rowOrder)));
//
//		// generating heatmap
//		//
//		jobStatus.addStatusMessage("60", "generating heatmap");
//		logger.debug("generating heatmap...\n");
//		int xHeatMap = 2;				
//		int yHeatMap = 2;
//		int rowDimension = dataset.getRowDimension();
//		int columnDimension = dataset.getColumnDimension();
//		int cellSide = 16;
//		int rowLabelsWidth = 70;
//		int colLabelsWidth = 70;
//		int infoWidth = 140;
//		double min = dataset.getDoubleMatrix().getMinValue();
//		double max = dataset.getDoubleMatrix().getMaxValue();
//		double mean = dataset.getDoubleMatrix().getMeanValue();
//		double deviation = dataset.getDoubleMatrix().getDeviationValue();
//		double standard;
//		System.out.println("heatmap dimensions: (rowDimension, columnDimension) = (" + rowDimension + ", " + columnDimension + ")(min, max) = (" + min + ", " + max + "), (mean, deviation) = (" + mean + ", " + deviation + ")");
//
//		Canvas canvas = new Canvas("");
//		canvas.setBorderWidth(0);
//		canvas.setBorderColor(Color.BLACK);
//		canvas.setBackGroundColor(Color.WHITE);
//
//		GridPanel gridPanel = new GridPanel("", (columnDimension * cellSide) + rowLabelsWidth + infoWidth + canvas.getBorderPadding() + canvas.getSpaceSeparator(), ((rowDimension+1) * cellSide) + colLabelsWidth + canvas.getBorderPadding() + canvas.getSpaceSeparator(), xHeatMap, yHeatMap);
//		GridTrack gridTrack = new GridTrack(rowDimension, columnDimension, cellSide, cellSide);
//		gridTrack.setRowLabels(ListUtils.ordered(dataset.getFeatureNames(), rowOrder));
//		List<String> columnLabels = new ArrayList<String>(dataset.getSampleNames().size());
//		for(int i=0 ; i <columnOrder.length ; i++) {
//			columnLabels.add(dataset.getSampleNames().get(columnOrder[i]) + " (" + dataset.getVariables().getByName(className).getValues().get(columnOrder[i]) + ")");
//		}
//		gridTrack.setColumnLabels(columnLabels);
//		//gridTrack.setName();
//		gridTrack.setTopRegion(colLabelsWidth);
//		gridTrack.setLeftRegion(rowLabelsWidth);
//		gridTrack.setRightRegion(infoWidth);
//
//		ScoreFeature feature;
//		for(int i=0 ; i<rowOrder.length ; i++) {
//			int row = rowOrder[i];
//			for(int j=0 ; j<columnOrder.length ; j++) {
//				int column = columnOrder[j];
//				//System.out.print("row, column = " + row + ", " + column + ": value = "); System.out.println(dataset.getDoubleMatrix().get(row, column));
//				//feature = new ScoreFeature("name (" + column + ", " + row + ")", "", 0, 0, (dataset.getDoubleMatrix().get(row, column)-min)/(max-min));
//				standard = (dataset.getDoubleMatrix().get(row, column)-mean)/(deviation);
//				System.out.println("(value, standard) = (" + dataset.getDoubleMatrix().get(row, column) + ", " + standard + ")");
//				feature = new ScoreFeature("name (" + column + ", " + row + ")", "", 0, 0, standard);
//				//feature.setJsFunction("http://www.cipf.es");
//				//					gridTrack.setFeature(row, column, feature);
//				gridTrack.setFeature(i, j, feature);
//			}
//		}
//
//		DecimalFormat df = new DecimalFormat("##0.000000000");
//		List<String> aux = new ArrayList<String>(res.getCorrelations().length);
//		for(int i=0 ; i<res.getCorrelations().length ; i++) {
//			aux.add(df.format(res.getCorrelations()[i]));
//		}
//		gridTrack.addInfo(new NamedArrayList("correlation", ListUtils.ordered(aux, rowOrder)));
//		aux.clear();
//		for(int i=0 ; i<res.getAdjPValues().length ; i++) {
//			aux.add(df.format(res.getAdjPValues()[i]));
//		}
//		gridTrack.addInfo(new NamedArrayList("adj. p-value", ListUtils.ordered(aux, rowOrder)));
//
//		gridPanel.add(gridTrack);
//		canvas.addPanel(gridPanel);		
//		canvas.render();		
//		canvas.save(getOutdir() + "/pearson_heatmap");
//
//		// saving data
//		//
//		jobStatus.addStatusMessage("80", "saving results");
//		logger.debug("saving results...");
//
//		//			PrintWriter writer = new PrintWriter(getOutdir() + "/pearson_correlation.txt");
//		//			writer.write(res.toString());
//		//			writer.close();
//
//		DataFrame dataFrame = new DataFrame(dataset.getFeatureNames().size(), 0);
//		dataFrame.setRowNames(ListUtils.ordered(dataset.getFeatureNames(), rowOrder));
//
//		//			System.out.println("columns ==> " + dataFrame.getColumnNames().size());
//		//			System.out.println("rows ==> "+dataFrame.getRowNames().size());
//		//			System.out.println("data size ==> "+dataFrame.getData().size());
//
//		dataFrame.addColumn("id", ListUtils.ordered(dataset.getFeatureNames(), rowOrder));
//		dataFrame.addColumn("correlation", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getCorrelations()), rowOrder)));
//		dataFrame.addColumn("statistic", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getStatistics()), rowOrder)));
//		dataFrame.addColumn("p-value", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getPValues()), rowOrder)));
//		dataFrame.addColumn("adj. p-value", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getAdjPValues()), rowOrder)));
//
//		//			System.out.println("dataframe : \n" + dataFrame.toString());
//
//
//		FeatureData featureData = new FeatureData(dataFrame);
//		//			System.out.println("--- id column = " + featureData.getDataFrame().getColumn("id").toString());
//		//			System.out.println("--- correlation column = " + featureData.getDataFrame().getColumn("correlation").toString());
//		//			System.out.println("--- statistic column = " + featureData.getDataFrame().getColumn("statistic").toString());
//		//
//		//			
//		featureData.write(new File(getOutdir() + "/pearson_correlation.txt"));
//
//		//			System.out.println("featuredata : \n" + featureData.toString());
//
//		result.addOutputItem(new Item("pearson_correlation_file", getOutdir() + "/pearson_correlation.txt", "The pearson correlation file is: ", TYPE.FILE));
//		result.addOutputItem(new Item("pearson_correlation_heatmap", getOutdir() + "/pearson_heatmap.png", "The pearson correlation heatmap is: ", TYPE.IMAGE));
//		//			result.addOutputItem(new Item("pearson_correlation_image","heatmap.png", "The pearson correlation image is: ", TYPE.IMAGE));
//
//		//			//generating boxplot
//		//			XYPlotChart xyplot = new XYPlotChart();
//		//			xyplot.setMainTitle("correlation Box plot");
//		//			CommandLine commandLine = parse(args);
//		//			//XYSeries xyseries = new XYSeries(dataset.getDoubleMatrix());
//		//			Dataset x = dataset.getSubDataset(commandLine.getOptionValue("class"), "1", "", "");
//		//			Dataset y = dataset.getSubDataset(commandLine.getOptionValue("class"), "0", "", "");		
//		//			
//		//						
//		//			double [] meanRowX = new double[dataset.getDoubleMatrix().getRowDimension()];
//		//			List<Double> lx  = new ArrayList<Double>(dataset.getDoubleMatrix().getRowDimension());
//		//			List<Double> ly  = new ArrayList<Double>(dataset.getDoubleMatrix().getRowDimension());
//		//			
//		//			for(int j=0; j<dataset.getDoubleMatrix().getRowDimension(); j++) {
//		//				lx.add(x.getDoubleMatrix().getRowMean(j));				
//		//				meanRowX[j] = x.getDoubleMatrix().getRowMean(j);
//		//			}
//		//			
//		//			double [] meanRowY = new double[dataset.getDoubleMatrix().getRowDimension()];
//		//			for(int j=0; j<dataset.getDoubleMatrix().getRowDimension(); j++) {
//		//				meanRowY[j] = 0;
//		//				ly.add(y.getDoubleMatrix().getRowMean(j));				
//		//				meanRowY[j] = y.getDoubleMatrix().getRowMean(j);
//		//			}
//		//			xyplot.addSeries(meanRowX,meanRowY, "male-female:");
//		//			xyplot.render();
//		//			File out = new File(outdir+"/correlationimage.png" );
//		//			logger.debug("saving pearson correlation Box Plot");			
//		//			jobStatus.addStatusMessage("90", "saving pearson Box PLot");			
//		//			ImageIO.write(xyplot.getBufferedImage(), "png", out);
//		//			//int[] order = ArrayUtils.order(lx);
//		////			List<Double> slx = ListUtils.ordered(lx, order);
//		////			List<Double> sly = ListUtils.ordered(ly, order);
//		//			
//		////			xyplot.addSeries(ListUtils.toArray(slx),ListUtils.toArray(sly), "male-female");
//		//		
//		//			
//		//			
//		////			for(int row=0 ; row < x.getDoubleMatrix().getRowDimension(); row++)
//		////			{		
//		////				rowListMeanX.add(x.getDoubleMatrix().getRowMean(row));
//		////				rowListMeanY.add(y.getDoubleMatrix().getRowMean(row));
//		////								
//		////			}
//		//			//			
//		////			xyplot.addSeries(dataset.getSubDataset("only_males", "1", "", ""), dataset.getSubDataset("only_males", "1", "", ""), "") ;
//
//		// done
//		//
//		jobStatus.addStatusMessage("100", "done");
//		logger.debug("pearson correlation done\n");
//		//			result.addOutputItem(new Item("spearman_correlation_file",outdir+"/pearsonCorrelation.txt", "The spearman correlation file is: ", TYPE.FILE));
//		//			result.addOutputItem(new Item("pearson_correlation_image","correlationimage.png", "The pearson correlation image is: ", TYPE.IMAGE));
//		//			result.addOutputItem(new Item("pearson_correlation_image","heatmap.png", "The pearson correlation image is: ", TYPE.IMAGE));
//		//			jobStatus.addStatusMessage("100", "done");
//
//
//
//		//	System.out.println("\n\n" + pearson.getMethod() + " results:\n" + pearson.compute().toString());
//
//	} catch (java.security.InvalidParameterException e) {
//		logger.error("not valid parameter: executePearson");
//		e.printStackTrace();
//	} catch (MatrixIndexException e) {
//		logger.error("MatrixIndexException: executePearson");
//		// TODO Auto-generated catch block
//		e.printStackTrace();
//	} catch (MathException e) {
//		// TODO Auto-generated catch block
//		logger.error("math exception: executePearson");
//		e.printStackTrace();
//	} catch (IOException e) {
//		// TODO Auto-generated catch block
//		logger.error("IOException: executePearson");
//		e.printStackTrace();
//	} catch (InvalidColumnIndexException e) {
//		// TODO Auto-generated catch block
//		e.printStackTrace();
//	}
//}
//
//private void executeSpearman(Dataset dataset, String className) {
//	logger.info("executing spearman");
//	List<String> vars = dataset.getVariables().getByName(className).getValues();
//	List<Double> doubleVars = new ArrayList<Double>(vars.size());
//	for(String str: vars) {
//		doubleVars.add(Double.parseDouble(str));
//	}
//
//	try {
//
//		// reading data
//		//
//		jobStatus.addStatusMessage("20", "reading data");
//		logger.debug("reading data...\n");
//
//		if ( dataset.getDoubleMatrix() == null ) { 
//			try {
//				dataset.load();
//			} catch (Exception e) {
//				logger.error("Error loading the dataset", e.toString());
//				return;
//			}
//			dataset.validate();
//		}
//
//		// spearman test
//		//
//		jobStatus.addStatusMessage("40", "computing spearman correlation");
//		logger.debug("computing spearman correlation...\n");
//
//		CorrelationTest spearman = new CorrelationTest(dataset.getDoubleMatrix(), doubleVars, "spearman");
//		TestResultList<CorrelationTestResult> res = spearman.compute();
//
//		int[] columnOrder = ListUtils.order(vars);
//		int[] rowOrder = ListUtils.order(ListUtils.toList(res.getCorrelations()), true);
//
//		//			System.out.println("corr = " + ListUtils.toString(corr));
//		//			System.out.println("sorted corr = " + ListUtils.toString(ListUtils.ordered(corr, rowOrder)));
//
//		// generating heatmap
//		//
//		jobStatus.addStatusMessage("60", "generating heatmap");
//		logger.debug("generating heatmap...\n");
//
//		Canvas heatmap = generateHeatmap(dataset, className, columnOrder, rowOrder, "correlation", res.getCorrelations(), "adj. p-value", res.getAdjPValues());
//		heatmap.save(getOutdir() + "/spearman_heatmap");
//		
//		
//		// saving data
//		//
//		jobStatus.addStatusMessage("80", "saving results");
//		logger.debug("saving results...");
//
//		DataFrame dataFrame = new DataFrame(dataset.getFeatureNames().size(), 0);
//		dataFrame.setRowNames(ListUtils.ordered(dataset.getFeatureNames(), rowOrder));
//
//		dataFrame.addColumn("id", ListUtils.ordered(dataset.getFeatureNames(), rowOrder));
//		dataFrame.addColumn("correlation", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getCorrelations()), rowOrder)));
//		dataFrame.addColumn("statistic", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getStatistics()), rowOrder)));
//		dataFrame.addColumn("p-value", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getPValues()), rowOrder)));
//		dataFrame.addColumn("adj. p-value", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getAdjPValues()), rowOrder)));
//
//		FeatureData featureData = new FeatureData(dataFrame);
//		featureData.write(new File(getOutdir() + "/spearman_correlation.txt"));
//
//		result.addOutputItem(new Item("spearman_correlation_file", getOutdir() + "/spearman_correlation.txt", "The spearman correlation file is: ", TYPE.FILE));
//		result.addOutputItem(new Item("spearman_correlation_heatmap", getOutdir() + "/spearman_heatmap.png", "The spearman correlation heatmap is: ", TYPE.IMAGE));
//
//
//		// done
//		//
//		jobStatus.addStatusMessage("100", "done");
//		logger.debug("spearman correlation done\n");
//		//			result.addOutputItem(new Item("spearman_correlation_file",outdir+"/pearsonCorrelation.txt", "The spearman correlation file is: ", TYPE.FILE));
//		//			result.addOutputItem(new Item("pearson_correlation_image","correlationimage.png", "The pearson correlation image is: ", TYPE.IMAGE));
//		//			result.addOutputItem(new Item("pearson_correlation_image","heatmap.png", "The pearson correlation image is: ", TYPE.IMAGE));
//		//			jobStatus.addStatusMessage("100", "done");
//
//
//
//		//	System.out.println("\n\n" + pearson.getMethod() + " results:\n" + pearson.compute().toString());
//
//	} catch (java.security.InvalidParameterException e) {
//		logger.error("not valid parameter: executeSpearman");
//		e.printStackTrace();
//	} catch (MatrixIndexException e) {
//		logger.error("MatrixIndexException: executeSpearman");
//		// TODO Auto-generated catch block
//		e.printStackTrace();
//	} catch (MathException e) {
//		// TODO Auto-generated catch block
//		logger.error("math exception: executeSpearman");
//		e.printStackTrace();
//	} catch (IOException e) {
//		// TODO Auto-generated catch block
//		logger.error("IOException: executeSpearman");
//		e.printStackTrace();
//	} catch (InvalidColumnIndexException e) {
//		// TODO Auto-generated catch block
//		e.printStackTrace();
//	}
//}
