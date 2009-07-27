package org.bioinfo.babelomics.tools.expression;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.security.InvalidParameterException;
import java.util.ArrayList;
import java.util.List;

import javax.imageio.ImageIO;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.ParseException;
import org.apache.commons.math.MathException;
import org.apache.commons.math.linear.MatrixIndexException;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.chart.XYPlotChart;
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
import org.bioinfo.math.stats.correlation.CorrelationTest;
import org.bioinfo.math.stats.inference.AnovaTest;
import org.bioinfo.math.stats.inference.MultipleTestCorrection;
import org.bioinfo.math.stats.inference.TTest;
import org.bioinfo.math.survival.CoxTest;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;
import org.bioinfo.utils.StringUtils;
import org.bioinfo.collections.matrix.DataFrame;




public class DifferentialAnalysis extends BabelomicsTool {


	public DifferentialAnalysis(String[] args) {
		super(args);
		initOptions();
	}

	@Override
	public void initOptions() {
		options.addOption(OptionFactory.createOption("dataset", "the data"));
		options.addOption(OptionFactory.createOption("test", "the test, possible values: t-test, bayes, sam, fold-change, anova, pearson, spearman, regression, cox"));
		options.addOption(OptionFactory.createOption("class", "class variable"));
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
		try {
			CommandLine cmd = parse(args);

			Dataset dataset = new Dataset(new File(cmd.getOptionValue("dataset")));
			String className = cmd.getOptionValue("class");
			String test = cmd.getOptionValue("test");

			String timeClass = cmd.getOptionValue("time-class", "");
			String censoredClass = cmd.getOptionValue("censor-class", "");
				
			if(cmd.hasOption("sample-filter") || cmd.hasOption("feature-filter")) {
				dataset = dataset.getSubDataset(cmd.getOptionValue("sample-filter"), "4", cmd.getOptionValue("feature-filter"), ""); 
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
				executePearson(dataset, className);
				return;
			}
			if(test.equals("spearman")) {
				executeSpearman(dataset, className);
				return;
			}
			if(test.equals("regression")) {
				executeRegression(dataset, className);
				return;
			}
			if(test.equals("cox")) {
				executeCox(dataset, className, timeClass, censoredClass);
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


	private void executeTTest(Dataset dataset, String className) {
		logger.info("executing t-test");
		try {

			if(dataset.getVariables().getByName(className).getLabels().size() == 2) {
				
				System.out.println("input dataset:\n" + dataset.toString());
				
				//loding data
				logger.info("2 classes t-test");
				jobStatus.addStatusMessage("10", "loading data");
				String label1 = dataset.getVariables().getByName(className).getLabels().get(0);
				String label2 = dataset.getVariables().getByName(className).getLabels().get(1);

				if ( dataset.getDoubleMatrix() == null ) { 
					dataset.load();
					dataset.validate();
				}
				System.out.println("input matrix:\n" + dataset.getDoubleMatrix().toString());

				System.out.println("input feature names:\n" + StringUtils.arrayToString(dataset.getFeatureNames(), "\t"));

				int[] colOrder = new int[dataset.getColumnDimension()];
				int[] cols = dataset.getColumnIndexesByVariableValue(className, label1);
				DoubleMatrix sample1 = dataset.getSubMatrixByColumns(cols);

				cols = dataset.getColumnIndexesByVariableValue(className, label2);
				DoubleMatrix sample2 = dataset.getSubMatrixByColumns(cols);

				System.out.println("sample 1, matrix:\n" + sample1.toString());
				System.out.println("sample 2, matrix:\n" + sample2.toString());
				
				//t-test
				TTest tTest = new TTest();
				TestResultList<TTestResult> tTestResultList = tTest.tTest(sample1, sample2);				
				MultipleTestCorrection.BHCorrection(tTestResultList);
				System.out.println("result\n" + tTestResultList.toString()+"\n");
				

				
				
				// creating output result (feature data format)
				//
				logger.info("creating output result");
				jobStatus.addStatusMessage("30", "creating output result");
				List<String> names = new ArrayList<String> (5);
				names.add("names");names.add("statistic");names.add("p-value");names.add("df");names.add("adj value");
				DataFrame dataFrame = new DataFrame(names, tTestResultList.size());
				for(int i=0 ; i<tTestResultList.size() ; i++) {
					dataFrame.setSingleValue(i, 0, dataset.getFeatureNames().get(i));
					dataFrame.setSingleValue(i, 1, ""+tTestResultList.get(i).getStatistic());
					dataFrame.setSingleValue(i, 2, ""+tTestResultList.get(i).getPValue());
					dataFrame.setSingleValue(i, 3, ""+tTestResultList.get(i).getDf());
					dataFrame.setSingleValue(i, 4, ""+tTestResultList.get(i).getAdjPValue());
				}
				FeatureData featureData = new FeatureData(dataFrame);
				System.out.println("result in feature format\n" + featureData.toString()+"\n");
				
				//saving data
				logger.info("saving output result");
				jobStatus.addStatusMessage("50", "saving output result");
				PrintWriter writer = new PrintWriter(outdir+"/stats.txt");
				writer.close();
				featureData.write(new File(outdir + "/stats.txt"));
				
				
				
				// creating output heatmap
				//
				logger.info("creating output heatmap");
				jobStatus.addStatusMessage("70", "creating output heatmap");
				int x = 2;				
				int y = 2;
				int rowDimension = dataset.getColumnDimension();
				int columnDimension = dataset.getRowDimension();
				int cellSide = 20;
				int rowLabelsWidth = 70;
				int colLabelsWidth = 70;
				int infoWidth = 0;
				double min = dataset.getDoubleMatrix().getMinValue();
				double max = dataset.getDoubleMatrix().getMaxValue();
				System.out.println("heatmap dimensions: (rowDimension, columnDimension) = (" + rowDimension + ", " + columnDimension + ")(min, max) = (" + min + ", " + max + ")");
				
				Canvas canvas = new Canvas("");
				canvas.setBorder(1);
				canvas.setBorderColor(Color.BLACK);
				canvas.setBgColor(Color.WHITE);
				
				GridPanel gridPanel = new GridPanel("", (rowDimension * cellSide) + rowLabelsWidth + infoWidth, (columnDimension * cellSide) + colLabelsWidth, x, y);
				GridTrack gridTrack = new GridTrack(rowDimension, columnDimension, cellSide, cellSide);
				gridTrack.setRowLabels(dataset.getFeatureNames());
				List<String> columnLabels = new ArrayList<String>(dataset.getSampleNames().size());
				for(int i=0 ; i <dataset.getSampleNames().size() ; i++) {
					columnLabels.add(dataset.getSampleNames().get(i) + " (" + dataset.getVariables().getByName(className).getValues().get(i) + ")");
				}
				gridTrack.setColumnLabels(columnLabels);
				//gridTrack.setName();
				gridTrack.setTopRegion(colLabelsWidth);
				gridTrack.setLeftRegion(rowLabelsWidth);
				ScoreFeature feature;
				for(int row=0 ; row<gridTrack.getColumnDimension() ; row++) {
					for(int column=0 ; column<gridTrack.getRowDimension() ; column++) {
//						System.out.print("row, column = " + row + ", " + column + ": value = "); System.out.println(dataset.getDoubleMatrix().get(row, column));
						feature = new ScoreFeature("name (" + column + ", " + row + ")", "", 0, 0, (dataset.getDoubleMatrix().get(row, column)-min)/(max-min));
						//feature.setJsFunction("http://www.cipf.es");
						gridTrack.setFeature(row, column, feature);
					}
				}
				gridPanel.add(gridTrack);
					
				canvas.addPanel(gridPanel);		
				canvas.render();
				
				
				//saving heatMap		
				logger.info("saving heatMap");
				jobStatus.addStatusMessage("70", "saving heatMap");
				canvas.save(outdir + "/heatmap");

				
				tTestResultList.setOrderCriteria(TTestResult.PVALUE_ORDER);
				TestResultList<TTestResult> ordered = tTestResultList.sort();
				System.out.println("order by p-value\n" + ordered.toString()+"\n");

				MultipleTestCorrection.FDRCorrection(ordered);
				System.out.println("order by fdr\n" + ordered.toString()+"\n");

				tTestResultList.setOrderCriteria(TTestResult.STATISTIC_ORDER);
				System.out.println("order by statistic\n" + ordered.sort().toString()+"\n");
				
				result.addOutputItem(new Item("t-test_file",outdir+"/stats.txt", "the t-test analisys file is: ", TYPE.FILE));				
				result.addOutputItem(new Item("heatMap_image",outdir+"heatmap.png", "The heat map image is: ", TYPE.IMAGE));				
				
			}else {
				logger.error("Number of labels distinct of 2");
			}
		} catch (InvalidParameterException e) {
			logger.error("Error opening the dataset", e.toString());
		} catch (MathException e) {
			logger.error("Error opening the dataset", e.toString());
		} catch (IOException e) {
			logger.error("Error opening the dataset", e.toString());
		} catch (org.bioinfo.math.exception.InvalidParameterException e) {
			// TODO Auto-generated catch block
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

	private void executeAnova(Dataset dataset, String className) {
		logger.info("executing anova");
		List<String> vars = dataset.getVariables().getByName(className).getValues();
		List<Double> doubleVars = new ArrayList<Double>(vars.size());
		for(String str: vars) {
			doubleVars.add(Double.parseDouble(str));
		}
		
		try {
			//dataset loadding
			if ( dataset.getDoubleMatrix() == null ) { 
				dataset.load();
				dataset.validate();
			}
			//job status+logger
			jobStatus.addStatusMessage("25", "reading data");
			logger.info("executing anova analisys...\n");
			
			//anova test
			AnovaTest anova = new AnovaTest(dataset.getDoubleMatrix(), vars);			
			TestResultList<AnovaTestResult> res = anova.compute();			
			
			//saving data			
			logger.info("saving anova analisys");			
			System.out.println("result = " + res.toString());
			jobStatus.addStatusMessage("80", "saving data");
			PrintWriter writer = new PrintWriter(outdir+"/anovaAnalisys.txt");
			writer.write(res.toString());
			writer.close();
			
			//generating heatmap
			//
			int xHeatMap = 2;				
			int yHeatMap = 2;
			int rowDimension = dataset.getColumnDimension();
			int columnDimension = dataset.getRowDimension();
			int cellSide = 20;
			int rowLabelsWidth = 70;
			int colLabelsWidth = 70;
			int infoWidth = 0;
			double min = dataset.getDoubleMatrix().getMinValue();
			double max = dataset.getDoubleMatrix().getMaxValue();
			System.out.println("heatmap dimensions: (rowDimension, columnDimension) = (" + rowDimension + ", " + columnDimension + ")(min, max) = (" + min + ", " + max + ")");
			
			Canvas canvas = new Canvas("");
			canvas.setBorder(1);
			canvas.setBorderColor(Color.BLACK);
			canvas.setBgColor(Color.WHITE);
			
			GridPanel gridPanel = new GridPanel("", (rowDimension * cellSide) + rowLabelsWidth + infoWidth, (columnDimension * cellSide) + colLabelsWidth, xHeatMap, yHeatMap);
			GridTrack gridTrack = new GridTrack(rowDimension, columnDimension, cellSide, cellSide);
			gridTrack.setRowLabels(dataset.getFeatureNames());
			List<String> columnLabels = new ArrayList<String>(dataset.getSampleNames().size());
			for(int i=0 ; i <dataset.getSampleNames().size() ; i++) {
				columnLabels.add(dataset.getSampleNames().get(i) + " (" + dataset.getVariables().getByName(className).getValues().get(i) + ")");
			}
			gridTrack.setColumnLabels(columnLabels);
			//gridTrack.setName();
			gridTrack.setTopRegion(colLabelsWidth);
			gridTrack.setLeftRegion(rowLabelsWidth);
			ScoreFeature feature;
			for(int row=0 ; row<gridTrack.getColumnDimension() ; row++) {
				for(int column=0 ; column<gridTrack.getRowDimension() ; column++) {
//					System.out.print("row, column = " + row + ", " + column + ": value = "); System.out.println(dataset.getDoubleMatrix().get(row, column));
					feature = new ScoreFeature("name (" + column + ", " + row + ")", "", 0, 0, (dataset.getDoubleMatrix().get(row, column)-min)/(max-min));
					//feature.setJsFunction("http://www.cipf.es");
					gridTrack.setFeature(row, column, feature);
				}
			}
			gridPanel.add(gridTrack);
			logger.info("saving heatmap");
			jobStatus.addStatusMessage("95", "saving heatmap");
			canvas.addPanel(gridPanel);		
			canvas.render();		
			canvas.save(outdir+"/heatmap");
			jobStatus.addStatusMessage("100", "done");
			
			
			//add outputResult file
			result.addOutputItem(new Item("anova_file",outdir+"/anovaAnalisys.txt", "The anova analisys file is: ", TYPE.FILE));			
			result.addOutputItem(new Item("heatMap_image","heatmap.png", "The pearson correlation image is: ", TYPE.IMAGE));
			jobStatus.addStatusMessage("100", "done");
			
			
			
		//	System.out.println("\n\n" + pearson.getMethod() + " results:\n" + pearson.compute().toString());
			
		} catch (java.security.InvalidParameterException e) {
			logger.error("not valid parametter: executePearson");
			e.printStackTrace();
		} catch (MatrixIndexException e) {
			logger.error("MatrixIndexException: executePearson");
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (MathException e) {
			// TODO Auto-generated catch block
			logger.error("math exception: executePearson");
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			logger.error("IOException: executePearson");
			e.printStackTrace();
		} 
}
		
	private void executePearson(Dataset dataset, String className) {
		logger.info("executing pearson");
		List<String> vars = dataset.getVariables().getByName(className).getValues();
		List<Double> doubleVars = new ArrayList<Double>(vars.size());
		for(String str: vars) {
			doubleVars.add(Double.parseDouble(str));
		}
		
		try {
			
			
			//dataset loadding
			if ( dataset.getDoubleMatrix() == null ) { 
				dataset.load();
				dataset.validate();
			}
			//job status+logger
			jobStatus.addStatusMessage("25", "reading data");
			logger.info("executing pearson correlation...\n");
			
			//pearson test
			CorrelationTest pearson = new CorrelationTest(dataset.getDoubleMatrix(), doubleVars, "pearson");
			TestResultList<CorrelationTestResult> res = pearson.compute();
			
			//saving data			
			logger.info("saving pearson analisys");			
			System.out.println("result = " + res.toString());			
			logger.info("saving pearson analisis\n");
			jobStatus.addStatusMessage("80", "saving data");
			PrintWriter writer = new PrintWriter(outdir+"/pearsonCorrelation.txt");
			writer.write(res.toString());
			writer.close();

			//generating boxplot
			XYPlotChart xyplot = new XYPlotChart();
			xyplot.setMainTitle("correlation Box plot");
			CommandLine cmd = parse(args);
			//XYSeries xyseries = new XYSeries(dataset.getDoubleMatrix());
			Dataset x = dataset.getSubDataset(cmd.getOptionValue("class"), "1", "", "");
			Dataset y = dataset.getSubDataset(cmd.getOptionValue("class"), "0", "", "");		
			
						
			double [] meanRowX = new double[dataset.getDoubleMatrix().getRowDimension()];
			List<Double> lx  = new ArrayList<Double>(dataset.getDoubleMatrix().getRowDimension());
			List<Double> ly  = new ArrayList<Double>(dataset.getDoubleMatrix().getRowDimension());
			
			for(int j=0; j<dataset.getDoubleMatrix().getRowDimension(); j++) {
				lx.add(x.getDoubleMatrix().getRowMean(j));				
				meanRowX[j] = x.getDoubleMatrix().getRowMean(j);
			}
			
			double [] meanRowY = new double[dataset.getDoubleMatrix().getRowDimension()];
			for(int j=0; j<dataset.getDoubleMatrix().getRowDimension(); j++) {
				meanRowY[j] = 0;
				ly.add(y.getDoubleMatrix().getRowMean(j));				
				meanRowY[j] = y.getDoubleMatrix().getRowMean(j);
			}
			xyplot.addSeries(meanRowX,meanRowY, "male-female:");
			xyplot.render();
			File out = new File(outdir+"/correlationimage.png" );
			logger.info("saving pearson correlation Box Plot");			
			jobStatus.addStatusMessage("90", "saving pearson Box PLot");			
			ImageIO.write(xyplot.getBufferedImage(), "png", out);
			//int[] order = ArrayUtils.order(lx);
//			List<Double> slx = ListUtils.ordered(lx, order);
//			List<Double> sly = ListUtils.ordered(ly, order);
			
//			xyplot.addSeries(ListUtils.toArray(slx),ListUtils.toArray(sly), "male-female");
		
			
			
//			for(int row=0 ; row < x.getDoubleMatrix().getRowDimension(); row++)
//			{		
//				rowListMeanX.add(x.getDoubleMatrix().getRowMean(row));
//				rowListMeanY.add(y.getDoubleMatrix().getRowMean(row));
//								
//			}
			//			
//			xyplot.addSeries(dataset.getSubDataset("only_males", "1", "", ""), dataset.getSubDataset("only_males", "1", "", ""), "") ;
			
			
			//generating heatmap
			//
			int xHeatMap = 2;				
			int yHeatMap = 2;
			int rowDimension = dataset.getColumnDimension();
			int columnDimension = dataset.getRowDimension();
			int cellSide = 20;
			int rowLabelsWidth = 70;
			int colLabelsWidth = 70;
			int infoWidth = 0;
			double min = dataset.getDoubleMatrix().getMinValue();
			double max = dataset.getDoubleMatrix().getMaxValue();
			System.out.println("heatmap dimensions: (rowDimension, columnDimension) = (" + rowDimension + ", " + columnDimension + ")(min, max) = (" + min + ", " + max + ")");
			
			Canvas canvas = new Canvas("");
			canvas.setBorder(1);
			canvas.setBorderColor(Color.BLACK);
			canvas.setBgColor(Color.WHITE);
			
			GridPanel gridPanel = new GridPanel("", (rowDimension * cellSide) + rowLabelsWidth + infoWidth, (columnDimension * cellSide) + colLabelsWidth, xHeatMap, yHeatMap);
			GridTrack gridTrack = new GridTrack(rowDimension, columnDimension, cellSide, cellSide);
			gridTrack.setRowLabels(dataset.getFeatureNames());
			List<String> columnLabels = new ArrayList<String>(dataset.getSampleNames().size());
			for(int i=0 ; i <dataset.getSampleNames().size() ; i++) {
				columnLabels.add(dataset.getSampleNames().get(i) + " (" + dataset.getVariables().getByName(className).getValues().get(i) + ")");
			}
			gridTrack.setColumnLabels(columnLabels);
			//gridTrack.setName();
			gridTrack.setTopRegion(colLabelsWidth);
			gridTrack.setLeftRegion(rowLabelsWidth);
			ScoreFeature feature;
			for(int row=0 ; row<gridTrack.getColumnDimension() ; row++) {
				for(int column=0 ; column<gridTrack.getRowDimension() ; column++) {
//					System.out.print("row, column = " + row + ", " + column + ": value = "); System.out.println(dataset.getDoubleMatrix().get(row, column));
					feature = new ScoreFeature("name (" + column + ", " + row + ")", "", 0, 0, (dataset.getDoubleMatrix().get(row, column)-min)/(max-min));
					//feature.setJsFunction("http://www.cipf.es");
					gridTrack.setFeature(row, column, feature);
				}
			}
			gridPanel.add(gridTrack);
			logger.info("saving heatmap");
			jobStatus.addStatusMessage("95", "saving heatmap");
			canvas.addPanel(gridPanel);		
			canvas.render();		
			canvas.save(outdir+"/heatmap");
			jobStatus.addStatusMessage("100", "done");
			result.addOutputItem(new Item("spearman_correlation_file",outdir+"/pearsonCorrelation.txt", "The spearman correlation file is: ", TYPE.FILE));
			result.addOutputItem(new Item("pearson_correlation_image","correlationimage.png", "The pearson correlation image is: ", TYPE.IMAGE));
			result.addOutputItem(new Item("pearson_heatMap_image","heatmap.png", "The pearson correlation image is: ", TYPE.IMAGE));
			jobStatus.addStatusMessage("100", "done");
			
			
			
		//	System.out.println("\n\n" + pearson.getMethod() + " results:\n" + pearson.compute().toString());
			
		} catch (java.security.InvalidParameterException e) {
			logger.error("not valid parametter: executePearson");
			e.printStackTrace();
		} catch (MatrixIndexException e) {
			logger.error("MatrixIndexException: executePearson");
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (MathException e) {
			// TODO Auto-generated catch block
			logger.error("math exception: executePearson");
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			logger.error("IOException: executePearson");
			e.printStackTrace();
		} catch (ParseException e) {
			logger.error("ParseException: executePearson");
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private void executeSpearman(Dataset dataset, String className) {
		logger.info("executing spearman");
		List<String> vars = dataset.getVariables().getByName(className).getValues();
		List<Double> doubleVars = new ArrayList<Double>(vars.size());
		for(String str: vars) {
			doubleVars.add(Double.parseDouble(str));
		}
		
		try {
			
			if ( dataset.getDoubleMatrix() == null ) { 
				dataset.load();
				dataset.validate();
			}
			
			jobStatus.addStatusMessage("25", "reading data");
			logger.info("executing pearson correlation...\n");

			CorrelationTest spearman = new CorrelationTest(dataset.getDoubleMatrix(), doubleVars, "spearman");
			TestResultList<CorrelationTestResult> res = spearman.compute();
			
			logger.info("saving..., res is null ? " + (res == null) + "\n");
			
			System.out.println("result = " + res.toString());

			
			logger.info("saving...\n");
			jobStatus.addStatusMessage("80", "saving data");
			PrintWriter writer = new PrintWriter(outdir+"/spearmanCorrelation.txt");
			writer.write(res.toString());
			writer.close();
			
			result.addOutputItem(new Item("spearman_correlation_file",outdir+"/spearmanCorrelation.txt", "The spearman correlation file is: ", TYPE.FILE));			
			jobStatus.addStatusMessage("100", "done");
				
			
			
			//generating boxplot
			XYPlotChart xyplot = new XYPlotChart();
			xyplot.setMainTitle("correlation Box plot");
			CommandLine cmd = parse(args);
			//XYSeries xyseries = new XYSeries(dataset.getDoubleMatrix());
			Dataset x = dataset.getSubDataset(cmd.getOptionValue("class"), "1", "", "");
			Dataset y = dataset.getSubDataset(cmd.getOptionValue("class"), "0", "", "");		
			
						
			double [] meanRowX = new double[dataset.getDoubleMatrix().getRowDimension()];
			List<Double> lx  = new ArrayList<Double>(dataset.getDoubleMatrix().getRowDimension());
			List<Double> ly  = new ArrayList<Double>(dataset.getDoubleMatrix().getRowDimension());
			
			for(int j=0; j<dataset.getDoubleMatrix().getRowDimension(); j++) {
				lx.add(x.getDoubleMatrix().getRowMean(j));				
				meanRowX[j] = x.getDoubleMatrix().getRowMean(j);
			}
			
			double [] meanRowY = new double[dataset.getDoubleMatrix().getRowDimension()];
			for(int j=0; j<dataset.getDoubleMatrix().getRowDimension(); j++) {
				meanRowY[j] = 0;
				ly.add(y.getDoubleMatrix().getRowMean(j));				
				meanRowY[j] = y.getDoubleMatrix().getRowMean(j);
			}
			xyplot.addSeries(meanRowX,meanRowY, "male-female:");
			xyplot.render();
			File out = new File(outdir+"/correlationimage.png" );
			logger.info("saving spearman correlation Box Plot");			
			jobStatus.addStatusMessage("90", "saving spearman Box PLot");			
			ImageIO.write(xyplot.getBufferedImage(), "png", out);
			//int[] order = ArrayUtils.order(lx);
//			List<Double> slx = ListUtils.ordered(lx, order);
//			List<Double> sly = ListUtils.ordered(ly, order);
			
//			xyplot.addSeries(ListUtils.toArray(slx),ListUtils.toArray(sly), "male-female");
		
			
			
//			for(int row=0 ; row < x.getDoubleMatrix().getRowDimension(); row++)
//			{		
//				rowListMeanX.add(x.getDoubleMatrix().getRowMean(row));
//				rowListMeanY.add(y.getDoubleMatrix().getRowMean(row));
//								
//			}
			//			
//			xyplot.addSeries(dataset.getSubDataset("only_males", "1", "", ""), dataset.getSubDataset("only_males", "1", "", ""), "") ;
			
			
			//generating heatmap
			//
			int xHeatMap = 2;				
			int yHeatMap = 2;
			int rowDimension = dataset.getColumnDimension();
			int columnDimension = dataset.getRowDimension();
			int cellSide = 20;
			int rowLabelsWidth = 70;
			int colLabelsWidth = 70;
			int infoWidth = 0;
			double min = dataset.getDoubleMatrix().getMinValue();
			double max = dataset.getDoubleMatrix().getMaxValue();
			System.out.println("heatmap dimensions: (rowDimension, columnDimension) = (" + rowDimension + ", " + columnDimension + ")(min, max) = (" + min + ", " + max + ")");
			
			Canvas canvas = new Canvas("");
			canvas.setBorder(1);
			canvas.setBorderColor(Color.BLACK);
			canvas.setBgColor(Color.WHITE);
			
			GridPanel gridPanel = new GridPanel("", (rowDimension * cellSide) + rowLabelsWidth + infoWidth, (columnDimension * cellSide) + colLabelsWidth, xHeatMap, yHeatMap);
			GridTrack gridTrack = new GridTrack(rowDimension, columnDimension, cellSide, cellSide);
			gridTrack.setRowLabels(dataset.getFeatureNames());
			List<String> columnLabels = new ArrayList<String>(dataset.getSampleNames().size());
			for(int i=0 ; i <dataset.getSampleNames().size() ; i++) {
				columnLabels.add(dataset.getSampleNames().get(i) + " (" + dataset.getVariables().getByName(className).getValues().get(i) + ")");
			}
			gridTrack.setColumnLabels(columnLabels);
			//gridTrack.setName();
			gridTrack.setTopRegion(colLabelsWidth);
			gridTrack.setLeftRegion(rowLabelsWidth);
			ScoreFeature feature;
			for(int row=0 ; row<gridTrack.getColumnDimension() ; row++) {
				for(int column=0 ; column<gridTrack.getRowDimension() ; column++) {
//					System.out.print("row, column = " + row + ", " + column + ": value = "); System.out.println(dataset.getDoubleMatrix().get(row, column));
					feature = new ScoreFeature("name (" + column + ", " + row + ")", "", 0, 0, (dataset.getDoubleMatrix().get(row, column)-min)/(max-min));
					//feature.setJsFunction("http://www.cipf.es");
					gridTrack.setFeature(row, column, feature);
				}
			}
			gridPanel.add(gridTrack);
			logger.info("saving heatmap");
			jobStatus.addStatusMessage("95", "saving heatmap");
			canvas.addPanel(gridPanel);		
			canvas.render();		
			canvas.save(outdir+"/heatmap");
			jobStatus.addStatusMessage("100", "done");
			result.addOutputItem(new Item("spearman_correlation_file",outdir+"/spearmanCorrelation.txt", "The spearman correlation file is: ", TYPE.FILE));
			result.addOutputItem(new Item("spearman_correlation_image","correlationimage.png", "The spearman correlation image is: ", TYPE.IMAGE));
			result.addOutputItem(new Item("spearman_heatMap_image","heatmap.png", "The spearman correlation image is: ", TYPE.IMAGE));
			jobStatus.addStatusMessage("100", "done");
			
			
		//	System.out.println("\n\n" + spearman.getMethod() + " results:\n" + spearman.compute().toString());
			
		} catch (java.security.InvalidParameterException e) {
			logger.error("not valid parametter: executePearson");
			e.printStackTrace();
		} catch (MatrixIndexException e) {
			logger.error("MatrixIndexException: executePearson");
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (MathException e) {
			// TODO Auto-generated catch block
			logger.error("math exception: executePearson");
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			logger.error("IOException: executePearson");
			e.printStackTrace();
		} catch (ParseException e) {
			logger.error("ParseException: executePearson");
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
			
			if ( dataset.getDoubleMatrix() == null ) { 
				dataset.load();
				dataset.validate();
			}
			
			jobStatus.addStatusMessage("25", "reading data");
			
						
			logger.info("executing regression analisis...\n");
			SimpleRegressionTest regression = new SimpleRegressionTest();
			TestResultList<SimpleRegressionTestResult> res = regression.compute(dataset.getDoubleMatrix(), doubleVars);
						
			logger.info("saving...\n");
			jobStatus.addStatusMessage("80", "saving data");
			PrintWriter writer = new PrintWriter(outdir+"/regression.txt");
			writer.write(res.toString());
			writer.close();
			
			
			
			result.addOutputItem(new Item("regressionAnalisis_file","regression.txt", "The regression file is: ", TYPE.FILE));
			jobStatus.addStatusMessage("100", "done");
			
			
						
		//	System.out.println("\n\n" + spearman.getMethod() + " results:\n" + spearman.compute().toString());
			
		} catch (java.security.InvalidParameterException e) {
			logger.error("not valid parametter: executePearson");
			e.printStackTrace();
		} catch (MatrixIndexException e) {
			logger.error("MatrixIndexException: executePearson");
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (MathException e) {
			// TODO Auto-generated catch block
			logger.error("math exception: executePearson");
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			logger.error("IOException: executePearson");
			e.printStackTrace();
		}
	}

	private void executeCox(Dataset dataset, String className, String timeClass, String censoredClass) {
		logger.info("executing Cox");
		List<String> vars = dataset.getVariables().getByName(className).getValues();
		List<Double> doubleVars = new ArrayList<Double>(vars.size());
		for(String str: vars) {
			doubleVars.add(Double.parseDouble(str));
		}
		
		try {
			
			if ( dataset.getDoubleMatrix() == null ) { 
				dataset.load();
				dataset.validate();
			}
			
			jobStatus.addStatusMessage("25", "reading data");
			
						
			logger.info("executing cox analisis...\n");
			
			List<String> strTimeVar = StringUtils.stringToList(",", "156,1040,59,421,329,769");
			List<String> strCensorVar = StringUtils.stringToList(",", "1,0,1,0,1,0");
			
			List<Double> timeVar = new ArrayList<Double>(strTimeVar.size());
			for(String str: strTimeVar) {
				timeVar.add(Double.parseDouble(str));
			}
			List<Double> censorVar = new ArrayList<Double>(strCensorVar.size());
			for(String str: strCensorVar) {
				censorVar.add(Double.parseDouble(str));
			}
			
			CoxTest coxtest = new CoxTest(dataset.getDoubleMatrix(), timeVar, censorVar);
			TestResultList<CoxTestResult> res ;
			jobStatus.addStatusMessage("50", "cox test running");
				res = coxtest.compute(dataset.getDoubleMatrix(), timeVar, censorVar);
				jobStatus.addStatusMessage("50", "cox test finished");
				System.out.println("Results cox test:\n" + coxtest.compute().toString());
				logger.info("saving...\n");
				jobStatus.addStatusMessage("80", "saving data");
				PrintWriter writer = new PrintWriter(outdir+"/cox.txt");
				writer.write(res.toString());
				writer.close();
			
			
			result.addOutputItem(new Item("coxAnalisis_file","cox.txt", "The cox file is: ", TYPE.FILE));
			jobStatus.addStatusMessage("100", "done");
			
			
						
		//	System.out.println("\n\n" + spearman.getMethod() + " results:\n" + spearman.compute().toString());
			
		
		} catch (java.security.InvalidParameterException e) {
			logger.error("not valid parametter: executePearson");
			e.printStackTrace();
		} catch (MatrixIndexException e) {
			logger.error("MatrixIndexException: executePearson");
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (MathException e) {
			// TODO Auto-generated catch block
			logger.error("math exception: executePearson");
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			logger.error("IOException: executePearson");
			e.printStackTrace();
		} 
	}

}



