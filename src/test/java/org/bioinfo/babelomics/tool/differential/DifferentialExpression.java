package org.bioinfo.babelomics.tool.differential;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.ParseException;
import org.apache.commons.math.MathException;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.collections.matrix.DataFrame;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.data.dataset.FeatureData;
import org.bioinfo.graphics.canvas.Canvas;
import org.bioinfo.graphics.canvas.feature.ScoreFeature;
import org.bioinfo.graphics.canvas.panel.GridPanel;
import org.bioinfo.graphics.canvas.track.GridTrack;
import org.bioinfo.math.data.DoubleMatrix;
import org.bioinfo.math.exception.InvalidParameterException;
import org.bioinfo.math.result.TTestResult;
import org.bioinfo.math.result.TestResultList;
import org.bioinfo.math.stats.inference.MultipleTestCorrection;
import org.bioinfo.math.stats.inference.TTest;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.utils.StringUtils;

public class DifferentialExpression extends BabelomicsTool {


	public DifferentialExpression(String[] args) {
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

	@Override
	public void execute() {
		try {
			CommandLine cmd = parse(args);

			Dataset dataset = new Dataset(new File(cmd.getOptionValue("dataset")));
			String className = cmd.getOptionValue("class");
			String test = cmd.getOptionValue("test");

			String timeClass = cmd.getOptionValue("time-class", null);
			String censoredClass = cmd.getOptionValue("time-class", null);
				
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
				
				logger.info("2 classes");
				String label1 = dataset.getVariables().getByName(className).getLabels().get(0);
				String label2 = dataset.getVariables().getByName(className).getLabels().get(1);

				dataset.load();
				System.out.println("input matrix:\n" + dataset.getDoubleMatrix().toString());

				System.out.println("input feature names:\n" + StringUtils.arrayToString(dataset.getFeatureNames(), "\t"));

				int[] colOrder = new int[dataset.getColumnDimension()];
				int[] cols = dataset.getColumnIndexesByVariableValue(className, label1);
				DoubleMatrix sample1 = dataset.getSubMatrixByColumns(cols);

				cols = dataset.getColumnIndexesByVariableValue(className, label2);
				DoubleMatrix sample2 = dataset.getSubMatrixByColumns(cols);

				System.out.println("sample 1, matrix:\n" + sample1.toString());
				System.out.println("sample 2, matrix:\n" + sample2.toString());
				
				TTest tTest = new TTest();
				TestResultList<TTestResult> tTestResultList = tTest.tTest(sample1, sample2);
				System.out.println("yeeeeeee");
				MultipleTestCorrection.BHCorrection(tTestResultList);
				System.out.println("result\n" + tTestResultList.toString()+"\n");
				

				String dir = "/mnt/commons/test/biodata/example/out/";
				
				// creating output result (feature data format)
				//
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
				featureData.write(new File(dir + "stats.txt"));
				
				// creating output heatmap
				//
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
				canvas.save(dir + "heatmap");

				
				
				
				

				tTestResultList.setOrderCriteria(TTestResult.PVALUE_ORDER);
				TestResultList<TTestResult> ordered = tTestResultList.sort();
				System.out.println("order by p-value\n" + ordered.toString()+"\n");

				MultipleTestCorrection.FDRCorrection(ordered);
				System.out.println("order by fdr\n" + ordered.toString()+"\n");

				tTestResultList.setOrderCriteria(TTestResult.STATISTIC_ORDER);
				System.out.println("order by statistic\n" + ordered.sort().toString()+"\n");
								
				
			}else {
				logger.error("Number of labels distinct of 2");
			}
		} catch (InvalidParameterException e) {
			logger.error("Error opening the dataset", e.toString());
		} catch (MathException e) {
			logger.error("Error opening the dataset", e.toString());
		} catch (IOException e) {
			logger.error("Error opening the dataset", e.toString());
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
		logger.info("executing anova, not implemented yet");
	}

	private void executePearson(Dataset dataset, String className) {
		logger.info("executing pearson, not implemented yet");
	}

	private void executeSpearman(Dataset dataset, String className) {
		logger.info("executing spearman, not implemented yet");
	}

	private void executeRegression(Dataset dataset, String className) {
		logger.info("executing regression, not implemented yet");
	}

	private void executeCox(Dataset dataset, String className, String timeClass, String censoredClass) {
		logger.info("executing cox, not implemented yet");
	}

}
