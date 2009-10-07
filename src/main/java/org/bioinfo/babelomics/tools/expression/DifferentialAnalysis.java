package org.bioinfo.babelomics.tools.expression;

import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.data.dataset.Dataset;

public class DifferentialAnalysis extends BabelomicsTool {

	public Dataset dataset = null;
	
	public DifferentialAnalysis() {
		initOptions();
	}

	@Override
	public void initOptions() {
		// TODO Auto-generated method stub
		
	}

	@Override
	protected void execute() {
		// TODO Auto-generated method stub
		
	}


	
	
	
	
	
	
	
	
	
//	//	
//	//	DEBUG_LEVEL = 1;
//	//	INFO_LEVEL = 2;;
//	//	WARNING_LEVEL = 3;
//	//	ERROR_LEVEL = 4;
//	//	FATAL_LEVEL = 5;
//	@Override
//	public void execute() {
//
//		try {
//			dataset = new Dataset(new File(commandLine.getOptionValue("dataset")));
//		} catch (Exception e) {
//			logger.error("Error opening the dataset", e.toString());
//			return;
//		}
//		String test = commandLine.getOptionValue("test");
//
//		String className = commandLine.getOptionValue("class", null);
//		String timeClass = commandLine.getOptionValue("time-class", null);
//		String censoredClass = commandLine.getOptionValue("censored-class", null);
//
//		if(commandLine.hasOption("sample-filter") || commandLine.hasOption("feature-filter")) {
//			try {
//				dataset = dataset.getSubDataset(commandLine.getOptionValue("sample-filter"), "4", commandLine.getOptionValue("feature-filter"), "");
//			} catch (Exception e) {
//				logger.error("Error filtering the dataset", e.toString());
//				return;
//			}
//		}
//
//
//		System.out.println(dataset.toString()+"\n");
//		if(test.equals("t")) {
//			executeTTest(dataset, className);
//			return;
//		}
//		if(test.equals("bayes")) {
//			executeBayes(dataset, className);
//			return;
//		}
//		if(test.equals("data-adaptive")) {
//			executeDataAdaptive(dataset, className);
//			return;
//		}
//		if(test.equals("sam")) {
//			executeSam(dataset, className);
//			return;
//		}
//		if(test.equals("anova")) {
//			executeAnova(dataset, className);
//			return;
//		}
//		if(test.equals("pearson")) {
//			executeCorrelation(dataset, className, "pearson");
//			return;
//		}
//		if(test.equals("spearman")) {
//			executeCorrelation(dataset, className, "spearman");
//			return;
//		}
//		if(test.equals("regression")) {
//			executeRegression(dataset, className);
//			return;
//		}
//		if(test.equals("cox")) {
//			executeCox(dataset, timeClass, censoredClass);
//			return;
//		}
//		if (test.equals("masigpro")) {
//			String continClass = commandLine.getOptionValue("contin-class", null);
//			String seriesClass = commandLine.getOptionValue("series-class", null);
//
//			executeMaSigPro(dataset, continClass, seriesClass);
//			return;
//		}
//
//		logger.warn("que raroo....");
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
//	
//	
//	
//	
//	
//	
//	
//	
//	
//	
//	
//	
//	
//	
//	
//	
//	
//	private void executeTTest(Dataset dataset, String className) {
//		logger.info("executing t-test");
//		
//		try {
//			// reading data
//			//
//			jobStatus.addStatusMessage("20", "reading data");
//			logger.debug("reading data...\n");
//
//			if(dataset.getVariables().getByName(className).getLabels().size() != 2) { 
//				throw new InvalidParameterException("Number of labels distinct of 2");
//			}
//
//			String label1 = dataset.getVariables().getByName(className).getLabels().get(0);
//			String label2 = dataset.getVariables().getByName(className).getLabels().get(1);
//
//			if ( dataset.getDoubleMatrix() == null ) { 
//				try {
//					dataset.load();
//				} catch (Exception e) {
//					logger.error("Error loading the dataset", e.toString());
//					return;
//				}
//				dataset.validate();
//			}			
//			
////			int[] colOrder = new int[dataset.getColumnDimension()];
//			int[] cols = dataset.getColumnIndexesByVariableValue(className, label1);
//			DoubleMatrix sample1 = dataset.getSubMatrixByColumns(cols);
//
//			cols = dataset.getColumnIndexesByVariableValue(className, label2);
//			DoubleMatrix sample2 = dataset.getSubMatrixByColumns(cols);
//
////			FileUtils.writeStringToFile(new File("/tmp/ttest/matrix1.txt"), sample1.toString());
////			FileUtils.writeStringToFile(new File("/tmp/ttest/matrix2.txt"), sample2.toString());
////			FileUtils.writeStringToFile(new File("/tmp/ttest/rownames.txt"), ListUtils.toString(dataset.getFeatureNames(), "\n"));
//
//			// t-test
//			//
//			jobStatus.addStatusMessage("40", "computing t-test");
//			logger.debug("computing t-test...\n");
//
//			TTest tTest = new TTest();
//			TestResultList<TTestResult> res = tTest.tTest(sample1, sample2);				
//			MultipleTestCorrection.BHCorrection(res);
//
//			int[] columnOrder = ListUtils.order(dataset.getVariables().getByName(className).getValues());
//			int[] rowOrder = ListUtils.order(ListUtils.toList(res.getStatistics()), true);
//
//			System.out.println("result size = " + res.size());
//			
//			// generating heatmap
//			//
//			jobStatus.addStatusMessage("60", "generating heatmap");
//			logger.debug("generating heatmap...\n");
//
//			Canvas heatmap = generateHeatmap(dataset, className, columnOrder, rowOrder, "statistic", res.getStatistics(), "adj. p-value", res.getAdjPValues());
//			heatmap.save(getOutdir() + "/t_heatmap");
//			
//			// saving data
//			//
//			jobStatus.addStatusMessage("80", "saving results");
//			logger.debug("saving results...");
//
//			DataFrame dataFrame = new DataFrame(dataset.getFeatureNames().size(), 0);
//			dataFrame.setRowNames(ListUtils.ordered(dataset.getFeatureNames(), rowOrder));
//
//			//dataFrame.addColumn("id", ListUtils.ordered(dataset.getFeatureNames(), rowOrder));
//			dataFrame.addColumn("statistic", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getStatistics()), rowOrder)));
//			dataFrame.addColumn("p-value", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getPValues()), rowOrder)));
//			dataFrame.addColumn("adj. p-value", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getAdjPValues()), rowOrder)));
//
//			FeatureData featureData = new FeatureData(dataFrame);
//			featureData.write(new File(getOutdir() + "/t.txt"));
//
//			BufferedWriter bw = new BufferedWriter(new FileWriter(getOutdir() + "/t_table.txt"));
//			bw.write(dataFrame.toString(true, true));
//			bw.close();
//
//
//			result.addOutputItem(new Item("t_file", "t.txt", "T-test file", TYPE.FILE));
//			
//			Item item = new Item("t_table", "t_table.txt", "T-test table", TYPE.FILE);
//			item.addTag("TABLE");
//			result.addOutputItem(item);
//
//			result.addOutputItem(new Item("t_heatmap", "t_heatmap.png", "T-test heatmap", TYPE.IMAGE));
//
//			
//			// done
//			//
//			jobStatus.addStatusMessage("100", "done");
//			logger.debug("t-test done\n");
//		} catch (java.security.InvalidParameterException e) {
//			logger.error("not valid parameter: execute t-test");
//			e.printStackTrace();
//		} catch (MatrixIndexException e) {
//			logger.error("MatrixIndexException: execute t-test");
//			e.printStackTrace();
//		} catch (MathException e) {
//			logger.error("math exception: execute t-test");
//			e.printStackTrace();
//		} catch (IOException e) {
//			logger.error("IOException: execute t-test");
//			e.printStackTrace();
//		} catch (InvalidColumnIndexException e) {
//			logger.error("IOException: execute t-test");
//			e.printStackTrace();
//		} catch (org.bioinfo.math.exception.InvalidParameterException e) {
//			logger.error("IOException: execute t-test");
//			e.printStackTrace();
//		}
//	}
//
//	private void executeBayes(Dataset dataset, String className) {
//		logger.info("executing bayes, not implemented yet");
//	}
//
//	private void executeDataAdaptive(Dataset dataset, String className) {
//		logger.info("executing data adaptive, not implemented yet");
//	}
//
//	private void executeSam(Dataset dataset, String className) {
//		logger.info("executing sam, not implemented yet");
//	}
//
//	// anova test
//    //
//	private void executeAnova(Dataset dataset, String className) {
//		logger.info("executing anova test");
//		List<String> vars = dataset.getVariables().getByName(className).getValues();
//
//		try {
//			// reading data
//			//
//			jobStatus.addStatusMessage("20", "reading data");
//			logger.debug("reading data...\n");
//
//			if ( dataset.getDoubleMatrix() == null ) { 
//				try {
//					dataset.load();
//				} catch (Exception e) {
//					logger.error("Error loading the dataset", e.toString());
//					return;
//				}
//				dataset.validate();
//			}
//
//			// anova test
//			//
//			jobStatus.addStatusMessage("40", "computing anova test");
//			logger.debug("computing anova test...\n");
//
//			AnovaTest anova = new AnovaTest(dataset.getDoubleMatrix(), vars);			
//			TestResultList<AnovaTestResult> res = anova.compute();			
//
//			int[] columnOrder = ListUtils.order(vars);
//			int[] rowOrder = ListUtils.order(ListUtils.toList(res.getStatistics()), true);
//
//			
//			// generating heatmap
//			//
//			jobStatus.addStatusMessage("60", "generating heatmap");
//			logger.debug("generating heatmap...\n");
//
//			Canvas heatmap = generateHeatmap(dataset, className, columnOrder, rowOrder, "statistic", res.getStatistics(), "adj. p-value", res.getAdjPValues());
//			heatmap.save(getOutdir() + "/anova_heatmap");
//			
//			// saving data
//			//
//			jobStatus.addStatusMessage("80", "saving results");
//			logger.debug("saving results...");
//
//			DataFrame dataFrame = new DataFrame(dataset.getFeatureNames().size(), 0);
//			dataFrame.setRowNames(ListUtils.ordered(dataset.getFeatureNames(), rowOrder));
//
//			dataFrame.addColumn("statistic", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getStatistics()), rowOrder)));
//			dataFrame.addColumn("p-value", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getPValues()), rowOrder)));
//			dataFrame.addColumn("adj. p-value", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getAdjPValues()), rowOrder)));
//
//			FeatureData featureData = new FeatureData(dataFrame);
//			featureData.write(new File(getOutdir() + "/anova.txt"));
//
//			BufferedWriter bw = new BufferedWriter(new FileWriter(getOutdir() + "/anova_table.txt"));
//			bw.write(dataFrame.toString(true, true));
//			bw.close();
//			
//			result.addOutputItem(new Item("anova_file", "anova.txt", "The anova file is: ", TYPE.FILE));
//			result.addOutputItem(new Item("anova_heatmap", "anova_heatmap.png", "The anova heatmap is: ", TYPE.IMAGE));
//
//			Item item = new Item("anova_table", "anova_table.txt", "The anova table is: ", TYPE.FILE);
//			item.addTag("TABLE");
//			result.addOutputItem(item);
//			
//			// done
//			//
//			jobStatus.addStatusMessage("100", "done");
//			logger.debug("anova test done\n");
//		} catch (java.security.InvalidParameterException e) {
//			logger.error("not valid parameter: execute anova");
//			e.printStackTrace();
//		} catch (MatrixIndexException e) {
//			logger.error("MatrixIndexException: execute anova");
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		} catch (MathException e) {
//			// TODO Auto-generated catch block
//			logger.error("math exception: execute anova");
//			e.printStackTrace();
//		} catch (IOException e) {
//			// TODO Auto-generated catch block
//			logger.error("IOException: execute anova");
//			e.printStackTrace();
//		} catch (InvalidColumnIndexException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
//	}
//
//
//
//
//
//	
//	
//	private void executeMaSigPro(Dataset dataset, String continClass, String seriesClass) {
//		logger.info("executing maSigPro");
//		
//		System.out.println("(continClass, seriesClass) = (" + continClass + ", " + seriesClass + ")");
//		
//		List<String> continVars = dataset.getVariables().getByName(continClass).getValues();
//		List<String> seriesVars = dataset.getVariables().getByName(seriesClass).getValues();
//
//		try {
//			// reading data
//			//
//			jobStatus.addStatusMessage("20", "reading data");
//			logger.debug("reading data...\n");
//
//			if ( dataset.getDoubleMatrix() == null ) { 
//				try {
//					dataset.load();
//				} catch (Exception e) {
//					logger.error("Error loading the dataset", e.toString());
//					return;
//				}
//				dataset.validate();
//			}
//
//			// masigpro test
//			//
//			jobStatus.addStatusMessage("40", "computing maSigPro");
//			logger.debug("computing maSigPro...\n");
//
//			String line, testPath = outdir + "/test/";
//			List<String> names = dataset.getFeatureNames();
//			TextFileWriter writer = new TextFileWriter(testPath + "input.txt");
//			writer.writeLine("#CONTIN\t" + ListUtils.toString(continVars, "\t"));
//			writer.writeLine("#SERIES\t" + ListUtils.toString(seriesVars, "\t"));
//			writer.writeLine("#NAMES\t" + ListUtils.toString(dataset.getSampleNames(), "\t"));
//			for(int i=0 ; i<names.size() ; i++) {
//				line = names.get(i) + "\t";
//				line = line + ListUtils.toString(ListUtils.toList(dataset.getDoubleMatrix().getRow(i)), "\t");
//				writer.writeLine(line);
//			}
//			writer.close();			
//		} catch (java.security.InvalidParameterException e) {
//			logger.error("not valid parameter: execute regression");
//			e.printStackTrace();
//		} catch (MatrixIndexException e) {
//			logger.error("MatrixIndexException: execute regression");
//			// TODO Auto-generated catch block
//			e.printStackTrace();
////		} catch (MathException e) {
////			// TODO Auto-generated catch block
////			logger.error("math exception: execute regression");
////			e.printStackTrace();
//		} catch (IOException e) {
//			// TODO Auto-generated catch block
//			logger.error("IOException: execute regression");
//			e.printStackTrace();
////		} catch (InvalidColumnIndexException e) {
////			// TODO Auto-generated catch block
////			e.printStackTrace();
//		}
//	}
//	
//	
//	

}




