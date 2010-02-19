package org.bioinfo.babelomics.tools.expression.differential;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.bioinfo.babelomics.methods.expression.differential.Limma;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ArrayUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.data.dataset.FeatureData;
import org.bioinfo.data.list.DataFrame;
import org.bioinfo.graphics.canvas.Canvas;
import org.bioinfo.math.data.DoubleMatrix;
import org.bioinfo.math.result.AnovaTestResult;
import org.bioinfo.math.result.LimmaTestResult;
import org.bioinfo.math.result.TTestResult;
import org.bioinfo.math.result.TestResultList;
import org.bioinfo.math.stats.inference.AnovaTest;
import org.bioinfo.math.stats.inference.FoldChangeTest;
import org.bioinfo.math.stats.inference.TTest;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;

public class ClassComparison extends BabelomicsTool {

	private Dataset dataset;
	private String test;
	private String className;
	private List<String> classValues;
	private String correction;

	public ClassComparison() {
		initOptions();
	}

	@Override
	public void initOptions() {
		options.addOption(OptionFactory.createOption("dataset", "the data"));
		options.addOption(OptionFactory.createOption("class-name", "class"));
		options.addOption(OptionFactory.createOption("class-values", "value"));
		options.addOption(OptionFactory.createOption("test", "test"));
		options.addOption(OptionFactory.createOption("correction", "Multiple-test correction: fdr, bh, by, bonferroni, hochberg, hold"));
		//options.addOption(OptionFactory.createOption("batch", "class class variable"));
	}

	@Override
	public void execute() {

		// reading dataset
		//
		dataset = initDataset(new File(commandLine.getOptionValue("dataset")));

		test = commandLine.getOptionValue("test", null);
		className = commandLine.getOptionValue("class-name", null);
		classValues = (commandLine.hasOption("class-values") ? StringUtils.toList(commandLine.getOptionValue("class-values", null), ",") : null);
		correction = commandLine.getOptionValue("correction", "fdr");

		if ( classValues == null ) {
			classValues = ListUtils.unique(dataset.getVariables().getByName(className).getValues());
		} else {
			classValues = ListUtils.unique(classValues);
			if ( classValues == null ) {
				abort("classvaluesmissing_execute_classcomparison", "class values missing", "class values missing", "class values missing");				
			}
		}

		// sanity check
		//
		if ( test == null ) {
			abort("testmissing_execute_classcomparison", "class comparison test missing", "class comparison test missing", "class comparison test missing");
		}

		// executing test
		//
		updateJobStatus("40", "computing " + test);
		if ( "limma".equalsIgnoreCase(test) ) {
			executeLimma();
		} else if ( "t".equalsIgnoreCase(test) ) {
			if ( classValues.size() == 2 ) {
				executeT();
			} else {
				abort("testmismatched_execute_classcomparison", "test " + test + " not supported for " + classValues.size() + "-class test", "test " + test + " not supported for " + classValues.size() + "-class test", "test " + test + " not supported for " + classValues.size() + "-class test");								
			}
		} else if ( "fold-change".equalsIgnoreCase(test) ) {
			if ( classValues.size() == 2 ) {
				executeFoldChange();
			} else {
				abort("testmismatched_execute_classcomparison", "test " + test + " not supported for " + classValues.size() + "-class test", "test " + test + " not supported for " + classValues.size() + "-class test", "test " + test + " not supported for " + classValues.size() + "-class test");								
			}
		} else if ( "anova".equalsIgnoreCase(test) ) {
			if ( classValues.size() > 2 ) {
				executeAnova();
			} else {
				abort("testmismatched_execute_classcomparison", "test " + test + " not supported for " + classValues.size() + "-class test", "test " + test + " not supported for " + classValues.size() + "-class test", "test " + test + " not supported for " + classValues.size() + "-class test");												
			}
		} else {
			abort("testunknown_execute_classcomparison", "unknown test " + test, "unknown test " + test, "unknown test " + test);															
		}
	}	

	/**
	 * 
	 */
	public void executeT() {

		int[] cols = dataset.getColumnIndexesByVariableValue(className, classValues.get(0));
		DoubleMatrix sample1 = dataset.getSubMatrixByColumns(cols);

		cols = dataset.getColumnIndexesByVariableValue(className, classValues.get(1));
		DoubleMatrix sample2 = dataset.getSubMatrixByColumns(cols);

		try {
			TTest tTest = new TTest();
			TestResultList<TTestResult> res = tTest.tTest(sample1, sample2);
			
			// apply multiple test correction according to input correction
			DiffExpressionUtils.multipleTestCorrection(res, correction);
			
			// generating heatmap
			//
			updateJobStatus("60", "generating heatmap");
			int[] columnOrder = ListUtils.order(dataset.getVariables().getByName(className).getValues());
			int[] rowOrder = ListUtils.order(ArrayUtils.toList(res.getStatistics()), true);
			Canvas heatmap = DiffExpressionUtils.generateHeatmap(dataset, className, columnOrder, rowOrder, "statistic", res.getStatistics(), "adj. p-value", res.getAdjPValues());
			String heatmapFilename = getOutdir() + "/" + test + "_heatmap.png";
			try {
				heatmap.save(heatmapFilename);
			} catch (IOException e) {
				printError("ioexception_executet_classcomparison", "error generating heatmap", e.toString(), e);
			}
			
			updateJobStatus("80", "saving results");
			DataFrame dataFrame = new DataFrame(dataset.getFeatureNames().size(), 0);

			dataFrame.addColumn("statistic", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getStatistics()), rowOrder)));
			dataFrame.addColumn("p-value", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getPValues()), rowOrder)));
			dataFrame.addColumn("adj. p-value", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getAdjPValues()), rowOrder)));

			dataFrame.setRowNames(ListUtils.ordered(dataset.getFeatureNames(), rowOrder));

			FeatureData featureData = new FeatureData(dataFrame);

			File file = new File(outdir + "/t.txt");
			featureData.save(file);
			if ( file.exists() ) {
				result.addOutputItem(new Item("tfile", file.getName(), "T-test output file", TYPE.FILE, new ArrayList<String>(2), new HashMap<String, String>(2), "T-test output files"));							
			}

			file = new File(outdir + "/t_table.txt");
			IOUtils.write(file, dataFrame.toString(true, true));
			if ( file.exists() ) {
				result.addOutputItem(new Item("ttable", file.getName(), "T-test output table", TYPE.FILE, StringUtils.toList("TABLE", ","), new HashMap<String, String>(2), "T-test output files"));											
			}
			
			if ( new File(heatmapFilename).exists() ) {
				result.addOutputItem(new Item(test + "_heatmap", test + "_heatmap.png", test.toUpperCase() + " heatmap", TYPE.IMAGE, new ArrayList<String>(2), new HashMap<String, String>(2), "Heatmap image"));
			}
			
			DiffExpressionUtils.addOutputLists(dataFrame, test, "statistic", result, outdir);
			
		} catch (Exception e) {
			e.printStackTrace();
			abort("exception_executet_classcomparison", "error running t-test", "error running t-test: " + e.getMessage(), "error running t-test: " + e.getMessage());
		}		
	}

	/**
	 * 
	 */
	public void executeFoldChange() {

		int[] cols = dataset.getColumnIndexesByVariableValue(className, classValues.get(0));
		DoubleMatrix sample1 = dataset.getSubMatrixByColumns(cols);

		cols = dataset.getColumnIndexesByVariableValue(className, classValues.get(1));
		DoubleMatrix sample2 = dataset.getSubMatrixByColumns(cols);

		try {
			FoldChangeTest foldChange = new FoldChangeTest();
			double[] logRes  = foldChange.logFoldChange(sample1, sample2);
			double[] diffRes = foldChange.diffFoldChange(sample1, sample2);

			updateJobStatus("80", "saving results");

			DataFrame dataFrame = new DataFrame(dataset.getFeatureNames().size(), 0);

			dataFrame.addColumn("log", ArrayUtils.toStringList(logRes));
			dataFrame.addColumn("diff", ArrayUtils.toStringList(diffRes));

			dataFrame.setRowNames(dataset.getFeatureNames());

			FeatureData featureData = new FeatureData(dataFrame);

			File file = new File(outdir + "/foldchange.txt");
			featureData.save(file);
			if ( file.exists() ) {
				result.addOutputItem(new Item("foldchangefile", file.getName(), "Fold-change output file", TYPE.FILE, new ArrayList<String>(2), new HashMap<String, String>(2), "Fold-change output files"));							
			}

			file = new File(outdir + "/foldchange_table.txt");
			IOUtils.write(file, dataFrame.toString(true, true));
			if ( file.exists() ) {
				result.addOutputItem(new Item("foldchangetable", file.getName(), "Fold-change output table", TYPE.FILE, StringUtils.toList("TABLE", ","), new HashMap<String, String>(2), "Fold-change output files"));											
			}						
			
			//DiffExpressionUtils.addOutputLists(dataFrame, test, "statistic", result, outdir);
			
		} catch (Exception e) {
			e.printStackTrace();
			abort("exception_executefoldchange_classcomparison", "error running fold-change", "error running fold-change: " + e.getMessage(), "error running fold-change: " + e.getMessage());
		}		
	}

	/**
	 * 
	 */
	public void executeAnova() {

		DoubleMatrix matrix = null;
		List<String> vars = new ArrayList<String>();
		List<Integer> indices = new ArrayList<Integer>();
		List<String> values = dataset.getVariables().getByName(className).getValues();

		if ( values.size() == classValues.size() ) {
			matrix = dataset.getDoubleMatrix();
			vars = values;
		} else {
			for(int i=0 ; i<values.size() ; i++) {
				if ( classValues.contains(values.get(i)) ) {
					indices.add(i);
					vars.add(values.get(i));
				}
			}
			matrix = dataset.getSubMatrixByColumns(ListUtils.toArray(indices));
		}

		try {
			AnovaTest anova = new AnovaTest(matrix, vars);			
			TestResultList<AnovaTestResult> res = anova.compute();
			
			// apply multiple test correction according to input correction
			DiffExpressionUtils.multipleTestCorrection(res, correction);

			// generating heatmap
			//
			updateJobStatus("60", "generating heatmap");
			int[] columnOrder = ListUtils.order(dataset.getVariables().getByName(className).getValues());
			int[] rowOrder = ListUtils.order(ArrayUtils.toList(res.getStatistics()), true);
			Canvas heatmap = DiffExpressionUtils.generateHeatmap(dataset, className, columnOrder, rowOrder, "statistic", res.getStatistics(), "adj. p-value", res.getAdjPValues());
			String heatmapFilename = getOutdir() + "/" + test + "_heatmap.png";
			try {
				heatmap.save(heatmapFilename);
			} catch (IOException e) {
				printError("ioexception_executeanova_classcomparison", "error generating heatmap", e.toString(), e);
			}
			
			updateJobStatus("80", "saving results");			
			DataFrame dataFrame = new DataFrame(dataset.getFeatureNames().size(), 0);

			dataFrame.addColumn("statistic", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getStatistics()), rowOrder)));
			dataFrame.addColumn("p-value", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getPValues()), rowOrder)));
			dataFrame.addColumn("adj. p-value", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getAdjPValues()), rowOrder)));

			dataFrame.setRowNames(ListUtils.ordered(dataset.getFeatureNames(), rowOrder));

			FeatureData featureData = new FeatureData(dataFrame);

			File file = new File(outdir + "/anova.txt");
			featureData.save(file);
			if ( file.exists() ) {
				result.addOutputItem(new Item("anovafile", file.getName(), "Anova output file", TYPE.FILE, new ArrayList<String>(2), new HashMap<String, String>(2), "Anova output files"));							
			}

			file = new File(outdir + "/anova_table.txt");
			IOUtils.write(file, dataFrame.toString(true, true));
			if ( file.exists() ) {
				result.addOutputItem(new Item("anovatable", file.getName(), "Anova output table", TYPE.FILE, StringUtils.toList("TABLE", ","), new HashMap<String, String>(2), "Anova output files"));											
			}
			
			if ( new File(heatmapFilename).exists() ) {
				result.addOutputItem(new Item(test + "_heatmap", test + "_heatmap.png", test.toUpperCase() + " heatmap", TYPE.IMAGE, new ArrayList<String>(2), new HashMap<String, String>(2), "Heatmap image"));
			}
			
			DiffExpressionUtils.addOutputLists(dataFrame, test, "statistic", result, outdir);
			
		} catch (Exception e) {
			e.printStackTrace();
			abort("exception_executeanova_classcomparison", "error running anova", "error running anova: " + e.getMessage(), "error running anova: " + e.getMessage());
		}		
	}

	/**
	 * @throws InvalidColumnIndexException 
	 * @throws IOException 
	 * 
	 */
	public void executeLimma() {

		Limma limma = null;
		if ( classValues.size() > 2 ) {
			limma = new Limma(babelomicsHomePath + "/bin/diffexp/limma_multiclasses.r");
		} else if ( classValues.size() == 2 ) {
			limma = new Limma(babelomicsHomePath + "/bin/diffexp/limma_twoclasses.r");			
		} else if ( classValues.size() == 1 ) {
			limma = new Limma(babelomicsHomePath + "/bin/diffexp/limma_oneclass.r");
		} else {
			abort("testmismatched_executelimma_classcomparison", "test " + test + " not supported for " + classValues.size() + "-class test", "test " + test + " not supported for " + classValues.size() + "-class test", "test " + test + " not supported for " + classValues.size() + "-class test");												
		}

		//System.out.println("dataset = " + dataset.toString());
		System.out.println("class name = " + className);

		limma.setInputFilename(dataset.getDatasetFile().getAbsolutePath());
		limma.setClasses(dataset.getVariables().getByName(className).getValues());
		limma.setContrast(classValues);

		try {
			TestResultList<LimmaTestResult> res = limma.compute();
			
			// apply multiple test correction according to input correction
			DiffExpressionUtils.multipleTestCorrection(res, correction);

			// generating heatmap
			//
			updateJobStatus("60", "generating heatmap");
			int[] columnOrder = ListUtils.order(dataset.getVariables().getByName(className).getValues());
			int[] rowOrder = ListUtils.order(ArrayUtils.toList(res.getStatistics()), true);
			Canvas heatmap = DiffExpressionUtils.generateHeatmap(dataset, className, columnOrder, rowOrder, "statistic", res.getStatistics(), "adj. p-value", res.getAdjPValues());
			String heatmapFilename = getOutdir() + "/" + test + "_heatmap.png";
			try {
				heatmap.save(heatmapFilename);
			} catch (IOException e) {
				printError("ioexception_executelimma_classcomparison", "error generating heatmap", e.toString(), e);
			}

			updateJobStatus("80", "saving results");
			DataFrame dataFrame = new DataFrame(dataset.getFeatureNames().size(), 0);

			dataFrame.addColumn("statistic", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getStatistics()), rowOrder)));
			//dataFrame.addColumn("lod", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getLods()), rowOrder)));
			dataFrame.addColumn("p-value", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getPValues()), rowOrder)));
			dataFrame.addColumn("adj. p-value", ListUtils.toStringList(ListUtils.ordered(ArrayUtils.toList(res.getAdjPValues()), rowOrder)));

			dataFrame.setRowNames(ListUtils.ordered(dataset.getFeatureNames(), rowOrder));

			FeatureData featureData = new FeatureData(dataFrame);

			File file = new File(outdir + "/limma.txt");
			featureData.save(file);
			if ( file.exists() ) {
				result.addOutputItem(new Item("limmafile", file.getName(), "Limma output file", TYPE.FILE, new ArrayList<String>(2), new HashMap<String, String>(2), "Lima output files"));							
			}

			file = new File(outdir + "/limma_table.txt");
			IOUtils.write(file, dataFrame.toString(true, true));
			if ( file.exists() ) {
				result.addOutputItem(new Item("limmatable", file.getName(), "Limma output table", TYPE.FILE, StringUtils.toList("TABLE", ","), new HashMap<String, String>(2), "Lima output files"));											
			}		
			
			if ( new File(heatmapFilename).exists() ) {
				result.addOutputItem(new Item(test + "_heatmap", test + "_heatmap.png", test.toUpperCase() + " heatmap", TYPE.IMAGE, new ArrayList<String>(2), new HashMap<String, String>(2), "Heatmap image"));
			}
			
			DiffExpressionUtils.addOutputLists(dataFrame, test, "statistic", result, outdir);
			
		} catch (Exception e) {
			e.printStackTrace();
			abort("exception_executelimma_classcomparison", "error running limma", "error running limma: " + e.toString(), "error running limma: " + e.toString());
		}		
	}	
}
