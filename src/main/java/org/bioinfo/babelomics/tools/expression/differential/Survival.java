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
import org.bioinfo.math.result.CoxTestResult;
import org.bioinfo.math.result.TestResultList;
import org.bioinfo.math.stats.survival.CoxTest;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;

public class Survival extends BabelomicsTool {

	public Survival() {
	}

	public void initOptions() {
		options.addOption(OptionFactory.createOption("dataset", "the data"));
		options.addOption(OptionFactory.createOption("test", "the test, possible values: cox"));
		options.addOption(OptionFactory.createOption("time-class", "class variable", false));
		options.addOption(OptionFactory.createOption("censored-class", "class class variable", false));
		options.addOption(OptionFactory.createOption("correction", "Multiple-test correction: fdr, bh, by, bonferroni, hochberg, hold"));
	}

	@Override
	public void execute() {
		
		// reading dataset
		//
		Dataset dataset = initDataset(new File(commandLine.getOptionValue("dataset")));
				
		String test = commandLine.getOptionValue("test", null);
		String timeClass = commandLine.getOptionValue("time-class", null);
		String censoredClass = commandLine.getOptionValue("censored-class", null);
		String correction = commandLine.getOptionValue("correction", "fdr");

		if ( ! "cox".equalsIgnoreCase(test) ) {
			abort("unknowntest_execute_survival", "unknown test (" + test + ")", "unknown test (" + test + ")", "unknown test (" + test + ")");
		}
		
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

		// cox test
		//
		updateJobStatus("40", "computing cox test");
		CoxTest coxtest = new CoxTest();
		TestResultList<CoxTestResult> res = null;
		try {
			res = coxtest.compute(dataset.getDoubleMatrix(), timeVars, censoredVars);
			
			// apply multiple test correction according to input correction
			DiffExpressionUtils.multipleTestCorrection(res, correction);						
		} catch (Exception e) {
			abort("exception_run_cox", "error running cox test", e.toString(), StringUtils.getStackTrace(e));
		}

		int[] columnOrder = ListUtils.order(vars);
		int[] rowOrder = ListUtils.order(ArrayUtils.toList(res.getStatistics()), true);

		// generating heatmap
		//
		updateJobStatus("60", "generating heatmap");
		Canvas heatmap = DiffExpressionUtils.generateHeatmap(dataset, timeClass, columnOrder, rowOrder, "coeff.", res.getCoefs(), "adj. p-value", res.getAdjPValues());
		String heatmapFilename = getOutdir() + "/cox_heatmap.png";
		try {
			heatmap.save(heatmapFilename);
		} catch (IOException e) {
			printError("ioexception_cox_cox", "error generating heatmap", e.toString(), e);
		}

		// saving data
		//
		updateJobStatus("80", "saving results");
		DataFrame dataFrame = new DataFrame(dataset.getFeatureNames().size(), 0);

		try {
			dataFrame.addColumn("statistic", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getStatistics()), rowOrder)));
			dataFrame.addColumn("coeff.", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getCoefs()), rowOrder)));
			dataFrame.addColumn("p-value", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getPValues()), rowOrder)));
			dataFrame.addColumn("adj. p-value", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getAdjPValues()), rowOrder)));

			dataFrame.setRowNames(ListUtils.ordered(dataset.getFeatureNames(), rowOrder));

			FeatureData featureData = new FeatureData(dataFrame);
			featureData.save(new File(getOutdir() + "/cox.txt"));
			result.addOutputItem(new Item("cox_file", "cox.txt", "Cox output file", TYPE.FILE, new ArrayList<String>(2), new HashMap<String, String>(2), "Output files"));
			
			IOUtils.write(new File(getOutdir() + "/cox_table.txt"), dataFrame.toString(true, true));			
			List<String> tags = new ArrayList<String>();
			tags.add("TABLE");
			result.addOutputItem(new Item("cox_table", "cox_table.txt", "Cox table", TYPE.FILE, tags, new HashMap<String, String>(2), "Output files"));
			
			DiffExpressionUtils.addOutputLists(dataFrame, test, "statistic", result, outdir);
			
		} catch (Exception e) {
			printError("ioexception_cox_cox", "error saving results", e.toString(), e);
		}

		if ( new File(heatmapFilename).exists() ) {
			result.addOutputItem(new Item("cox_heatmap", "cox_heatmap.png", "Cox heatmap", TYPE.IMAGE, new ArrayList<String>(2), new HashMap<String, String>(2), "Heatmap image"));
		}
	}
}
