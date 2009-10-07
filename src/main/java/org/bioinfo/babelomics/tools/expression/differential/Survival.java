package org.bioinfo.babelomics.tools.expression.differential;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.security.InvalidParameterException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.apache.commons.math.MathException;
import org.apache.commons.math.linear.MatrixIndexException;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.babelomics.tools.expression.DifferentialAnalysis;
import org.bioinfo.collections.exceptions.InvalidColumnIndexException;
import org.bioinfo.collections.matrix.DataFrame;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.data.dataset.FeatureData;
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
	}

	@Override
	public void execute() {
		
		// reading dataset
		//
		Dataset dataset = initDataset(new File(commandLine.getOptionValue("dataset")));
				
		String test = commandLine.getOptionValue("test", null);
		String timeClass = commandLine.getOptionValue("time-class", null);
		String censoredClass = commandLine.getOptionValue("censored-class", null);

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
		} catch (Exception e) {
			abort("exception_run_cox", "error running cox test", e.toString(), StringUtils.getStackTrace(e));
		}

		int[] columnOrder = ListUtils.order(vars);
		int[] rowOrder = ListUtils.order(ListUtils.toList(res.getStatistics()), true);

		// generating heatmap
		//
		updateJobStatus("60", "generating heatmap");
		Canvas heatmap = DiffExpressionUtils.generateHeatmap(dataset, timeClass, columnOrder, rowOrder, "coeff.", res.getCoefs(), "adj. p-value", res.getAdjPValues());
		String heatmapFilename = getOutdir() + "/cox_heatmap";
		try {
			heatmap.save(heatmapFilename);
		} catch (IOException e) {
			printError("ioexception_cox_cox", "error generating heatmap", e.toString(), e);
		}

		// saving data
		//
		updateJobStatus("80", "saving results");
		DataFrame dataFrame = new DataFrame(dataset.getFeatureNames().size(), 0);
		dataFrame.setRowNames(ListUtils.ordered(dataset.getFeatureNames(), rowOrder));

		try {
			dataFrame.addColumn("statistic", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getStatistics()), rowOrder)));
			dataFrame.addColumn("coeff.", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getCoefs()), rowOrder)));
			dataFrame.addColumn("p-value", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getPValues()), rowOrder)));
			dataFrame.addColumn("adj. p-value", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getAdjPValues()), rowOrder)));

			FeatureData featureData = new FeatureData(dataFrame);
			featureData.write(new File(getOutdir() + "/cox.txt"));
			result.addOutputItem(new Item("cox_file", "cox.txt", "Cox output file", TYPE.FILE));
			
			IOUtils.write(new File(getOutdir() + "/cox_table.txt"), dataFrame.toString(true, true));			
			List<String> tags = new ArrayList<String>();
			tags.add("TABLE");
			result.addOutputItem(new Item("cox_table", "cox_table.txt", "Cox table", TYPE.FILE, tags, new HashMap<String, String>(2)));
		} catch (Exception e) {
			printError("ioexception_cox_cox", "error saving results", e.toString(), e);
		}

		if ( new File(heatmapFilename + ".png").exists() ) {
			result.addOutputItem(new Item("cox_heatmap", "cox_heatmap.png", "Cox heatmap", TYPE.IMAGE));
		}
	}
}