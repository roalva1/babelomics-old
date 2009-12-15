package org.bioinfo.babelomics.tools;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import org.bioinfo.chart.HistogramChart;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.math.data.DoubleMatrix;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;
import org.jfree.chart.ChartUtilities;
public class Histogram extends BabelomicsTool {
	Dataset dataset;
	String test;
	String className;
	List<String> values;
	List<Double> doubleVars;
	String correction;
	
	@Override
	public void initOptions() {
		options.addOption(OptionFactory.createOption("datalist", "the feature data containig the ranked list"));		
		options.addOption(OptionFactory.createOption("column", "",false));
		options.addOption(OptionFactory.createOption("class", "",false,true));
		options.addOption(OptionFactory.createOption("width", "",false,true));
		options.addOption(OptionFactory.createOption("height", "",false,true));
	}

	
	@Override
	protected void execute() {
		dataset = null;
		className = commandLine.getOptionValue("class", null);
		HistogramChart hc = new HistogramChart("Histogram chart", "label", "values");
		try {
			if (className!=null){
				dataset = new Dataset(new File(commandLine.getOptionValue("datalist")));
				values = dataset.getVariables().getByName(className).getLabels();
				for (String str: values){
					int[] colIndexByVariableValue = dataset.getColumnIndexesByVariableValue(className, str);
					doubleVars = new ArrayList<Double>(colIndexByVariableValue.length);
					DoubleMatrix matrixByVal =dataset.getSubMatrixByColumns(colIndexByVariableValue);
					addSeries(matrixByVal, hc, str);
				}
			
			}
			else{
				dataset = new Dataset(new File(commandLine.getOptionValue("datalist")));
				addSeries(dataset.getDoubleMatrix(),hc,"PUT HERE THE COLUMN_NAME/NUMBER");
			}
			
		} catch (IOException e4) {
			// TODO Auto-generated catch block
			e4.printStackTrace();
		}
		int progress = 1;
		int finalProgress = 3;
		
		try {
			jobStatus.addStatusMessage("" + (progress*100/finalProgress), "reading ranked list");
		} catch (FileNotFoundException e) {
			abort("filenotfoundexception_execute_preprocessing", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
		}

		progress ++;
		updateJobStatus(""+progress,"reading ok");
		
		
		// Generate histogram
		String imgFilename = this.getOutdir() + "/histogram.png";
		//int columnDimension = dataset.getColumnDimension();		
		
		
//		if (className!=null){
//			hc.addSeries(dataset.getDoubleMatrix().getColumn(0),"PUT HERE THE COLUMN_NAME/NUMBER");
//		
//		}
		
//		else{
//			hc.addSeries(dataset.getDoubleMatrix().getColumn(0),"PUT HERE THE COLUMN_NAME/NUMBER");
//		}
			
		// Save histogram
		
		progress++;
		updateJobStatus(""+(progress*100/finalProgress),"saving graph");

		try {
			ChartUtilities.saveChartAsPNG(new File(imgFilename), hc, 700, 500);
			
			result.addOutputItem(new Item("histogram_image","histogram.png","histogram image (png format)",TYPE.IMAGE, Arrays.asList("IMAGE","HISTOGRAM_IMAGE"),new HashMap<String,String>(),"histogram image"));
			//result.addOutputItem(item)
		} catch (IOException e) {
			e.printStackTrace();			
		}
		
		
	}



	private void addSeries(DoubleMatrix matrixByVal, HistogramChart hc, String str) {
		hc.addSeries(matrixByVal.getColumn(0),str);
	}
	
}