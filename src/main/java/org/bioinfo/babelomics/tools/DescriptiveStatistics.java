package org.bioinfo.babelomics.tools;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import org.bioinfo.babelomics.methods.expression.clustering.ClusteringUtils;
import org.bioinfo.chart.BoxPlotChart;
import org.bioinfo.chart.HistogramChart;
import org.bioinfo.collections.tree.multiple.MultipleTree;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.data.format.io.exception.InvalidFileFormatException;
import org.bioinfo.data.format.io.parser.NewickParser;
import org.bioinfo.math.data.DoubleMatrix;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;
import org.jfree.chart.ChartUtilities;


public class DescriptiveStatistics extends BabelomicsTool {
	Dataset dataset;
	String test;
	String className;
	List<String> values;
	List<Double> doubleVars;
	String correction;
	String type;
	 
	@Override
	public void initOptions() {
		options.addOption(OptionFactory.createOption("datalist", "the feature data containig the ranked list"));		
		options.addOption(OptionFactory.createOption("column", "",false));
		options.addOption(OptionFactory.createOption("tree", "",false));
		options.addOption(OptionFactory.createOption("histogram", "",false));
		options.addOption(OptionFactory.createOption("boxplot", "",false));
		//options.addOption(OptionFactory.createOption("histogramboxplot", "",false,false));
		options.addOption(OptionFactory.createOption("class", "",false,true));
		options.addOption(OptionFactory.createOption("width", "",false,true));
		options.addOption(OptionFactory.createOption("height", "",false,true));
	}

	
	@Override
	protected void execute() {
		dataset = null;
		className = commandLine.getOptionValue("class", null);
		
		if(commandLine.getOptionValue("histogram",null)!=null){
			BoxPlotChart bp = (commandLine.getOptionValue("boxplot", null) != null) ? new BoxPlotChart("Box-plot", "", ""):null;
			HistogramChart hc = new HistogramChart("Histogram chart", "label", "frec");
			preloadSeries(hc,bp);
			finishGraph(hc, bp);
		}
		
		else if (commandLine.getOptionValue("boxplot",null)!=null){
			 BoxPlotChart bp = new BoxPlotChart("Box-plot", "", "");
			 preloadSeries(null,bp);
			 finishGraph(null, bp);
		}
		
		
		else if(commandLine.getOptionValue("tree", null) != null){
			 doTree();
		}
		

	}


	private void preloadSeries(HistogramChart hc,BoxPlotChart bp) {
		try {
			if (className!=null){
				dataset = new Dataset(new File(commandLine.getOptionValue("datalist")));
				values = dataset.getVariables().getByName(className).getLabels();
				for (String str: values){
					int[] colIndexByVariableValue = dataset.getColumnIndexesByVariableValue(className, str);
					doubleVars = new ArrayList<Double>(colIndexByVariableValue.length);
					DoubleMatrix matrixByVal =dataset.getSubMatrixByColumns(colIndexByVariableValue);
					addSeries(matrixByVal, hc,bp, str);
				
					
				}
			}
			else{
				dataset = new Dataset(new File(commandLine.getOptionValue("datalist")));
				addSeries(dataset.getDoubleMatrix(),hc,bp, "PUT HERE THE COLUMN_NAME/NUMBER");
			}
			
		} catch (IOException e4) {
			// TODO Auto-generated catch block
			e4.printStackTrace();
		}
	}


	private void addSeries(DoubleMatrix matrixByVal, HistogramChart hc, BoxPlotChart bp, String str) {
		if (hc != null){
			hc.addSeries(matrixByVal.getColumn(0),str);
		}
		if(bp != null ){
			bp.addSeries(matrixByVal.getColumn(0),str, "hola");
		}
	}
	

	private void finishGraph(HistogramChart hc, BoxPlotChart bp){
		String imgFilename ="";
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
		if (hc != null){
			imgFilename = this.getOutdir() + "/histogram.png";
			try {
				ChartUtilities.saveChartAsPNG(new File(imgFilename), hc, 700, 500);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			result.addOutputItem(new Item("histogram_image","histogram.png","histogram image (png format)",TYPE.IMAGE, Arrays.asList("IMAGE","HISTOGRAM_IMAGE"),new HashMap<String,String>(),"histogram image"));
		}
		
		if (bp != null){
			imgFilename = this.getOutdir() + "/boxplot.png";
			try {
				ChartUtilities.saveChartAsPNG(new File(imgFilename), bp, 700, 500);
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}
			result.addOutputItem(new Item("boxplot_image","boxplot.png","boxplot image (png format)",TYPE.IMAGE, Arrays.asList("IMAGE","BOXPLOT_IMAGE"),new HashMap<String,String>(),"boxplot image"));
		}
		
		
		// Save histogram
		progress++;
		updateJobStatus(""+(progress*100/finalProgress),"saving graph");
		
	}
	
	private void doTree() {
		NewickParser nwParser = new NewickParser();
		MultipleTree tree;
		try {
			tree = nwParser.parse(new File(commandLine.getOptionValue("datalist")));
			String imgFilename = this.getOutdir() + "/tree.png";
			ClusteringUtils.saveImageTree(tree, "Tree",  imgFilename, true);
			result.addOutputItem(new Item("tree_image","tree.png","tree image (png format)",TYPE.IMAGE, Arrays.asList("IMAGE","TREE_IMAGE"),new HashMap<String,String>(),"tree image"));
		} catch (IOException e) {
			e.printStackTrace();
		} catch (InvalidFileFormatException e) {
			e.printStackTrace();
		}

	}

}
