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
import org.bioinfo.commons.exec.Command;
import org.bioinfo.commons.exec.SingleProcess;
import org.bioinfo.commons.io.utils.FileUtils;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ArrayUtils;
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

public class MaSigPro extends BabelomicsTool {

	public MaSigPro() {
		initOptions();
	}

	public void initOptions() {
		options.addOption(OptionFactory.createOption("dataset", "the data"));
		options.addOption(OptionFactory.createOption("contin-class", "class variable"));
		options.addOption(OptionFactory.createOption("series-class", "class class variable"));
	}

	@Override
	public void execute() {
		
		// reading dataset
		//
		Dataset dataset = initDataset(new File(commandLine.getOptionValue("dataset")));
				
		String continClass = commandLine.getOptionValue("contin-class", null);
		String seriesClass = commandLine.getOptionValue("series-class", null);
		
		List<String> continVars = dataset.getVariables().getByName(continClass).getValues();
		List<String> seriesVars = dataset.getVariables().getByName(seriesClass).getValues();

		// masigpro
		//
		updateJobStatus("40", "computing maSigPro test");
		List<String> lines = new ArrayList<String>();
		lines.add("#names\t" + ListUtils.toString(dataset.getSampleNames(), "\t"));
		lines.add("#contin\t" + ListUtils.toString(continVars, "\t"));
		lines.add("#series\t" + ListUtils.toString(seriesVars, "\t"));
		for(int i=0 ; i<dataset.getRowDimension() ; i++) {
			lines.add(dataset.getFeatureNames().get(i) + "\t" + ListUtils.toString(ListUtils.toStringList(dataset.getDoubleMatrix().getRow(i)), "\t"));
		}
		File inputFile = new File(outdir + "/input.maSigPro.txt");
		try {
			IOUtils.write(inputFile, lines);
		} catch (IOException e) {
			abort("ioexception_execute_masigpro", "error writting intermediate file", e.toString(), StringUtils.getStackTrace(e));
		}
		
		List<String> env = new ArrayList<String>();
		env.add("data=" + inputFile.getAbsolutePath());
		env.add("degree=2");
		env.add("Q=0.05");
		env.add("adjust=BH");
		env.add("alfa=0.05");
		env.add("clustermethod=hclust");
		env.add("k=9");
		env.add("main=");
		env.add("outdir=" + outdir);
	
		Command cmd = new Command("/usr/local/R-2.9.2/bin/R CMD BATCH --no-save --no-restore " + System.getenv("BABELOMICS_HOME") + "/bin/masigpro/masigpro.R " + outdir + "/r.log", env);
		SingleProcess sp = new SingleProcess(cmd);
		sp.runSync();
		
		System.out.println("cmd line = " + cmd.getCommandLine());
//		System.out.println("cmd env = " + ListUtils.toString(cmd.getEnvironment(), "\t"));
//		System.out.println("cmd exit value = " + cmd.getExitValue());
		System.out.println("cmd output = " + cmd.getOutput());
		System.out.println("cmd error = " + cmd.getError());
		
		// saving data
		//
		updateJobStatus("80", "saving results");
		File outDirFile = new File(outdir);
		
		File outFile = null;
		File[] outFiles = null;
		
		// general info
		outFiles = FileUtils.listFiles(outDirFile, ".*summary.*");
		for (File f: outFiles) {
			result.addOutputItem(new Item("summaryfile", f.getName(), "Significant genes", TYPE.FILE, new ArrayList<String>(2), new HashMap<String, String>(2), "List of significant genes"));			
		}
		if ( (outFile = new File(outdir + "/pvalues.txt")).exists() ) {
			result.addOutputItem(new Item("pvaluesfile", outFile.getName(), "p-values file", TYPE.FILE, new ArrayList<String>(2), new HashMap<String, String>(2), "p-values and adjusted p-values of global model for all genes"));						
		}
		if ( (outFile = new File(outdir + "/influ_info.txt")).exists() ) {
			result.addOutputItem(new Item("influinfofile", outFile.getName(), "Influence data file", TYPE.FILE, new ArrayList<String>(2), new HashMap<String, String>(2), "Influence data (genes with influential values, possible outliers)"));						
		}
		if ( (outFile = new File(outdir + "/influ_data.png")).exists() ) {
			result.addOutputItem(new Item("infludataimg", outFile.getName(), "Influence data plot", TYPE.IMAGE, new ArrayList<String>(2), new HashMap<String, String>(2), "Influence data (genes with influential values, possible outliers)"));						
		}
		
		// model data
		outFiles = FileUtils.listFiles(outDirFile, ".*coefficients\\.txt");
		for (File f: outFiles) {
			result.addOutputItem(new Item("sigcoeffile", f.getName(), "Significant coefficients", TYPE.FILE, new ArrayList<String>(2), new HashMap<String, String>(2), "Significant coefficients"));			
		}	
		outFiles = FileUtils.listFiles(outDirFile, ".*sig\\.pvalues\\.txt");
		for (File f: outFiles) {
			result.addOutputItem(new Item("sigpvaluefile", f.getName(), "Significant p-values", TYPE.FILE, new ArrayList<String>(2), new HashMap<String, String>(2), "Significant p-values"));			
		}
		outFiles = FileUtils.listFiles(outDirFile, ".*profiles\\.txt");
		for (File f: outFiles) {
			result.addOutputItem(new Item("sigprofilefile", f.getName(), "Significant profiles", TYPE.FILE, new ArrayList<String>(2), new HashMap<String, String>(2), "Significant profiles"));			
		}
		
		
		// visualization
		outFiles = FileUtils.listFiles(outDirFile, ".*heatmap.*");
		for (File f: outFiles) {
			result.addOutputItem(new Item("heatmapimg", f.getName(), "Heatmap", TYPE.IMAGE, new ArrayList<String>(2), new HashMap<String, String>(2), "Heatmap plot"));			
		}
		outFiles = FileUtils.listFiles(outDirFile, ".*Groups\\.png");
		for (File f: outFiles) {
			result.addOutputItem(new Item("groupimg", f.getName(), "Group", TYPE.IMAGE, new ArrayList<String>(2), new HashMap<String, String>(2), "Group plot"));			
		}
		outFiles = FileUtils.listFiles(outDirFile, ".*heatmap.*");
		for (File f: outFiles) {
			result.addOutputItem(new Item("profileimg", f.getName(), "Profile", TYPE.IMAGE, new ArrayList<String>(2), new HashMap<String, String>(2), "Profiles plot"));			
		}
				
		// Summary: list of significant genes (*summary*)
		// P-values and adjusted p-values of global model for all genes (*pvalues.txt) 
		// Influence data (genes with influential values, possible outliers) (*influ_info.txt, *influ_data.png)
		// Significant coefficients 
		
		
		
//
//		// saving data
//		//
//		updateJobStatus("80", "saving results");
//		DataFrame dataFrame = new DataFrame(dataset.getFeatureNames().size(), 0);
//		dataFrame.setRowNames(ListUtils.ordered(dataset.getFeatureNames(), rowOrder));
//
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
	}
}
