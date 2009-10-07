package org.bioinfo.babelomics.tools.expression.differential;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.commons.exec.Command;
import org.bioinfo.commons.exec.SingleProcess;
import org.bioinfo.commons.io.utils.FileUtils;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.dataset.Dataset;
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
	
		//Command cmd = new Command("/usr/local/R-2.9.2/bin/R CMD BATCH --no-save --no-restore " + System.getenv("BABELOMICS_HOME") + "/bin/masigpro/masigpro.R " + outdir + "/r.log", env);
		Command cmd = new Command("R CMD BATCH --no-save --no-restore " + System.getenv("BABELOMICS_HOME") + "/bin/masigpro/masigpro.R " + outdir + "/r.log", env);
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
			result.addOutputItem(new Item("summaryfile", f.getName(), "Significant genes for '" + getCleanName(f) + "'", TYPE.FILE, new ArrayList<String>(2), new HashMap<String, String>(2), "List of significant genes"));			
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
			result.addOutputItem(new Item("sigcoeffile", f.getName(), "Significant coefficients for '" + getCleanName(f) + "'", TYPE.FILE, new ArrayList<String>(2), new HashMap<String, String>(2), "Significant coefficients"));			
		}	
		outFiles = FileUtils.listFiles(outDirFile, ".*sig\\.pvalues\\.txt");
		for (File f: outFiles) {
			result.addOutputItem(new Item("sigpvaluefile", f.getName(), "Significant p-values for '" + getCleanName(f) + "'", TYPE.FILE, new ArrayList<String>(2), new HashMap<String, String>(2), "Significant p-values"));			
		}
		outFiles = FileUtils.listFiles(outDirFile, ".*profiles\\.txt");
		for (File f: outFiles) {
			result.addOutputItem(new Item("sigprofilefile", f.getName(), "Significant profiles for '" + getCleanName(f) + "'", TYPE.FILE, new ArrayList<String>(2), new HashMap<String, String>(2), "Significant profiles"));			
		}
		
		
		// visualization
		outFiles = FileUtils.listFiles(outDirFile, ".*heatmap.*");
		for (File f: outFiles) {
			result.addOutputItem(new Item("heatmapimg", f.getName(),  "'" + getCleanName(f) + "' heatmap", TYPE.IMAGE, new ArrayList<String>(2), new HashMap<String, String>(2), "Heatmap plot"));			
		}
		outFiles = FileUtils.listFiles(outDirFile, ".*Groups\\.png");
		for (File f: outFiles) {
			result.addOutputItem(new Item("groupimg", f.getName(), "'" + getCleanName(f) + "' group", TYPE.IMAGE, new ArrayList<String>(2), new HashMap<String, String>(2), "Group plot"));			
		}
		outFiles = FileUtils.listFiles(outDirFile, ".*heatmap.*");
		for (File f: outFiles) {
			result.addOutputItem(new Item("profileimg", f.getName(), "'" + getCleanName(f) + "' profile", TYPE.IMAGE, new ArrayList<String>(2), new HashMap<String, String>(2), "Profiles plot"));			
		}				
	}
	
	private String getCleanName(File file) {
		String res = file.getName().replace(".txt", "").replace(".png", "");
		res = res.replace("groups_", "").replace("summary", "").replace("_sig.profiles", "").replace("_sig.pvalues", "").replace("_coefficients", "");
		res = res.replace("_Groups", "").replace("_heatmap", "").replace("_Profiles", "").replace("_sig.pvalues", "_heatmap");
		res = res.replace("_", " ");
		return res;
	}
	
}
