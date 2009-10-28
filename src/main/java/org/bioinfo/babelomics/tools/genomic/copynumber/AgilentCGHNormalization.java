package org.bioinfo.babelomics.tools.genomic.copynumber;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.commons.io.utils.FileUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.preprocess.microarray.AgilentCghUtils;
import org.bioinfo.io.file.compress.CompressFactory;
import org.bioinfo.io.file.compress.GenericCompressManager;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;

public class AgilentCGHNormalization extends BabelomicsTool {

	public AgilentCGHNormalization() {
		initOptions();
	}
	
	@Override
	public void initOptions() {
		getOptions().addOption(OptionFactory.createOption("compressed-file", "Compressed file containning the raw files", false));
		getOptions().addOption(OptionFactory.createOption("raw-dir", "Directory where the raw files are located", false));
				
		getOptions().addOption(OptionFactory.createOption("bg-correction", "Background correction: none, substract, half, minimum, movingmin, edwards, normexp. Default background correction: minimum", false));
		getOptions().addOption(OptionFactory.createOption("wa-normalization", "Within arrays normalization: none, median, loess, composite, robustline. Default within arrays normalization: median", false));
		getOptions().addOption(OptionFactory.createOption("ba-normalization", "Between arrays normalization: none, scale, quantile. Default between arrays normalization: none", false));
		getOptions().addOption(OptionFactory.createOption("design", "Design: 1 when Cy3 is the reference, -1 when Cy5 is the reference. Default design: 1", false));
		
		getOptions().addOption(OptionFactory.createOption("sample-filter", "class variable", false));
		getOptions().addOption(OptionFactory.createOption("feature-filter", "class variable", false));
	}

	@Override
	protected void execute() {
		File tmpDir = new File(outdir + "/tmp");
		String compressedFilename = commandLine.getOptionValue("compressed-file", null);
		String rawDirname = commandLine.getOptionValue("raw-dir", null);
		String bgCorrection = commandLine.getOptionValue("bg-correction", "minimum");
		String waNormalization = commandLine.getOptionValue("wa-normalization", "median");
		String baNormalization = commandLine.getOptionValue("ba-normalization", "none");
		String design = commandLine.getOptionValue("design", "1");

		// checking analysis methods and input data (compressed file or directory name)
		//
		if ( compressedFilename == null && rawDirname == null ) {
			abort("missingdata_execute_agilentcghnormalization", "missing input data", "missing input data", "missing input data");						
		}

		try {
			jobStatus.addStatusMessage("10", "reading dataset");
		} catch (FileNotFoundException e) {
			abort("filenotfoundexception_execute_agilentcghnormalization", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
		}

		List<String> rawFilenames = null;
		if ( compressedFilename != null ) {
			// getting raw files from compressed file
			//
			try {
				jobStatus.addStatusMessage("30", "decomprising dataset");
			} catch (FileNotFoundException e) {
				abort("filenotfoundexception_execute_agilentcghnormalization", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
			}

			System.out.println("input dataset = " + compressedFilename);
			System.out.println("tmp dir = " + tmpDir.getAbsolutePath());
			try {
				GenericCompressManager compresor = CompressFactory.getCompressManager(new File(compressedFilename));
				rawFilenames = compresor.decompress(compressedFilename, tmpDir.getAbsolutePath());
			} catch (Exception e) {
				abort("exception_execute_agilentcghnormalization", "error decompressing dataset", e.toString(), StringUtils.getStackTrace(e));
			}

		} else {
			// getting raw files from directory
			//
			tmpDir = new File(rawDirname); 
			File[] rawFiles = FileUtils.listFiles(new File(rawDirname), ".+");
			rawFilenames = new ArrayList<String>(rawFiles.length);
			for(File file: rawFiles) {
				rawFilenames.add(file.getAbsolutePath());
			}
		}

		// sanity check
		//
		if ( rawFilenames == null || rawFilenames.size() == 0 ) {
			abort("missingrawfiles_execute_agilentcghnormalization", "missing raw files", "missing raw files", "missing raw files");
		}

		// creating the cel_files containning the cel files to normalize
		//
//		File celFiles = new File(outdir + "/cel_files.txt");
//		try {
//			IOUtils.write(celFiles, "cel_files\n" + ListUtils.toString(rawFilenames, "\n"));
//		} catch (IOException e) {
//			abort("ioexception_execute_affynormalization", "error writting cel files", e.toString(), StringUtils.getStackTrace(e));
//		}

		System.out.println("raw files = " + ListUtils.toString(rawFilenames, "\n"));

		// normalizing data
		//
		try {
			jobStatus.addStatusMessage("50", "normalizing data");
		} catch (FileNotFoundException e) {
			abort("filenotfoundexception_execute_agilentcghnormalization", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
		}

		AgilentCghUtils.normalization(System.getenv("BABELOMICS_HOME") + "/bin/copynumber/normalization.agilent.r", bgCorrection, waNormalization, baNormalization, design, tmpDir.getAbsolutePath(), outdir);
		
		// saving normalization results
		//
		try {
			jobStatus.addStatusMessage("90", "saving normalization results");
		} catch (FileNotFoundException e) {
			abort("filenotfoundexception_execute_agilentcghnormalization", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
		}

		File file;
		List<String> tags = StringUtils.toList("data,datamatrix,cgh", ",");

		file = new File(outdir + "/cgh_normalization.txt"); 
		if ( file.exists() ) {
				result.addOutputItem(new Item("cghnormalization", file.getName(), "CGH normalization ", TYPE.FILE, tags, new HashMap<String, String>(2), "CGH normalization files"));
		} else {
			printError("error cgh normalization", "error cgh normalization", "error cgh normalization");
		}

		file = new File(outdir + "/cgh_positions.txt"); 
		if ( file.exists() ) {
				result.addOutputItem(new Item("cghpositions", file.getName(), "Positions ", TYPE.FILE, new ArrayList<String>(2), new HashMap<String, String>(2), "CGH normalization files"));
		}

		file = new File(outdir + "/boxplot.png");
		if ( file.exists() ) {
			result.addOutputItem(new Item("boxplot", file.getName(), "Normalized arrays box-plot", TYPE.IMAGE, new ArrayList<String>(2), new HashMap<String, String>(2), "CGH normalization images"));
		}
		
		file = new File(outdir + "/densityplot.png");
		if ( file.exists() ) {
			result.addOutputItem(new Item("densityplot", file.getName(), "Normalized arrays density plot", TYPE.IMAGE, new ArrayList<String>(2), new HashMap<String, String>(2), "CGH normalization images"));
		}
		
		for(String filename: rawFilenames) {
			file = new File(filename + "densityplot.png");
			if ( file.exists() ) {
				result.addOutputItem(new Item("densityplotfor " + filename, file.getName(), "Density plot for " + file.getAbsolutePath(), TYPE.IMAGE, new ArrayList<String>(2), new HashMap<String, String>(2), "CGH arrays density plots"));
			}			
			file = new File(filename + "MAplot.png");
			if ( file.exists() ) {
				result.addOutputItem(new Item("maplotfor " + filename, file.getName(), "MA plot for " + file.getAbsolutePath(), TYPE.IMAGE, new ArrayList<String>(2), new HashMap<String, String>(2), "CGH arrays MA plots"));
			}
		}
	}
}
