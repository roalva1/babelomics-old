package org.bioinfo.babelomics.tools.expression.normalization;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.bioinfo.babelomics.methods.expression.AffyUtils;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.chart.BoxPlotChart;
import org.bioinfo.collections.exceptions.InvalidColumnIndexException;
import org.bioinfo.commons.Config;
import org.bioinfo.commons.io.utils.FileUtils;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.commons.utils.MapUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.io.file.compress.CompressFactory;
import org.bioinfo.io.file.compress.GenericCompressManager;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;
import org.jfree.chart.ChartUtilities;


public class AffyExpressionNormalization extends BabelomicsTool {

	private String aptBinPath = System.getenv("BABELOMICS_HOME") + "/bin/apt";

	public AffyExpressionNormalization() {
		initOptions();
	}

	@Override
	public void initOptions() {
		getOptions().addOption(OptionFactory.createOption("compressed-file", "Compressed file containning the raw files", false));
		getOptions().addOption(OptionFactory.createOption("raw-dir", "Directory where the raw files are located", false));
		getOptions().addOption(OptionFactory.createOption("no-cel-convert", "Disable to convert CEL files to GCOS text file format (by default, the conversion is enabled)", false, false));
		getOptions().addOption(OptionFactory.createOption("rma", "RMA analysis", false, false));
		getOptions().addOption(OptionFactory.createOption("plier", "Plier analysis", false, false));
		getOptions().addOption(OptionFactory.createOption("calls", "Present-absent calls analysis", false, false));
		getOptions().addOption(OptionFactory.createOption("sample-filter", "class variable", false));
		getOptions().addOption(OptionFactory.createOption("feature-filter", "class variable", false));

	}

	@Override
	public void execute() {
		File tmpDir = new File(outdir + "/tmp");
		String compressedFilename = commandLine.getOptionValue("compressed-file", null);
		String rawDirname = commandLine.getOptionValue("raw-dir", null);
		boolean celConvert = !commandLine.hasOption("no-cel-convert");
		boolean rma = commandLine.hasOption("rma");
		boolean plier = commandLine.hasOption("plier");
		boolean calls = commandLine.hasOption("calls");

		// checking analysis methods and input data (compressed file or directory name)
		//
		if ( !rma && !plier && !calls ) {
			abort("missinganalysis_execute_affynormalization", "missing analysis", "missing analysis, valid values are rma, plier and calls", "missing analysis, valid values are rma, plier and calls");			
		}
		if ( compressedFilename == null && rawDirname == null ) {
			abort("missingdata_execute_affynormalization", "missing input data", "missing input data", "missing input data");						
		}

		try {
			jobStatus.addStatusMessage("10", "reading dataset");
		} catch (FileNotFoundException e) {
			abort("filenotfoundexception_execute_affynormalization", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
		}

		List<String> rawFilenames = null;
		if ( compressedFilename != null ) {
			// getting raw files from compressed file
			//
			try {
				jobStatus.addStatusMessage("30", "decomprising dataset");
			} catch (FileNotFoundException e) {
				abort("filenotfoundexception_execute_affynormalization", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
			}

			System.out.println("input dataset = " + compressedFilename);
			System.out.println("tmp dir = " + tmpDir.getAbsolutePath());
			try {
				GenericCompressManager compresor = CompressFactory.getCompressManager(new File(compressedFilename));
				rawFilenames = compresor.decompress(compressedFilename, tmpDir.getAbsolutePath());
			} catch (Exception e) {
				abort("exception_execute_affynormalization", "error decompressing dataset", e.toString(), StringUtils.getStackTrace(e));
			}

		} else {
			// getting raw files from directory
			//
			File[] rawFiles = FileUtils.listFiles(new File(rawDirname), ".+.CEL");
			rawFilenames = new ArrayList<String>(rawFiles.length);
			for(File file: rawFiles) {
				rawFilenames.add(file.getAbsolutePath());
			}
		}

		// sanity check
		//
		if ( rawFilenames == null || rawFilenames.size() == 0 ) {
			abort("missingrawfiles_execute_affynormalization", "missing raw files", "missing raw files", "missing raw files");
		}

		// creating the cel_files containning the cel files to normalize
		//
		File celFiles = new File(outdir + "/cel_files.txt");
		try {
			IOUtils.write(celFiles, "cel_files\n" + ListUtils.toString(rawFilenames, "\n"));
		} catch (IOException e) {
			abort("ioexception_execute_affynormalization", "error writting cel files", e.toString(), StringUtils.getStackTrace(e));
		}

		System.out.println("raw files = " + ListUtils.toString(rawFilenames, "\n"));

		// converting to CEL text
		//
		if ( celConvert ) {			
			try {
				jobStatus.addStatusMessage("40", "converting CEL to GCOS text file format");
			} catch (FileNotFoundException e) {
				abort("filenotfoundexception_execute_affynormalization", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
			}

			AffyUtils.aptCelConvert(aptBinPath + "/apt-cel-convert", tmpDir.getAbsolutePath(), celFiles);			

			File[] rawFiles = FileUtils.listFiles(tmpDir, ".+.CEL", true);
			rawFilenames = ListUtils.toStringList(rawFiles);
			System.out.println("-----------> raw file names = " + ListUtils.toString(rawFilenames, ","));

			try {
				IOUtils.write(celFiles, "cel_files\n" + ListUtils.toString(rawFilenames, "\n"));
			} catch (IOException e) {
				abort("ioexception_execute_affynormalization", "error writting cel files", e.toString(), StringUtils.getStackTrace(e));
			}
		}		

		Map<String, String> chipInfo = null;
		try {
			Config config = new Config(System.getenv("BABELOMICS_HOME") + "/conf/apt.conf");		
			String chipName = getChipName(rawFilenames, config.getKeys());
			String infoStr = config.getProperty(chipName);
			chipInfo = MapUtils.stringToMap(infoStr);
			chipInfo.put("name", chipName);

			System.out.println(" chip info = " + chipInfo.toString());
		} catch (Exception e) {
			abort("exception_execute_affynormalization", "error getting chip type from raw files", e.toString(), StringUtils.getStackTrace(e));
		}

		String chipType = chipInfo.get("type");
		if ( chipType == null ) {
			abort("exception_execute_affynormalization", "could not find out the chip type", "could not find out the chip type", "could not find out the chip type");			
		}

		if ( !"3-prime".equalsIgnoreCase(chipType) && !"whole-transcript".equalsIgnoreCase(chipType) ) {
			abort("exception_execute_affynormalization", "array type (" + chipType + ") not supported", "array type (" + chipType + ") not supported", "array type (" + chipType + ") not supported");			
		}
		// normalizing data
		//
		try {
			jobStatus.addStatusMessage("50", "normalizing data");
		} catch (FileNotFoundException e) {
			abort("filenotfoundexception_execute_affynormalization", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
		}

		List<String> analysis = new ArrayList<String>();
		if ( "3-prime".equalsIgnoreCase(chipInfo.get("type")) ) {

			if ( rma ) analysis.add("rma");
			if ( calls ) analysis.add("pm-mm,mas5-detect.calls=1.pairs=1");
			if ( plier ) analysis.add("plier-mm");
			
			AffyUtils.aptProbesetSummarize(aptBinPath + "/apt-probeset-summarize", outdir, analysis, new File(System.getenv("BABELOMICS_HOME") + chipInfo.get("cdf")), celFiles);
			
		} else {
			
			if ( rma ) analysis.add("rma");
			if ( calls ) analysis.add("dabg");
			if ( plier ) analysis.add("plier-gcbg");
			
			AffyUtils.aptProbesetSummarize(aptBinPath + "/apt-probeset-summarize", outdir, analysis, new File(System.getenv("BABELOMICS_HOME") + chipInfo.get("cdf")), celFiles);
		}

		//		System.err.println("cmd output: " + sp.getRunnableProcess().getOutput());
		//		System.err.println("cmd error: " + sp.getRunnableProcess().getError());

		// saving normalization results
		//
		try {
			jobStatus.addStatusMessage("90", "saving normalization results");
		} catch (FileNotFoundException e) {
			abort("filenotfoundexception_execute_affynormalization", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
		}

		File file;
		List<String> tags = StringUtils.toList("data,datamatrix,expression", ",");

		file = new File(outdir + "/rma.summary.txt"); 
		if ( file.exists() ) {
			try {
				IOUtils.write(file, cleanLines(IOUtils.readLines(file)));
				saveAsDataset(file);
				
				result.addOutputItem(new Item("rma.summary", file.getName(), "Summary ", TYPE.FILE, tags, new HashMap<String, String>(2), "RMA analysis"));
				saveBoxPlot(file, "RMA box-plot", "rmaimg", "RMA analysis");				
			} catch (Exception e) {
				e.printStackTrace();
				printError("error saving rma results", "error saving rma results", "error saving rma results");
			}
		}

		file = new File(outdir + "/plier-mm.summary.txt"); 
		if ( file.exists() ) {
			try {
				IOUtils.write(file, cleanLines(IOUtils.readLines(file)));
				saveAsDataset(file);
				
				result.addOutputItem(new Item("plier-mm.summary", file.getName(), "MM summary ", TYPE.FILE, tags, new HashMap<String, String>(2), "Plier analysis"));								
				saveBoxPlot(file, "Plier MM box-plot", "plierimg", "Plier analysis");				
			} catch (Exception e) {
				e.printStackTrace();
			}
		}

		file = new File(outdir + "/plier-gcbg.summary.txt"); 
		if ( file.exists() ) {
			try {
				IOUtils.write(file, cleanLines(IOUtils.readLines(file)));
				saveAsDataset(file);
				
				result.addOutputItem(new Item("plier-gcbg.summary", file.getName(), "GCBG summary ", TYPE.FILE, tags, new HashMap<String, String>(2), "Plier analysis"));								
				saveBoxPlot(file, "Plier GCBG box-plot", "plierimg", "Plier analysis");				
			} catch (Exception e) {
				e.printStackTrace();
			}
		}

		file = new File(outdir + "/pm-mm.mas5-detect.summary.txt"); 
		if ( file.exists() ) {
			try {
				IOUtils.write(file, cleanLines(IOUtils.readLines(file)));
				saveAsDataset(file);
				
				result.addOutputItem(new Item("pm-mm.summary", file.getName(), "PM-MM summary ", TYPE.FILE, tags, new HashMap<String, String>(2), "Present-absent calls"));								
				saveBoxPlot(file, "PM-MM box-plot", "pmmmimg", "Present-absent calls");				
			} catch (Exception e) {
				e.printStackTrace();
			}
		}

		file = new File(outdir + "/pm-mm.mas5-detect.calls.txt"); 
		if ( file.exists() ) {
			try {
				IOUtils.write(file, cleanLines(IOUtils.readLines(file)));
				result.addOutputItem(new Item("pm-mm.calls", file.getName(), "Calls ", TYPE.FILE, new ArrayList<String>(2), new HashMap<String, String>(2), "Present-absent calls"));								
			} catch (IOException e) {
				e.printStackTrace();
			}
		}

		file = new File(outdir + "/dabg.summary.txt"); 
		if ( file.exists() ) {
			try {
				IOUtils.write(file, cleanLines(IOUtils.readLines(file)));
				saveAsDataset(file);
				
				result.addOutputItem(new Item("dabg.summary", file.getName(), "DABG summary ", TYPE.FILE, tags, new HashMap<String, String>(2), "Present-absent calls"));								
				saveBoxPlot(file, "DABG box-plot", "dabgimg", "Present-absent calls");				
			} catch (Exception e) {
				e.printStackTrace();
			}
		}
	}

	/**
	 * 
	 * @param rawFilenames
	 * @param chipNames
	 * @return
	 * @throws Exception
	 */
	private String getChipName(List<String> rawFilenames, List<String> chipNames) throws Exception {
		List<String> lines = null;
		List<String> results = new ArrayList<String>(rawFilenames.size());

		System.out.println("----------> getChipPath, chip names = " + ListUtils.toString(chipNames, ", "));

		for(int i=0 ; i<rawFilenames.size() ; i++) {
			results.add(i, null);
			lines = IOUtils.head(new File(rawFilenames.get(i)), 20);
			for(String chipName: chipNames) {
				for(String line: lines) {
					if ( line.contains(chipName) ) {
						System.out.println("found " + chipName + " in file " + rawFilenames.get(i));
						results.add(i, chipName);
						break;
					}
				}
				if ( results.get(i) != null ) break;
			}
		}

		String chipName = results.get(0);
		for (int i=0 ; i<rawFilenames.size() ; i++) {
			if ( ! chipName.equalsIgnoreCase(results.get(i)) ) {
				String msg = "mismatch CEL files corresponding to different chips.\n";
				for(int j=0 ; j<rawFilenames.size() ; j++) {
					msg = msg + "Raw file '" + new File(rawFilenames.get(j)).getName() + "' is a '" + results.get(j) + "' array\n";
				}
				throw new Exception(msg);
			}
		}
		return chipName;
	}

	/**
	 * 
	 * @param lines
	 * @return
	 */
	private List<String> cleanLines(List<String> lines) {
		List<String> result = new ArrayList<String>();
		for(String line: lines) {
			if ( line.startsWith("#") ) {
			} else if ( line.startsWith("probeset_id") ) {
				//line = line.replace("\t", ",");
				//result.add(line.replace("probeset_id,", "#NAMES\t"));
				result.add(line.replace("probeset_id", "#NAMES"));
			} else {
				result.add(line);
			}
		}
		return result;
	}
	
	/**
	 * 
	 * @param file
	 * @throws IOException
	 * @throws InvalidColumnIndexException
	 */
	private void saveAsDataset(File file) throws IOException, InvalidColumnIndexException {
		Dataset dataset = new Dataset(file);
		dataset.load();
		dataset.save();
	}
	
	public void saveBoxPlot(File file, String title, String resultId, String group) throws IOException, InvalidColumnIndexException {
		File imgFile = new File(file.getAbsolutePath().replace(".txt", ".png"));
		Dataset dataset = new Dataset(file, true);
		BoxPlotChart bpc = new BoxPlotChart(title, "", "");
		bpc.getLegend().setVisible(false);
		for(int i=0; i<dataset.getColumnDimension(); i++) {
			bpc.addSeries(dataset.getDoubleMatrix().getColumn(i), "samples", dataset.getSampleNames().get(i));
		}
		try {
			ChartUtilities.saveChartAsPNG(imgFile, bpc, 400, 256);
			if ( imgFile.exists() ) {
				result.addOutputItem(new Item(resultId, imgFile.getName(), title, TYPE.IMAGE, new ArrayList<String>(2), new HashMap<String, String>(2), group));							
			} else {
				printError("error saving " + title, "error saving " + title, "error saving " + title);
			}
		} catch (IOException e) {
			e.printStackTrace();
			printError("error generating " + title, "error generating " + title, "error generating " + title);
		}
	}

}
