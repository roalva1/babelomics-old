package org.bioinfo.babelomics.tools.expression.normalization;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.bioinfo.babelomics.exception.InvalidParameterException;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.chart.BoxPlotChart;
import org.bioinfo.commons.io.utils.FileUtils;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ArrayUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.commons.utils.MapUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.data.list.exception.InvalidIndexException;
import org.bioinfo.io.file.compress.CompressFactory;
import org.bioinfo.io.file.compress.GenericCompressManager;
import org.bioinfo.microarray.AffymetrixExpresionUtils;
import org.bioinfo.microarray.ExpressionUtils;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;
import org.jfree.chart.ChartUtilities;

public class ExpressionNormalizationTool extends BabelomicsTool {

	int nbChannels = 0;
	String technology = "";
	List<String> rawFileNames = null;
	File tmpDir = null;

	public ExpressionNormalizationTool() {
		initOptions();
	}

	@Override
	public void initOptions() {
		getOptions().addOption(OptionFactory.createOption("compressed-file", "Compressed file containning the raw files", false));
		getOptions().addOption(OptionFactory.createOption("raw-dir", "Directory where the raw files are located", false));
		getOptions().addOption(OptionFactory.createOption("compressed-file-tags", "Data tags separated by commas, valid values are: microarray, expression, affymetrix, agilent, genepix, one-channel, two-channels", false));

		// affy normalization
		//
		getOptions().addOption(OptionFactory.createOption("no-cel-convert", "Disable to convert CEL files to GCOS text file format (by default, the conversion is enabled). Only for Affymetrix normalization.", false, false));
		getOptions().addOption(OptionFactory.createOption("rma", "RMA analysis. Only for Affymetrix normalization.", false, false));
		getOptions().addOption(OptionFactory.createOption("plier", "Plier analysis. Only for Affymetrix normalization.", false, false));
		getOptions().addOption(OptionFactory.createOption("calls", "Present-absent calls analysis. Only for Affymetrix normalization.", false, false));


		// agilent or genepix normalization
		//
		getOptions().addOption(OptionFactory.createOption("sample-names", "Sample names. Only for agilent or genepix normalization.", false));
		getOptions().addOption(OptionFactory.createOption("bg-correction", "Background correction: none, normexp, half, subtract, minimum, movingmin, edwards, rma. Default background correction: none. Only for Agilent or GenePix normalization.", false));
		getOptions().addOption(OptionFactory.createOption("wa-normalization", "Within arrays normalization: none, loess, printtiploess, median, composite, control, robustspline. Default within arrays normalization: loess. Only for two-colors Agilent or GenePix normalization.", false));
		getOptions().addOption(OptionFactory.createOption("ba-normalization", "Between arrays normalization: none, quantiles, scale, vsn. Default between arrays normalization: scale. Only for Agilent or GenePix normalization.", false));
		getOptions().addOption(OptionFactory.createOption("flags-no-fitted", "If this option is set then spots will not be used in the fitting of the parameters of the normalization steps. Only for Agilent or GenePix normalization.", false, false));
		getOptions().addOption(OptionFactory.createOption("flags-as-missing", "If this option is set then spots will have a missing (NA) normalized value and A-value as well. Only for Agilent or GenePix normalization.", false, false));		
	}

	/**
	 * 
	 * @throws IOException
	 */
	public void prepare() throws IOException {
		outdir = outdir + "/";
		tmpDir = new File(outdir + "tmp");
		String compressedFileName = commandLine.getOptionValue("compressed-file", null);
		String rawDirName = commandLine.getOptionValue("raw-dir", null);
		List<String> tags = StringUtils.toList(commandLine.getOptionValue("compressed-file-tags", ""), ",");

		
		System.out.println("----------> " + ListUtils.toString(tags, ","));
		
		// checking analysis methods and input data (compressed file or directory name)
		//
		if ( compressedFileName == null && rawDirName == null ) {
			abort("missingdata_execute_expressionnormalization", "missing input data", "missing input data", "missing input data");						
		}

		if ( tags == null || tags.size() == 0 ) {
			abort("unidentifieddata_execute_expressionnormalization", "unidentified input data", "unidentified input data", "unidentified input data");						
		}


		//		String technology = commandLine.getOptionValue("technology", null);
		//		String channels = commandLine.getOptionValue("channels", null);
		//		int nbChannels = 0;

		if ( tags.contains("one-channel") ) {
			nbChannels = 1;
		} else if ( tags.contains("two-channels") ) {
			nbChannels = 2;
		}

		if ( tags.contains("agilent") ) {
			technology = "agilent";
		} else if ( tags.contains("genepix") ) {
			technology = "genepix";
		} else if ( tags.contains("affymetrix") ) {
			technology = "affy";
			nbChannels = 1;
		}

		System.out.println("nbChannels = " + nbChannels + ", technology = " + technology);

		if ( technology == null ) {
			abort("missingtechnology_execute_expressionnormalization", "missing technology tag", "missing technology tag", "missing technology tag");						
		}

		if ( !technology.equalsIgnoreCase("agilent") && !technology.equalsIgnoreCase("genepix") && !technology.equalsIgnoreCase("affy") ) {
			abort("invalidtechnology_execute_expressionnormalization", "invalid technology tag. Valid values are: agilent, genepix, affymetrix", "invalid technology type " + technology + ". Valid values are: agilent, genepix", "invalid technology type " + technology + ". Valid values are: agilent, genepix, affymetrix.");									
		}

		if ( nbChannels == 0 ) {
			abort("missingnbchannels_execute_expressionnormalization", "missing number of channels", "missing number of channels", "missing number of channels");						
		}


		jobStatus.addStatusMessage("10", "reading dataset");

		if ( compressedFileName != null ) {
			// getting raw files from compressed file
			//
			jobStatus.addStatusMessage("30", "decomprising dataset");

			System.out.println("input dataset = " + compressedFileName);
			System.out.println("tmp dir = " + tmpDir.getAbsolutePath());
			GenericCompressManager compresor = CompressFactory.getCompressManager(new File(compressedFileName));
			rawFileNames = compresor.decompress(compressedFileName, tmpDir.getAbsolutePath());
		} else {
			// getting raw files from directory
			//
			tmpDir = new File(rawDirName); 
		}
		
		File[] rawFiles = FileUtils.listFiles(tmpDir, "affy".equalsIgnoreCase(technology) ? ".+.cel" : ".+", true);
		rawFileNames = new ArrayList<String>(rawFiles.length);
		for(File file: rawFiles) {
			rawFileNames.add(file.getAbsolutePath());
		}

		// sanity check
		//
		if ( rawFileNames == null || rawFileNames.size() == 0 ) {
			abort("missingrawfiles_execute_expressionnormalization", "missing raw files", "missing raw files", "missing raw files");
		}		
	}

	/**
	 * 
	 */
	@Override
	protected void execute() {
		try {
			prepare();

			if ( "affy".equalsIgnoreCase(technology) ) {
				affyNormalization();
			} else {
				normalization();
			}
		} catch (Exception e) {
			printError("exception_execute_genepixexpression2cnormalization", e.toString(), e.getMessage(), e);
		}
	}

	private void normalization() throws IOException {
		
		String sampleNames = commandLine.getOptionValue("sample-names", null);
		String bgCorrection = commandLine.getOptionValue("bg-correction", "minimum");
		String waNormalization = commandLine.getOptionValue("wa-normalization", "median");
		String baNormalization = commandLine.getOptionValue("ba-normalization", "none");
		boolean flagsNotFitted = commandLine.hasOption("flags-not-fitted");
		boolean flagsAsMissing = commandLine.hasOption("flags-as-missing");

		String readingScript = System.getenv("BABELOMICS_HOME") + "/bin/normalizexp/" + (nbChannels == 1 ? "onecolor" : "twocolor");
		String normalizationScript = System.getenv("BABELOMICS_HOME") + "/bin/normalizexp/" + (nbChannels == 1 ? "onecolor" : "twocolor");

		if ( "genepix".equalsIgnoreCase(technology) ) {
			readingScript += "_genepix_reading.r";
			normalizationScript += "_genepix_normalizing.r";
		} else if ( "agilent".equalsIgnoreCase(technology) ) {
			readingScript += "_agilent_reading.r";
			normalizationScript += "_agilent_normalizing.r";			
		} else {
			abort("invalidtechnology_execute_expressionnormalization", "invalid technology tag. Valid values are: agilent, genepix, affymetrix", "invalid technology type " + technology + ". Valid values are: agilent, genepix", "invalid technology type " + technology + ". Valid values are: agilent, genepix, affymetrix");
		}

		// normalizing data
		//
		jobStatus.addStatusMessage("50", "normalizing data");

		if ( nbChannels == 1 ) {
			ExpressionUtils.OneColorNormalization(readingScript, normalizationScript, rawFileNames, (sampleNames != null ? StringUtils.toList(sampleNames, ","): getSamples(rawFileNames)), bgCorrection, baNormalization, flagsNotFitted, flagsAsMissing, outdir);			
		} else if ( nbChannels == 2 ){
			ExpressionUtils.TwoColorsNormalization(readingScript, normalizationScript, rawFileNames, (sampleNames != null ? StringUtils.toList(sampleNames, ","): getSamples(rawFileNames)), bgCorrection, waNormalization, baNormalization, flagsNotFitted, flagsAsMissing, outdir);			
		} else {
			abort("invalidchannels_execute_expressionnormalization", "invalid number of channels (" + nbChannels + ")", "invalid number of channels (" + nbChannels + ")", "invalid number of channels (" + nbChannels + ")");			
		}

		// saving normalization results
		//
		jobStatus.addStatusMessage("90", "saving normalization results");


		File file = null;
		if ( new File(outdir + "/" + ExpressionUtils.getNormalizedFileName()).exists() && 
			 new File(outdir + "/" + ExpressionUtils.getFeatureDataFileName()).exists() ) {

			file = new File(outdir + "/normalized_dataset.txt"); 			
			ExpressionUtils.createDataset(outdir + "/" + ExpressionUtils.getNormalizedFileName(), outdir + "/" + ExpressionUtils.getFeatureDataFileName(), 5, file.getAbsolutePath());

			if ( file.exists() ) {				
				String tags = "DATA,DATAMATRIX,EXPRESSION";
				File redirectionFile = new File(outdir + "/normalized.redirection");
				createPreprocessingRedirectionFile(redirectionFile, file);
				if ( redirectionFile.exists() ) {
					tags = tags + ",REDIRECTION(" + redirectionFile.getName() + ":Send to Preprocessing tool...)";
				}
				result.addOutputItem(new Item("normalized", file.getName(), "Normalized dataset ", TYPE.FILE, StringUtils.toList(tags, ","), new HashMap<String, String>(2), "Two-colors GenePix normalization files"));
			}

			file = new File(outdir + "/normalized_dataset.featdata"); 			
			if ( file.exists() ) {				
				result.addOutputItem(new Item("normalized", file.getName(), "Feature data ", TYPE.FILE, StringUtils.toList("idlist", ","), new HashMap<String, String>(2), "Two-colors GenePix normalization files"));
			}
		}


		file = new File(outdir + "/" + ExpressionUtils.getaValuesFileName()); 
		if ( file.exists() ) {
			result.addOutputItem(new Item("avalues", file.getName(), "A-values", TYPE.FILE, new ArrayList<String>(2), new HashMap<String, String>(2), "Two-colors GenePix normalization files"));
		}
	}


	/**
	 * 
	 * @throws IOException
	 * @throws InvalidParameterException
	 * @throws InvalidIndexException
	 */
	private void affyNormalization() throws IOException, InvalidParameterException, InvalidIndexException {

		String aptBinPath = System.getenv("BABELOMICS_HOME") + "/bin/apt";

		boolean celConvert = !commandLine.hasOption("no-cel-convert");
		boolean rma = commandLine.hasOption("rma");
		boolean plier = commandLine.hasOption("plier");
		boolean calls = commandLine.hasOption("calls");

		// checking analysis methods and input data (compressed file or directory name)
		//
		if ( !rma && !plier && !calls ) {
			abort("missinganalysis_execute_affynormalization", "missing analysis", "missing analysis, valid values are rma, plier and calls", "missing analysis, valid values are rma, plier and calls");			
		}

		// creating the cel_files containning the cel files to normalize
		//
		File celFiles = new File(outdir + "/cel_files.txt");
		IOUtils.write(celFiles, "cel_files\n" + ListUtils.toString(rawFileNames, "\n"));

		System.out.println("raw files = " + ListUtils.toString(rawFileNames, "\n"));
		
		// converting to CEL text
		//
		if ( celConvert ) {			
			jobStatus.addStatusMessage("40", "converting CEL to GCOS text file format");

			File convertDir = new File(tmpDir.getAbsolutePath() + "/convert");
			convertDir.mkdir();
			AffymetrixExpresionUtils.aptCelConvert(aptBinPath + "/apt-cel-convert", celFiles.getAbsolutePath(), convertDir.getAbsolutePath());
			File[] rawFiles = FileUtils.listFiles(convertDir, ".+.CEL", true);
			rawFileNames = ArrayUtils.toStringList(rawFiles);
			//System.out.println("-----------> converting to gcos text file format, raw file names = " + ListUtils.toString(rawFileNames, ","));
			IOUtils.write(celFiles, "cel_files\n" + ListUtils.toString(rawFileNames, "\n"));

//			AffymetrixExpresionUtils.aptCelConvert(aptBinPath + "/apt-cel-convert", celFiles.getAbsolutePath(), tmpDir.getAbsolutePath());
//
//			File[] rawFiles = FileUtils.listFiles(tmpDir, ".+.CEL", true);
//			rawFileNames = ArrayUtils.toStringList(rawFiles);
//			System.out.println("-----------> converting to gcos text file format, raw file names = " + ListUtils.toString(rawFileNames, ","));
//
//			IOUtils.write(celFiles, "cel_files\n" + ListUtils.toString(rawFileNames, "\n"));
		}		
		
		//Config config = new Config();
		config.append(new File(System.getenv("BABELOMICS_HOME") + "/conf/apt.conf"));
		String chipName = getChipName(rawFileNames, config.getKeys());

		if ( chipName == null ) {
			abort("exception_execute_affynormalization", "array type not supported", "array type not supported", "array type not supported");			
		}

		String infoStr = config.getProperty(chipName);
		Map<String, String> chipInfo = MapUtils.stringToMap(infoStr);
		chipInfo.put("name", chipName);

		System.out.println(" chip info = " + chipInfo.toString());

		
		String chipType = chipInfo.get("type");
		if ( chipType == null ) {
			abort("exception_execute_affynormalization", "could not find out the chip type", "could not find out the chip type", "could not find out the chip type");			
		}

		if ( !"3-prime".equalsIgnoreCase(chipType) && !"whole-transcript".equalsIgnoreCase(chipType) ) {
			abort("exception_execute_affynormalization", "array type (" + chipType + ") not supported", "array type (" + chipType + ") not supported", "array type (" + chipType + ") not supported");			
		}
		// normalizing data
		//
		jobStatus.addStatusMessage("50", "normalizing data");

		List<String> analysis = new ArrayList<String>();
		if ( "3-prime".equalsIgnoreCase(chipInfo.get("type")) ) {

			if ( rma ) analysis.add("rma");
			if ( calls ) analysis.add("pm-mm,mas5-detect.calls=1.pairs=1");
			if ( plier ) analysis.add("plier-mm");

			AffymetrixExpresionUtils.aptProbesetSummarize(aptBinPath + "/apt-probeset-summarize", analysis, config.getProperty("BABELOMICS_DATA_HOME") + chipInfo.get("cdf"), celFiles.getAbsolutePath(), outdir);

		} else {

			if ( rma ) analysis.add("rma");
			if ( calls ) analysis.add("dabg");
			if ( plier ) analysis.add("plier-gcbg");

			AffymetrixExpresionUtils.aptProbesetSummarize(aptBinPath + "/apt-probeset-summarize", analysis, config.getProperty("BABELOMICS_DATA_HOME") + chipInfo.get("cdf"), celFiles.getAbsolutePath(), outdir);
		}

		//		System.err.println("cmd output: " + sp.getRunnableProcess().getOutput());
		//		System.err.println("cmd error: " + sp.getRunnableProcess().getError());

		// saving normalization results
		//
		jobStatus.addStatusMessage("90", "saving normalization results");

		File file;
		List<String> tags = StringUtils.toList("data,datamatrix,expression", ",");

		file = new File(outdir + "/rma.summary.txt"); 
		if ( file.exists() ) {
			IOUtils.write(file, cleanLines(IOUtils.readLines(file)));
			saveAsDataset(file);

			result.addOutputItem(new Item("rma.summary", file.getName(), "Summary ", TYPE.FILE, tags, new HashMap<String, String>(2), "RMA analysis"));
			saveBoxPlot(file, "RMA box-plot", "rmaimg", "RMA analysis");				
		}

		file = new File(outdir + "/plier-mm.summary.txt"); 
		if ( file.exists() ) {
			IOUtils.write(file, cleanLines(IOUtils.readLines(file)));
			saveAsDataset(file);

			result.addOutputItem(new Item("plier-mm.summary", file.getName(), "MM summary ", TYPE.FILE, tags, new HashMap<String, String>(2), "Plier analysis"));								
			saveBoxPlot(file, "Plier MM box-plot", "plierimg", "Plier analysis");				
		}

		file = new File(outdir + "/plier-gcbg.summary.txt"); 
		if ( file.exists() ) {
			IOUtils.write(file, cleanLines(IOUtils.readLines(file)));
			saveAsDataset(file);

			result.addOutputItem(new Item("plier-gcbg.summary", file.getName(), "GCBG summary ", TYPE.FILE, tags, new HashMap<String, String>(2), "Plier analysis"));								
			saveBoxPlot(file, "Plier GCBG box-plot", "plierimg", "Plier analysis");				
		}

		file = new File(outdir + "/pm-mm.mas5-detect.summary.txt"); 
		if ( file.exists() ) {
			IOUtils.write(file, cleanLines(IOUtils.readLines(file)));
			saveAsDataset(file);

			result.addOutputItem(new Item("pm-mm.summary", file.getName(), "PM-MM summary ", TYPE.FILE, tags, new HashMap<String, String>(2), "Present-absent calls"));								
			saveBoxPlot(file, "PM-MM box-plot", "pmmmimg", "Present-absent calls");				
		}

		file = new File(outdir + "/pm-mm.mas5-detect.calls.txt"); 
		if ( file.exists() ) {
			IOUtils.write(file, cleanLines(IOUtils.readLines(file)));
			result.addOutputItem(new Item("pm-mm.calls", file.getName(), "Calls ", TYPE.FILE, new ArrayList<String>(2), new HashMap<String, String>(2), "Present-absent calls"));								
		}

		file = new File(outdir + "/dabg.summary.txt"); 
		if ( file.exists() ) {
			IOUtils.write(file, cleanLines(IOUtils.readLines(file)));
			saveAsDataset(file);

			result.addOutputItem(new Item("dabg.summary", file.getName(), "DABG summary ", TYPE.FILE, tags, new HashMap<String, String>(2), "Present-absent calls"));								
			saveBoxPlot(file, "DABG box-plot", "dabgimg", "Present-absent calls");				
		}
	}
	
	
	/**
	 * 
	 * @param rawFilenames
	 * @param chipNames
	 * @return
	 * @throws InvalidParameterException 
	 * @throws IOException 
	 * @throws Exception
	 */
	private String getChipName(List<String> rawFilenames, List<String> chipNames) throws InvalidParameterException, IOException {
		List<String> lines = null;
		List<String> results = new ArrayList<String>(rawFilenames.size());

		//System.out.println("----------> getChipName, chip names = " + ListUtils.toString(chipNames, ", "));
		
		for(int i=0 ; i<rawFilenames.size() ; i++) {
			
			//System.out.println("----------> getChipName, file name = " + rawFilenames.get(i));
			
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

		//System.out.println("***** results = " + ListUtils.toString(results));
		String chipName = results.get(0);
		for (int i=0 ; i<rawFilenames.size() ; i++) {
			if ( ! chipName.equalsIgnoreCase(results.get(i)) ) {
				String msg = "mismatch CEL files corresponding to different chips.\n";
				for(int j=0 ; j<rawFilenames.size() ; j++) {
					msg = msg + "Raw file '" + new File(rawFilenames.get(j)).getName() + "' is a '" + results.get(j) + "' array\n";
				}
				throw new InvalidParameterException(msg);
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
	 * @throws InvalidIndexException
	 */
	private void saveAsDataset(File file) throws IOException, InvalidIndexException {
		Dataset dataset = new Dataset(file);
		dataset.load();
		dataset.save();
	}

	public void saveBoxPlot(File file, String title, String resultId, String group) throws IOException, InvalidIndexException {
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

	/**
	 * 
	 * @param rawFileNames
	 * @return
	 */
	private List<String> getSamples(List<String> rawFileNames) {
		int index;
		File file;
		List<String> samples = new ArrayList<String>(rawFileNames.size());
		for(int i=0 ; i<rawFileNames.size() ; i++) {
			file = new File(rawFileNames.get(i));
			index = file.getName().lastIndexOf('.');
			if ( (index > 0) && (index <= file.getName().length() - 2) ) {
				samples.add(file.getName().substring(0, index));
			} else {
				samples.add("Sample_" + i);
			}
		}
		return samples;
	}

	public void createPreprocessingRedirectionFile(File redirectionFile, File fileToRedirect) {
		List<String> redirectionInputs = new ArrayList<String>();
		redirectionInputs.add("tool=preprocessing");
		redirectionInputs.add("jobname=preprocessing");
		redirectionInputs.add("jobdescription=redirected from job $JOB_NAME");
		redirectionInputs.add("dataset_databox=" + fileToRedirect.getName() + " (from job $JOB_NAME)");
		redirectionInputs.add("dataset=$JOB_FOLDER/" + fileToRedirect.getName());
		redirectionInputs.add("dataset_wum_data=true");
		redirectionInputs.add("merge_replicates=mean");
		redirectionInputs.add("impute_missing=mean");
		try {
			IOUtils.write(redirectionFile.getAbsolutePath(), redirectionInputs);
		} catch (IOException e) {
			e.printStackTrace();
		}
	}	
}
