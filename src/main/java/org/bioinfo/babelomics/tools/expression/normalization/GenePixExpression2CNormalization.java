package org.bioinfo.babelomics.tools.expression.normalization;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.commons.io.utils.FileUtils;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.data.dataset.FeatureData;
import org.bioinfo.io.file.compress.CompressFactory;
import org.bioinfo.io.file.compress.GenericCompressManager;
import org.bioinfo.microarray.AgilentExpressionUtils;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;


public class GenePixExpression2CNormalization extends BabelomicsTool {

	public GenePixExpression2CNormalization() {
		initOptions();
	}

	@Override
	public void initOptions() {

		getOptions().addOption(OptionFactory.createOption("compressed-file", "Compressed file containning the raw files", false));
		getOptions().addOption(OptionFactory.createOption("raw-dir", "Directory where the raw files are located", false));

		getOptions().addOption(OptionFactory.createOption("sample-names", "Sample names", false));
		getOptions().addOption(OptionFactory.createOption("bg-correction", "Background correction: none, normexp, half, subtract, minimum, movingmin, edwards, rma. Default background correction: none", false));
		getOptions().addOption(OptionFactory.createOption("wa-normalization", "Within arrays normalization: none, loess, printtiploess, median, composite, control, robustspline. Default within arrays normalization: loess", false));
		getOptions().addOption(OptionFactory.createOption("ba-normalization", "Between arrays normalization: none, quantiles, scale, vsn. Default between arrays normalization: scale", false));
		getOptions().addOption(OptionFactory.createOption("flags-no-fitted", "If this option is set then spots will not be used in the fitting of the parameters of the normalization steps", false, false));
		getOptions().addOption(OptionFactory.createOption("flags-as-missing", "If this option is set then spots will have a missing (NA) normalized value and A-value as well", false, false));

		getOptions().addOption(OptionFactory.createOption("sample-filter", "class variable", false));
		getOptions().addOption(OptionFactory.createOption("feature-filter", "class variable", false));
	}

	@Override
	public void execute() {
		File tmpDir = new File(outdir + "/tmp");
		String compressedFileName = commandLine.getOptionValue("compressed-file", null);
		String rawDirName = commandLine.getOptionValue("raw-dir", null);
		String sampleNames = commandLine.getOptionValue("sample-names", null);
		String bgCorrection = commandLine.getOptionValue("bg-correction", "minimum");
		String waNormalization = commandLine.getOptionValue("wa-normalization", "median");
		String baNormalization = commandLine.getOptionValue("ba-normalization", "none");
		boolean flagsNotFitted = commandLine.hasOption("flags-not-fitted");
		boolean flagsAsMissing = commandLine.hasOption("flags-as-missing");

		// checking analysis methods and input data (compressed file or directory name)
		//
		if ( compressedFileName == null && rawDirName == null ) {
			abort("missingdata_execute_genepixexpression2cnormalization", "missing input data", "missing input data", "missing input data");						
		}

		try {
			jobStatus.addStatusMessage("10", "reading dataset");

			List<String> rawFileNames = null;
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
				File[] rawFiles = FileUtils.listFiles(new File(rawDirName), ".+");
				rawFileNames = new ArrayList<String>(rawFiles.length);
				for(File file: rawFiles) {
					rawFileNames.add(file.getAbsolutePath());
				}
			}

			// sanity check
			//
			if ( rawFileNames == null || rawFileNames.size() == 0 ) {
				abort("missingrawfiles_execute_genepixexpression2cnormalization", "missing raw files", "missing raw files", "missing raw files");
			}

			System.out.println("raw files = " + ListUtils.toString(rawFileNames, "\n"));

			// normalizing data
			//
			jobStatus.addStatusMessage("50", "normalizing data");

			AgilentExpressionUtils.TwoColorsNormalization(System.getenv("BABELOMICS_HOME") + "/bin/normalizexp/twocolor_genepix_reading.r", System.getenv("BABELOMICS_HOME") + "/bin/normalizexp/twocolor_genepix_normalizing.r", rawFileNames, (sampleNames != null ? StringUtils.toList(sampleNames, ","): getSamples(rawFileNames)), bgCorrection, waNormalization, baNormalization, flagsNotFitted, flagsAsMissing, outdir);

			// saving normalization results
			//
			jobStatus.addStatusMessage("90", "saving normalization results");

			File file;

//			file = new File(outdir + "/" + AgilentExpressionUtils.getNormalizedFileName()); 
//			if ( file.exists() ) {
//				result.addOutputItem(new Item("normalized", file.getName(), "Two-colors GenePix normalization ", TYPE.FILE, tags, new HashMap<String, String>(2), "Two-colors GenePix normalization files"));
//			} else {
//				printError("error two-colors genepix normalization", "error two-colors genepix normalization", "error two-colors genepix normalization");
//			}			
//			
//			file = new File(outdir + "/" + AgilentExpressionUtils.getFeatureDataFileName()); 
//			if ( file.exists() ) {
//				result.addOutputItem(new Item("normalized", file.getName(), "Feature data", TYPE.FILE, tags, new HashMap<String, String>(2), "Two-colors GenePix normalization files"));
//			}

			if ( new File(outdir + "/" + AgilentExpressionUtils.getNormalizedFileName()).exists() && 
				 new File(outdir + "/" + AgilentExpressionUtils.getFeatureDataFileName()).exists() ) {
				
				file = new File(outdir + "/normalized_dataset.txt"); 			
				AgilentExpressionUtils.createDataset(outdir + "/" + AgilentExpressionUtils.getNormalizedFileName(), outdir + "/" + AgilentExpressionUtils.getFeatureDataFileName(), 5, file.getAbsolutePath());
				
				if ( file.exists() ) {				
					String tags = "data,datamatrix,expression";
					File redirectionFile = new File(outdir + "/normalized.redirection");
					NormalizationUtils.createPreprocessingRedirectionFile(redirectionFile, file);
					if ( redirectionFile.exists() ) {
						tags = tags + ",redirection(" + redirectionFile.getName() + ":Send to Preprocessing tool...)";
					}
					result.addOutputItem(new Item("normalized", file.getName(), "Normalized dataset ", TYPE.FILE, StringUtils.toList(tags, ","), new HashMap<String, String>(2), "Two-colors GenePix normalization files"));
				}
				
				file = new File(outdir + "/normalized_dataset.featdata"); 			
				if ( file.exists() ) {				
					result.addOutputItem(new Item("normalized", file.getName(), "Feature data ", TYPE.FILE, StringUtils.toList("idlist", ","), new HashMap<String, String>(2), "Two-colors GenePix normalization files"));
				}
			}
			
			
			file = new File(outdir + "/" + AgilentExpressionUtils.getaValuesFileName()); 
			if ( file.exists() ) {
				result.addOutputItem(new Item("avalues", file.getName(), "A-values", TYPE.FILE, new ArrayList<String>(2), new HashMap<String, String>(2), "Two-colors GenePix normalization files"));
			}
			
		} catch (FileNotFoundException e) {
			printError("filenotfoundexception_execute_genepixexpression2cnormalization", e.toString(), e.getMessage(), e);
		} catch (IOException e) {
			printError("ioexception_execute_genepixexpression2cnormalization", e.toString(), e.getMessage(), e);
		}

	}


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

}
