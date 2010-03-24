package org.bioinfo.babelomics.tools.expression;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.commons.io.utils.FileUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.io.file.compress.CompressFactory;
import org.bioinfo.io.file.compress.GenericCompressManager;
import org.bioinfo.microarray.AffymetrixExpressionUtils;
import org.bioinfo.microarray.AgilentExpressionUtils;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;

public class RawExpressionViewer extends BabelomicsTool {

	public RawExpressionViewer() {
		initOptions();
	}

	@Override
	public void initOptions() {

		getOptions().addOption(OptionFactory.createOption("compressed-file", "Compressed file containning the raw files", false));
		getOptions().addOption(OptionFactory.createOption("raw-dir", "Directory where the raw files are located", false));
		getOptions().addOption(OptionFactory.createOption("compressed-file-tags", "Data tags separated by commas, valid values are: microarray, expression, affymetrix, agilent, genepix, one-channel, two-channels", false));
		//		getOptions().addOption(OptionFactory.createOption("technology", "Techonology type: agilent, genepix", false));
		//		getOptions().addOption(OptionFactory.createOption("channels", "Number of channels: 1 or 2", false));

		getOptions().addOption(OptionFactory.createOption("sample-filter", "class variable", false));
		getOptions().addOption(OptionFactory.createOption("feature-filter", "class variable", false));
	}

	@Override
	public void execute() {
		outdir = outdir + "/";
		File tmpDir = new File(outdir + "tmp");
		String compressedFileName = commandLine.getOptionValue("compressed-file", null);
		String rawDirName = commandLine.getOptionValue("raw-dir", null);
		List<String> tags = StringUtils.toList(commandLine.getOptionValue("compressed-file-tags", ""), ",");

		// checking analysis methods and input data (compressed file or directory name)
		//
		if ( compressedFileName == null && rawDirName == null ) {
			abort("missingdata_execute_rawexpressionviewer", "missing input data", "missing input data", "missing input data");						
		}

		if ( tags == null || tags.size() == 0 ) {
			abort("unidentifieddata_execute_rawexpressionviewer", "unidentified input data", "unidentified input data", "unidentified input data");						
		}


		//		String technology = commandLine.getOptionValue("technology", null);
		//		String channels = commandLine.getOptionValue("channels", null);
		//		int nbChannels = 0;

		int nbChannels = 0;
		if ( tags.contains("one-channel") ) {
			nbChannels = 1;
		} else if ( tags.contains("two-channels") ) {
			nbChannels = 2;
		}

		String technology = null;
		if ( tags.contains("agilent") ) {
			technology = "agilent";
		} else if ( tags.contains("genepix") ) {
			technology = "genepix";
		} else if ( tags.contains("affymetrix") ) {
			technology = "affy";
			nbChannels = 1;
		}


		if ( technology == null ) {
			abort("missingtechnology_execute_rawexpressionviewer", "missing technology typee", "missing technology type", "missing technology type");						
		}

		if ( !technology.equalsIgnoreCase("agilent") && !technology.equalsIgnoreCase("genepix") && !technology.equalsIgnoreCase("affy") ) {
			abort("invalidtechnology_execute_rawexpressionviewer", "invalid technology type " + technology + ". Valid values are: agilent, genepix, affymetrix", "invalid technology type " + technology + ". Valid values are: agilent, genepix", "invalid technology type " + technology + ". Valid values are: agilent, genepix");									
		}

		if ( nbChannels == 0 ) {
			abort("missingdata_execute_rawexpressionviewer", "missing number of channels", "missing number of channels", "missing number of channels");						
		}

		//		if ( channels == null ) {
		//			abort("missingdata_execute_rawexpressionviewer", "missing number of channels", "missing number of channels", "missing number of channles");						
		//		}
		//		
		//		try {
		//			nbChannels = Integer.parseInt(channels);
		//			if ( nbChannels <= 0 || nbChannels >=3 ) {
		//				abort("invalidnbchannels_execute_rawexpressionviewer", "invalid number of channels " + channels + ". Valid values are: 1, 2", "invalid number of channels " + channels + ". Valid values are: 1, 2", "invalid number of channels " + channels + ". Valid values are: 1, 2");													
		//			}
		//		} catch (Exception e) {
		//			abort("invalidnbchannels_execute_rawexpressionviewer", "invalid number of channels " + channels + ". Valid values are: 1, 2", "invalid number of channels " + channels + ". Valid values are: 1, 2", "invalid number of channels " + channels + ". Valid values are: 1, 2");									
		//		}


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
				File[] rawFiles = FileUtils.listFiles(new File(rawDirName), tags.contains("affymetrix") ? ".+.CEL" : ".+");
				rawFileNames = new ArrayList<String>(rawFiles.length);
				for(File file: rawFiles) {
					rawFileNames.add(file.getAbsolutePath());
				}
			}

			// sanity check
			//
			if ( rawFileNames == null || rawFileNames.size() == 0 ) {
				abort("missingrawfiles_execute_rawexpressionviewer", "missing raw files", "missing raw files", "missing raw files");
			}

			System.out.println("raw files = " + ListUtils.toString(rawFileNames, "\n"));

			// normalizing data
			//
			jobStatus.addStatusMessage("50", "normalizing data");

			String readingBinPath, getRawValuesBinPath, intensityPlotBinPath, maPlotBinPath;
			if ( "affy".equalsIgnoreCase(technology) ) {
				readingBinPath = System.getenv("BABELOMICS_HOME") + "/bin/normalizexp/" + technology.toLowerCase() + "_reading.r";
				getRawValuesBinPath = System.getenv("BABELOMICS_HOME") + "/bin/normalizexp/" + technology.toLowerCase() + "_get_raw_values.r";
			} else {
				readingBinPath = System.getenv("BABELOMICS_HOME") + "/bin/normalizexp/" + (nbChannels == 1 ? "onecolor" : "twocolor") + "_" + technology.toLowerCase() + "_reading.r";
				getRawValuesBinPath = System.getenv("BABELOMICS_HOME") + "/bin/normalizexp/" + (nbChannels == 1 ? "onecolor" : "twocolor") + "_" + technology.toLowerCase() + "_get_raw_values.r";
			}
			intensityPlotBinPath = System.getenv("BABELOMICS_HOME") + "/bin/plots/plot_image_" + technology.toLowerCase() + ".r";
			maPlotBinPath = System.getenv("BABELOMICS_HOME") + "/bin/plots/plotMA_from_" + (nbChannels == 1 ? "single_matrix" : "MA_matrices") + ".r";

			List<String> values = (nbChannels == 1) ? StringUtils.toList("F,B,features", ",") : StringUtils.toList("M,A,G,R,features", ",");

			File rawDataRObjectFile = new File(outdir + "/raw_data.RData");		
			File featuresFile = new File(outdir + "/features.txt");		
			File fFile = new File(outdir + "/F.txt");		
			File bFile = new File(outdir + "/B.txt");		
			File mFile = new File(outdir + "/M.txt");		
			File aFile = new File(outdir + "/A.txt");		
			File gFile = new File(outdir + "/G.txt");		
			File rFile = new File(outdir + "/R.txt");		

			// indir=input/ files=array1.gpr,array2.gpr samplenames=array1,array2 outfile=genepix_onecolor_rawdata.RData R CMD BATCH --no-save --no-restore ~/appl/bioinfo-installer/babelomics/bin/normalizexp/onecolor_genepix_reading.r reading.Rout
			//
			AgilentExpressionUtils.createRawDataRObject(readingBinPath, rawFileNames, getSamples(rawFileNames), rawDataRObjectFile.getAbsolutePath(), "raw_data.Rout", outdir);			
			if ( ! rawDataRObjectFile.exists() ) {
				abort("errorcreatingreadingrobject_execute_rawexpressionviewer", "impossible to create R object from your raw files", "impossible to create R object from your raw files", "impossible to create R object from your raw files");				
			}

			if ( "affy".equalsIgnoreCase(technology) ) {
				AffymetrixExpressionUtils.getRawValues(getRawValuesBinPath, rawDataRObjectFile.getAbsolutePath(), fFile.getName(), "get_raw_values.Rout", outdir);
			} else {
				AgilentExpressionUtils.getRawValues(getRawValuesBinPath, rawDataRObjectFile.getAbsolutePath(), values, "get_raw_values.Rout", outdir);
			}

			// intensity plots
			//
			if ( "affy".equalsIgnoreCase(technology) ) {
				AffymetrixExpressionUtils.createIntensityPlot(intensityPlotBinPath, rawDataRObjectFile.getAbsolutePath(), "intensity_", "intensity_plot.Rout", outdir);
				addOutputItems(outdir, "intensity_", "png", "intensity_image", "Intensity image");					
			} else {
				if ( featuresFile.exists() ) {
					if ( fFile.exists() ) {
						AgilentExpressionUtils.createIntensityPlot(intensityPlotBinPath, fFile.getAbsolutePath(), featuresFile.getAbsolutePath(), "fg_", "fg_intensity_plot.Rout", outdir);
						addOutputItems(outdir, "fg_", "png", "fg_intensity_image", "Foreground intensity image");					
					}

					if ( bFile.exists() ) {
						AgilentExpressionUtils.createIntensityPlot(intensityPlotBinPath, bFile.getAbsolutePath(), featuresFile.getAbsolutePath(), "bg_", "bg_intensity_plot.Rout", outdir);
						addOutputItems(outdir, "bg_", "png", "bg_intensity_image", "Background intensity image");					
					}

					if ( mFile.exists() ) {
						AgilentExpressionUtils.createIntensityPlot(intensityPlotBinPath, mFile.getAbsolutePath(), featuresFile.getAbsolutePath(), "m_", "m_intensity_plot.Rout", outdir);
						addOutputItems(outdir, "m_", "png", "m_intensity_image", "M-values intensity image");					
					}

					if ( gFile.exists() ) {
						AgilentExpressionUtils.createIntensityPlot(intensityPlotBinPath, gFile.getAbsolutePath(), featuresFile.getAbsolutePath(), "g_", "g_intensity_plot.Rout", outdir);
						addOutputItems(outdir, "g_", "png", "g_intensity_image", "Green-values intensity image");					
					}

					if ( rFile.exists() ) {
						AgilentExpressionUtils.createIntensityPlot(intensityPlotBinPath, gFile.getAbsolutePath(), featuresFile.getAbsolutePath(), "r_", "r_intensity_plot.Rout", outdir);
						addOutputItems(outdir, "r_", "png", "r_intensity_image", "Red-values intensity image");					
					}
				}
			}

			// ma plots
			//
			if ( nbChannels == 1 ) {
				if ( fFile.exists() ) {
					AgilentExpressionUtils.createMAPlot(maPlotBinPath, fFile.getAbsolutePath(), "MA_", "ma_plot.Rout", outdir);
				}
			} else {
				if ( mFile.exists() && aFile.exists() ) {
					AgilentExpressionUtils.createMAPlot(maPlotBinPath, mFile.getAbsolutePath(), aFile.getAbsolutePath(), "MA_", "ma_plot.Rout", outdir);
				}
			}
			if ( "affy".equalsIgnoreCase(technology) ) {
				addOutputItems(outdir, "MA_", "", "ma_plot", "MA plot");				
			} else {
				addOutputItems(outdir, "MA_", "png", "ma_plot", "MA plot");				
			}

		} catch (FileNotFoundException e) {
			printError("filenotfoundexception_execute_rawexpressionviewer", e.toString(), e.getMessage(), e);
		} catch (IOException e) {
			printError("ioexception_execute_rawexpressionviewer", e.toString(), e.getMessage(), e);
		}

	}


	private void addOutputItems(String dir, String prefix, String suffix, String id, String groupName) {
		String sampleName;
		String pattern = prefix + ".*" + suffix;
		File [] files = FileUtils.listFiles(new File(dir), pattern);

		System.out.println("addOutputItems...");
		for (int i=0 ; i<files.length ; i++) {
			System.out.println("file " + i + " -> " + files[i].getName());
			
			sampleName = files[i].getName().replace(prefix, "").replace(".png", "");
			result.addOutputItem(new Item(id + "_" + i, files[i].getName(), groupName + " for " + sampleName, TYPE.IMAGE, new ArrayList<String>(2), new HashMap<String, String>(2), groupName + "s"));
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
