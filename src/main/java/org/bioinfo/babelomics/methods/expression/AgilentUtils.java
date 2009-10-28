package org.bioinfo.babelomics.methods.expression;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.bioinfo.commons.exec.Command;
import org.bioinfo.commons.exec.SingleProcess;
import org.bioinfo.commons.utils.ListUtils;

public class AgilentUtils {

	/*
	 * 
	 */
	private static String normalizedFileName = "normalized.txt";
	private static String aValuesFileName = "Avalues.txt";
	private static String geneIDsFileName = "geneIDs.txt";
	
	public static void TwoColorsNormalization(String readingBinPath, String normalizingBinPath, List<String> fileNames, List<String> sampleNames, String bgCorrection, String waNormalization, String baNormalization, boolean flagsNotFitted, boolean flagsAsMissing, String outDirName) {
		
		File rFile = new File(outDirName + "/agilent_rawdata.RData");
		List<String> env = new ArrayList<String>();
		env.add("files=" + ListUtils.toString(fileNames, ","));
		if ( sampleNames != null && sampleNames.size() > 0 ) {
			env.add("samplenames=" + ListUtils.toString(sampleNames, ","));
		}
		env.add("outfile=" + rFile.getAbsolutePath());
	
		Command cmd = new Command("R CMD BATCH --no-save --no-restore " + readingBinPath + " " + outDirName + "/twocolor_agilent_reading.Rout", env);
		System.out.println("cmd line = " + cmd.getCommandLine());
		System.out.println("cmd env = " + ListUtils.toString(env, " "));		
		SingleProcess sp = new SingleProcess(cmd);
		sp.runSync();

		if ( rFile.exists() ) {
			env = new ArrayList<String>();
			
			env.add("infile=" + rFile.getAbsolutePath());
			env.add("normdatafile=" + outDirName + "/" + normalizedFileName);
			env.add("Avaluesfile=" + outDirName + "/" + aValuesFileName);
			env.add("featuresfile=" + outDirName + "/" +  geneIDsFileName);
			env.add("BGmethod=" + bgCorrection);
			env.add("WAmethod=" + waNormalization);
			env.add("BAmethod=" + baNormalization);
			env.add("flagsnotfitted=" + (flagsNotFitted ? "TRUE" : "FALSE"));
			env.add("flagsasmissing=" + (flagsAsMissing ? "TRUE" : "FALSE"));
						
			cmd = new Command("R CMD BATCH --no-save --no-restore " +  normalizingBinPath + " " + outDirName + "/twocolor_agilent_normalizing.Rout", env);			
			System.out.println("cmd line = " + cmd.getCommandLine());
			System.out.println("cmd env = " + ListUtils.toString(env, " "));		
			sp = new SingleProcess(cmd);
			sp.runSync();			
		}
	}
	
	public static String getNormalizedFileName() {
		return normalizedFileName;
	}
	
	public static String getaValuesFileName() {
		return aValuesFileName;
	}

	public static String getGeneIDsFileName() {
		return geneIDsFileName;
	}
}
