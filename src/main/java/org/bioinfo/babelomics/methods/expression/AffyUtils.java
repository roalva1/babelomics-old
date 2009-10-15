package org.bioinfo.babelomics.methods.expression;

import java.io.File;
import java.util.List;

import org.bioinfo.commons.exec.Command;
import org.bioinfo.commons.exec.SingleProcess;

public class AffyUtils {

	/*
	 * 
	 */
	public static void aptCelConvert(String binPath, String dirname, File celFiles) {
		StringBuilder aptCmd =  new StringBuilder(binPath);
		aptCmd.append(" -f text -o ").append(dirname).append(" --cel-files ").append(celFiles.getAbsolutePath());
		System.out.println("apt cmd = " + aptCmd.toString());
		
		Command cmd = new Command(aptCmd.toString());
		SingleProcess sp = new SingleProcess(cmd);
		sp.runSync();
	}
	
	/*
	 * 
	 */
	public static void aptProbesetSummarize(String binPath, String outdir, List<String> analysis, File cdfFile, File celFiles) {
		StringBuilder aptCmd = new StringBuilder(binPath);
		aptCmd.append(" -o ").append(outdir).append(" -d ").append(cdfFile.getAbsolutePath());
		for (String a: analysis) {
			aptCmd.append(" -a ").append(a);
		}
		aptCmd.append(" --cel-files ").append(celFiles.getAbsolutePath());
		//cmdStr = cmdStr + " " + tmpDirname + "/*.CEL";
		System.out.println("apt cmd = " + aptCmd.toString());
		
		Command cmd = new Command(aptCmd.toString());
		SingleProcess sp = new SingleProcess(cmd);
		sp.runSync();
	}

}
