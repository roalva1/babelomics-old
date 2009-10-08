package org.bioinfo.babelomics.methods.genomic.copynumber;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.bioinfo.babelomics.exception.InvalidParameterException;
import org.bioinfo.commons.exec.Command;
import org.bioinfo.commons.exec.SingleProcess;

public class CopyNumberAnalysisExecutor {

	private File normalizedFile;
	private String segmentationBinPath;
	private String cghMcrBinPath;
	private String outdir;
	private int gapAllowed = 500;
	private double alteredLow = 0.20;
	private double alteredHigh = 0.80;
	private int recurrence = 50;
	

	public CopyNumberAnalysisExecutor (String normalizedFile) {
		this.normalizedFile = new File(normalizedFile);
	}
	
	public CopyNumberAnalysisExecutor (File normalizedFile) {
		this.normalizedFile = normalizedFile;
	}
	
	
	public void run() throws InvalidParameterException {
		
		if ( segmentationBinPath == null ) {
			throw new InvalidParameterException("copy number binary path missing");
		}
		
		if ( outdir == null ) {
			throw new InvalidParameterException("copy number out directory path missing");
		}
		
		List<String> env = new ArrayList<String>();
		env.add("infile=" + normalizedFile.getAbsolutePath());		
		env.add("outdir=" + outdir);
	
		Command cmd = new Command("R CMD BATCH --no-save --no-restore " + segmentationBinPath + " " + outdir + "/r.segmentation.log", env);
		
		System.out.println("cmd = " + cmd.getCommandLine());
		
		SingleProcess sp = new SingleProcess(cmd);
		sp.runSync();		

		if ( cghMcrBinPath != null ) {		
			File segmentedFile = new File(outdir + "/segmented.txt");
			if ( segmentedFile.exists() ) {
				env = new ArrayList<String>();
				env.add("normalized.file=" + normalizedFile.getAbsolutePath());		
				env.add("segmented.file=" + segmentedFile.getAbsolutePath());		
				env.add("outdir=" + outdir);

				env.add("gapAllowed=" + gapAllowed);
				env.add("alteredLow=" + alteredLow);
				env.add("alteredHigh=" + alteredHigh);
				env.add("recurrence=" + recurrence);

				cmd = new Command("R CMD BATCH --no-save --no-restore " + cghMcrBinPath + " " + outdir + "/r.cgh.log", env);
				
				System.out.println("cmd = " + cmd.getCommandLine());
				
				sp = new SingleProcess(cmd);
				sp.runSync();
			} else {
				System.err.println("erro creating segmented file");
			}
		}
	}

	public String getOutdir() {
		return outdir;
	}

	public void setOutdir(String outdir) {
		this.outdir = outdir;
	}

	public int getGapAllowed() {
		return gapAllowed;
	}

	public void setGapAllowed(int gapAllowed) {
		this.gapAllowed = gapAllowed;
	}

	public double getAlteredLow() {
		return alteredLow;
	}

	public void setAlteredLow(double alteredLow) {
		this.alteredLow = alteredLow;
	}

	public double getAlteredHigh() {
		return alteredHigh;
	}

	public void setAlteredHigh(double alteredHigh) {
		this.alteredHigh = alteredHigh;
	}

	public int getRecurrence() {
		return recurrence;
	}

	public void setRecurrence(int recurrence) {
		this.recurrence = recurrence;
	}

	public void setNormalizedFile(File normalizedFile) {
		this.normalizedFile = normalizedFile;
	}

	public File getNormalizedFile() {
		return normalizedFile;
	}

	public void setCghMcrBinPath(String cghMcrBinPath) {
		this.cghMcrBinPath = cghMcrBinPath;
	}

	public String getCghMcrBinPath() {
		return cghMcrBinPath;
	}
	public String getSegmentationBinPath() {
		return segmentationBinPath;
	}

	public void setSegmentationBinPath(String segmentationBinPath) {
		this.segmentationBinPath = segmentationBinPath;
	}
}
