package org.bioinfo.babelomics.methods.genomic.copynumber;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.bioinfo.babelomics.exception.InvalidParameterException;
import org.bioinfo.commons.exec.Command;
import org.bioinfo.commons.exec.SingleProcess;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.data.dataset.Dataset;

public class CopyNumberAnalysisExecutor {
	
	private String normalizedFilename;
	private String segmentedFilename;
	private String cghFilename;
	private String segmentationBinPath;
	private String cghMcrBinPath;
	private int gapAllowed = 500;
	private double alteredLow = 0.20;
	private double alteredHigh = 0.80;
	private int recurrence = 50;
	

	public CopyNumberAnalysisExecutor (String segmentationBinPath) {
		this.segmentationBinPath = segmentationBinPath;
	}	
	
	public void run() throws InvalidParameterException, IOException {
		String normNames = null;
		String segNames = null;
		
		if ( segmentationBinPath == null ) {
			throw new InvalidParameterException("copy number binary path missing");
		}
		
		if ( normalizedFilename == null ) {
			throw new InvalidParameterException("copy number input normalized file name missing");
		}

		if ( segmentedFilename == null ) {
			throw new InvalidParameterException("copy number output segmented file name missing");
		}
		
		if ( cghMcrBinPath != null && cghFilename == null ) {
			throw new InvalidParameterException("copy number output cgh file name missing");			
		}
		
		List<String> env = new ArrayList<String>();
		env.add("infile=" + normalizedFilename);		
		env.add("outfile=" + segmentedFilename);
	
		normNames = getColumnNames(normalizedFilename);
		env.add("colnames=" + normNames);
						
		Command cmd = new Command("R CMD BATCH --no-save --no-restore " + segmentationBinPath + " " + segmentedFilename + ".log", env);
		
		System.out.println("cmd = " + cmd.getCommandLine());
		System.out.println("env = " + ListUtils.toString(env, " "));
		
		SingleProcess sp = new SingleProcess(cmd);
		sp.runSync();		

		if ( cghMcrBinPath != null ) {		
			File segmentedFile = new File(segmentedFilename);
			if ( segmentedFile.exists() ) {
				env = new ArrayList<String>();
				env.add("normalizedfile=" + normalizedFilename);		
				env.add("normnames=" + normNames);
				env.add("segmentedfile=" + segmentedFilename);		
				env.add("segnames=" + getColumnNames(segmentedFilename));
				env.add("outfile=" + cghFilename);

				env.add("gapAllowed=" + gapAllowed);
				env.add("alteredLow=" + alteredLow);
				env.add("alteredHigh=" + alteredHigh);
				env.add("recurrence=" + recurrence);

				cmd = new Command("R CMD BATCH --no-save --no-restore " + cghMcrBinPath + " " + cghFilename + ".log", env);
				
				System.out.println("cmd = " + cmd.getCommandLine());
				System.out.println("env = " + ListUtils.toString(env, " "));
				
				sp = new SingleProcess(cmd);
				sp.runSync();
				
				File cghFile = new File(cghFilename);
				if ( ! cghFile.exists() ) {
					System.err.println("error creating cgh file");
				}
			} else {
				System.err.println("error creating segmented file");
			}
		}
	}

	
	private String getColumnNames(String filename) throws IOException {
		String names = null;
		List<String> lines = IOUtils.grep(filename, "#NAMES.*");
		if ( lines != null ) {
			names = lines.get(0).replace("#", "");
			names = names.replace("\t", ",");
		} else {
			throw new IOException("File " + filename + " does not contain #NAMES label");
		}
		return names;
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

	public void setNormalizedFile(String normalizedFilename) {
		this.normalizedFilename = normalizedFilename;
	}

	public String getNormalizedFilename() {
		return normalizedFilename;
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

	public String getSegmentatedFilename() {
		return segmentedFilename;
	}
	
	public void setSegmentatedFilename(String segmentedFilename) {
		this.segmentedFilename = segmentedFilename;
	}

	public String getCghFilename() {
		return cghFilename;
	}
	
	public void setCghFilename(String cghFilename) {
		this.cghFilename = cghFilename;
	}
}
