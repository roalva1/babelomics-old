package org.bioinfo.babelomics.methods.genomic.genotype;

import java.io.File;
import java.io.IOException;
import java.security.InvalidParameterException;

import org.bioinfo.collections.exceptions.InvalidColumnIndexException;
import org.bioinfo.commons.exec.Command;
import org.bioinfo.commons.exec.SingleProcess;
import org.bioinfo.commons.io.utils.FileUtils;
import org.bioinfo.data.dataset.FeatureData;

public class GenotypeAnalysis {

	private File pedFile;
	private File mapFile;
	private String plinkPath;
	private String outdir;

	public GenotypeAnalysis() {
		this.pedFile = null;
		this.mapFile = null;
	}

	public GenotypeAnalysis(String pedFilePath, String mapFilePath) {
		this.pedFile = new File(pedFilePath);
		this.mapFile = new File(mapFilePath);
	}

	public GenotypeAnalysis(File pedFile, File mapFile) {
		this.pedFile = pedFile;
		this.mapFile = mapFile;
	}

	/*
	 * 
	 * PLINK TESTS
	 * 
	 */
	public void association(String assocTest) throws IOException {
		association(assocTest, 0.01);
	}

	public void association(String assocTest, double maf) throws IOException {
		checkParameters();
		if(assocTest == null) {
			throw new InvalidParameterException("association test is null");
		}
		if(!assocTest.equalsIgnoreCase("assoc") || !assocTest.equalsIgnoreCase("fisher") || !assocTest.equalsIgnoreCase("linear") || !assocTest.equalsIgnoreCase("logistic")) {
			throw new InvalidParameterException("association test is not valid, valid options are: 'assoc', 'fisher', 'linear' or 'logistic', parameter: " + assocTest);
		}
		StringBuilder plinkCommandLine = createBasicPlinkCommand();
		plinkCommandLine.append(" --" + assocTest.toLowerCase() + " --maf " + maf + " ");
		executePlinkCommand(plinkCommandLine.toString());

		FeatureData featureData = new FeatureData();
		if(assocTest.equalsIgnoreCase("assoc")) {
			// columns:  CHR, SNP, BP, A1, F_A, F_U, A2, CHISQ, P, OR

		}
		if(assocTest.equalsIgnoreCase("fisher")) {
			// columns:  CHR, SNP, BP, A1, F_A, F_U, A2, P, OR
			
		}
		if(assocTest.equalsIgnoreCase("linear")) {

		}
		if(assocTest.equalsIgnoreCase("logistic")) {
			// columns: CHR, SNP, BP, A1, TEST, NMISS, OR, STAT, P
			
		}
	}


	public void stratification() throws IOException {
		checkParameters();
		StringBuilder plinkCommandLine = createBasicPlinkCommand();
		plinkCommandLine.append(" --cluster ");
		executePlinkCommand(plinkCommandLine.toString());
	}


	/*
	 * 
	 * PRIVATE MOETHODS
	 * 
	 */

	private StringBuilder createBasicPlinkCommand() {
		StringBuilder plinkCommandLine = new StringBuilder(plinkPath);
		plinkCommandLine.append(" --ped " + pedFile.getAbsolutePath());
		plinkCommandLine.append(" --map " + mapFile.getAbsolutePath());
		plinkCommandLine.append(" --out " + outdir + "/plink ");
		return plinkCommandLine;
	}

	private void executePlinkCommand(String command) {
		Command plinkCommand = new Command(command);
		SingleProcess sp = new SingleProcess(plinkCommand);
		sp.runSync();
	}

	private void checkParameters() throws IOException {
		FileUtils.checkFile(pedFile);
		FileUtils.checkFile(mapFile);
		FileUtils.checkFile(plinkPath);
		FileUtils.checkFile(outdir);
	}


	/**
	 * @param pedFile the pedFile to set
	 */
	public void setPedFile(File pedFile) {
		this.pedFile = pedFile;
	}


	/**
	 * @return the pedFile
	 */
	public File getPedFile() {
		return pedFile;
	}


	/**
	 * @param mapFile the mapFile to set
	 */
	public void setMapFile(File mapFile) {
		this.mapFile = mapFile;
	}


	/**
	 * @return the mapFile
	 */
	public File getMapFile() {
		return mapFile;
	}

	/**
	 * @param plinkPath the plinkPath to set
	 */
	public void setPlinkPath(String plinkPath) {
		this.plinkPath = plinkPath;
	}

	/**
	 * @return the plinkPath
	 */
	public String getPlinkPath() {
		return plinkPath;
	}

	/**
	 * @param outdir the outdir to set
	 */
	public void setOutdir(String outdir) {
		this.outdir = outdir;
	}

	/**
	 * @return the outdir
	 */
	public String getOutdir() {
		return outdir;
	}

}
