package org.bioinfo.babelomics.methods.genomic.genotype;

import java.io.File;
import java.security.InvalidParameterException;

import org.bioinfo.commons.exec.Command;
import org.bioinfo.commons.exec.SingleProcess;

public class PlinkTestsExecutor {

	private File pedFile;
	private File mapFile;
	private String plinkPath;
	private String outdir;
	
	public PlinkTestsExecutor() {
		this.pedFile = null;
		this.mapFile = null;
	}
	
	public PlinkTestsExecutor(String pedFilePath, String mapFilePath) {
		this.pedFile = new File(pedFilePath);
		this.mapFile = new File(mapFilePath);
	}
	
	public PlinkTestsExecutor(File pedFile, File mapFile) {
		this.pedFile = pedFile;
		this.mapFile = mapFile;
	}
	
	/*
	 * 
	 * PLINK TESTS
	 * 
	 */
	
	public void association(String assocTest) throws InvalidParameterException {
		checkParameters();
		if(assocTest == null) {
			throw new InvalidParameterException("association test is null");
		}
		if(!assocTest.equalsIgnoreCase("assoc") || !assocTest.equalsIgnoreCase("fisher")) {
			throw new InvalidParameterException("association test is not valid, valid options are: 'assoc' or 'fisher', parameter: " + assocTest);
		}
		StringBuilder plinkCommandLine = createBasicPlinkCommand();
		plinkCommandLine.append(" --" +  assocTest.toLowerCase());
		executePlinkCommand(plinkCommandLine.toString());
	}

	
	public void stratification() throws InvalidParameterException {
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
	
	private void checkParameters() throws InvalidParameterException {
		if(pedFile == null) {
			throw new InvalidParameterException("pedFile is null");
		}
		if(mapFile == null) {
			throw new InvalidParameterException("mapFile is null");
		}
		if(!pedFile.exists()) {
			throw new InvalidParameterException("pedFile file does not exist: " + getPedFile().getAbsolutePath());
		}
		if(!mapFile.exists()) {
			throw new InvalidParameterException("mapFile file does not exist: " + getMapFile().getAbsolutePath());
		}
		if(plinkPath == null) {
			throw new InvalidParameterException("plinkPath is null");
		}
		if(!new File(plinkPath).exists()) {
			throw new InvalidParameterException("plinkPath does not exist: " + getPedFile().getAbsolutePath());
		}
		if(outdir == null) {
			throw new InvalidParameterException("outdir is null");
		}
		if(!new File(outdir).exists()) {
			throw new InvalidParameterException("outdir does not exist: " + getPedFile().getAbsolutePath());
		}
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
