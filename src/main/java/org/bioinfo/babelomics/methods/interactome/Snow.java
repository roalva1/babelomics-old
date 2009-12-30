package org.bioinfo.babelomics.methods.interactome;

import org.bioinfo.babelomics.exception.InvalidParameterException;
import org.bioinfo.commons.exec.Command;
import org.bioinfo.commons.exec.SingleProcess;

public class Snow {

	String mode = "unknown";
	String interactome = "unknown";
	String list1FileName = null;
	String list2FileName = null;
	String interactionsFileName = null;
	boolean checkInteractions = true;
	String idNature = "protein";
	int interactionsNumber = 1;
	
	String snowBinPath = null;
	String outDirName = null;
	
	public Snow(String snowBinPath, String list1FileName, String interactome, String idNature, int interactionsNumber, String outdirName) {
		this.snowBinPath = snowBinPath;
		this.list1FileName = list1FileName;
		this.interactome = interactome;
		this.idNature = idNature;
		this.interactionsNumber = interactionsNumber;
		this.outDirName = outDirName;

		this.mode = "one-list";
	}
	
	public Snow(String snowBinPath, String list1FileName, String interactionsFileName, boolean checkInteractions, String idNature, int interactionsNumber, String outdirName) {
		this.snowBinPath = snowBinPath;
		this.list1FileName = list1FileName;
		this.interactome = "own";
		this.interactionsFileName = interactionsFileName;
		this.checkInteractions = checkInteractions;
		this.idNature = idNature;
		this.interactionsNumber = interactionsNumber;
		this.outDirName = outDirName;

		this.mode = "one-list";
	}

	public Snow(String snowBinPath, String list1FileName, String list2FileName, String interactome, String idNature, int interactionsNumber, String outdirName) {
		this.snowBinPath = snowBinPath;
		this.list1FileName = list1FileName;
		this.list2FileName = list2FileName;
		this.interactome = interactome;
		this.idNature = idNature;
		this.interactionsNumber = interactionsNumber;
		this.outDirName = outDirName;
		
		this.mode = "two-list";
	}

	public Snow(String snowBinPath, String list1FileName, String list2FileName, String interactionsFileName, boolean checkInteractions, String idNature, int interactionsNumber, String outdirName) {		
		this.snowBinPath = snowBinPath;
		this.list1FileName = list1FileName;
		this.list2FileName = list2FileName;
		this.interactome = "own";
		this.interactionsFileName = interactionsFileName;
		this.checkInteractions = checkInteractions;
		this.idNature = idNature;
		this.interactionsNumber = interactionsNumber;
		this.outDirName = outDirName;
		
		this.mode = "two-list";
	}


	public void run() throws InvalidParameterException {
		if ( snowBinPath == null ) {
			throw new InvalidParameterException("Missing binary path to snow script");
		}
		
		String cmd = createSnowCommand();
		executeSnowCommand(cmd);
		
	}
	
	private String createSnowCommand() {
		StringBuilder cmd = new StringBuilder(snowBinPath);
		if ( )
		cmd.append(" ").append(list1FileName);
		cmd.append(" ").append(outDirName);
		return cmd.toString();
	}

	private void executeSnowCommand(String cmd) {
		Command snowCommand = new Command(cmd);
		SingleProcess sp = new SingleProcess(snowCommand);
		sp.runSync();
	}
}
