package org.bioinfo.babelomics.methods.genomic.genotype;

import java.io.File;
import java.io.IOException;
import java.security.InvalidParameterException;

import org.bioinfo.commons.exec.Command;
import org.bioinfo.commons.exec.SingleProcess;
import org.bioinfo.commons.io.TextFileWriter;
import org.bioinfo.commons.io.utils.FileUtils;


public class AffyGenotypeUtils {

	public static void affyGenotypeNormalization(String aptProbesetPath, String chipType, String cdfPath, String chrXPath) {
		StringBuilder aptCommandLine = new StringBuilder();
		aptCommandLine.append(aptProbesetPath);
		Command aptCommand = new Command(aptCommandLine.toString());
		SingleProcess sp = new SingleProcess(aptCommand);
		sp.runSync();
	}
	
	public static void affyToPedAndMap(String chpDirPath, String pedigreFilePath, String outdir) throws InvalidParameterException, IOException {
		File pedigreeFile = new File(chpDirPath);
		File chpDir = new  File(pedigreFilePath);
		if(!pedigreeFile.exists() || !chpDir.exists()) {
			throw new InvalidParameterException("Some parameters are not valid, chp-dir: " +chpDirPath+",  pedigree: "+pedigreFilePath);
		}
		String affyArrayType = "affymetrix-array-type";
		
		TextFileWriter pedFile = new TextFileWriter(outdir+".ped");
		TextFileWriter mapFile = new TextFileWriter(outdir+".map");
		File[] files = FileUtils.listFiles(chpDir, ".+.chp.*", true);
		
		StringBuilder pedLine;
		for(File file: files) {
			pedLine = new StringBuilder();
			
			
			
			
			pedFile.writeLine(pedLine.toString());
		}
		
		
		
		pedFile.close();
		mapFile.close();
	}

}
