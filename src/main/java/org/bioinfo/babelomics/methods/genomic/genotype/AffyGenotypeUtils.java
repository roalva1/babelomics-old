package org.bioinfo.babelomics.methods.genomic.genotype;

import java.io.File;
import java.io.IOException;

import org.bioinfo.commons.exec.Command;
import org.bioinfo.commons.exec.SingleProcess;
import org.bioinfo.commons.io.TextFileWriter;
import org.bioinfo.commons.io.utils.FileUtils;


public class AffyGenotypeUtils {

	public static void affyGenotypeNormalization(String aptProbesetPath, String chipType, String cdfPath, String chrXPath, String outdir, String dataDir) throws IOException {
		StringBuilder aptCommandLine = new StringBuilder();
		aptCommandLine.append(aptProbesetPath);
		
		FileUtils.checkFile(cdfPath);
		FileUtils.checkFile(chrXPath);
		FileUtils.checkDirectory(outdir, true);
		FileUtils.checkDirectory(dataDir, false);
		
		aptCommandLine.append(" -c " + cdfPath);
		aptCommandLine.append(" --chrX-snps " + chrXPath);

		if(chipType.equalsIgnoreCase("6.0")) {
			aptCommandLine.append(" -a birdseed ");
		}		
		aptCommandLine.append(" --cc-chp-output ");
		aptCommandLine.append(" -o " + outdir);
		
		File celFiles = new File(dataDir);
		File[] files = FileUtils.listFiles(celFiles, ".+.cel", true);
		for(int i=0; i<files.length; i++) {
			aptCommandLine.append(" " + files[i].getAbsolutePath());
		}
		
		Command aptCommand = new Command(aptCommandLine.toString());
		SingleProcess sp = new SingleProcess(aptCommand);
		sp.runSync();
	}
	
	public static void affyToPedAndMap(String pedigreFilePath, String outdir, String chpDirPath) throws IOException {
		FileUtils.checkFile(pedigreFilePath);
		FileUtils.checkFile(chpDirPath);
		String affyArrayType = "affymetrix-array-type";
		
		TextFileWriter pedFile = new TextFileWriter(outdir+".ped");
		TextFileWriter mapFile = new TextFileWriter(outdir+".map");
		File[] files = FileUtils.listFiles(new File(chpDirPath), ".+.chp.*", true);
		
		StringBuilder pedLine;
		for(File file: files) {
			pedLine = new StringBuilder();
			
			
			
			
			pedFile.writeLine(pedLine.toString());
		}
		
		
		
		pedFile.close();
		mapFile.close();
	}

}
