package org.bioinfo.babelomics.methods.genomic.genotype;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import org.bioinfo.commons.exec.Command;
import org.bioinfo.commons.exec.SingleProcess;
import org.bioinfo.commons.io.TextFileReader;
import org.bioinfo.commons.io.TextFileWriter;
import org.bioinfo.commons.io.utils.FileUtils;
import org.bioinfo.commons.io.utils.IOUtils;


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

	public static void affyGenomeWide6ChpCallsToPedAndMap(String genomeWideSnpAnnotPath, String pedigreFilePath, String outdir, String chpDirPath) throws IOException {
		FileUtils.checkFile(genomeWideSnpAnnotPath);
		FileUtils.checkFile(pedigreFilePath);
		FileUtils.checkFile(chpDirPath);

		HashMap<String, String> affyIdToDbSnpId = getAffyGenomeWide6Annotation(genomeWideSnpAnnotPath);

		TextFileWriter pedFile = new TextFileWriter(outdir+".ped");
		TextFileWriter mapFile = new TextFileWriter(outdir+".map");
		TextFileReader chpFile;

		File[] files = FileUtils.listFiles(new File(chpDirPath), ".+.chp.*", true);
		StringBuilder pedStringBuilder;
		StringBuilder mapStringBuilder = null;
		String line;
		String[] chpFields;
		String[] snpFields;
		boolean createMapFile = true;
		for(File file: files) {
			pedStringBuilder = new StringBuilder();
			chpFile = new TextFileReader(file.getAbsolutePath());

			//			pedStringBuilder.append("$numLine " + file.getName() + " 0 0 0");
			//			if(file.name.startsWith("C") || file.name.startsWith("c")) {
			//				pedStringBuilder.append(" 1 ");
			//			}else {
			//				pedStringBuilder.append(" 2 ");
			//			}
			while((line = chpFile.readLine().trim()) != null) {
				if(!line.startsWith("#") && !line.startsWith("Probe") && !line.equals("")) {
					chpFields =  line.split("\t");

					if(affyIdToDbSnpId.get(chpFields[0]) != null) {
						snpFields = affyIdToDbSnpId.get(chpFields[0]).split("\t");
						if(chpFields[1].equalsIgnoreCase("AA")) {
							pedStringBuilder.append(snpFields[4] + " " + snpFields[4] + "  ");
						}else {
							if(chpFields[1].equalsIgnoreCase("AB")) {
								pedStringBuilder.append(snpFields[4] + " " + snpFields[5] + "  ");
							}else {
								if(chpFields[1].equalsIgnoreCase("BB")) {
									pedStringBuilder.append(snpFields[5] + " " + snpFields[5] + "  ");
								}else {
									pedStringBuilder.append("0 0  ");
								}
							}
						}
						// with the first file, we create the map file
						if(createMapFile) {
							mapStringBuilder = new StringBuilder();
							if(snpFields[1].equals("---")) {
								mapStringBuilder.append("0 " + snpFields[0] + " 0 " + "-1\n");
							}else{
								mapStringBuilder.append(snpFields[1] + " " + snpFields[0] + " 0 " + snpFields[2] + "\n");
							}
						}
					}
				}
			}
			chpFile.close();
			pedFile.writeLine(pedStringBuilder.toString());
			// just the first time
			if(createMapFile && mapStringBuilder != null) {
				mapFile.writeStringToFile(mapStringBuilder.toString());
				mapFile.close();
				createMapFile = false;
			}
		}
		pedFile.close();
	}

	public static void affyHumanMapping500kChpCallsToPedAndMap(String mapping250kSty, String mapping250kNsp, String pedigreFilePath, String outdir, String chpDirPath) throws IOException {
		FileUtils.checkFile(mapping250kSty);
		FileUtils.checkFile(mapping250kNsp);
		FileUtils.checkFile(pedigreFilePath);
		FileUtils.checkDirectory(chpDirPath);

		System.out.println("Getting the annotation...");
		HashMap<String, String> affyIdToDbSnpId = getAffyHumanMapping500kAnnotation(mapping250kNsp, mapping250kSty);
		File chpDirPathFile = new File(chpDirPath);

		StringBuilder pedStringBuilder;
		String line;
		String[] chpFields;
		String[] pedigreeFields;
		String[] snpFields;

		File[] styFiles;
		File[] nspFiles;
		List<String> styIds = null;
		List<String> nspIds = null;
		List<String> pedegreeLines = IOUtils.readLines(pedigreFilePath);
		List<String> pedHeaderData = new ArrayList<String>(pedegreeLines.size());
		List<String> styPedData = new ArrayList<String>(pedegreeLines.size());
		List<String> nspPedData = new ArrayList<String>(pedegreeLines.size());
		boolean styAddData = true;
		boolean nspAddData = true;
		TextFileReader tfr = null;

		for(String pedLine: pedegreeLines) {
			if(!pedLine.startsWith("#")) {
				pedigreeFields = pedLine.split("\t");
				
				pedHeaderData.add(pedigreeFields[0]+" "+pedigreeFields[1]+" "+pedigreeFields[2]+" "+pedigreeFields[3]+" "+pedigreeFields[4]+" "+pedigreeFields[5]+" ");
				// working with STY file
				if(!pedigreeFields[6].equalsIgnoreCase("NA")) {
					styFiles = FileUtils.listFiles(chpDirPathFile, pedigreeFields[6]+".+");
					System.out.println("Processing "+chpDirPathFile+"/"+pedigreeFields[6] + " files matched: "+Arrays.toString(styFiles));
					if(styFiles.length == 1) {
						tfr = new TextFileReader(styFiles[0].getAbsolutePath());
						if(styIds == null) {
							styIds = new ArrayList<String>(240000);
						}
						pedStringBuilder = new StringBuilder();
						while((line = tfr.readLine()) != null) {
							line = line.trim();
							if(!line.startsWith("#") && !line.startsWith("Probe") && !line.equals("")) {
								chpFields =  line.split("\t");
								if(affyIdToDbSnpId.get(chpFields[0]) != null) {
									snpFields = affyIdToDbSnpId.get(chpFields[0]).split("\t");
									if(chpFields[1].equalsIgnoreCase("AA")) {
										pedStringBuilder.append(snpFields[4] + " " + snpFields[4] + "  ");
									}else {
										if(chpFields[1].equalsIgnoreCase("AB")) {
											pedStringBuilder.append(snpFields[4] + " " + snpFields[5] + "  ");
										}else {
											if(chpFields[1].equalsIgnoreCase("BB")) {
												pedStringBuilder.append(snpFields[5] + " " + snpFields[5] + "  ");
											}else {
												pedStringBuilder.append("0 0  ");
											}
										}
									}
									// with the first STY file, we fiil it with the IDs
									if(styAddData) {
										styIds.add(chpFields[0]);
									}
								}
							}
						}
						styPedData.add(pedStringBuilder.toString());
						styAddData = false;
						tfr.close();
					}else {
						
					}
				}else {
					styPedData.add("");
				}
				
				// working with NSP file
				if(!pedigreeFields[7].equalsIgnoreCase("NA")) {
					nspFiles = FileUtils.listFiles(chpDirPathFile, pedigreeFields[7]+".+");
					System.out.println("Processing "+chpDirPathFile+"/"+pedigreeFields[7] + ", files matched: "+Arrays.toString(nspFiles));
					if(nspFiles.length == 1) {
						tfr = new TextFileReader(nspFiles[0].getAbsolutePath());					
						if(nspIds == null) {
							nspIds = new ArrayList<String>(270000);
						}
						pedStringBuilder = new StringBuilder();
						while((line = tfr.readLine()) != null) {
							line = line.trim();
							if(!line.startsWith("#") && !line.startsWith("Probe") && !line.equals("")) {
								chpFields =  line.split("\t");
								if(affyIdToDbSnpId.get(chpFields[0]) != null) {
									snpFields = affyIdToDbSnpId.get(chpFields[0]).split("\t");
									if(chpFields[1].equalsIgnoreCase("AA")) {
										pedStringBuilder.append(snpFields[4]).append(" ").append(snpFields[4]).append("  ");
									}else {
										if(chpFields[1].equalsIgnoreCase("AB")) {
											pedStringBuilder.append(snpFields[4]).append(" ").append(snpFields[5]).append("  ");
										}else {
											if(chpFields[1].equalsIgnoreCase("BB")) {
												pedStringBuilder.append(snpFields[5]).append(" ").append(snpFields[5]).append("  ");
											}else {
												pedStringBuilder.append("0 0  ");
											}
										}
									}
									// with the first STY file, we fiil it with the IDs
									if(nspAddData) {
										nspIds.add(chpFields[0]);
									}
								}
							}
						}
						nspPedData.add(pedStringBuilder.toString());
						nspAddData = false;
						tfr.close();
					}
					else {
						
					}
				}else {
					nspPedData.add("");
				}
			}
		}

		System.out.println("header:"+pedHeaderData.size() + " lines, sty:"+ styPedData.size()+" lines, nsp: "+nspPedData.size()+" lines");
		////////////////////////////////////////////////
		// create the PED file
		//////////////////////////////////////////////
		TextFileWriter pedFile = new TextFileWriter(outdir+"job.ped");
		StringBuilder aux;
		for(int i = 0; i < pedHeaderData.size(); i++) {
			aux = new StringBuilder();
			aux.append(pedHeaderData.get(i));
			if(styPedData.get(i).equals("")) {
				for(int j=0; j<styIds.size(); j++) {
					aux.append("0 0  ");
				}
			}else {
				aux.append(styPedData.get(i));
			}
			if(nspPedData.get(i).equals("")) {
				for(int j=0; j<nspIds.size(); j++) {
					aux.append("0 0  ");
				}
			}else {
				aux.append(nspPedData.get(i));
			}
			pedFile.writeLine(aux.toString());
		}
		pedFile.close();

		////////////////////////////////////////////////
		// create the MAP file
		//////////////////////////////////////////////
		TextFileWriter mapFile = new TextFileWriter(outdir+"job.map");
		StringBuilder mapStringBuilder = new StringBuilder();
		for(String s: styIds) {
			snpFields = affyIdToDbSnpId.get(s).split("\t");
			if(snpFields[1].equals("---")) {
				mapStringBuilder.append("0 " + snpFields[0] + " 0 " + "-1\n");
			}else{
				mapStringBuilder.append(snpFields[1] + " " + snpFields[0] + " 0 " + snpFields[2] + "\n");
			}
		}
		for(String s: nspIds) {
			snpFields = affyIdToDbSnpId.get(s).split("\t");
			if(snpFields[1].equals("---")) {
				mapStringBuilder.append("0 " + snpFields[0] + " 0 " + "-1\n");
			}else{
				mapStringBuilder.append(snpFields[1] + " " + snpFields[0] + " 0 " + snpFields[2] + "\n");
			}
		}
		mapFile.writeStringToFile(mapStringBuilder.toString());
		mapFile.close();
	}

	public static void affyHumanMapping500kMapFileCreator(String mapping250kNsp, String mapping250kSty, String outfile) throws IOException {
		FileUtils.checkFile(mapping250kNsp);
		FileUtils.checkFile(mapping250kSty);
		FileUtils.checkDirectory(new File(outfile).getParentFile());
		String line;
		String[] fields;
		// NSP chip
		TextFileWriter twr = new TextFileWriter(outfile);
		TextFileReader tfr = new TextFileReader(mapping250kNsp);
		while((line = tfr.readLine()) != null) {
			if(!line.startsWith("#") && !line.startsWith("\"Probe")) {
				fields = line.replace("\"", "").split(",");
				twr.writeLine(fields[2]+"\t"+fields[1]+"\t0\t"+fields[3]);
			}
		}
		tfr.close();
		// STY chip
		tfr = new TextFileReader(mapping250kSty);
		while((line = tfr.readLine()) != null) {
			if(!line.startsWith("#") && !line.startsWith("\"Probe")) {
				fields = line.replace("\"", "").split(",");
				twr.writeLine(fields[2]+"\t"+fields[1]+"\t0\t"+fields[3]);
			}
		}
		tfr.close();
		twr.close();
	}

	public static void affyGenomeWide6MapFileCreator(String genomeWideSnpAnnotPath, String outfile) throws IOException {
		FileUtils.checkFile(genomeWideSnpAnnotPath);
		FileUtils.checkDirectory(new File(outfile).getParentFile());
		String line;
		String[] fields;
		TextFileWriter twr = new TextFileWriter(outfile);
		TextFileReader tfr = new TextFileReader(genomeWideSnpAnnotPath);
		while((line = tfr.readLine()) != null) {
			if(!line.startsWith("#") && !line.startsWith("\"Probe")) {
				fields = line.replace("\"", "").split(",");
				twr.writeLine(fields[2]+"\t"+fields[1]+"\t0\t"+fields[3]);
			}
		}
		tfr.close();
		twr.close();
	}

	public static void affyHumanMapping500kParsedAnnotationFileCreator(String mapping250kNsp, String mapping250kSty, String outfile) throws IOException {
		FileUtils.checkFile(mapping250kNsp);
		FileUtils.checkFile(mapping250kSty);
		FileUtils.checkDirectory(new File(outfile).getParentFile());
		String line;
		String[] fields;
		// NSP chip
		TextFileWriter twr = new TextFileWriter(outfile);
		TextFileReader tfr = new TextFileReader(mapping250kNsp);
		while((line = tfr.readLine()) != null) {
			if(!line.startsWith("#") && !line.startsWith("\"Probe")) {
				fields = line.replace("\"", "").split(",");
				twr.writeLine(fields[0]+"\t"+fields[1]+"\t"+fields[2]+"\t"+fields[3]+"\t"+fields[4]+"\t"+fields[8]+"\t"+fields[9]+"\t"+fields[7]);
			}
		}
		tfr.close();
		// STY chip
		tfr = new TextFileReader(mapping250kSty);
		while((line = tfr.readLine()) != null) {
			if(!line.startsWith("#") && !line.startsWith("\"Probe")) {
				fields = line.replace("\"", "").split(",");
				twr.writeLine(fields[0]+"\t"+fields[1]+"\t"+fields[2]+"\t"+fields[3]+"\t"+fields[4]+"\t"+fields[8]+"\t"+fields[9]+"\t"+fields[7]);
			}
		}
		tfr.close();
		twr.close();
	}

	public static HashMap<String, String> getAffyHumanMapping500kAnnotation(String mapping250kNsp, String mapping250kSty) throws IOException {
		FileUtils.checkFile(mapping250kNsp);
		FileUtils.checkFile(mapping250kSty);
		HashMap<String, String> affyIdToDbSnpId = new HashMap<String, String>(1000000);
		String line;
		String[] fields;
		// NSP chip
		TextFileReader tfr = new TextFileReader(mapping250kNsp);
		while((line = tfr.readLine()) != null) {
			if(!line.startsWith("#") && !line.startsWith("\"Probe")) {
				fields = line.replace("\"", "").split(",");
				affyIdToDbSnpId.put(fields[0], fields[1]+"\t"+fields[2]+"\t"+fields[3]+"\t"+fields[4]+"\t"+fields[8]+"\t"+fields[9]+"\t"+fields[7]);
			}
		}
		tfr.close();
		// STY chip
		tfr = new TextFileReader(mapping250kSty);
		while((line = tfr.readLine()) != null) {
			if(!line.startsWith("#") && !line.startsWith("\"Probe")) {
				fields = line.replace("\"", "").split(",");
				affyIdToDbSnpId.put(fields[0], fields[1]+"\t"+fields[2]+"\t"+fields[3]+"\t"+fields[4]+"\t"+fields[8]+"\t"+fields[9]+"\t"+fields[7]);
			}
		}
		tfr.close();
		return affyIdToDbSnpId;
	}

	public static void affyGenomeWide6ParsedAnnotationFileCreator(String genomeWideSnpAnnotPath, String outfile) throws IOException {
		FileUtils.checkFile(genomeWideSnpAnnotPath);
		FileUtils.checkDirectory(new File(outfile).getParentFile());
		String line;
		String[] fields;
		TextFileWriter twr = new TextFileWriter(outfile);
		TextFileReader tfr = new TextFileReader(genomeWideSnpAnnotPath);
		while((line = tfr.readLine()) != null) {
			if(!line.startsWith("#") && !line.startsWith("\"Probe")) {
				fields = line.replace("\"", "").split(",");
				twr.writeLine(fields[0]+"\t"+fields[1]+"\t"+fields[2]+"\t"+fields[3]+"\t"+fields[4]+"\t"+fields[8]+"\t"+fields[9]+"\t"+fields[7]);
			}
		}
		tfr.close();
		twr.close();
	}

	public static HashMap<String, String> getAffyGenomeWide6Annotation(String genomeWideSnpAnnotPath) throws IOException {
		FileUtils.checkFile(genomeWideSnpAnnotPath);
		HashMap<String, String> affyIdToDbSnpId = new HashMap<String, String>(1500000);
		String line;
		String[] fields;
		TextFileReader tfr = new TextFileReader(genomeWideSnpAnnotPath);
		int cont = 0;
		while((line = tfr.readLine()) != null) {
			if(!line.startsWith("#") && !line.startsWith("\"Probe")) {
				fields = line.replace("\"", "").split(",");
				affyIdToDbSnpId.put(fields[0], fields[1]+"\t"+fields[2]+"\t"+fields[3]+"\t"+fields[4]+"\t"+fields[8]+"\t"+fields[9]+"\t"+fields[7]);
				if(cont++ % 50000 == 0) {
					System.out.println(cont+"  "+affyIdToDbSnpId.get(fields[0]));
				}
			}
		}
		tfr.close();
		return affyIdToDbSnpId;
	}

}
