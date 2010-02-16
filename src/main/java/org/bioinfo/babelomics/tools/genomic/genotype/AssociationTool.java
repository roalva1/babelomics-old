package org.bioinfo.babelomics.tools.genomic.genotype;

import java.io.IOException;
import java.util.Arrays;
import java.util.HashMap;

import org.bioinfo.babelomics.methods.genomic.genotype.GenotypeAnalysis;
import org.bioinfo.commons.io.utils.FileUtils;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;

public class AssociationTool extends GenotypeAnalysisTool {

	private static final String PLINK_OUTFILE_NAME = "plink";
	
	public AssociationTool() {

	}

	@Override
	public void initOptions() {
		options.addOption(OptionFactory.createOption("test", "Valid values: assoc, fisher, linear, logistic, tdt", true, true));
		options.addOption(OptionFactory.createOption("maf", "Minor allele frequency, default value: 0.02", false));
		options.addOption(OptionFactory.createOption("log", "Wether to apply Odd ratio logarithm to the file result", false, false));
	}

	@Override
	public void execute() {
		logger.debug("executing association test");
		try {
			// update status
			jobStatus.addStatusMessage("10", "Parsing parameters");
			
			parseGenotypeCommonOptions();
			// specific options
			String test = commandLine.getOptionValue("test", "assoc");
			double maf = Double.parseDouble(commandLine.getOptionValue("maf", "0.02"));

			// prepare the GenotypeAnalysis object for execution
			genotypeAnalysis = new GenotypeAnalysis(pedFilePath, mapFilePath);
			genotypeAnalysis.setPlinkPath(plinkPath);
			genotypeAnalysis.setOutdir(outdir);
			
			logger.debug("executing: "+plinkPath+" --ped "+pedFilePath+" --map "+mapFilePath+" --out "+outdir+"/"+PLINK_OUTFILE_NAME+" --maf "+maf + "--"+test);
			jobStatus.addStatusMessage("40", "Executing association test");
			
			genotypeAnalysis.association(test, maf);
			
			jobStatus.addStatusMessage("90", "Saving association results");
			// check if everything has gone well
			if(test != null) {
				saveResultXml(test, PLINK_OUTFILE_NAME);
			}
			
		} catch (IOException e) {
			e.printStackTrace();
			logger.error("Error opening the dataset", e.toString());
			result.addOutputItem(new Item("ioexcepcion_error", e.toString(), "An error occured", Item.TYPE.MESSAGE, Arrays.asList("ERROR"), new HashMap<String, String>(2), "Error message"));
		} catch (Exception e) {
			e.printStackTrace();
			logger.error("Error opening the dataset", e.toString());
			result.addOutputItem(new Item("excepcion_error", e.toString(), "An error occured", Item.TYPE.MESSAGE, Arrays.asList("ERROR"), new HashMap<String, String>(2), "Error message"));
		} 				
	}

	private void saveResultXml(String test, String filename) throws IOException {
		if(test.equalsIgnoreCase("assoc")) {
			FileUtils.checkFile(outdir+"/"+filename+".assoc");
			result.addOutputItem(new Item(filename+".assoc_file", filename+".assoc", "Association result file form PLINK (test: chi-square)", Item.TYPE.FILE, Arrays.asList(""), new HashMap<String, String>(2), "Association results for test "+test));
		}
		if(test.equalsIgnoreCase("fisher")) {
			FileUtils.checkFile(outdir+"/"+filename+".assoc."+test);
			result.addOutputItem(new Item(filename+".assoc."+test+"_file", filename+".assoc."+test, "Association result file form PLINK (test: "+test+")", Item.TYPE.FILE, Arrays.asList(""), new HashMap<String, String>(2), "Association results for test "+test));
		}
		if(test.equalsIgnoreCase("linear") || test.equalsIgnoreCase("logistic")) {
			FileUtils.checkFile(outdir+"/"+filename+".assoc.logistic");
			result.addOutputItem(new Item(filename+".assoc.logistic_file", filename+".assoc.logistic", "Association result file form PLINK (test: "+test+")", Item.TYPE.FILE, Arrays.asList(""), new HashMap<String, String>(2), "Association results for test "+test));
		}
		if(test.equalsIgnoreCase("tdt")) {
			FileUtils.checkFile(outdir+"/"+filename+".tdt");
			result.addOutputItem(new Item(filename+"."+test+"_file", filename+"."+test, "Association result file form PLINK (test: "+test+")", Item.TYPE.FILE, Arrays.asList(""), new HashMap<String, String>(2), "Association results for test "+test));
		}
		result.addOutputItem(new Item(filename+".hh_file", filename+".hh", "List of heterozygous haploid genotypes (SNPs/individuals)", Item.TYPE.FILE, Arrays.asList(""), new HashMap<String, String>(2), "Association results for test "+test));
		result.addOutputItem(new Item(filename+".log_file", filename+".log", "Log file from PLINK", Item.TYPE.FILE, Arrays.asList(""), new HashMap<String, String>(2), "Association results for test "+test));
	}
}
