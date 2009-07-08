package org.bioinfo.babelomics.tools.functional;

import java.io.File;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.ParseException;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.tool.OptionFactory;

public class Blast2Go extends BabelomicsTool {


	public Blast2Go(String[] args) {
		super(args);
		initOptions();
	}

	@Override
	public void initOptions() {
		options.addOption(OptionFactory.createOption("xmlfile", "Blast file (xml format)"));
		options.addOption(OptionFactory.createOption("xmlVersion", ""));
		options.addOption(OptionFactory.createOption("blastNumber", "Minimum number of genes with a score (0-10000)"));
		options.addOption(OptionFactory.createOption("blastMin-NumberFilter", "Number of bio-entities in results (0-10000)"));
		options.addOption(OptionFactory.createOption("blastHit-DescriptFilter", "gene name in list"));
		options.addOption(OptionFactory.createOption("blastHit-DescriptPosition", "Number of partitions "));
		options.addOption(OptionFactory.createOption("blastId", "Add ID to Blast definition", false));
		options.addOption(OptionFactory.createOption("descriptAnno", "Sort list", false));
		
		options.addOption(OptionFactory.createOption("eValue-filter", "E-value hit filter ((1e-)"));
		options.addOption(OptionFactory.createOption("annCutOff-filter", "Annotation cut-off (0-100)"));
		options.addOption(OptionFactory.createOption("goWeight-filter", "GO weight (0-30)"));
		options.addOption(OptionFactory.createOption("hsp-filter", "Hsp hit coverage cut-off (0-100)"));
		
		//EC codes
		options.addOption(OptionFactory.createOption("ida", "ida", false));
		options.addOption(OptionFactory.createOption("imp", "imp", false));
		options.addOption(OptionFactory.createOption("igi", "igi", false));
		options.addOption(OptionFactory.createOption("ipi", "ipi", false));
		options.addOption(OptionFactory.createOption("iep", "iep", false));
		options.addOption(OptionFactory.createOption("tas", "tas", false));
		options.addOption(OptionFactory.createOption("nac", "nac", false));
		options.addOption(OptionFactory.createOption("ic", "ic", false));
		options.addOption(OptionFactory.createOption("iss", "iss", false));
		options.addOption(OptionFactory.createOption("rca", "rca", false));		
		options.addOption(OptionFactory.createOption("iea", "iea", false));
		options.addOption(OptionFactory.createOption("nd", "nd", false));
		options.addOption(OptionFactory.createOption("nr", "nr", false));
		
		
		
		
	}

	@Override
	public void execute() {
		try {
			CommandLine cmd = parse(args);

			File xmlfile = new File(cmd.getOptionValue("xmlfile"));
			String xmlVersion = cmd.getOptionValue("xmlVersion");
			String blastMinNum = cmd.getOptionValue("blastMin-NumberFilter");
			String blastHit = cmd.getOptionValue("blastHit-DescriptFilter");
			String blastHitDescPosition = cmd.getOptionValue("blastHit-DescriptPosition");	
			String blastId = cmd.getOptionValue("blastId", null);
			String descriptAnno = cmd.getOptionValue("descriptAnno", null);			
			
			File eValue = new File(cmd.getOptionValue("eValue-filter"));
			File annCutOff = new File(cmd.getOptionValue("annCutOff-filter"));
			File goWeight = new File(cmd.getOptionValue("goWeight-filter"));
			File hsp = new File(cmd.getOptionValue("hsp-filter"));
			
			String ida = cmd.getOptionValue("ida", null);
			String ima = cmd.getOptionValue("ima", null);
			String igi = cmd.getOptionValue("igi", null);
			String ipi = cmd.getOptionValue("ipi", null);
			String iep = cmd.getOptionValue("iep", null);
			String tas = cmd.getOptionValue("tas", null);
			String nac = cmd.getOptionValue("nac", null);
			String ic = cmd.getOptionValue("ic", null);
			String iss = cmd.getOptionValue("iss", null);
			String rca = cmd.getOptionValue("rca", null);
			String iea = cmd.getOptionValue("iea", null);
			String nd = cmd.getOptionValue("nd", null);
			String nr = cmd.getOptionValue("nr", null);
			
			
			
			//executeBlast2Go(xmlfile, xmlVersion, blastMinNum, blastHit,blastHitDescPosition, blastId, descriptAnno,  eValue,annCutOff,goWeight,hsp,ida,ima);
			
		} catch (ParseException e) {
			logger.error("Error parsing command line", e.toString());
			System.out.println("\n");
			printUsage();
		} 
//		catch (IOException e) {
//			logger.error("Error opening the dataset", e.toString());
//		}
		
		
		
	}
	
	private void executeBlast2Go(File xmlfile, String xmlVersion, String blastMinNum, String bioentityNumberFilter,String geneNameList, String partitionNumber, String significance, String sort) {
		logger.info("executing svm, not implemented yet");
	}
	
}
