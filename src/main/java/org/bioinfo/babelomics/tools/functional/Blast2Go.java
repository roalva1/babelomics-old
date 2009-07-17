package org.bioinfo.babelomics.tools.functional;

import java.io.File;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.ParseException;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.tool.OptionFactory;

public class Blast2Go extends BabelomicsTool {


	public Blast2Go(String[] args) {
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
		//			CommandLine cmd = parse(args);

		File xmlfile = new File(commandLine.getOptionValue("xmlfile"));
		String xmlVersion = commandLine.getOptionValue("xmlVersion");
		String blastMinNum = commandLine.getOptionValue("blastMin-NumberFilter");
		String blastHit = commandLine.getOptionValue("blastHit-DescriptFilter");
		String blastHitDescPosition = commandLine.getOptionValue("blastHit-DescriptPosition");	
		String blastId = commandLine.getOptionValue("blastId", null);
		String descriptAnno = commandLine.getOptionValue("descriptAnno", null);			

		File eValue = new File(commandLine.getOptionValue("eValue-filter"));
		File annCutOff = new File(commandLine.getOptionValue("annCutOff-filter"));
		File goWeight = new File(commandLine.getOptionValue("goWeight-filter"));
		File hsp = new File(commandLine.getOptionValue("hsp-filter"));

		String ida = commandLine.getOptionValue("ida", null);
		String ima = commandLine.getOptionValue("ima", null);
		String igi = commandLine.getOptionValue("igi", null);
		String ipi = commandLine.getOptionValue("ipi", null);
		String iep = commandLine.getOptionValue("iep", null);
		String tas = commandLine.getOptionValue("tas", null);
		String nac = commandLine.getOptionValue("nac", null);
		String ic = commandLine.getOptionValue("ic", null);
		String iss = commandLine.getOptionValue("iss", null);
		String rca = commandLine.getOptionValue("rca", null);
		String iea = commandLine.getOptionValue("iea", null);
		String nd = commandLine.getOptionValue("nd", null);
		String nr = commandLine.getOptionValue("nr", null);



		//executeBlast2Go(xmlfile, xmlVersion, blastMinNum, blastHit,blastHitDescPosition, blastId, descriptAnno,  eValue,annCutOff,goWeight,hsp,ida,ima);

		//		catch (IOException e) {
		//			logger.error("Error opening the dataset", e.toString());
		//		}



	}

	private void executeBlast2Go(File xmlfile, String xmlVersion, String blastMinNum, String bioentityNumberFilter,String geneNameList, String partitionNumber, String significance, String sort) {
		logger.info("executing svm, not implemented yet");
	}

}
