package org.bioinfo.babelomics.tools.functional;

import java.io.File;

import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.tool.OptionFactory;

public class Blast2Go extends BabelomicsTool {


	public Blast2Go() {
		initOptions();
	}

	@Override
	public void initOptions() {
		getOptions().addOption(OptionFactory.createOption("xmlfile", "Blast file (xml format)"));
		getOptions().addOption(OptionFactory.createOption("xmlVersion", ""));
		getOptions().addOption(OptionFactory.createOption("blastNumber", "Minimum number of genes with a score (0-10000)"));
		getOptions().addOption(OptionFactory.createOption("blastMin-NumberFilter", "Number of bio-entities in results (0-10000)"));
		getOptions().addOption(OptionFactory.createOption("blastHit-DescriptFilter", "gene name in list"));
		getOptions().addOption(OptionFactory.createOption("blastHit-DescriptPosition", "Number of partitions "));
		getOptions().addOption(OptionFactory.createOption("blastId", "Add ID to Blast definition", false));
		getOptions().addOption(OptionFactory.createOption("descriptAnno", "Sort list", false));

		getOptions().addOption(OptionFactory.createOption("eValue-filter", "E-value hit filter ((1e-)"));
		getOptions().addOption(OptionFactory.createOption("annCutOff-filter", "Annotation cut-off (0-100)"));
		getOptions().addOption(OptionFactory.createOption("goWeight-filter", "GO weight (0-30)"));
		getOptions().addOption(OptionFactory.createOption("hsp-filter", "Hsp hit coverage cut-off (0-100)"));

		//EC codes
		getOptions().addOption(OptionFactory.createOption("ida", "ida", false));
		getOptions().addOption(OptionFactory.createOption("imp", "imp", false));
		getOptions().addOption(OptionFactory.createOption("igi", "igi", false));
		getOptions().addOption(OptionFactory.createOption("ipi", "ipi", false));
		getOptions().addOption(OptionFactory.createOption("iep", "iep", false));
		getOptions().addOption(OptionFactory.createOption("tas", "tas", false));
		getOptions().addOption(OptionFactory.createOption("nac", "nac", false));
		getOptions().addOption(OptionFactory.createOption("ic", "ic", false));
		getOptions().addOption(OptionFactory.createOption("iss", "iss", false));
		getOptions().addOption(OptionFactory.createOption("rca", "rca", false));		
		getOptions().addOption(OptionFactory.createOption("iea", "iea", false));
		getOptions().addOption(OptionFactory.createOption("nd", "nd", false));
		getOptions().addOption(OptionFactory.createOption("nr", "nr", false));




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
