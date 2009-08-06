package org.bioinfo.babelomics.tools.preprocessing;

	import java.io.File;
import java.io.IOException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.ParseException;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.tool.OptionFactory;



public class IdConverter extends BabelomicsTool  {

	
		public IdConverter(String[] args) {
			initOptions();
		}

		@Override
		public void initOptions() {
			options.addOption(OptionFactory.createOption("dataset", "the data"));
			options.addOption(OptionFactory.createOption("microarray", "Agilent microarray G2518A", false));
			options.addOption(OptionFactory.createOption("go", "GO converter", false,false));			
			options.addOption(OptionFactory.createOption("pdb", "PDB converter", false,false));
			options.addOption(OptionFactory.createOption("refseq", "RefSeq DNA predicted", false,false));
			

		}

		@Override
		public void execute() {
			try {
				
				
//				CommandLine cmd = parse(args, true);
				
				Dataset dataset = new Dataset(new File(commandLine.getOptionValue("dataset")));
				String microarray = commandLine.getOptionValue("microarray", null);
				String go = commandLine.getOptionValue("go", null);
				String pdb = commandLine.getOptionValue("pdb", null);
				String refseq = commandLine.getOptionValue("refseq", null);
				
				System.out.println(dataset.toString()+"\n");
				
				
				if ( microarray != null ) {
					logger.info("Agilent microarray G2518A converter, not yet implemented");
					}

				if ( go != null ) {
					logger.info("go converter, not yet implemented");
					}
				
				if ( pdb != null ) {
					logger.info("pdb converter, not yet implemented");
					}
				
				if ( refseq != null ) {
					logger.info("refseq converter, not yet implemented");
					}
			} catch (Exception e) {
				logger.error("Error opening the dataset", e.toString());
			} 
		}

				
}
