package org.bioinfo.babelomics.tools.variation;

import java.io.File;
import java.io.IOException;
import java.lang.reflect.Constructor;
import java.lang.reflect.Type;
import java.util.List;

import org.apache.commons.cli.CommandLine;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.commons.io.utils.FileUtils;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.data.graph.edge.AbstractEdge;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.variation.data.converter.Converter;
import org.bioinfo.variation.data.core.DataSet;
import org.bioinfo.variation.data.core.MapFile;
import org.bioinfo.variation.data.core.PedFile;
import org.bioinfo.variation.data.core.PedMap;
import org.bioinfo.variation.data.core.filter.Filter;
import org.bioinfo.variation.data.core.filter.Snp.BasePairPosSnpFilter;
import org.bioinfo.variation.data.core.filter.Snp.ChromosomeSnpFilter;
import org.bioinfo.variation.data.core.filter.Snp.GeneticDistanceSnpFilter;
import org.bioinfo.variation.data.core.filter.Snp.RegionSnpFilter;
import org.bioinfo.variation.data.core.filter.Snp.SnpFilter;
import org.bioinfo.variation.data.core.filter.Snp.SnpIdSnpFilter;
import org.bioinfo.variation.data.core.filter.individual.FamilyIndividualFilter;
import org.bioinfo.variation.data.core.filter.individual.FounderIndividualFilter;
import org.bioinfo.variation.data.core.filter.individual.IDIndividualFilter;
import org.bioinfo.variation.data.core.filter.individual.MaternalIndividualFilter;
import org.bioinfo.variation.data.core.filter.individual.PaternalIndividualFilter;
import org.bioinfo.variation.data.core.filter.individual.PhenotypeIndividualFilter;
import org.bioinfo.variation.data.core.filter.individual.SexIndividualFilter;

public class Ped extends BabelomicsTool {

	private File inputPedFile, inputMapFile, inputSnpFile;
	private File outputPedFile, outputMapFile, outputFile;
	private PedMap pedmap;
	private DataSet ds;
//	private List<String> snpList;
//	private PedParser parser;
//	private String format;
	
	@Override
	public void initOptions() {
		options.addOption(OptionFactory.createOption("input-ped-file", "An input file containing a PED file", true));
		options.addOption(OptionFactory.createOption("input-map-file", "An input file containing a MAP file", true));
		options.addOption(OptionFactory.createOption("output-ped-file", "An output file containing a PED file", false));
		options.addOption(OptionFactory.createOption("output-map-file", "An output file containing a MAP file", false));
		options.addOption(OptionFactory.createOption("output-file", "Name of an output file", false));
		
		//Inputs for filters
		options.addOption(OptionFactory.createOption("input-snp-file", "An input file with snps separeted by line", false, true));
		options.addOption(OptionFactory.createOption("input-base-file", "An input file with bases separeted by line", false, true));
		options.addOption(OptionFactory.createOption("input-ch-file", "An input file with chromosomes separeted by line", false, true));
		options.addOption(OptionFactory.createOption("input-gd-file", "An input file with genetic ditances separeted by line", false, true));
		options.addOption(OptionFactory.createOption("input-gd-file", "An input file with genetic ditances separeted by line", false, true));
		options.addOption(OptionFactory.createOption("ch", "A chromosome for region filter", false, true));
		options.addOption(OptionFactory.createOption("from-base", "Region filter: from base", false, true));
		options.addOption(OptionFactory.createOption("to-base", "Region filter: to base", false, true));
		options.addOption(OptionFactory.createOption("family", "Filter by family ID", false, true));
		options.addOption(OptionFactory.createOption("founder", "Filter by family ID", false, false));
		options.addOption(OptionFactory.createOption("individual-id", "Filter by individual ID", false, true));
		options.addOption(OptionFactory.createOption("maternal-id", "Filter by maternal ID", false, true));
		options.addOption(OptionFactory.createOption("paternal-id", "Filter by paternal ID", false, true));		
		options.addOption(OptionFactory.createOption("phenotype", "Filter by phenotype", false, true));
		options.addOption(OptionFactory.createOption("sex", "Filter by sex", false, true));

		options.addOption(OptionFactory.createOption("format", "Format to a tped or a proper ped file. Args are ped or tped.", false, true));
		
	}

	@Override
	protected void execute() {
		
		logger.debug("reading input-file...");	
		inputPedFile = new File(commandLine.getOptionValue("input-ped-file"));
		inputMapFile = new File(commandLine.getOptionValue("input-map-file"));
//		outputPedFile = new File(commandLine.getOptionValue("output-ped-file"));
//		outputMapFile = new File(commandLine.getOptionValue("output-map-file"));
		logger.debug("checking if inputFile exist...");
		try {
			FileUtils.checkFile(inputPedFile);
			FileUtils.checkFile(inputMapFile);
			
			PedFile pedFile = new PedFile(inputPedFile.getAbsolutePath());
			MapFile mapFile = new MapFile(inputMapFile.getAbsolutePath());
			ds = new DataSet(pedFile, mapFile);
			
			if(commandLine.hasOption("input-snp-file")){
				FileUtils.checkFile(commandLine.getOptionValue("input-snp-file"));
				List<String> list = IOUtils.readLines(commandLine.getOptionValue("input-snp-file"));
				SnpIdSnpFilter snpIdFilter = new SnpIdSnpFilter(list);
				ds.filter(snpIdFilter);
			}
			if(commandLine.hasOption("input-base-file")){
				FileUtils.checkFile(commandLine.getOptionValue("input-base-file"));
				List<String> list = IOUtils.readLines(commandLine.getOptionValue("input-base-file"));
				BasePairPosSnpFilter filter = new BasePairPosSnpFilter(list);
				ds.filter(filter);
			}
			if(commandLine.hasOption("input-ch-file")){
				FileUtils.checkFile(commandLine.getOptionValue("input-ch-file"));
				List<String> list = IOUtils.readLines(commandLine.getOptionValue("input-ch-file"));
				ChromosomeSnpFilter filter = new ChromosomeSnpFilter(list);
				ds.filter(filter);
			}
			if(commandLine.hasOption("input-gd-file")){
				FileUtils.checkFile(commandLine.getOptionValue("input-gd-file"));
				List<String> list = IOUtils.readLines(commandLine.getOptionValue("input-gd-file"));
				GeneticDistanceSnpFilter filter = new GeneticDistanceSnpFilter(list);
				ds.filter(filter);
			}
			if(commandLine.hasOption("ch") && commandLine.hasOption("from-base") && commandLine.hasOption("to-base")){
				RegionSnpFilter filter = new RegionSnpFilter(commandLine.getOptionValue("ch"), Integer.parseInt(commandLine.getOptionValue("from-base")),Integer.parseInt(commandLine.getOptionValue("to-base")));
				ds.filter(filter);
			}
			if(commandLine.hasOption("family")){
				FamilyIndividualFilter filter = new FamilyIndividualFilter(commandLine.getOptionValue("family"));
				ds.filter(filter);
			}
			if(commandLine.hasOption("founder")){
				FounderIndividualFilter filter = new FounderIndividualFilter();
				ds.filter(filter);
			}
			if(commandLine.hasOption("individual-id")){
				IDIndividualFilter filter = new IDIndividualFilter(commandLine.getOptionValue("individual-id"));
				ds.filter(filter);
			}
			if(commandLine.hasOption("maternal-id")){
				MaternalIndividualFilter filter = new MaternalIndividualFilter(commandLine.getOptionValue("maternal-id"));
				ds.filter(filter);
			}
			if(commandLine.hasOption("paternal-id")){
				PaternalIndividualFilter filter = new PaternalIndividualFilter(commandLine.getOptionValue("paternal-id"));
				ds.filter(filter);
			}
			if(commandLine.hasOption("phenotype")){
				PhenotypeIndividualFilter filter = new PhenotypeIndividualFilter(commandLine.getOptionValue("phenotype"));
				ds.filter(filter);
			}
			if(commandLine.hasOption("sex")){
				SexIndividualFilter filter = new SexIndividualFilter(commandLine.getOptionValue("sex"));
				ds.filter(filter);
			}
			Converter c = new Converter(ds);
c.toPed();
//			System.out.println(ds.getPedfile();
//			System.out.println(ds.snpIndexesToString(ds.getSnpIndex()));
			
		} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}


