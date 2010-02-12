package org.bioinfo.babelomics.tools.preprocessing;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.cli.CommandLine;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.commons.utils.MapUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.infrared.common.dbsql.DBConnector;
import org.bioinfo.infrared.common.feature.FeatureList;
import org.bioinfo.infrared.core.dbsql.GeneDBManager;
import org.bioinfo.infrared.funcannot.AnnotationItem;
import org.bioinfo.infrared.funcannot.dbsql.AnnotationDBManager;
import org.bioinfo.infrared.funcannot.filter.BiocartaFilter;
import org.bioinfo.infrared.funcannot.filter.Filter;
import org.bioinfo.infrared.funcannot.filter.GOFilter;
import org.bioinfo.infrared.funcannot.filter.KeggFilter;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;

public class CreateAnnotation extends BabelomicsTool {


	List<Filter> filterList = new ArrayList<Filter>();

	@Override
	public void initOptions() {

		getOptions().addOption(OptionFactory.createOption("listfile", "File containning the IDs to convert", false));
		getOptions().addOption(OptionFactory.createOption("list", "IDs to convert separated by commas", false));
		getOptions().addOption(OptionFactory.createOption("all-genome", "All genome",false));
		
		// GO biological process options
		addGOOptions("bp");
		addGOOptions("cc");
		addGOOptions("mf");

		// kegg
		getOptions().addOption(OptionFactory.createOption("kegg", "Kegg database",false,false));
		getOptions().addOption(OptionFactory.createOption("kegg-min-num-genes", "Kegg min number of genes",false));
		getOptions().addOption(OptionFactory.createOption("kegg-max-num-genes", "Kegg max number of genes",false));
		getOptions().addOption(OptionFactory.createOption("kegg-count-genes-from-genome", "computes the number of annotated genes from all genome, otherwise from you input list",false,false));

		// biocarta
		getOptions().addOption(OptionFactory.createOption("biocarta", "Biocarta database",false,false));
		getOptions().addOption(OptionFactory.createOption("biocarta-min-num-genes", "Biocarta min number of genes",false));
		getOptions().addOption(OptionFactory.createOption("biocarta-max-num-genes", "Biocarta max number of genes",false));
		getOptions().addOption(OptionFactory.createOption("biocarta-count-genes-from-genome", "computes the number of annotated genes from all genome, otherwise from you input list",false,false));

		// output format
		getOptions().addOption(OptionFactory.createOption("output-format", "Output format: compact or extended. Default: compact", false));
	}


	private void addGOOptions(String namespace){
		String namespaceTitle = "biological process";		
		if(namespace.equalsIgnoreCase("cc")) namespaceTitle = "cellular component";
		if(namespace.equalsIgnoreCase("mf")) namespaceTitle = "molecular function";
		getOptions().addOption(OptionFactory.createOption("go-" + namespace, "GO " + namespaceTitle + " database",false,false));
		getOptions().addOption(OptionFactory.createOption("go-" + namespace + "-inclusive", "GO " + namespaceTitle + ", inclusive analysis (one per GO level), otherwise joins GO levels",false,false));
		getOptions().addOption(OptionFactory.createOption("go-" + namespace + "-min-level", "GO " + namespaceTitle + ", min go level to take into account, default 3",false));
		getOptions().addOption(OptionFactory.createOption("go-" + namespace + "-max-level", "GO " + namespaceTitle + ", max GO level to take into account, default 15",false));
		getOptions().addOption(OptionFactory.createOption("go-" + namespace + "-min-num-genes", "GO " + namespaceTitle + ", min number of genes filter",false));
		getOptions().addOption(OptionFactory.createOption("go-" + namespace + "-max-num-genes", "GO " + namespaceTitle + ", max number of genes filter",false));
		getOptions().addOption(OptionFactory.createOption("go-" + namespace + "-all-genome", "GO " + namespaceTitle + ", computes the number of annotated genes from all genome, otherwise from you input list",false,false));
		getOptions().addOption(OptionFactory.createOption("go-" + namespace + "-keywords", "GO " + namespaceTitle + ", keywords filter",false));
		getOptions().addOption(OptionFactory.createOption("go-" + namespace + "-keywords-logic", "GO " + namespaceTitle + ", keywords filter logic: all or any",false));
	}


	@Override
	protected void execute() {
		try {	
			
			DBConnector dbConnector = new DBConnector(species, new File(System.getenv("BABELOMICS_HOME") + "/conf/infrared.conf"));	

			System.out.println("species = " + species + ", db connector = " + dbConnector.toString());

			List<String> ids = null;

			// check input IDs
			//
			String inputIds = commandLine.getOptionValue("list", null);			
			if ( inputIds == null ) {
				String fileName = commandLine.getOptionValue("listfile", null);
				if ( fileName != null ) {
					try {
						ids = IOUtils.column(new File(fileName), 0);
					} catch (IOException e) {
						abort("ioexception_execute_newannotation", e.getMessage(), e.getMessage(), e.getMessage());
					}
				} else {
					if ( commandLine.hasOption("all-genome") ) {
						GeneDBManager geneDBManager = new GeneDBManager(dbConnector);
						ids = geneDBManager.getAllEnsemblIds();
					} else {
						abort("error_execute_idconverter", "Missing input file", "Missing input file", "Missing input file");
					}
				}
			} else {
				ids = StringUtils.toList(inputIds, ",");
			}

			if ( ids == null || ids.size() == 0) {
				abort("error_execute_idconverter", "No IDs found", "No IDs found", "No IDs found");
			}

			// go bp
			if(commandLine.hasOption("go-bp")) {
				parseGODb(commandLine,"bp");
			}
			// go cc
			if(commandLine.hasOption("go-cc")) {
				parseGODb(commandLine,"cc");
			}
			// go mf
			if(commandLine.hasOption("go-mf")) {
				parseGODb(commandLine,"mf");
			}
			// kegg
			if(commandLine.hasOption("kegg")) {			
				KeggFilter keggFilter = new KeggFilter(Integer.parseInt(commandLine.getOptionValue("kegg-min-num-genes","1")),Integer.parseInt(commandLine.getOptionValue("kegg-max-num-genes","500")));
				filterList.add(keggFilter);
			}
			// biocarta
			if(commandLine.hasOption("biocarta")) {
				BiocartaFilter biocartaFilter = new BiocartaFilter(Integer.parseInt(commandLine.getOptionValue("biocarta-min-num-genes","1")),Integer.parseInt(commandLine.getOptionValue("biocarta-max-num-genes","500")));
				filterList.add(biocartaFilter);
			}

			String outputFormat = commandLine.getOptionValue("output-format", "compact");			

			FeatureList<AnnotationItem> al = null;
			AnnotationDBManager af = new AnnotationDBManager(dbConnector);

			//System.out.println("dbConnector = " + dbConnector.toString());
			//System.out.println("ids = " + ListUtils.toString(ids, ","));

			String fileName = null;

			for(Filter filter: filterList) {
				al = null;
				if ( filter instanceof GOFilter ) {
					//al = af.getGOAnnotation(ids);
					al = af.getGOAnnotation(ids, (GOFilter) filter);
					System.out.println("GO Filter: " + ((GOFilter) filter).getNamespace());
					fileName = ((GOFilter) filter).getNamespace();
					// ((GOFilter) filter).getNamespace()
				} else if ( filter instanceof KeggFilter ) {
					al = af.getKeggAnnotation(ids, (KeggFilter) filter);
					System.out.println("Kegg Filter");
					fileName = "kegg";
					// ((GOFilter) filter).getNamespace()
				} else if ( filter instanceof BiocartaFilter ) {
					al = af.getBiocartaAnnotation(ids, (BiocartaFilter) filter);
					System.out.println("Biocarta Filter");
					fileName = "biocarta";
					// ((GOFilter) filter).getNamespace()
				}
				if ( al != null ) {
					System.out.println("size of list : " + al.size() + ", saving in file: " + fileName);
					if ( outputFormat.equalsIgnoreCase("extended") ) {
						IOUtils.write(new File(outdir + "/" + fileName + ".txt"), al.toString());
						System.out.println(al.toString());
					} else {
						Map<String, List<String>> map = new HashMap<String, List<String>>();
						for(int i=0 ; i<al.size() ; i++) {
							if ( ! map.containsKey(al.get(i).getId()) ) {
								map.put(al.get(i).getId(), new ArrayList<String>()); 
							}
							map.get(al.get(i).getId()).add(al.get(i).getFunctionalTermId());
						}
						StringBuilder sb = new StringBuilder();
						for(String key: MapUtils.getKeys(map)) {
							sb.append(key).append("\t").append(ListUtils.toString(map.get(key), ",")).append("\n");
						}
						IOUtils.write(new File(outdir + "/" + fileName + ".txt"), sb.toString());
						System.out.println(sb.toString());
					}
					
					result.addOutputItem(new Item(fileName, fileName + ".txt", "Database name: filename", Item.TYPE.FILE, new ArrayList<String>(), new HashMap<String,String>(), "Results"));					
				}
			}
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	public void parseGODb(CommandLine cmdLine, String namespace){
		if(cmdLine.hasOption("go-" + namespace + "")) {
			if(cmdLine.hasOption("go-" + namespace + "-inclusive")) {				
				int min = Integer.parseInt(cmdLine.getOptionValue("go-" + namespace + "-min-level","5"));
				int max = Integer.parseInt(cmdLine.getOptionValue("go-" + namespace + "-min-level","12"));
				for(int i=min; i<=max; i++){
					GOFilter goBpFilter = parseGOFilter(cmdLine, namespace);
					goBpFilter.setMinLevel(i);
					goBpFilter.setMaxLevel(i);					
					filterList.add(goBpFilter);	
				}
			} else {
				GOFilter goBpFilter = parseGOFilter(cmdLine, namespace);				
				filterList.add(goBpFilter);
			}			
		}		
	}

	public GOFilter parseGOFilter(CommandLine cmdLine, String namespace){
		String infraredNamespace = "biological_process";		
		if(namespace.equalsIgnoreCase("cc")) infraredNamespace = "cellular_component";
		if(namespace.equalsIgnoreCase("mf")) infraredNamespace = "molecular_function";
		GOFilter goFilter = new GOFilter(infraredNamespace);
		goFilter.setMinLevel(Integer.parseInt(cmdLine.getOptionValue("go-" + namespace + "-min-level","5")));
		goFilter.setMaxLevel(Integer.parseInt(cmdLine.getOptionValue("go-" + namespace + "-max-level","12")));
		goFilter.setMinNumberGenes(Integer.parseInt(cmdLine.getOptionValue("go-" + namespace + "-min-num-genes","1")));	
		//goFilter.setMinNumberGenes(Integer.parseInt(cmdLine.getOptionValue("go-" + namespace + "-min-num-genes","5")));	
		goFilter.setMaxNumberGenes(Integer.parseInt(cmdLine.getOptionValue("go-" + namespace + "-max-num-genes","500")));
		//goFilter.addKeywords(arg0);
		goFilter.setLogicalOperator(cmdLine.getOptionValue("go-" + namespace + "-keywords-logic","AND"));
		return goFilter;
	}


	protected String getDBName(Filter filter){
		String name = StringUtils.randomString();
		if(filter instanceof GOFilter) {						
			GOFilter goFilter = (GOFilter) filter.clone();
			name = "go_" + goFilter.getNamespace() + "_" + goFilter.getMinLevel() + "_" + goFilter.getMaxLevel();
		} else if(filter instanceof KeggFilter) {
			name = "kegg";
		}
		else if(filter instanceof BiocartaFilter) {
			name = "biocarta";
		}
		return name;
	}


	protected String getDBTitle(Filter filter){
		String title = "Untitled",levels;
		if(filter instanceof GOFilter) {						
			GOFilter goFilter = (GOFilter) filter.clone();
			if(goFilter.getMinLevel()==goFilter.getMaxLevel()) {
				levels = "(level " + goFilter.getMinLevel() + ")";
			} else{
				levels = "(levels from " + goFilter.getMinLevel() + " to " + goFilter.getMaxLevel() + ")"; 
			}						
			title = "GO " + goFilter.getNamespace().replace("_", " ") + " " + levels;
		} else if(filter instanceof KeggFilter) {
			title = "Kegg";
		} else if(filter instanceof BiocartaFilter) {
			title = "Biocarta";
		}
		return title;
	}


}
