package org.bioinfo.babelomics.tools.preprocessing;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.commons.utils.MapUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.infrared.common.DBConnector;
import org.bioinfo.infrared.core.XRefDBManager;
import org.bioinfo.infrared.core.common.FeatureList;
import org.bioinfo.infrared.core.feature.XRef;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;

public class IDConverter extends BabelomicsTool {

	public IDConverter() {
		initOptions();
	}

	@Override
	public void initOptions() {
		getOptions().addOption(OptionFactory.createOption("listfile", "File containning the IDs to convert", false));
		getOptions().addOption(OptionFactory.createOption("list", "IDs to convert separated by commas", false));

		getOptions().addOption(OptionFactory.createOption("db-names", "DB names to convert", false));
		getOptions().addOption(OptionFactory.createOption("output-format", "Output format: compact or extended. Default: compact", false));		
	}

	@Override
	protected void execute() {
		//try {	
		List<String> ids = null;
		List<String> dbNames = new ArrayList<String> ();

		// check input IDs
		//
		String inputIds = commandLine.getOptionValue("list", null);			
		if ( inputIds == null ) {
			String fileName = commandLine.getOptionValue("listfile", null);
			if ( fileName != null ) {
				try {
					ids = IOUtils.column(new File(fileName), 0);
				} catch (IOException e) {
					abort("ioexception_execute_idconverter", e.getMessage(), e.getMessage(), e.getMessage());
				}
			} else {
				abort("error_execute_idconverter", "Missing input file", "Missing input file", "Missing input file");					
			}
		} else {
			ids = StringUtils.toList(inputIds, ",");
		}

		if ( ids == null || ids.size() == 0) {
			abort("error_execute_idconverter", "No IDs found", "No IDs found", "No IDs found");
		}

		// check output refs
		//
		if ( commandLine.hasOption("db-names") ) {
			dbNames = StringUtils.toList(commandLine.getOptionValue("db-names"), ","); 
		}

		String outputFormat = commandLine.getOptionValue("output-format", "compact");			

		
		if ( dbNames == null || dbNames.size() == 0 ) {
			abort("error_execute_idconverter", "Missing DB names", "Missing DB names", "Missing DB names");
		}

		try {
			DBConnector dbConnector = new DBConnector(species, new File(babelomicsHomePath + "/conf/infrared.properties"));
			
			StringBuilder line = new StringBuilder();
			List<String> alls = new ArrayList<String>();
			List<String> outIds = new ArrayList<String>(); 
			Map<String, List<String>> idsMap = new HashMap<String, List<String>>(); 
			
			System.out.println("-> infrared config = " + babelomicsHomePath + "/conf/infrared.conf");
			
			System.out.println("-> db connector = " + dbConnector.toString());
			System.out.println("-> db names = " + ListUtils.toString(dbNames, ","));
			List<Map<String, FeatureList<XRef>>> list = new XRefDBManager(dbConnector).getListByDBNames(ids, dbNames);

			if ( list != null && list.size() > 0 ) {
				
				line.append("#NAMES\t").append(ListUtils.toString(MapUtils.getKeys(list.get(0)), "\t"));
				alls.add(line.toString());
				
				String idFileName;
				List<String> idLines = new ArrayList<String>();
				for(int i=0 ; i<list.size() ; i++) {
										
					line.delete(0, line.length());
					line.append(ids.get(i));
					//System.out.println("Converting " + ids.get(i));
					for(String key: MapUtils.getKeys(list.get(i))) {
						
						
						//System.out.print("to " + key + " ---> ");
						outIds.clear();
						for(XRef xref: list.get(i).get(key)) {
							outIds.add(xref.getId());
						}
						
						if ( ! idsMap.containsKey(key) ) {
							idsMap.put(key, new ArrayList<String>());
							idsMap.get(key).add("#NAMES\t" + key);
						}
						idsMap.get(key).add(ids.get(i) + "\t" + ListUtils.toString(outIds,","));
						line.append("\t").append(ListUtils.toString(outIds,","));
					}
					alls.add(line.toString());
				}
				
				// save results
				//
				String fileName = "ids_summary.txt";
				IOUtils.write(new File(outdir + "/" + fileName), alls);
				result.addOutputItem(new Item("all", fileName, "Conversion (summary file)", Item.TYPE.FILE, new ArrayList<String>(), new HashMap<String,String>(), "Output files"));
				
				String[] arr;
				List<String> items = null;
				List<String> lines = new ArrayList<String>();
				
				for(String key: MapUtils.getKeys(idsMap)) {
										
					outIds.clear();
					idFileName = key +  ".txt";		
					for(int i=1 ; i<idsMap.get(key).size() ; i++) {
						arr = idsMap.get(key).get(i).split("\t");
						if ( arr.length > 1 ) {
							items = StringUtils.toList(arr[1], ",");
							for(String item: items) {
								outIds.add(item);
							}
						}
					}
					IOUtils.write(new File(outdir + "/" + idFileName), outIds);
					result.addOutputItem(new Item(key + "_ids", idFileName, key + " IDs", Item.TYPE.DATA, StringUtils.toList("idlist", ","), new HashMap<String,String>(), "Output files"));
					result.addOutputItem(new Item(key + "_file_ids", idFileName, key + " IDs", Item.TYPE.FILE, new ArrayList<String>(), new HashMap<String,String>(), "Output files"));

					fileName = key + "_conversion.txt";
															
					if ( outputFormat.equalsIgnoreCase("extended") ) {
						lines.add(idsMap.get(key).get(0));
						for(int i=1 ; i<idsMap.get(key).size() ; i++) {
							arr = idsMap.get(key).get(i).split("\t");
							if ( arr.length > 1 ) {
								items = StringUtils.toList(arr[1], ",");
								for(String item: items) {
									lines.add(arr[0] + "\t" + item);
								}
							} else {
								lines.add(arr[0] + "\t");
							}
						}
						IOUtils.write(new File(outdir + "/" + fileName), ListUtils.toString(lines, "\n"));
					} else {
						IOUtils.write(new File(outdir + "/" + fileName), idsMap.get(key));											
					}

					//result.addOutputItem(new Item(key + " IDs", fileName, "ID conversion", Item.TYPE.FILE, new ArrayList<String>(), new HashMap<String,String>(), "Results per file"));
					
					result.addOutputItem(new Item(key, fileName, "ID conversion table", Item.TYPE.FILE, StringUtils.toList("TABLE,IDCONVERTER_TABLE", ","), new HashMap<String, String>(), "ID conversion tables (" + outputFormat + " format)"));
				}
			}
		} catch (Exception e) {
			this.printError("exception_execute_idconverter", e.toString(), e.toString(), e);
		}
	}

}
