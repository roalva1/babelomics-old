package org.bioinfo.babelomics.tools.functional.tissues;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.babelomics.tools.functional.FunctionalProfilingTool;
import org.bioinfo.chart.BoxPlotChart;
import org.bioinfo.collections.matrix.DataFrame;
import org.bioinfo.data.dataset.FeatureData;
import org.bioinfo.db.DBConnection;
import org.bioinfo.db.Query;
import org.bioinfo.db.handler.BeanArrayListHandler;
import org.bioinfo.db.handler.MatrixHandler;
import org.bioinfo.db.handler.ResultSetHandler;
import org.bioinfo.infrared.common.dbsql.DBConnector;
import org.bioinfo.io.FileSystemUtils;
import org.bioinfo.math.MathUtils;
import org.bioinfo.math.data.DoubleMatrix;
import org.bioinfo.math.result.TTestResult;
import org.bioinfo.math.result.TestResultList;
import org.bioinfo.math.stats.inference.MultipleTestCorrection;
import org.bioinfo.math.stats.inference.TTest;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;
import org.bioinfo.utils.ListUtils;
import org.bioinfo.utils.StringUtils;
import org.jfree.chart.plot.PlotOrientation;

public class TMT  extends BabelomicsTool {

	private String dbDriver = "mysql";
	private String dbHost = "mem20";
	private String dbPort = "3306"; 
	private String dbName;
	private String dbUser = "ensembl_user";
	private String dbPass = "ensembl_user";
	
	private Map<String, String> ensemblIdtoGene1 = new HashMap<String, String>(), ensemblIdtoGene2 = new HashMap<String, String>();
	private int replicated1 = 0, replicated2 = 0;
	private List<String> genesWithoutEnsemblIds1 = new ArrayList<String>(), genesWithoutEnsemblIds2 = new ArrayList<String>();
	private Map<String, List<String>> genesWithMultipleEnsemblIds1 = new HashMap<String, List<String>>(), genesWithMultipleEnsemblIds2 = new HashMap<String, List<String>>();
	private Map<String, List<String>> probesMap1 = null, probesMap2 = null;
	private List<String> genesWithoutProbes1 = null, genesWithoutProbes2 = null;

	public TMT() {
		initOptions();
	}

	@Override
	public void initOptions() {
		options.addOption(OptionFactory.createOption("organism", "Organism, valid values are 'human' and 'mouse'"));
		options.addOption(OptionFactory.createOption("platform", "Valid values are 'sage' (for SAGE Tags) and 'affy' (for Microarray Affymetrix expression data)"));
		options.addOption(OptionFactory.createOption("tag-type", "Type of tags for SAGE tags, valid values are 'short' and 'long'", false));
		options.addOption(OptionFactory.createOption("list1", "the feature data containig the list #1 of genes, or the feature data file"));
		options.addOption(OptionFactory.createOption("list2", "the feature data containig the list #2 of genes, or the feature data file", false));
		options.addOption(OptionFactory.createOption("tissues", "the list of tissues separated by commas. Enter 'all tissues' to take into account all available tissues"));
		options.addOption(OptionFactory.createOption("normalization-method", "Normalization method, valid values are 'mas5' and 'gcrma'. Defalut value: 'mas5'", false));
		options.addOption(OptionFactory.createOption("multiple-probes", "Multiple probes expression value, valid values are 'mean', 'greatest', 'lowest', 'perc25', 'perc50' and 'perc75'. Default value: 'mean'", false));
	}

	@Override
	public void execute() {
		//				try {
		//					FeatureData fd1 = new FeatureData(new File(commandLine.getOptionValue("list1")));
		//					FeatureData fd2 = new FeatureData(new File(commandLine.getOptionValue("list2")));

		File f1 = new File(commandLine.getOptionValue("list1"));
		File f2 = commandLine.hasOption("list2") ? new File(commandLine.getOptionValue("list2")) :  null;

		String platform = commandLine.getOptionValue("platform");
		String organism = commandLine.getOptionValue("organism");
		String tagType = commandLine.getOptionValue("tag-type", "short");
		String normMethod = commandLine.getOptionValue("normalization-method", "mas5");
		String multipleProbes = commandLine.getOptionValue("multiple-probes", "mean");

		List<String> tissues = StringUtils.stringToList(",", commandLine.getOptionValue("tissues"));
		if ( tissues.contains("all tissues") ) {
			tissues = "affy".equalsIgnoreCase(platform) ? get_affy_tissues(organism) : null;
		}

		System.out.println("tissues:\n" + ListUtils.toString(tissues));

		executeTMT(organism, f1, f2, platform, tagType, tissues, normMethod, multipleProbes);

		//				} catch (IOException e) {
		//					logger.error("Error opening the dataset", e.toString());
		//				}
	}

	private void executeTMT(String organism, File f1, File f2, String platform,
			String tagType, List<String> tissues, 
			String normMethod, String multipleProbes) {

		try {
			DBConnector dbConnector = new DBConnector("mouse".equalsIgnoreCase(organism) ? "mus" : "hsa");
			System.out.println("db connector = " + dbConnector.toString());

			// reading data
			//
			jobStatus.addStatusMessage("20", "reading data");
			logger.debug("reading data...\n");

			List<String> geneList1 = new ArrayList<String>(), geneList2 = new ArrayList<String>();

//			System.out.println("gene list 1:\n" + ListUtils.toString(StringUtils.stringToList(FileSystemUtils.fileToString(f1))));
//			System.out.println("gene list 2:\n" + ListUtils.toString(StringUtils.stringToList(FileSystemUtils.fileToString(f2))));

			Map<String, List<String>> geneMap1, geneMap2;
			
			geneMap1 = FunctionalProfilingTool.getEnsemblMap(dbConnector, StringUtils.stringToList(FileSystemUtils.fileToString(f1)));
			geneMap1 = cleanGeneMap(geneMap1, true);
			geneList1 = createListFromMap(geneMap1);

			for(String key: geneMap1.keySet()) {
				ensemblIdtoGene1.put(geneMap1.get(key).get(0), key);
			}
			
			replicated1 = geneList1.size();
			geneList1 = ListUtils.unique(geneList1);
			replicated1 -= geneList1.size();

			probesMap1 = getProbes(organism, platform, geneList1);
			
			if ( f2 == null ) {
				geneList2 = getAllGenes(organism, platform);
				geneMap2 = new HashMap<String, List<String>>();
				for(String gene: geneList2) {
					geneMap2.put(gene, new ArrayList<String>());
					geneMap2.get(gene).add(gene);
				}
			} else {
				geneMap2 = FunctionalProfilingTool.getEnsemblMap(dbConnector, StringUtils.stringToList(FileSystemUtils.fileToString(f2)));
				geneMap2 = cleanGeneMap(geneMap2, false);
				geneList2 = createListFromMap(geneMap2);
			}
				replicated2 = geneList1.size();
				geneList2 = ListUtils.unique(geneList2);
				replicated2 -= geneList1.size();

				for(String key: geneMap2.keySet()) {
					ensemblIdtoGene2.put(geneMap2.get(key).get(0), key);
				}

				probesMap2 = getProbes(organism, platform, geneList2);

			
			
			System.out.println("ensembl gene list 1 (" + geneList1.size() + "):\n" + ListUtils.toString(geneList1));
			System.out.println("ensembl gene list 2 (" + geneList2.size() + "):\n" + ListUtils.toString(geneList2));

			Map<String, String> libraryMap = getTissueLibraries(organism, tissues, platform);
			// in order to guarantee the order
			//
			List<String> libraryNames = new ArrayList<String> (libraryMap.size());
			for(String key: libraryMap.keySet()) {
				libraryNames.add(key);
			}
			//			System.out.println("library map");
			//			for(String key: libraryMap.keySet()) {
			//				System.out.println(key + " -> " + libraryMap.get(key));
			//			}

			// getting frequency matrixes
			//
			jobStatus.addStatusMessage("40", "getting frequency matrixes");
			logger.debug("getting frequency matrixes...\n");

			DoubleMatrix matrix1 = getFreqMatrix(geneList1, libraryNames, libraryMap, probesMap1, organism, normMethod, multipleProbes);
			DoubleMatrix matrix2 = getFreqMatrix(geneList2, libraryNames, libraryMap, probesMap2, organism, normMethod, multipleProbes);

			System.out.println("**************** matrix 1");
			System.out.print("#names\t");
			for(String gene: geneList1) {
				System.out.print(gene + "\t");
			}
			System.out.println("");
			for(int i=0 ; i<matrix1.getRowDimension() ; i++) {
				System.out.print(libraryNames.get(i) + "\t");
				for(int j=0 ; j<matrix1.getColumnDimension() ; j++) {
					System.out.print(matrix1.get(i, j) + "\t");
				}
				System.out.println("");
			}

			System.out.println("****************** matrix 2");
			System.out.print("#names\t");
			for(String gene: geneList2) {
				System.out.print(gene + "\t");
			}
			System.out.println("");
			for(int i=0 ; i<matrix2.getRowDimension() ; i++) {
				System.out.print(libraryNames.get(i) + "\t");
				for(int j=0 ; j<matrix2.getColumnDimension() ; j++) {
					System.out.print(matrix2.get(i, j) + "\t");
				}
				System.out.println("");
			}

			//			System.out.println("libraries = " + libraryNames.toString());
			//			System.out.println("matrix1\n" + matrix1.toString());
			//			System.out.println("matrix2\n" + matrix2.toString());

			// computing t-test
			//
			jobStatus.addStatusMessage("60", "computing t-test");
			logger.debug("computing t-test...\n");

			TTest tTest = new TTest();
			TestResultList<TTestResult> res = tTest.tTest(matrix1, matrix2);				
			MultipleTestCorrection.BHCorrection(res);

			int[] rowOrder = ListUtils.order(ListUtils.toList(res.getStatistics()), true);

			// saving data
			//
			jobStatus.addStatusMessage("80", "saving results");
			logger.debug("saving results...");

			DataFrame dataFrame = new DataFrame(libraryNames.size(), 0);
			dataFrame.setRowNames(ListUtils.ordered(libraryNames, rowOrder));

			//dataFrame.addColumn("id", ListUtils.ordered(dataset.getFeatureNames(), rowOrder));
			dataFrame.addColumn("statistic", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getStatistics()), rowOrder)));
			dataFrame.addColumn("p-value", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getPValues()), rowOrder)));
			dataFrame.addColumn("adj. p-value", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getAdjPValues()), rowOrder)));

			FeatureData featureData = new FeatureData(dataFrame);
			featureData.write(new File(getOutdir() + "/tmt.txt"));

//			FileUtils.writeStringToFile(new File(getOutdir() + "/replicates.txt"), "Removed " + replicated1 + " genes from list 1\nRemoved " + replicated2 + " genes from list 2");

			writeProbesMap(new File(getOutdir() + "/gene_probes_1.txt"), probesMap1, ensemblIdtoGene1);
			writeProbesMap(new File(getOutdir() + "/gene_probes_2.txt"), probesMap2, ensemblIdtoGene1);
		
			if ( genesWithoutEnsemblIds1 != null ) writeGeneList(new File(getOutdir() + "/genes_without_ensemblid_1.txt"), genesWithoutEnsemblIds1, "#Genes from list 1 without ensembl ID\n");
			if ( genesWithoutEnsemblIds2 != null ) writeGeneList(new File(getOutdir() + "/genes_without_ensemblid_2.txt"), genesWithoutEnsemblIds2, "#Genes from list 2 without ensembl ID\n");

			if ( genesWithoutProbes1 != null ) writeGeneList(new File(getOutdir() + "/genes_without_probes_1.txt"), genesWithoutProbes1, "#Genes from list 1 without probes\n");
			if ( genesWithoutProbes2 != null ) writeGeneList(new File(getOutdir() + "/genes_without_probes_2.txt"), genesWithoutProbes2, "#Genes from list 2 without probes\n");
			
			if ( genesWithMultipleEnsemblIds1 != null ) writeGeneList(new File(getOutdir() + "/genes_with_multiple_ensemblid_1.txt"), genesWithMultipleEnsemblIds1, "#Genes from list 1 with multiple Ensembl IDs, they have been removed from the analysis\n");
			if ( genesWithMultipleEnsemblIds1 != null ) writeGeneList(new File(getOutdir() + "/genes_with_multiple_ensemblid_2.txt"), genesWithMultipleEnsemblIds2, "#Genes from list 2 with multiple Ensembl IDs, they have been removed from the analysis\n");

			result.addOutputItem(new Item("tmt_file", getOutdir() + "/tmt.txt", "The TMT file is: ", TYPE.FILE));

			// done
			//
			jobStatus.addStatusMessage("100", "done");
			logger.debug("tmt done\n");

		} catch (Exception e) {
			printError("exception_execute_tmt", "execute tmt error", e.toString(), e);
		}
	}

	private List<String> getAllGenes(String organism, String platform) {
		List<String> genes;

		dbName = "mouse".equalsIgnoreCase(organism) ? "gnf_mouse" : "gnf_human";

		DBConnection dbConn = new DBConnection(dbDriver, dbHost, dbPort, dbName, dbUser, dbPass);

		ResultSetHandler rsh = new BeanArrayListHandler(String.class);
		String q = "select distinct ensemblGene_id from ensemblGene"; 
		Query query;
		try {
			query = dbConn.createSQLQuery(q);
			genes = (ArrayList<String>) query.execute(rsh);
		} catch (Exception e) {
			genes = null;
			printError("exception_get_all_genes_tmt", "get all genes tmt error", e.toString(), e);
		}

		//		System.out.println("dbName = " + dbName + ", query = " + q);

		return genes;
	}

	private void writeProbesMap(File file, Map<String, List<String>> probes, Map<String, String> genes) throws IOException {
		StringBuilder sb = new StringBuilder();
		sb.append("#gene").append("\tensembl ID\tprobes\n");
		for(String key: probes.keySet()) {
			sb.append(genes.get(key)).append("\t").append(key).append("\t").append(ListUtils.toString(probes.get(key), " ")).append("\n");
		}
//		FileUtils.writeStringToFile(file, sb.toString());
	}

	private void writeGeneList(File file, List<String> list, String msg) throws IOException {
		if ( list != null && list.size() > 0 ) {
			StringBuilder sb = new StringBuilder();
			sb.append(msg);
			for(String item: list) {
				sb.append(item).append("\n");
			}
//			FileUtils.writeStringToFile(file, sb.toString());
		}
	}
	
	private void writeGeneList(File file, Map<String, List<String>> map, String msg) throws IOException {
		if ( map != null && map.size() > 0 ) {
			StringBuilder sb = new StringBuilder();
			sb.append(msg);
			for(String key: map.keySet()) {
				sb.append(key).append("\t").append(ListUtils.toString(map.get(key), " ")).append("\n");
			}
//			FileUtils.writeStringToFile(file, sb.toString());
		}
	}

	
	private List<String> createListFromMap(Map<String, List<String>> map) {
		List<String> list = new ArrayList<String> ();
		for(String key: map.keySet()) {
			list.add(map.get(key).get(0));
		}
		return list;
	}

	private Map<String, List<String>> cleanGeneMap(Map<String, List<String>> geneMap, boolean firstList) {
		List<String> list;
		Map<String, List<String>> map = new HashMap<String, List<String>> ();
		for(String key: geneMap.keySet()) {
			list = geneMap.get(key);
			if ( list == null ) {
				if ( firstList ) {
					genesWithoutEnsemblIds1.add(key);
				} else {
					genesWithoutEnsemblIds2.add(key);					
				}
			} else if ( list.size() > 1 ) {
				if ( firstList ) {
					genesWithMultipleEnsemblIds1.put(key, list);					
				} else {
					genesWithMultipleEnsemblIds2.put(key, list);					
				}
			} else {
				map.put(key, geneMap.get(key));
			}
		}
		return map;
	}

	private DoubleMatrix getFreqMatrix(List<String> geneList, List<String> libraryNames, Map<String, String> libraryMap, Map<String, List<String>> probeMap, String organism, String normMethod, String freqMethod) {

		List<String> genes = new ArrayList<String>();
		for(String gene: geneList) {
			if ( probeMap.get(gene) != null ) {
				genes.add(gene);
			}
		}

		DoubleMatrix matrix = new DoubleMatrix(libraryNames.size(), genes.size());

		dbName = "mouse".equalsIgnoreCase(organism) ? "gnf_mouse" : "gnf_human";

		DBConnection dbConn = new DBConnection(dbDriver, dbHost, dbPort, dbName, dbUser, dbPass);

		double freq;
		List<Double> values = new ArrayList<Double>();

		Query query;
		Object[][] expressions;
		ResultSetHandler rsh = new MatrixHandler();
		String q, prefix = "select expression from " + ("gcrma".equalsIgnoreCase(normMethod) ? "expression_gcRMA" : "expression_MAS5") + " where";

		int total = 0, row = 0, column = 0;
		for(String libraryKey: libraryNames) {
			column = 0;
			for(String gene: genes) {

				//				System.out.println("(row, column) = (" + row + ", " + column + "), (library, gene) = (" + libraryKey + ", " + gene + "), freqMethod = " + freqMethod); 

				values.clear();
				for(String probe: probeMap.get(gene)) {
					q = prefix + " probeset_id = \"" + probe + "\" and tissue = \"" + libraryMap.get(libraryKey) + "\""; 
					//						System.out.println("dbName = " + dbName + ", query = " + q);
					try {
						query = dbConn.createSQLQuery(q);
						expressions = (Object[][]) query.execute(rsh);
						if ( expressions != null && expressions.length > 0 ) {
							total = 0;
							for(int i=0 ; i<expressions.length ; i++) {
								//									System.out.print(expressions[i][0] + " ");
								total += (Integer) expressions[i][0];
							}
							//								System.out.println("mean = " + (1.0 * total / expressions.length));
							values.add(1.0 * total / expressions.length);
						}
					} catch (Exception e) {
						probeMap = null;
						printError("exception_get_probes_tmt", "get probes tmt error", e.toString(), e);
					}
				}
				freq = getFrequency(values, freqMethod); 
				//System.out.println("probes for gene " +  gene + " values = " + ListUtils.toString(values) + ", freq = " + freq);					

				matrix.set(row, column, freq);
				column++;
			}
			row++;
		}			
		return matrix;
	}

	private double getFrequency(List<Double> values, String freqMethod) {

		double[] array = ListUtils.toDoubleArray(values);


		if ( "mean".equalsIgnoreCase(freqMethod) ) {
			return MathUtils.mean(array);
		} else if ( "max".equalsIgnoreCase(freqMethod) ) {
			return MathUtils.max(array);
		} else if ( "min".equalsIgnoreCase(freqMethod) ) {
			return MathUtils.min(array);
		} else if ( "25".equalsIgnoreCase(freqMethod) ) {
			return MathUtils.percentile(array, 25);
		} else if ( "50".equalsIgnoreCase(freqMethod) ) {
			return MathUtils.percentile(array, 50);
		} else if ( "75".equalsIgnoreCase(freqMethod) ) {			
			return MathUtils.percentile(array, 75);
		}

		return Double.NaN;
	}

	private Map<String, List<String>> getProbes(String organism, String library, List<String> geneList) {
		Map<String, List<String>> probeMap = new HashMap<String, List<String>>();

		dbName = "mouse".equalsIgnoreCase(organism) ? "gnf_mouse" : "gnf_human";

		DBConnection dbConn = new DBConnection(dbDriver, dbHost, dbPort, dbName, dbUser, dbPass);

		String q;
		Query query;
		List<String> probes;
		ResultSetHandler rsh = new BeanArrayListHandler(String.class);
		for(String gene: geneList) {
			q = "select probeset_id from ensemblGene where ensemblGene_id = \"" + gene + "\""; 
			//			System.out.println("dbName = " + dbName + ", query = " + q);
			try {
				query = dbConn.createSQLQuery(q);
				probes = (ArrayList<String>) query.execute(rsh);
				if ( probes != null && probes.size() > 0 ) {
					probeMap.put(gene, probes);
				}
			} catch (Exception e) {
				probeMap = null;
				printError("exception_get_probes_tmt", "get probes tmt error", e.toString(), e);
			}
		}

		return probeMap;
	}

	private Map<String, String> getTissueLibraries(String organism, List<String> tissues, String libraryName) {
		Map<String, String> libraryMap = new HashMap<String, String>();

		dbName = "mouse".equalsIgnoreCase(organism) ? "gnf_mouse" : "gnf_human";

		DBConnection dbConn = new DBConnection(dbDriver, dbHost, dbPort, dbName, dbUser, dbPass);

		String q;
		Query query;
		List<Library> libraries;
		ResultSetHandler rsh = new BeanArrayListHandler(Library.class);
		for(String tissue: tissues) {
			q = "select distinct tissue_id, tissue_name from tissue_info where tissue_name = \"" + tissue + "\""; 
			//			System.out.println("dbName = " + dbName + ", query = " + q);
			try {
				query = dbConn.createSQLQuery(q);
				libraries = (ArrayList<Library>) query.execute(rsh);
				for(Library library: libraries) {
					libraryMap.put(library.name, library.id);
				}
			} catch (Exception e) {
				libraryMap = null;
				printError("exception_get_tissue_libraries_tmt", "get tissue libraries tmt error", e.toString(), e);
			}
		}

		return libraryMap;
	}

	/**
	 * 
	 * @param organism
	 * @return
	 */
	private List<String> get_affy_tissues(String organism) {

		List<String> tissues;

		dbName = "mouse".equalsIgnoreCase(organism) ? "gnf_mouse" : "gnf_human";

		DBConnection dbConn = new DBConnection(dbDriver, dbHost, dbPort, dbName, dbUser, dbPass);

		ResultSetHandler rsh = new BeanArrayListHandler(String.class);
		String q = "select distinct tissue_name from tissue_info order by tissue_name asc"; 
		Query query;
		try {
			query = dbConn.createSQLQuery(q);
			tissues = (ArrayList<String>) query.execute(rsh);
		} catch (Exception e) {
			tissues = null;
			printError("exception_get_affy_tissues_tmt", "get affy tissues tmt error", e.toString(), e);
		}

		//		System.out.println("dbName = " + dbName + ", query = " + q);

		return tissues;
	}

	//	private List<String> get_sage_tissues(String organism, String tagType) {
	//		
	//		List<String> tissues;
	//		
	//		dbName = "mouse".equalsIgnoreCase(organism) ? "gnf_mouse" : "gnf_human";
	//		
	//		DBConnection dbConn = new DBConnection(dbDriver, dbHost, dbPort, dbName, dbUser, dbPass);
	//		
	//		ResultSetHandler rsh = new BeanArrayListHandler(String.class);
	//		String q = "select distinct tissue_name from tissue_info order by tissue_name asc"; 
	//		Query query;
	//		try {
	//			query = dbConn.createSQLQuery(q);
	//			tissues = (ArrayList<String>) query.execute(rsh);
	//		} catch (Exception e) {
	//			tissues = null;
	//			printError("exception_get_affy_tissues_tmt", "get affy tissues tmt error", e.toString(), e);
	//		}
	//		
	//		System.out.println("dbName = " + dbName + ", query = " + q);
	//
	//		return tissues;
	//	}

	public BoxPlotChart createBoxplot(String entity, double pValue, List<Double> list1, List<Double> list2) {

		DecimalFormat df = new DecimalFormat("##0.00000");
		BoxPlotChart bpc = new BoxPlotChart(entity.toUpperCase() + "\nadj. p-value = " + df.format(pValue), "", "score");
		bpc.getLegend().setVisible(false);
		bpc.addSeries(list1, "list 1", "list 1");
		bpc.addSeries(list2, "list 2", "list 2");
		bpc.setOrientation(PlotOrientation.HORIZONTAL);

		return bpc;
	}
}
