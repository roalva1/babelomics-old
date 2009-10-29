package org.bioinfo.babelomics.tools.functional.tissues;

import java.io.File;
import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.bioinfo.babelomics.methods.functional.InfraredUtils;
import org.bioinfo.collections.matrix.DataFrame;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ArrayUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.commons.utils.MapUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.dataset.FeatureData;
import org.bioinfo.db.DBConnection;
import org.bioinfo.db.api.PreparedQuery;
import org.bioinfo.db.api.Query;
import org.bioinfo.db.handler.BeanArrayListHandler;
import org.bioinfo.db.handler.MatrixHandler;
import org.bioinfo.db.handler.ResultSetHandler;
import org.bioinfo.infrared.common.dbsql.DBConnector;
import org.bioinfo.math.data.DoubleMatrix;
import org.bioinfo.math.result.TTestResult;
import org.bioinfo.math.result.TestResultList;
import org.bioinfo.math.stats.MultipleTestCorrection;
import org.bioinfo.math.stats.inference.TTest;
import org.bioinfo.math.util.MathUtils;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;

public class SageTmt extends Tmt {

	// best place for these paramteres, a config file
	//
	protected String dbDriver = "mysql";
	protected String dbHost = "mem20";
	protected String dbPort = "3306"; 
	protected String dbUser = "ensembl_user";
	protected String dbPass = "ensembl_user";
	protected String dbName;

	//	private Map<String, String> ensemblIdtoGene1 = new HashMap<String, String>(), ensemblIdtoGene2 = new HashMap<String, String>();
	//	private int replicated1 = 0, replicated2 = 0;
	//	private List<String> genesWithoutEnsemblIds1 = new ArrayList<String>(), genesWithoutEnsemblIds2 = new ArrayList<String>();
	//	private Map<String, List<String>> genesWithMultipleEnsemblIds1 = new HashMap<String, List<String>>(), genesWithMultipleEnsemblIds2 = new HashMap<String, List<String>>();
	//	
	//	private Map<String, List<String>> null, probesMap2 = null;
	//	private List<String> genesWithoutProbes1 = null, genesWithoutProbes2 = null;

	public SageTmt() {
	}

	public void initOptions() {
		super.initOptions();
		options.addOption(OptionFactory.createOption("tag-type", "Type of SAGE tags, valid values are 'short' and 'long'", false));
		options.addOption(OptionFactory.createOption("exclude-cell-lines", "Exclude cell lines", false,false));
		options.addOption(OptionFactory.createOption("min-tags", "Minimum number of tags. Defalut value: 5000", false));
		options.addOption(OptionFactory.createOption("perc-null-libraries", "Percentage for null values accepted in libraries (ranging from 0 to 100). Default value: 80", false));
		options.addOption(OptionFactory.createOption("perc-null-genes", "Percentage for null values accepted in genes (ranging from 0 to 100). Default value: 80", false));
	}

	public void execute() {
		File f1 = new File(commandLine.getOptionValue("list1"));
		File f2 = commandLine.hasOption("list2") ? new File(commandLine.getOptionValue("list2")) :  null;
 
		boolean areEnsemblIds = commandLine.hasOption("ensembl-ids");
//		String organism = commandLine.getOptionValue("organism");
		String normMethod = commandLine.getOptionValue("normalization-method", "mas5");
		String multipleProbes = commandLine.getOptionValue("multiple-probes", "mean");

		List<String> tissues = StringUtils.toList(commandLine.getOptionValue("tissues"), ",");
		
		dbName = "mouse".equalsIgnoreCase(species) ? "gnf_mouse" : "gnf_human";

		try {			
			if ( tissues.contains("all tissues") ) {
				tissues = getAllTissues(species);
			}

			System.out.println("tissues:\n" + ListUtils.toString(tissues));

			DBConnector dbConnector = new DBConnector("mouse".equalsIgnoreCase(species) ? "mmu" : "hsa");
			System.out.println("db connector = " + dbConnector.toString());

			// reading data
			//
			jobStatus.addStatusMessage("20", "reading data");
			logger.debug("reading data...\n");

			// handling list #1
			//
			List<String> geneList1, uniqueGeneList1, dupGeneList1 = null;
			List<String> ensemblList1, noConverted1 = null;
			Map<String, List<String>> ensemblMap1 = null;
			geneList1 = IOUtils.readLines(f1);
			uniqueGeneList1 = ListUtils.unique(geneList1);
			if ( geneList1.size() != uniqueGeneList1.size() ) {
				dupGeneList1 = ListUtils.duplicated(geneList1);
				logger.debug("removing " + dupGeneList1.size() + " duplicated genes from List #1: " + ListUtils.toString(dupGeneList1, ","));
			}
			if ( ! areEnsemblIds ) {
				ensemblMap1 = InfraredUtils.getEnsemblMap(dbConnector, uniqueGeneList1);
				if ( ensemblMap1 == null || ensemblMap1.size() == 0 ) {
					throw new Exception("No Ensembl IDs found when converting your gene list #1 to Ensembl IDs");
				}
				noConverted1 = new ArrayList<String> ();
				ensemblList1 = new ArrayList<String> ();
				for(String key: ensemblMap1.keySet()) {
					if ( ensemblMap1.get(key) != null && ensemblMap1.get(key).size() > 0 ) {
						ensemblList1.addAll(ensemblMap1.get(key));
					} else {
						noConverted1.add(key);
					}
				}
				uniqueGeneList1 = ListUtils.unique(ensemblList1);
			}

			Map<String, List<String>> probesMap1 = getProbes(species, uniqueGeneList1);
			if ( probesMap1 == null || probesMap1.size() == 0 ) {
				throw new Exception("No Affymetrix probes found for gene list #1");
			}

			//			Map<String, List<String>> geneMap1, geneMap2;
			//			
			//			geneMap1 = InfraredUtils.getEnsemblMap(dbConnector, StringUtils.stringToList(IOUtils.toString(f1)));
			//			geneMap1 = cleanGeneMap(geneMap1, true);
			//			geneList1 = createListFromMap(geneMap1);
			//
			//			for(String key: geneMap1.keySet()) {
			//				ensemblIdtoGene1.put(geneMap1.get(key).get(0), key);
			//			}
			//			
			//			replicated1 = geneList1.size();
			//			geneList1 = ListUtils.unique(geneList1);
			//			replicated1 -= geneList1.size();


			// handling list #2
			//
			List<String> geneList2, uniqueGeneList2, dupGeneList2 = null;
			List<String> ensemblList2, noConverted2 = null;
			Map<String, List<String>> ensemblMap2 = null;

			if ( f2 == null ) {
				List<String> genes = getAllGenes(species);
				
				uniqueGeneList2 = new ArrayList<String>();
				for (String gene: genes) {
					if ( !uniqueGeneList1.contains(gene) ) {
						uniqueGeneList2.add(gene);
					}
				}
			} else {
				geneList2 = IOUtils.readLines(f2);
				uniqueGeneList2 = ListUtils.unique(geneList2);
				if ( geneList2.size() != uniqueGeneList2.size() ) {
					dupGeneList2 = ListUtils.duplicated(geneList2);
					logger.debug("removing " + dupGeneList1.size() + " duplicated genes from List #2: " + ListUtils.toString(dupGeneList2, ","));
				}
				if ( ! areEnsemblIds ) {
					ensemblMap2 = InfraredUtils.getEnsemblMap(dbConnector, uniqueGeneList2);
					if ( ensemblMap2 == null || ensemblMap2.size() == 0 ) {
						throw new Exception("No Ensembl IDs found when converting your gene list #2 to Ensembl IDs");
					}
					noConverted2 = new ArrayList<String> ();
					ensemblList2 = new ArrayList<String> ();
					for(String key: ensemblMap2.keySet()) {
						if ( ensemblMap2.get(key) != null && ensemblMap2.get(key).size() > 0 ) {
							ensemblList2.addAll(ensemblMap2.get(key));
						} else {
							noConverted2.add(key);
						}
					}
					uniqueGeneList2 = ListUtils.unique(ensemblList2);
				}
			}
			
			Map<String, List<String>> probesMap2 = getProbes(species, uniqueGeneList2);
			if ( probesMap2 == null || probesMap2.size() == 0 ) {
				throw new Exception("No Affymetrix probes found for gene list #2");
			}

			// getting libraries
			//
			Map<String, String> libraryMap = getTissueLibraries(species, tissues);
			List<String> libraryNames = MapUtils.getKeys(libraryMap);

			// getting frequency matrixes
			//
			jobStatus.addStatusMessage("40", "getting frequency matrixes");
			logger.debug("getting frequency matrixes...\n");

			DoubleMatrix matrix1 = getFreqMatrix(uniqueGeneList1, libraryNames, libraryMap, probesMap1, normMethod, multipleProbes);
			DoubleMatrix matrix2 = getFreqMatrix(uniqueGeneList2, libraryNames, libraryMap, probesMap2, normMethod, multipleProbes);

			// computing t-test
			//
			jobStatus.addStatusMessage("60", "computing t-test");
			logger.debug("computing t-test...\n");

			TTest tTest = new TTest();
			TestResultList<TTestResult> res = tTest.tTest(matrix1, matrix2);				
			MultipleTestCorrection.BHCorrection(res);

//			int[] rowOrder = ListUtils.order(ListUtils.toList(res.getStatistics()), true);
			int[] rowOrder = ArrayUtils.order(res.getStatistics(), true);

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
			featureData.save(new File(outdir + "/tmt.txt"));
			result.addOutputItem(new Item("tmt_file", "tmt.txt", "TMT output file:", TYPE.FILE));
			
			if ( dupGeneList1 != null && dupGeneList1.size() > 0 ) {
				IOUtils.write(new File(outdir + "/duplicated_list1.txt"), dupGeneList1);
				result.addOutputItem(new Item("duplicated_list1_file", "duplicated_list1.txt", "Duplicated genes from list #1 (they were not taken into accout)", TYPE.FILE));
			}
			if ( dupGeneList2 != null && dupGeneList2.size() > 0 ) {
				IOUtils.write(new File(outdir + "/duplicated_list2.txt"), dupGeneList2);
				result.addOutputItem(new Item("duplicated_list2_file", "duplicated_list2.txt", "Duplicated genes from list #2 (they were not taken into accout)", TYPE.FILE));
			}

			if ( noConverted1 != null && noConverted1.size() > 0 ) {
				IOUtils.write(new File(outdir + "/no_converted_list1.txt"), noConverted1);
				result.addOutputItem(new Item("no_converted_list1_file", "no_converted_list1.txt", "Not found Ensembl IDs for the following genes of list #1", TYPE.FILE));
			}
			if ( noConverted2 != null && noConverted2.size() > 0 ) {
				IOUtils.write(new File(outdir + "/no_converted_list2.txt"), noConverted2);
				result.addOutputItem(new Item("no_converted_list2_file", "no_converted_list2.txt", "Not found Ensembl IDs for the following genes of list #2", TYPE.FILE));
			}
			
			
			//			FileUtils.writeStringToFile(new File(getOutdir() + "/replicates.txt"), "Removed " + replicated1 + " genes from list 1\nRemoved " + replicated2 + " genes from list 2");

//			writeProbesMap(new File(getOutdir() + "/gene_probes_1.txt"), probesMap1, ensemblIdtoGene1);
//			writeProbesMap(new File(getOutdir() + "/gene_probes_2.txt"), probesMap2, ensemblIdtoGene1);
//
//			if ( genesWithoutEnsemblIds1 != null ) writeGeneList(new File(getOutdir() + "/genes_without_ensemblid_1.txt"), genesWithoutEnsemblIds1, "#Genes from list 1 without ensembl ID\n");
//			if ( genesWithoutEnsemblIds2 != null ) writeGeneList(new File(getOutdir() + "/genes_without_ensemblid_2.txt"), genesWithoutEnsemblIds2, "#Genes from list 2 without ensembl ID\n");
//
//			if ( genesWithoutProbes1 != null ) writeGeneList(new File(getOutdir() + "/genes_without_probes_1.txt"), genesWithoutProbes1, "#Genes from list 1 without probes\n");
//			if ( genesWithoutProbes2 != null ) writeGeneList(new File(getOutdir() + "/genes_without_probes_2.txt"), genesWithoutProbes2, "#Genes from list 2 without probes\n");
//
//			if ( genesWithMultipleEnsemblIds1 != null ) writeGeneList(new File(getOutdir() + "/genes_with_multiple_ensemblid_1.txt"), genesWithMultipleEnsemblIds1, "#Genes from list 1 with multiple Ensembl IDs, they have been removed from the analysis\n");
//			if ( genesWithMultipleEnsemblIds1 != null ) writeGeneList(new File(getOutdir() + "/genes_with_multiple_ensemblid_2.txt"), genesWithMultipleEnsemblIds2, "#Genes from list 2 with multiple Ensembl IDs, they have been removed from the analysis\n");


		} catch (Exception e) {
			abort("exception_execute_tmt", "execute tmt error", e.getMessage(), StringUtils.getStackTrace(e));
		}
	}


//	private void writeProbesMap(File file, Map<String, List<String>> probes, Map<String, String> genes) throws IOException {
//		StringBuilder sb = new StringBuilder();
//		sb.append("#gene").append("\tensembl ID\tprobes\n");
//		for(String key: probes.keySet()) {
//			sb.append(genes.get(key)).append("\t").append(key).append("\t").append(ListUtils.toString(probes.get(key), " ")).append("\n");
//		}
//		//		FileUtils.writeStringToFile(file, sb.toString());
//	}

//	private void writeGeneList(File file, List<String> list, String msg) throws IOException {
//		if ( list != null && list.size() > 0 ) {
//			StringBuilder sb = new StringBuilder();
//			sb.append(msg);
//			for(String item: list) {
//				sb.append(item).append("\n");
//			}
//			//			FileUtils.writeStringToFile(file, sb.toString());
//		}
//	}

//	private void writeGeneList(File file, Map<String, List<String>> map, String msg) throws IOException {
//		if ( map != null && map.size() > 0 ) {
//			StringBuilder sb = new StringBuilder();
//			sb.append(msg);
//			for(String key: map.keySet()) {
//				sb.append(key).append("\t").append(ListUtils.toString(map.get(key), " ")).append("\n");
//			}
//			//			FileUtils.writeStringToFile(file, sb.toString());
//		}
//	}


//	private List<String> createListFromMap(Map<String, List<String>> map) {
//		List<String> list = new ArrayList<String> ();
//		for(String key: map.keySet()) {
//			list.add(map.get(key).get(0));
//		}
//		return list;
//	}

//	private Map<String, List<String>> cleanGeneMap(Map<String, List<String>> geneMap, boolean firstList) {
//		List<String> list;
//		Map<String, List<String>> map = new HashMap<String, List<String>> ();
//		for(String key: geneMap.keySet()) {
//			list = geneMap.get(key);
//			if ( list == null ) {
//				if ( firstList ) {
//					genesWithoutEnsemblIds1.add(key);
//				} else {
//					genesWithoutEnsemblIds2.add(key);					
//				}
//			} else if ( list.size() > 1 ) {
//				if ( firstList ) {
//					genesWithMultipleEnsemblIds1.put(key, list);					
//				} else {
//					genesWithMultipleEnsemblIds2.put(key, list);					
//				}
//			} else {
//				map.put(key, geneMap.get(key));
//			}
//		}
//		return map;
//	}

	private DoubleMatrix getFreqMatrix(List<String> geneList, List<String> libraryNames, Map<String, String> libraryMap, Map<String, List<String>> probeMap, String normMethod, String freqMethod) {

		List<String> genes = new ArrayList<String>();
		for(String gene: geneList) {
			if ( probeMap.get(gene) != null ) {
				genes.add(gene);
			}
		}

		DoubleMatrix matrix = new DoubleMatrix(libraryNames.size(), genes.size());

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

		double[] array = ArrayUtils.toDoubleArray(values);

		if ( "mean".equalsIgnoreCase(freqMethod) ) {
			return MathUtils.mean(array);
		} else if ( "max".equalsIgnoreCase(freqMethod) ) {
			return ArrayUtils.max(array);
		} else if ( "min".equalsIgnoreCase(freqMethod) ) {
			return ArrayUtils.min(array);
		} else if ( "25".equalsIgnoreCase(freqMethod) ) {
			return MathUtils.percentile(array, 25);
		} else if ( "50".equalsIgnoreCase(freqMethod) ) {
			return MathUtils.percentile(array, 50);
		} else if ( "75".equalsIgnoreCase(freqMethod) ) {			
			return MathUtils.percentile(array, 75);
		}

		return Double.NaN;
	}

	
	public List<String> getAllTissues(String organism) throws SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException {
		List<String> tissues;


		DBConnection dbConn = new DBConnection(dbDriver, dbHost, dbPort, dbName, dbUser, dbPass);

		ResultSetHandler rsh = new BeanArrayListHandler(String.class);
		String q = "select distinct tissue_name from tissue_info order by tissue_name asc"; 
		Query query;
		query = dbConn.createSQLQuery(q);
		tissues = (ArrayList<String>) query.execute(rsh);

		return tissues;
	}

	public Map<String, List<String>> getProbes(String organism, List<String> geneList) throws SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException {
		Map<String, List<String>> probeMap = new HashMap<String, List<String>>();

		DBConnection dbConn = new DBConnection(dbDriver, dbHost, dbPort, dbName, dbUser, dbPass);

		String q;
		PreparedQuery query;
		List<String> probes;
		ResultSetHandler rsh = new BeanArrayListHandler(String.class);
		for(String gene: geneList) {
			q = "select probeset_id from ensemblGene where ensemblGene_id = ?"; 
			//			System.out.println("dbName = " + dbName + ", query = " + q);
			query = dbConn.createSQLPrepQuery(q);
			query.setParams(gene);
			probes = (ArrayList<String>) query.execute(rsh);
			if ( probes != null && probes.size() > 0 ) {
				probeMap.put(gene, probes);
			}
		}

		return probeMap;
	}

	public List<String> getAllGenes(String organism) throws SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException {
		List<String> genes;

		String dbDriver = "mysql";
		String dbHost = "mem20";
		String dbPort = "3306"; 
		String dbUser = "ensembl_user";
		String dbPass = "ensembl_user";
		String dbName = "mouse".equalsIgnoreCase(organism) ? "gnf_mouse" : "gnf_human";

		DBConnection dbConn = new DBConnection(dbDriver, dbHost, dbPort, dbName, dbUser, dbPass);

		ResultSetHandler rsh = new BeanArrayListHandler(String.class);
		String q = "select distinct ensemblGene_id from ensemblGene"; 
		Query query = dbConn.createSQLQuery(q);
		genes = (ArrayList<String>) query.execute(rsh);

		return genes;
	}

	public Map<String, String> getTissueLibraries(String organism, List<String> tissues) throws SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException {
		Map<String, String> libraryMap = new HashMap<String, String>();

		DBConnection dbConn = new DBConnection(dbDriver, dbHost, dbPort, dbName, dbUser, dbPass);

		String q;
		PreparedQuery query;
		List<Library> libraries;
		ResultSetHandler rsh = new BeanArrayListHandler(Library.class);
		for(String tissue: tissues) {
			q = "select distinct tissue_id, tissue_name from tissue_info where tissue_name = ?"; 
			query = dbConn.createSQLPrepQuery(q);
			query.setParams(tissue);
			libraries = (ArrayList<Library>) query.execute(rsh);
			for(Library library: libraries) {
				libraryMap.put(library.name, library.id);
			}
		}

		return libraryMap;
	}

}
