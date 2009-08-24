package org.bioinfo.babelomics.tools.functional.textmining;



import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.math.linear.MatrixIndexException;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.babelomics.tools.functional.textmining.Score;
import org.bioinfo.chart.BoxPlotChart;
import org.bioinfo.collections.exceptions.InvalidColumnIndexException;
import org.bioinfo.collections.matrix.DataFrame;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.commons.utils.MapUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.dataset.FeatureData;
import org.bioinfo.db.DBConnection;
import org.bioinfo.db.api.Query;
import org.bioinfo.db.handler.BeanArrayListHandler;
import org.bioinfo.db.handler.ResultSetHandler;
import org.bioinfo.math.exception.InvalidParameterException;
import org.bioinfo.math.result.KolmogorovSmirnovTestResult;
import org.bioinfo.math.result.TestResultList;
import org.bioinfo.math.stats.inference.KolmogorovSmirnovTest;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.plot.PlotOrientation;

public class Marmite extends BabelomicsTool {

	public Marmite() {
		initOptions();
	}

	@Override
	public void initOptions() {

		options.addOption(OptionFactory.createOption("list1", "gene list #1"));
		options.addOption(OptionFactory.createOption("list2", "gene list #2"));
		options.addOption(OptionFactory.createOption("bioentity-name", "Valid values: wordroot, disease, chemical"));
		options.addOption(OptionFactory.createOption("bioentity-score-filter", "Minimum number of genes with a score (0-10000)", false));
		options.addOption(OptionFactory.createOption("bioentity-number-filter", "Number of bio-entities in results (0-10000)", false));
	}

	@Override
	public void execute() {
		//		try {
		//			FeatureData fd1 = new FeatureData(new File(commandLine.getOptionValue("list1")));
		//			FeatureData fd2 = new FeatureData(new File(commandLine.getOptionValue("list2")));

		File f1 = new File(commandLine.getOptionValue("list1"));
		File f2 = new File(commandLine.getOptionValue("list2"));

		String bioentity = commandLine.getOptionValue("bioentity-name");
		int scoreFilter = Integer.parseInt(commandLine.getOptionValue("bioentity-score-filter", "5"));
		int numberFilter = Integer.parseInt(commandLine.getOptionValue("bioentity-number-filter", "50"));

		executeMarmite(f1, f2, bioentity, scoreFilter, numberFilter);

		//		} catch (IOException e) {
		//			logger.error("Error opening the dataset", e.toString());
		//		}
	}



	private void executeMarmite(File f1, File f2, String bioentity, int scoreFilter, int numberFilter) {

		long time;

		try {
			// reading data
			//
			jobStatus.addStatusMessage("20", "reading data");
			logger.debug("reading data...\n");

			List<String> geneList1 = ListUtils.unique(StringUtils.stringToList(IOUtils.toString(f1)));
			List<String> geneList2 = ListUtils.unique(StringUtils.stringToList(IOUtils.toString(f2)));

			jobStatus.addStatusMessage("#", "get Entity Map1");
			time = System.currentTimeMillis();
			Map<String, List<Score>> entityMap1 = getEntityMap(geneList1, bioentity);
			jobStatus.addStatusMessage("#", "getEntityMap1, elapsed time = " + ((System.currentTimeMillis()-time)/1000F) + " s");

			jobStatus.addStatusMessage("#", "get Entity Map2");
			time = System.currentTimeMillis();
			Map<String, List<Score>> entityMap2 = getEntityMap(geneList2, bioentity);	
			jobStatus.addStatusMessage("#", "getEntityMap2, elapsed time = " + ((System.currentTimeMillis()-time)/1000F) + " s");

			jobStatus.addStatusMessage("#", "creating the entities list");
			time = System.currentTimeMillis();
			List<String> entities = MapUtils.getKeys(entityMap1);
			entities.addAll(MapUtils.getKeys(entityMap2));
			entities = ListUtils.unique(entities);

			String entity;
			List<Score> scoreList1, scoreList2;
			List<String> validEntities = new ArrayList<String>();
			for(int i=0 ; i<entities.size() ; i++) {
				entity = entities.get(i);
				scoreList1 = entityMap1.get(entity); 
				scoreList2 = entityMap2.get(entity); 
				if ( scoreList1 != null && scoreList1.size() >= scoreFilter && scoreList2 != null && scoreList2.size() >= scoreFilter ) {
					validEntities.add(entity);
				}
			}
			entities.clear();
			jobStatus.addStatusMessage("#", "entity list, elapsed time = " + ((System.currentTimeMillis()-time)/1000F) + " s");

			jobStatus.addStatusMessage("#", "creating score matrixes");
			time = System.currentTimeMillis();
			double[][] scoreMatrix1 = new double[validEntities.size()][];
			double[][] scoreMatrix2 = new double[validEntities.size()][];
			for(int i=0 ; i<validEntities.size() ; i++) {
				entity = validEntities.get(i);
				scoreList1 = entityMap1.get(entity); 
				scoreList2 = entityMap2.get(entity); 

				scoreMatrix1[i] = new double[scoreList1.size()];
				for (int j=0 ; j<scoreList1.size(); j++) {
					scoreMatrix1[i][j] = scoreList1.get(j).getValue();
				}
				scoreMatrix2[i] = new double[scoreList2.size()];
				for (int j=0 ; j<scoreList2.size(); j++) {
					scoreMatrix2[i][j] = scoreList2.get(j).getValue();
				}	
			}								
			jobStatus.addStatusMessage("#", "score matrixes, elapsed time = " + ((System.currentTimeMillis()-time)/1000F) + " s");


			//			System.out.println("scoreMatrix1 size = " + scoreMatrix1.length + ", " + scoreMatrix1[0].length);
			//			System.out.println("scoreMatrix2 size = " + scoreMatrix2.length + ", " + scoreMatrix2[0].length);

			// marmite functional enrichment
			//
			jobStatus.addStatusMessage("40", "computing marmite functional enrichment");
			logger.debug("computing marmite functional enrichment...\n");

			jobStatus.addStatusMessage("#", "computing ks test");
			time = System.currentTimeMillis();
			TestResultList<KolmogorovSmirnovTestResult> res = new KolmogorovSmirnovTest().compute(scoreMatrix1, scoreMatrix2);
			//			System.out.println("kolmogorov result :\n" +  res.toString());
			jobStatus.addStatusMessage("#", "ks test, elapsed time = " + ((System.currentTimeMillis()-time)/1000F) + " s");

			int[] rowOrder = ListUtils.order(ListUtils.toList(res.getAdjPValues()));
			//System.out.println("order = " + ListUtils.toString(ListUtils.toList(rowOrder)));

			DataFrame dataFrame = new DataFrame(validEntities.size(), 0);
			dataFrame.setRowNames(ListUtils.ordered(validEntities, rowOrder));

			//			System.out.println("********* sides = " + res.getSides());

			dataFrame.addColumn("statistic", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getStatistics()), rowOrder)));
			dataFrame.addColumn("side", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getSides()), rowOrder)));
			dataFrame.addColumn("p-value", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getPValues()), rowOrder)));
			dataFrame.addColumn("adj. p-value", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getAdjPValues()), rowOrder)));

			// output file
			String filename = "marmite_output.txt";
			FeatureData featureData = new FeatureData(dataFrame);
			featureData.write(new File(getOutdir() + "/" + filename));
			result.addOutputItem(new Item("marmite_file", filename, "The MARMITE output file", TYPE.FILE));

			// output table
			filename = "marmite_table.txt";
			BufferedWriter bw = new BufferedWriter(new FileWriter(getOutdir() + "/" + filename));
			//System.out.println("marmite table:\n" + dataFrame.toString(true, true));
			bw.write(dataFrame.toString(true, true));
			bw.close();
			Item item = new Item("marmite_table", filename, "The MARMITE output table", TYPE.FILE);
			item.addTag("TABLE");
			result.addOutputItem(item);

			System.out.println("marmite output:\n" +  featureData.toString());


			// generating boxplots
			//
			jobStatus.addStatusMessage("60", "generating boxplots");
			logger.debug("generating boxplots...\n");
			BoxPlotChart boxplot;
			List<Score> scoreList;
			List<Double> list1, list2;
			//			System.out.println("row name size = " + entities.size());
			jobStatus.addStatusMessage("#", "generating boxplots");
			time = System.currentTimeMillis();
			for (int i=0 ; i<validEntities.size() && i < numberFilter ; i++) {
				try {

					entity = dataFrame.getRowNames().get(i);

					//					System.out.println("ploting entity = " + entity);

					scoreList = entityMap1.get(entity);

					list1 = new ArrayList<Double> (scoreList.size());
					for (int j=0 ; j<scoreList.size(); j++) {
						list1.add((double) scoreList.get(j).getValue());
					}						

					scoreList = entityMap2.get(entity); 
					list2 = new ArrayList<Double> (scoreList.size());
					for (int j=0 ; j<scoreList.size(); j++) {
						list2.add((double) scoreList.get(j).getValue());
					}						

					//					System.out.println("list1, entity = " + entity + ": " + ListUtils.toString(list1));
					//					System.out.println("list2, entity = " + entity + ": " + ListUtils.toString(list2));

					boxplot = createBoxplot(entity, Double.parseDouble(dataFrame.getColumn("adj. p-value").get(i)), list1, list2);
					filename = entity.toLowerCase().replace(" ", "_") + ".png";
					ChartUtilities.saveChartAsPNG(new File(getOutdir() + "/" + filename), boxplot, 400, 200);
					result.addOutputItem(new Item(entity.toLowerCase().replace(" ", "_") + "_boxplot", filename, "Boxplot for " + entity, TYPE.IMAGE));

				} catch (IOException e) {
					e.printStackTrace();
				}				
			}

			jobStatus.addStatusMessage("#", "boxplots, elapsed time = " + ((System.currentTimeMillis()-time)/1000F) + " s");

			// saving data
			//
			jobStatus.addStatusMessage("80", "saving results");
			logger.debug("saving results...");


			// done
			//
			jobStatus.addStatusMessage("100", "done");
			logger.debug("marmite funcitonal enrichment done\n");

		} catch (java.security.InvalidParameterException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (MatrixIndexException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
			logger.error("Error opening the dataset", e.toString());
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (InvalidParameterException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (InvalidColumnIndexException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}


	private Map<String, List<Score>> getEntityMap(List<String> genes, String entity) {	
		ArrayList<Score> list;
		Map<String, List<Score>> map = new HashMap<String, List<Score>>();

		Query query;
		ResultSetHandler rsh = new BeanArrayListHandler(Score.class);
		DBConnection dbConn = new DBConnection("mysql", "mem20", "3306", "Babelomics", "ensembl_user", "ensembl_user");

		String prefix = "select " + (entity.equalsIgnoreCase("roots") ? "word_root" : "entity") + ", hgnc_name, score from";
		if ( entity.equalsIgnoreCase("diseases") ) prefix = prefix + " hsa_bioalmaDis";
		if ( entity.equalsIgnoreCase("products") ) prefix = prefix + " hsa_bioalmaChem";
		if ( entity.equalsIgnoreCase("roots") )    prefix = prefix + " hsa_bioalmaWordRoot";
		prefix = prefix + " where hgnc_name ="; 
		for(String gene: genes) {
			try {
				//System.out.println("query = " + prefix + " '" + gene.toUpperCase() + "'");
				query = dbConn.createSQLQuery(prefix + " \"" + gene.toUpperCase() + "\"");
				list = (ArrayList<Score>) query.execute(rsh);
				if ( list != null ) {
					for(Score score: list) {
						if ( ! map.containsKey(score.getEntity()) ) {
							map.put(score.getEntity(), new ArrayList<Score>());
						}
						map.get(score.getEntity()).add(score);
					}
				}
			} catch (Exception e) {
				System.out.println("marmite, getEntityMap(), error accessing db: " + prefix + " \"" + gene.toUpperCase() + "\", error = " + e.toString());
			}				
		}

		//System.out.println("getEntityMap, genes (" + genes.size() + ") :\n" + genes.toString());
		return map;
	}

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
