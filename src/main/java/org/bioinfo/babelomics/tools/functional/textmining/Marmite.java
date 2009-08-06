package org.bioinfo.babelomics.tools.functional.textmining;

import static org.junit.Assert.fail;

import java.awt.Color;
import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.ParseException;
import org.apache.commons.math.MathException;
import org.apache.commons.math.linear.MatrixIndexException;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.chart.BoxPlotChart;
import org.bioinfo.collections.exceptions.InvalidColumnIndexException;
import org.bioinfo.collections.list.NamedArrayList;
import org.bioinfo.collections.matrix.DataFrame;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.data.dataset.FeatureData;
import org.bioinfo.db.DBConnection;
import org.bioinfo.db.Query;
import org.bioinfo.db.handler.BeanArrayListHandler;
import org.bioinfo.db.handler.ResultSetHandler;
import org.bioinfo.graphics.canvas.Canvas;
import org.bioinfo.graphics.canvas.feature.ScoreFeature;
import org.bioinfo.graphics.canvas.panel.GridPanel;
import org.bioinfo.graphics.canvas.track.GridTrack;
import org.bioinfo.io.FileSystemUtils;
import org.bioinfo.math.comparison.KolmogorovSmirnovTest;
import org.bioinfo.math.exception.InvalidParameterException;
import org.bioinfo.math.result.CorrelationTestResult;
import org.bioinfo.math.result.KolmogorovSmirnovTestResult;
import org.bioinfo.math.result.TestResultList;
import org.bioinfo.math.stats.correlation.CorrelationTest;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;
import org.bioinfo.utils.ListUtils;
import org.bioinfo.utils.StringUtils;
import org.jfree.chart.ChartUtilities;
import org.jfree.chart.plot.PlotOrientation;

public class Marmite extends BabelomicsTool {

	public Marmite(String[] args) {
		initOptions();
	}

	@Override
	public void initOptions() {
		options.addOption(OptionFactory.createOption("list1", "the feature data containig the list #1 of genes, or the feature data file"));
		options.addOption(OptionFactory.createOption("list2", "the feature data containig the list #2 of genes, or the feature data file"));
		options.addOption(OptionFactory.createOption("bioentity-name", "Valid values are: diseases (for disease associated words), products (for chemical products) or roots (for word roots)"));
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
		try {

			// reading data
			//
			jobStatus.addStatusMessage("20", "reading data");
			logger.debug("reading data...\n");

			Map<String, List<Score>> entityMap1 = getEntityMap(StringUtils.stringToList(FileSystemUtils.fileToString(f1)), bioentity);	
			Map<String, List<Score>> entityMap2 = getEntityMap(StringUtils.stringToList(FileSystemUtils.fileToString(f2)), bioentity);	

			List<String> entities = new ArrayList<String> ();
			for(String entity: entityMap1.keySet()) {
				if ( entityMap1.get(entity) != null && entityMap1.get(entity).size() >= scoreFilter && 
						entityMap2.get(entity) != null && entityMap2.get(entity).size() >= scoreFilter ) {
					entities.add(entity);
					
//					System.out.println("map1, entity = " + entity + ", score: " + ListUtils.toString(entityMap1.get(entity)));
//					System.out.println("map2, entity = " + entity + ", score: " + ListUtils.toString(entityMap2.get(entity)));
//					System.out.println("-----");
				}				
			}

			String entity;
			List<Score> scoreList1, scoreList2;
			double[][] scoreMatrix1 = new double[entities.size()][];
			double[][] scoreMatrix2 = new double[entities.size()][];
			for(int i=0 ; i<entities.size() ; i++) {
				entity = entities.get(i);
				scoreList1 = entityMap1.get(entity); 
				scoreList2 = entityMap2.get(entity); 
				if ( scoreList1 != null && scoreList1.size() >= scoreFilter && scoreList2 != null && scoreList2.size() >= scoreFilter ) {

					scoreMatrix1[i] = new double[scoreList1.size()];
					for (int j=0 ; j<scoreList1.size(); j++) {
						scoreMatrix1[i][j] = scoreList1.get(j).getValue();
					}

					scoreMatrix2[i] = new double[scoreList2.size()];
					for (int j=0 ; j<scoreList2.size(); j++) {
						scoreMatrix2[i][j] = scoreList2.get(j).getValue();
					}						
				}								
			}

//			System.out.println("matrix1 :");
//			for (int i=0 ; i<scoreMatrix1.length; i++) {
//				for (int j=0 ; j<scoreMatrix1[i].length; j++) {
//					System.out.print(scoreMatrix1[i][j] + "\t");
//				}
//				System.out.println("");
//			}						
//			
//			System.out.println("matrix2 :");
//			for (int i=0 ; i<scoreMatrix2.length; i++) {
//				for (int j=0 ; j<scoreMatrix2[i].length; j++) {
//					System.out.print(scoreMatrix2[i][j] + "\t");
//				}
//				System.out.println("");
//			}						


			System.out.println("scoreMatrix1 size = " + scoreMatrix1.length + ", " + scoreMatrix1[0].length);
			System.out.println("scoreMatrix2 size = " + scoreMatrix2.length + ", " + scoreMatrix2[0].length);

			// marmite functional enrichment
			//
			jobStatus.addStatusMessage("40", "computing marmite functional enrichment");
			logger.debug("computing marmite functional enrichment...\n");

			TestResultList<KolmogorovSmirnovTestResult> res = new KolmogorovSmirnovTest().compute(scoreMatrix1, scoreMatrix2);
			System.out.println("kolmogorov result :\n" +  res.toString());

			int[] rowOrder = ListUtils.order(ListUtils.toList(res.getAdjPValues()));
			System.out.println("order = " + ListUtils.toString(ListUtils.toList(rowOrder)));

			DataFrame dataFrame = new DataFrame(entities.size(), 0);
			dataFrame.setRowNames(ListUtils.ordered(entities, rowOrder));
								
			System.out.println("********* sides = " + res.getSides());
			
			dataFrame.addColumn("statistic", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getStatistics()), rowOrder)));
			dataFrame.addColumn("side", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getSides()), rowOrder)));
			dataFrame.addColumn("p-value", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getPValues()), rowOrder)));
			dataFrame.addColumn("adj. p-value", ListUtils.toStringList(ListUtils.ordered(ListUtils.toList(res.getAdjPValues()), rowOrder)));
			
			// generating boxplots
			//
			jobStatus.addStatusMessage("60", "generating boxplots");
			logger.debug("generating boxplots...\n");
			String path;
			BoxPlotChart boxplot;
			List<Score> scoreList;
			List<Double> list1, list2;
			System.out.println("row name size = " + entities.size());
			for (int i=0 ; i<entities.size() && i < numberFilter ; i++) {
				try {
					
					entity = dataFrame.getRowNames().get(i);
					
					System.out.println("ploting entity = " + entity);
					
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
					
					System.out.println("list1, entity = " + entity + ": " + ListUtils.toString(list1));
					System.out.println("list2, entity = " + entity + ": " + ListUtils.toString(list2));
					
					boxplot = createBoxplot(entity, Double.parseDouble(dataFrame.getColumn("adj. p-value").get(i)), list1, list2);
					path = getOutdir() + "/" + entity.toLowerCase().replace(" ", "_") + ".png";
					ChartUtilities.saveChartAsPNG(new File(path), boxplot, 400, 200);
					result.addOutputItem(new Item(entity.toLowerCase().replace(" ", "_") + "_boxplot", path, "Boxplot for " + entity + ": ", TYPE.IMAGE));
					
				} catch (IOException e) {
					e.printStackTrace();
				}				
			}


			// saving data
			//
			jobStatus.addStatusMessage("80", "saving results");
			logger.debug("saving results...");

								
			FeatureData featureData = new FeatureData(dataFrame);
			path = getOutdir() + "/marmite_output.txt";
			featureData.write(new File(path));
			result.addOutputItem(new Item("marmite_file", path, "The marmite output is: ", TYPE.FILE));
			
			System.out.println("marmite output:\n" +  featureData.toString());

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
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (InvalidParameterException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		} catch (InvalidColumnIndexException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private Map<String, List<Score>> getEntityMap(List<String> genes, String entity) {
		
		
//		String q = "select s.id, s.name from symptom s, symptom_syndrome ss where ss.symptom_id = s.id and ss.syndrome_id = " + syndromeId;
//		//System.out.println("getSyndromes, q = " + q);
//		try {
//			ResultSetHandler rsh = new BeanArrayListHandler(Symptom.class);
//			Query query = dbConn.createSQLQuery(q);
//			list = (ArrayList<Symptom>) query.execute(rsh);
//		} catch (Exception e) {
//			list = null;
//		}
//		return list;
		
		
		ArrayList<Score> list;
		Map<String, List<Score>> map = new HashMap<String, List<Score>>();
		
		Query query;
		ResultSetHandler rsh = new BeanArrayListHandler(Score.class);
		DBConnection dbConn = new DBConnection("mysql", "mem20", "3306", "Babelomics", "ensembl_user", "ensembl_user");

		String prefix = "select " + (entity.equalsIgnoreCase("roots") ? "word_root" : "entity") + ", score from";
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
						if ( ! map.containsKey(score.getName()) ) {
							map.put(score.getName(), new ArrayList<Score>());
						}
						map.get(score.getName()).add(new Score(gene, score.getValue()));
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
