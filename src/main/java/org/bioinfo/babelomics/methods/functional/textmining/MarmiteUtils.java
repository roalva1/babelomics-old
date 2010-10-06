package org.bioinfo.babelomics.methods.functional.textmining;

import java.io.File;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.bioinfo.chart.BoxPlotChart;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.db.DBConnection;
import org.bioinfo.db.api.PreparedQuery;
import org.bioinfo.db.handler.BeanArrayListHandler;
import org.bioinfo.db.handler.ResultSetHandler;
import org.jfree.chart.plot.PlotOrientation;

public class MarmiteUtils {

	public static int NO_SORT = 0;
	public static int ASCENDING_SORT = 1;
	public static int DESCENDING_SORT = 2;

	public static int REMOVE_NEVER = 1;
	public static int REMOVE_EACH = 2;
	public static int REMOVE_REF = 3;
	public static int REMOVE_ALL = 4;	
	
	public static Map<String, List<Score>> getEntityMap(List<String> genes, String entity) {	
		ArrayList<Score> list;
		Map<String, List<Score>> map = new HashMap<String, List<Score>>();

		PreparedQuery query;
		//Query query;
		ResultSetHandler rsh = new BeanArrayListHandler(Score.class);
		
		DBConnection dbConn = new DBConnection("mysql", "mem20", "3306", "Babelomics", "ensembl_user", "ensembl_user");

		String prefix = "select " + (entity.equalsIgnoreCase("roots") ? "word_root" : "entity") + ", hgnc_name, score from";
		if ( entity.equalsIgnoreCase("diseases") ) prefix = prefix + " hsa_bioalmaDis";
		if ( entity.equalsIgnoreCase("products") ) prefix = prefix + " hsa_bioalmaChem";
		if ( entity.equalsIgnoreCase("roots") )    prefix = prefix + " hsa_bioalmaWordRoot";
		prefix = prefix + " where hgnc_name = ?"; 
		//prefix = prefix + " where hgnc_name ="; 
		for(String gene: genes) {
			try {
				//System.out.println("query = " + prefix + " '" + gene.toUpperCase() + "'");
				//query = dbConn.createSQLPrepQuery(prefix + " \"" + gene.toUpperCase() + "\"");
				query = dbConn.createSQLPrepQuery(prefix);
				query.setParams(gene.toUpperCase());
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
				System.out.println("marmite, getEntityMap(), error accessing db: " + prefix + ", " + gene.toUpperCase() + ", error = " + e.toString());
			}				
		}

		//System.out.println("getEntityMap, genes (" + genes.size() + ") :\n" + genes.toString());
		return map;
	}

	public static BoxPlotChart createBoxplot(String entity, double pValue, List<Double> list1, List<Double> list2) {

		DecimalFormat df = new DecimalFormat("##0.00000");
		//BoxPlotChart bpc = new BoxPlotChart(entity.toUpperCase() + "\nadj. p-value = " + df.format(pValue), "", "score");
		BoxPlotChart bpc = new BoxPlotChart("", "", "");
		bpc.getLegend().setVisible(false);
		bpc.addSeries(list1, "list 1", "list 1");
		bpc.addSeries(list2, "list 2", "list 2");
		bpc.setOrientation(PlotOrientation.HORIZONTAL);

		return bpc;
	}

	public static int getThresholdPosition(double acum, List<Double> statistic, int order){
		int position = 0;
		for(int i=0; i<statistic.size(); i++){
			if( (order==ASCENDING_SORT && statistic.get(i)>=acum) || (order==DESCENDING_SORT && statistic.get(i)<acum) ) {
				position = i;
				break;
			}			
		}
		return position;
	}
	
	public static List<List<String>> removeDuplicates(List<String> list1, List<String> list2, int duplicatesMode){
		List<List<String>> res = new ArrayList<List<String>>(2);
		
		if(duplicatesMode==REMOVE_REF){
			
			//complementary
			res.add(ListUtils.unique(list1));
			for (String id:list1) {
				if(list2.contains(id)){
					list2.remove(id);
				}
			}
			res.add(list2);
			
		} else if (duplicatesMode==REMOVE_ALL){
			
			// all
			list1 = ListUtils.unique(list1);
			list2 = ListUtils.unique(list2);
			for (String id:list1) {
				if(list2.contains(id)) {
					list1.remove(id);
					list2.remove(id);
				}	
			}
			res.add(list1);
			res.add(list2);
		} else {		
			// each list, by default			
			res.add(ListUtils.unique(list1));
			res.add(ListUtils.unique(list2));
		}
		
		return res;
	}

	public static String getDuplicatesModeString(int duplicatesMode) {
		if (duplicatesMode == REMOVE_REF) {
			return "Removed common IDs in the reference list (list #2)";			
		} else if (duplicatesMode == REMOVE_ALL) {
			return "Removed common IDs in both list";						
		} else {
			return "No removed commons IDs but removed duplicates in each list independently";
		}
	}

	public static void writeExcludedEntitiesFile(File file, List<String> list, Map<String, List<Score>> map) throws IOException {
		String line;
		List<String> lines = new ArrayList<String> ();
		lines.add("#NAMES\tgenes");
		for(String entity: list) {
			lines.add(entity + "\t" + getEntityGenes(map.get(entity)));
		}
		IOUtils.write(file, lines);
	}
	
	public static String getEntityGenes(List<Score> list) {
		String genes = "";
		if ( list != null && list.size() > 0 ) {
			int i;
			for(i=0 ; i<list.size()-1 ; i++) {
				genes = genes + list.get(i).getGene() + ",";
			}
			genes = genes + list.get(i).getGene();
		}
		return genes;
	}
	
}
