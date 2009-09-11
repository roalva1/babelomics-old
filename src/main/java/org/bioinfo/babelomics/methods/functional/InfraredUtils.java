package org.bioinfo.babelomics.methods.functional;

import java.sql.SQLException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.infrared.common.dbsql.DBConnector;
import org.bioinfo.infrared.common.feature.FeatureList;
import org.bioinfo.infrared.core.Gene;
import org.bioinfo.infrared.core.dbsql.GeneDBManager;
import org.bioinfo.infrared.core.dbsql.XRefDBManager;
import org.bioinfo.infrared.funcannot.AnnotationItem;
import org.bioinfo.infrared.funcannot.dbsql.AnnotationDBManager;
import org.bioinfo.infrared.funcannot.filter.Filter;
import org.bioinfo.infrared.funcannot.filter.GOFilter;
import org.bioinfo.infrared.funcannot.filter.KeggFilter;

public class InfraredUtils {

	// get annotations from a filter
	public static FeatureList<AnnotationItem> getAnnotations(DBConnector dbConnector, List<String> ids, Filter filter) throws SQLException, IllegalAccessException, ClassNotFoundException, InstantiationException{
		AnnotationDBManager annotationMng = new AnnotationDBManager(dbConnector);
		if(filter instanceof GOFilter){
			return annotationMng.getGOAnnotation(ids, (GOFilter) filter);
		} else if(filter instanceof KeggFilter){
			return annotationMng.getKeggAnnotation(ids, (KeggFilter) filter);
		}
		return null;
	}	
	
	// get genome (ensembl gene format)
	public static List<String> getGenome(DBConnector dbConnector){
		try {
			return new GeneDBManager(dbConnector).getAllEnsemblIds();
		} catch (Exception e) {
			e.printStackTrace();
		}
		return null;
	}

	// get chromosome region genes (ensembl gene format)
	public static List<String> getChromosomeRegionGenes(DBConnector dbConnector, String chromosome, int start, int end) {
		try {
			FeatureList<Gene> genes = new GeneDBManager(dbConnector).getAllByLocation(chromosome, start, end);
			List<String> list = new ArrayList<String> (genes.size());
			for (Gene gene: genes) {
				list.add(gene.getId());
			}
			return list;
		} catch (Exception e) {
			e.printStackTrace();
			return null;
		}
	}

	// translate a list of ids to ensembl format
	public static List<String> toEnsemblId(DBConnector dbConnector, List<String> ids) {
		try {	
			List<List<String>> list = new XRefDBManager(dbConnector).getIdsByDBName(ids, "ensembl_gene");			
			List<String> ensemblIds = new ArrayList<String>();
			if(list==null) {				
				ensemblIds = null;
			} else {
				
				for(List<String> stringList: list) {					
					if(stringList!=null){
						for(String ensemblId: stringList) {						
							ensemblIds.add(ensemblId);
						}
					}
				}
				ensemblIds = ListUtils.unique(ensemblIds);
			}
			return ensemblIds;
		} catch (Exception e) {
			e.printStackTrace();
			return null;
		}
	}
	
	// translate a list of ids to ensembl format
	public static Map<String, List<String>> getEnsemblMap(DBConnector dbConnector, List<String> ids) {
		try {
			Map<String, List<String>> map = new HashMap<String, List<String>>();			
			List<List<String>> lists = new XRefDBManager(dbConnector).getIdsByDBName(ids, "ensembl_gene");

			List<String> ensemblIds = new ArrayList<String>();
			if ( lists == null ) {
				map = null;
			} else {
				for(int i=0 ; i<lists.size() ; i++) {
					if ( lists.get(i) == null || lists.get(i).size() <= 0) {
						map.put(ids.get(i), null);
					} else {
						if ( !map.containsKey(ids.get(i)) ) {
							map.put(ids.get(i), new ArrayList<String>());
						}
						map.get(ids.get(i)).addAll(lists.get(i));
					}
				}
			}
			return map;
		} catch (Exception e) {
			e.printStackTrace();
			return null;
		}
	}

	
}
