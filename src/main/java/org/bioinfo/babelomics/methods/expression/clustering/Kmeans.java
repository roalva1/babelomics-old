package org.bioinfo.babelomics.methods.expression.clustering;

import java.io.File;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.bioinfo.commons.exec.Command;
import org.bioinfo.commons.exec.SingleProcess;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ArrayUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.commons.utils.MapUtils;
import org.bioinfo.data.format.core.newick.NewickTree;
import org.bioinfo.data.format.io.NewickParser;
import org.bioinfo.math.data.DoubleMatrix;

public class Kmeans extends Cluster {
	private int kvalue;
	
	public Kmeans(DoubleMatrix matrix, List<String> rowNames, List<String> colNames, String distance, int kvalue) {
		super(matrix, rowNames, colNames, distance);
		this.kvalue = kvalue;
	}

	@Override
	public NewickTree run() throws Exception {
		NewickTree nw = null;
		
		File tmpDir = File.createTempFile("input", ".dir");
		tmpDir.delete();
		tmpDir.mkdir();
		
		File inputFile = new File(tmpDir.getAbsoluteFile() + "/in.txt");
		File outputFile = new File(tmpDir.getAbsoluteFile() + "/in_K_G" + kvalue + ".kgg");

//		System.out.println("(infile, outfile) = (" + inputFile.getAbsolutePath() + ", " + outputFile.getAbsolutePath() + ")");
		
		List<String> lines = new ArrayList<String>(rowNames.size() + 1);
		lines.add("NAMES\t" + ListUtils.toString(colNames, "\t"));
		for(int i=0 ; i<rowNames.size() ; i++) {
			lines.add(rowNames.get(i) + "\t" + ArrayUtils.toString(matrix.getRow(i), "\t"));
		}
		IOUtils.write(inputFile, lines);
		
		String cmdStr = System.getenv("BABELOMICS_HOME") + "/bin/clustering/kmeans -f " + inputFile.getAbsolutePath() + " -g " + getDistance(distance) + " -k " + kvalue;
		System.out.println("clustering command KMEANS: " + cmdStr);		
		
		Command cmd = new Command(cmdStr); 
		SingleProcess sp = new SingleProcess(cmd);
		sp.runSync();

		if ( outputFile.exists() && outputFile.getTotalSpace() > 0 ) {
			lines = IOUtils.readLines(outputFile);
			lines.remove(0);
			Map<String, List<String>> map = new HashMap<String, List<String>> ();
			String[] values = null;
			for(String line: lines) {
				values = line.split("\t");
				if ( !map.containsKey(values[1]) ) {
					map.put(values[1], new ArrayList<String>());
				}
				map.get(values[1]).add(values[0]);
			}
			StringBuilder sb = new StringBuilder();
			sb.append("(");
			List<String> keys = MapUtils.getKeys(map);
			int i=0;
			for(i=0 ; i<keys.size()-1 ; i++) {
				sb.append("(").append(ListUtils.toString(map.get(keys.get(i)), ",")).append("),");
			}
			sb.append("(").append(ListUtils.toString(map.get(keys.get(i)), ",")).append("));");
//			System.out.println("newick " + sb.toString());
			nw = new NewickParser().parse(sb.toString());
		} else {
			throw new Exception("Impossible to generate newick");
		}
		
		return nw;
	}
	
	private int getDistance(String distance) {
        if ( "none".equalsIgnoreCase(distance) ) {
			return 1;
		} else if ( "uncentered".equalsIgnoreCase(distance) ) {
			return 1;
		} else if ( "pearson".equalsIgnoreCase(distance) ) {
			return 2;
		} else if ( "spearman".equalsIgnoreCase(distance) ) {
			return 5;
		} else if ( "kendall".equalsIgnoreCase(distance) ) {
			return 6;
		} else if ( "euclidean".equalsIgnoreCase(distance) ) {
			return 7;
		}
        return 7;
	}

}
