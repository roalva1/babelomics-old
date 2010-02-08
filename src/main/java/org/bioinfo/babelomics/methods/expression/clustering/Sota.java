package org.bioinfo.babelomics.methods.expression.clustering;

import java.io.File;
import java.util.ArrayList;
import java.util.List;

import org.bioinfo.commons.exec.Command;
import org.bioinfo.commons.exec.SingleProcess;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ArrayUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.data.format.io.parser.NewickParser;
import org.bioinfo.data.tree.multiway.MultipleTree;
import org.bioinfo.math.data.DoubleMatrix;

public class Sota extends Cluster {

	public Sota(DoubleMatrix matrix, List<String> rowNames, List<String> colNames, String distance) {
		super(matrix, rowNames, colNames, distance);
	}
	
	@Override
	public MultipleTree run() throws Exception {
		MultipleTree nw = null;
		
		File inputFile = File.createTempFile("input", null);
		File outputFile = File.createTempFile("output", null);

		System.out.println("(infile, outfile) = (" + inputFile.getAbsolutePath() + ", " + outputFile.getAbsolutePath() + ")");
		
		List<String> lines = new ArrayList<String>(rowNames.size() + 1);
		lines.add("#NAMES\t" + ListUtils.toString(colNames, "\t"));
		for(int i=0 ; i<rowNames.size() ; i++) {
			lines.add(rowNames.get(i) + "\t" + ArrayUtils.toString(matrix.getRow(i), "\t"));
		}
		IOUtils.write(inputFile, lines);
		
		String cmdStr = System.getenv("BABELOMICS_HOME") + "/bin/clustering/sota " + inputFile.getAbsolutePath() + " " + outputFile.getAbsolutePath() + " " + distance + " -newick";
		System.out.println("clustering command SOTA: " + cmdStr);
		Command cmd = new Command(cmdStr); 
		SingleProcess sp = new SingleProcess(cmd);
		sp.runSync();

		if ( outputFile.exists() && outputFile.getTotalSpace() > 0 ) {
			nw = new NewickParser().parse(IOUtils.toString(outputFile));
		} else {
			throw new Exception("Impossible to generate newick");
		}
		
		return nw;
	}
}
