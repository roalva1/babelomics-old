package org.bioinfo.babelomics.tools.expression.clustering;

import java.util.List;

import org.bioinfo.data.format.core.newick.NewickTree;
import org.bioinfo.math.data.DoubleMatrix;

public abstract class Cluster {
	
	protected DoubleMatrix matrix;
	protected List<String> rowNames;
	protected List<String> colNames;
	protected String distance;

	public Cluster(DoubleMatrix matrix, List<String> rowNames, List<String> colNames, String distance) {
		this.matrix = matrix;
		this.rowNames = rowNames;
		this.colNames = colNames;
		this.distance = distance;
	}

	public abstract NewickTree run() throws Exception;
}
