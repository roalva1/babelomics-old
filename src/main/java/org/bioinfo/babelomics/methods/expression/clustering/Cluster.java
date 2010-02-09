package org.bioinfo.babelomics.methods.expression.clustering;

import java.util.List;

import org.bioinfo.data.tree.multiway.MultipleTree;
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

	public abstract MultipleTree run() throws Exception;
}
