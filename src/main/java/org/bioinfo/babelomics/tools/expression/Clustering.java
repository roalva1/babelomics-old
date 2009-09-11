package org.bioinfo.babelomics.tools.expression;

import java.io.File;
import java.io.FileNotFoundException; 
import java.io.IOException;
import java.util.List;

import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.babelomics.tools.expression.clustering.ClusteringUtils;
import org.bioinfo.babelomics.tools.expression.clustering.Kmeans;
import org.bioinfo.babelomics.tools.expression.clustering.Sota;
import org.bioinfo.babelomics.tools.expression.clustering.Upgma;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.data.format.core.newick.NewickTree;
import org.bioinfo.math.data.DoubleMatrix;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;


public class Clustering extends BabelomicsTool {


	public Clustering() {
	}

	@Override
	public void initOptions() {
		getOptions().addOption(OptionFactory.createOption("dataset", "the data"));
		getOptions().addOption(OptionFactory.createOption("method", "the method, possible values: upgma, sota, som, kmeans"));
		getOptions().addOption(OptionFactory.createOption("distance", "the distance, possible values: euclidean, spearman, pearson. Default value: euclidean, false"));
		getOptions().addOption(OptionFactory.createOption("kvalue", "k-value for kmeans clustering. Default value: 15", false));
		getOptions().addOption(OptionFactory.createOption("sample-filter", "class variable", false));
		getOptions().addOption(OptionFactory.createOption("feature-filter", "class variable", false));

	}

	@Override
	public void execute() {
		int kvalue = 15;
		Dataset dataset = null;
		String method = commandLine.getOptionValue("method");

		String distance = commandLine.getOptionValue("distance", "euclidean");
		
		try {
			kvalue = Integer.parseInt(commandLine.getOptionValue("kvalue", "15"));
		} catch (NumberFormatException e ) {
			if ( "kmeans".equalsIgnoreCase(method) ) {
				abort("invalidkvalue_execute_clustering", "Invalid k-value", e.toString(), StringUtils.getStackTrace(e));
			}
		}

		String datasetPath = commandLine.getOptionValue("dataset");
		if ( datasetPath == null ) {
			abort("missingdataset_execute_clustering", "Missing dataset", "Missing dataset", "Missing dataset");
		}

		if ( method == null ) {
			abort("missingclusteringmethod_execute_clustering", "Missing clustering method", "Missing clustering method", "Missing clustering method");
		}

		if ( !"upgma".equalsIgnoreCase(method) && !"sota".equalsIgnoreCase(method) && !"som".equalsIgnoreCase(method) && !"kmeans".equalsIgnoreCase(method) ) {
			abort("unknownclusteringmethod_execute_clustering", "Unknown clustering method", "Unknown clustering method '" + method + "'", "Unknown clustering method '" + method + "'");			
		}
		
		try {
			jobStatus.addStatusMessage("20", "reading dataset");
		} catch (FileNotFoundException e) {
			abort("filenotfoundexception_execute_clustering", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
		}

		File datasetFile = new File(commandLine.getOptionValue("dataset"));
		try {
			dataset = new Dataset(datasetFile);
		} catch (Exception e) {
			abort("exception_execute_clustering", "error reading dataset '" + datasetFile.getName() + "'", e.toString(), StringUtils.getStackTrace(e));
		}		

		if ( dataset.getDoubleMatrix() == null ) { 
			try {
				dataset.load();
			} catch (Exception e) {
				abort("exception_execute_clustering", "Error loading dataset '" + datasetFile.getName() + "'", e.toString(), StringUtils.getStackTrace(e));
			}
			dataset.validate();
		}

//		if(commandLine.hasOption("sample-filter") || commandLine.hasOption("feature-filter")) {
//			dataset = dataset.getSubDataset(commandLine.getOptionValue("sample-filter"), "4", commandLine.getOptionValue("feature-filter"), ""); 
//		}

		
		if ( !"kmeans".equalsIgnoreCase(method) && !"upgma".equalsIgnoreCase(method)  &&
			 !"sota".equalsIgnoreCase(method)   && !"som".equalsIgnoreCase(method) ) {
			abort("unknownclusteringmethod_execute_clustering", "Unknown clustering method '" + method + "'", "Unknown clustering method '" + method + "'", "Unknown clustering method '" + method + "'");
		}
		
		NewickTree nwGenes = null, nwSamples = null;
		
		try {
			jobStatus.addStatusMessage("40", "generating genes clusters");
		} catch (FileNotFoundException e) {
			printError("filenotfoundexception_execute_clustering", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
		}

		try {
			nwGenes = runClustering(dataset.getDoubleMatrix(), dataset.getFeatureNames(), dataset.getSampleNames(), method, distance, kvalue);
		} catch (Exception e) {
			printError("exception_executesota_clustering", "error running " + method + " algorithm for genes", e.toString(), e);
		}
		
		if ( nwGenes != null ) {
			try {
				IOUtils.write(new File(this.getOutdir() + "/genes.nw"), nwGenes.toString());
				result.addOutputItem(new Item("gene_newick_file", "genes.nw", "Clusters of genes (newick format file)", TYPE.FILE));
			} catch (IOException e) {
				printError("ioexception_executesota_clustering", "error saving genes newick", e.toString(), e);
				nwGenes = null;
			}			
		}
		
		try {
			jobStatus.addStatusMessage("60", "generating samples clusters");
		} catch (FileNotFoundException e) {
			abort("filenotfoundexception_execute_clustering", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
		}

		try {
			nwSamples = runClustering(new DoubleMatrix(dataset.getDoubleMatrix().transpose().getData()), dataset.getSampleNames(), dataset.getFeatureNames(), method, distance, kvalue);
		} catch (Exception e) {
			printError("exception_execute" + method + "_clustering", "error running " + method + " algorithm for samples", e.toString(), e);
		}
		
		if ( nwSamples != null ) {
			try {
				IOUtils.write(new File(this.getOutdir() + "/samples.nw"), nwSamples.toString());
				result.addOutputItem(new Item("sample_newick_file", "samples.nw", "Clusters of samples (newick format file)", TYPE.FILE));
			} catch (IOException e) {
				printError("ioexception_execute" + method + "_clustering", "error saving samples newick", e.toString(), e);
				nwSamples = null;
			}			
		}
		
		try {
			jobStatus.addStatusMessage("80", "generating clustering image");
		} catch (FileNotFoundException e) {
			abort("filenotfoundexception_execute_clustering", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
		}

		if ( nwGenes != null && nwSamples != null ) {
			try {
				String imgFilename = this.getOutdir() + "/" + method;
				
				int rowOrder[] = getOrder(nwGenes.getLabels(), dataset.getFeatureNames());
				int columnOrder[] = getOrder(nwSamples.getLabels(), dataset.getSampleNames());
				DoubleMatrix matrix = orderMatrix(dataset.getDoubleMatrix(), rowOrder, columnOrder);
				
				ClusteringUtils.saveImageTree(matrix, nwGenes, nwSamples, imgFilename);
				File imgFile = new File(imgFilename + ".png");
				if ( imgFile.exists() ) {
					result.addOutputItem(new Item(method + "_clustering_image", method + ".png", method.toUpperCase() + " clustering image (png format)", TYPE.IMAGE));					
				} else {
					printError("execute" + method + "_clustering", "error saving clustering image", "error saving clustering image");					
				}
			} catch (IOException e) {
				printError("ioexception_execute" + method + "_clustering", "error saving clustering image", e.toString(), e);
			}
		}
	}

	/**
	 * 
	 * @param matrix
	 * @param rowNames
	 * @param colNames
	 * @param method
	 * @param distance
	 * @param kvalue
	 * @return
	 * @throws Exception
	 */
	private NewickTree runClustering(DoubleMatrix matrix, List<String> rowNames, List<String> colNames, String method, String distance, int kvalue) throws Exception {
		NewickTree tree = null;

		
		if ( "sota".equalsIgnoreCase(method) ) {
			Sota sota = new Sota(matrix, rowNames, colNames, distance);
			tree = sota.run();
		} else if ( "som".equalsIgnoreCase(method) ) {
			throw new Exception("SOM algorithm is not implemented yet");
		} else if ( "upgma".equalsIgnoreCase(method) ) {
			Upgma upgma = new Upgma(matrix, rowNames, colNames, distance);
			tree = upgma.run();
		} else if ( "kmeans".equalsIgnoreCase(method) ) {
			Kmeans kmeans = new Kmeans(matrix, rowNames, colNames, distance, kvalue);
			tree = kmeans.run();			
		}
		return tree;
	}
		
		
	private int[] getOrder(List<String> src, List<String> dest) {
		int order[] = new int[src.size()];
		int i=0;
		for(String name: src) {
			order[i++] = dest.indexOf(name);
		}
		return order;
	}

	private DoubleMatrix orderMatrix(DoubleMatrix inputMatrix, int rowOrder[], int columnOrder[]) {
		DoubleMatrix matrix = new DoubleMatrix(inputMatrix.getRowDimension(), inputMatrix.getColumnDimension());
		for(int row=0 ; row<matrix.getRowDimension() ; row++) {
			for(int col=0 ; col<matrix.getColumnDimension() ; col++) {
				matrix.set(row, col, inputMatrix.get(rowOrder[row], columnOrder[col]));
			}
		}				
		return matrix;
	}

}
