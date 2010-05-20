package org.bioinfo.babelomics.tools.expression;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import org.bioinfo.babelomics.methods.expression.clustering.ClusteringUtils;
import org.bioinfo.babelomics.methods.expression.clustering.Kmeans;
import org.bioinfo.babelomics.methods.expression.clustering.Som;
import org.bioinfo.babelomics.methods.expression.clustering.Sota;
import org.bioinfo.babelomics.methods.expression.clustering.Upgma;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.data.tree.multiway.MultipleTree;
import org.bioinfo.data.tree.multiway.MultipleTreeUtils;
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

		if ( "kmeans".equalsIgnoreCase(method) ) {
			try {
				kvalue = Integer.parseInt(commandLine.getOptionValue("kvalue", "15"));
			} catch (NumberFormatException e ) {
				abort("invalidkvalue_execute_clustering", "Invalid k-value", "Invalid value (" + commandLine.getOptionValue("method") + ") for k-value", "Invalid value (" + commandLine.getOptionValue("method") + ") for k-value");
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
			dataset.load();
			dataset.validate();
		} catch (Exception e) {
			abort("exception_execute_clustering", "error reading dataset '" + datasetFile.getName() + "'", e.toString(), StringUtils.getStackTrace(e));
		}		

		//		if(commandLine.hasOption("sample-filter") || commandLine.hasOption("feature-filter")) {
		//			dataset = dataset.getSubDataset(commandLine.getOptionValue("sample-filter"), "4", commandLine.getOptionValue("feature-filter"), ""); 
		//		}


		if ( !"kmeans".equalsIgnoreCase(method) && !"upgma".equalsIgnoreCase(method)  &&
				!"sota".equalsIgnoreCase(method)   && !"som".equalsIgnoreCase(method) ) {
			abort("unknownclusteringmethod_execute_clustering", "Unknown clustering method '" + method + "'", "Unknown clustering method '" + method + "'", "Unknown clustering method '" + method + "'");
		}

		MultipleTree nwGenes = null, nwSamples = null;

		try {
			jobStatus.addStatusMessage("40", "generating genes clusters");
		} catch (FileNotFoundException e) {
			printError("filenotfoundexception_execute_clustering", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
		}

		try {
			nwGenes = runClustering(dataset.getDoubleMatrix(), dataset.getFeatureNames(), dataset.getSampleNames(), method, distance, kvalue, true);
		} catch (Exception e) {
			printError("exception_executesota_clustering", "error running " + method + " algorithm for genes", e.toString(), e);
		}

		List<String> tags = StringUtils.toList("data,newick", ",");
		if ( nwGenes != null ) {
			try {
				String clusterFolder = outdir + "/clusters/";
				new  File(clusterFolder).mkdir();
				MultipleTreeUtils.saveClusters(nwGenes, "", clusterFolder);
				IOUtils.write(new File(this.getOutdir() + "/genes.nw"), nwGenes.toString());		
				result.addOutputItem(new Item("gene_newick_file", "genes.nw", "Clusters of genes", TYPE.FILE, tags, new HashMap<String, String>(2), "Clusters in newick format"));
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
				result.addOutputItem(new Item("sample_newick_file", "samples.nw", "Clusters of samples", TYPE.FILE, tags, new HashMap<String, String>(2), "Clusters in newick format"));
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

		File imgFile;
		String imgFilename;
//		if ( nwSamples != null) {
//			try {
//
//				imgFilename = this.getOutdir() + "/samples." + method + ".png";
//				ClusteringUtils.saveImageTree(nwSamples, "Clusters of samples",  imgFilename, false, false);
//				imgFile = new File(imgFilename);
//				if ( imgFile.exists() ) {
//					result.addOutputItem(new Item(method + "_clustering_image", imgFile.getName(), method.toUpperCase() + " sample clustering image (png format)", TYPE.IMAGE, new ArrayList<String>(2), new HashMap<String, String>(2), "Cluster images"));					
//				} else {
//					printError("execute" + method + "_clustering", "error saving sample clustering image", "error saving sample clustering image");										
//				}
//			} catch (IOException e) {
//				printError("ioexception_execute" + method + "_clustering", "error saving clustering image", e.toString(), StringUtils.getStackTrace(e));
//			}
//		}

//		if ( nwGenes != null && nwSamples != null) {
//			try {
//
//				imgFilename = this.getOutdir() + "/genes." + method + ".png";
//				ClusteringUtils.saveImageTree(nwGenes, "Clusters of genes",  imgFilename, true, true);
//				imgFile = new File(imgFilename);
//				if ( imgFile.exists() ) {
//					result.addOutputItem(new Item(method + "_clustering_image", imgFile.getName(), method.toUpperCase() + " gene clustering image (png format)", TYPE.IMAGE, new ArrayList<String>(2), new HashMap<String, String>(2), "Cluster images"));
//					
//					String mapFilename = imgFilename + ".html.map";
//					File mapFile = new File(mapFilename);
//					if ( mapFile.exists() ) {
//						String[] values;
//						String cleanLine;
//						String clusterDir = outdir + "/clusters/";
//						List<String> lines = IOUtils.readLines(mapFile);
//						new File(clusterDir).mkdir();
//						for(String line: lines) {
//							if ( line.startsWith("<!--cluster_") ) {
//								cleanLine = line.replace("<!--", "").replace("-->", "");
//								values = cleanLine.split(":");
//								IOUtils.write(new File(clusterDir + values[0] + ".txt"), StringUtils.toList(values[1]));
//							}
//						}
//					}
//				} else {
//					printError("execute" + method + "_clustering", "error saving gene clustering image", "error saving gene clustering image");										
//				}
//			} catch (IOException e) {
//				printError("ioexception_execute" + method + "_clustering", "error saving clustering image", e.toString(), StringUtils.getStackTrace(e));
//			}
//		}

		if ( nwGenes != null && nwSamples != null) {
			try {

				imgFilename = this.getOutdir() + "/" + method + ".png";

				int rowOrder[] = getOrder(nwGenes.getLabels(), dataset.getFeatureNames());
				int columnOrder[] = getOrder(nwSamples.getLabels(), dataset.getSampleNames());

				//				System.out.println("row order = \n" + ArrayUtils.toString(rowOrder));
				//				System.out.println("column order = \n" + ArrayUtils.toString(columnOrder));
				//				System.out.println("nw samples labels = " + ListUtils.toString(nwSamples.getLabels(), ",") + "\ndataset sample names = " + ListUtils.toString(dataset.getSampleNames(), ","));
				DoubleMatrix matrix = orderMatrix(dataset.getDoubleMatrix(), rowOrder, columnOrder);

				ClusteringUtils.saveImageTree(matrix, nwGenes, nwSamples, imgFilename, true);
				imgFile = new File(imgFilename);
				if ( imgFile.exists() ) {
					result.addOutputItem(new Item(method + "_clustering_image", imgFile.getName(), method.toUpperCase() + " heatmap image (png format)", TYPE.IMAGE, Arrays.asList("CLUSTER"), new HashMap<String, String>(2), "Cluster images"));					
				} else {
					printError("execute" + method + "_clustering", "error saving clustering image", "error saving clustering image");					
				}
			} catch (IOException e) {
				printError("ioexception_execute" + method + "_clustering", "error saving clustering image", e.toString(), StringUtils.getStackTrace(e));
			}
		}



		if ( nwGenes != null ) {

			if ( "kmeans".equalsIgnoreCase(method) ) {
				String redirectionTags = null; 
				List<String> redirectionInputs = null;
				File redirectionFile = null, list1 = null;
				File list2 = new File(outdir + "/rownames.txt");
				for(int i=1 ; i<=kvalue ; i++) {
					// preparing significant list (top and bottom)
					//
					list1 = new File(outdir + "/cluster_" + i + ".txt");
					if ( list1.exists() || list2.exists() ) {
						redirectionFile = new File(outdir + "/cluster_" + i + "_to_fatigo.redirection");

						redirectionInputs = new ArrayList<String>();

						redirectionInputs.add("comparison=list2list");

						redirectionInputs.add("list1_wum_data=true");
						redirectionInputs.add("list1_databox=" + list1.getName() + " (cluster " + i + " from job $JOB_NAME)");
						redirectionInputs.add("list1=$JOB_FOLDER/" + list1.getName());

						redirectionInputs.add("list2_wum_data=true");
						redirectionInputs.add("list2_databox=" + list2.getName() + " (all genes from job $JOB_NAME)");
						redirectionInputs.add("list2=$JOB_FOLDER/" + list2.getName());

						redirectionInputs.add("duplicates=ref");

						redirectionInputs.add("tool=fatigo");
						redirectionInputs.add("jobname=fatigo from kmeans cluster " + 1);						
						redirectionInputs.add("jobdescription=redirected from job $JOB_NAME");

						try {
							IOUtils.write(redirectionFile.getAbsolutePath(), redirectionInputs);
						} catch (IOException e) {
							// TODO Auto-generated catch block
							e.printStackTrace();
						}								

						if ( redirectionFile.exists() ) {
							redirectionTags = "REDIRECTION(" + redirectionFile.getName() + ":Send to FatiGO tool...)";
							result.addOutputItem(new Item("cluster_" + i + "_to_fatigo", "", "Send cluster " + i + " to FatiGO tool", TYPE.TEXT, StringUtils.toList(redirectionTags, ","), new HashMap<String, String>(2), "Continue processing"));
						}
					}
				}
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
	private MultipleTree runClustering(DoubleMatrix matrix, List<String> rowNames, List<String> colNames, String method, String distance, int kvalue) throws Exception {
		return runClustering(matrix, rowNames, colNames, method, distance, kvalue, false);
	}

	private MultipleTree runClustering(DoubleMatrix matrix, List<String> rowNames, List<String> colNames, String method, String distance, int kvalue, boolean createClusterFiles) throws Exception {
		MultipleTree tree = null;

		if ( "sota".equalsIgnoreCase(method) ) {
			Sota sota = new Sota(matrix, rowNames, colNames, distance, babelomicsHomePath);
			tree = sota.run(createClusterFiles);
		} else if ( "som".equalsIgnoreCase(method) ) {
			Som som = new Som(matrix, rowNames, colNames, distance, babelomicsHomePath);
			tree = som.run(createClusterFiles);
		} else if ( "upgma".equalsIgnoreCase(method) ) {
			Upgma upgma = new Upgma(matrix, rowNames, colNames, distance, babelomicsHomePath);
			tree = upgma.run(createClusterFiles);
		} else if ( "kmeans".equalsIgnoreCase(method) ) {
			int k = (kvalue > rowNames.size() ? rowNames.size() : kvalue);

			Kmeans kmeans = new Kmeans(matrix, rowNames, colNames, distance, k, outdir, babelomicsHomePath);
			tree = kmeans.run(createClusterFiles);			
		}
		return tree;
	}


	private int[] getOrder(List<String> src, List<String> dest) {
		int order[] = new int[src.size()];
		int i=0;
		for(String name: src) {
			name = name.trim();
			//System.out.println("getOrder, name = (" + name + ")");
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
