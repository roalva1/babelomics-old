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
import org.bioinfo.babelomics.tools.expression.differential.DiffExpressionUtils;
import org.bioinfo.commons.io.utils.FileUtils;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ListUtils;
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
		getOptions().addOption(OptionFactory.createOption("sample-clustering", "clustering of samples", false));
		getOptions().addOption(OptionFactory.createOption("gene-clustering", "clustering of genes", false));
		getOptions().addOption(OptionFactory.createOption("method", "the method, possible values: upgma, sota, som, kmeans"));
		getOptions().addOption(OptionFactory.createOption("distance", "the distance, possible values: euclidean, spearman, pearson. Default value: euclidean", false));
		getOptions().addOption(OptionFactory.createOption("sample-kvalue", "k-value for kmeans sample-clustering. Default value: 5", false));
		getOptions().addOption(OptionFactory.createOption("gene-kvalue", "k-value for kmeans gene-clustering. Default value: 15", false));
		getOptions().addOption(OptionFactory.createOption("sample-filter", "class variable", false));
		getOptions().addOption(OptionFactory.createOption("feature-filter", "class variable", false));
	}

	@Override
	public void execute() {
		boolean sampleClustering = true;
		boolean geneClustering = false;
		int sampleKvalue = 5;
		int geneKvalue = 15;
		Dataset dataset = null;

		String sampleClusteringParam = commandLine.getOptionValue("sample-clustering",  null);
		String geneClusteringParam = commandLine.getOptionValue("gene-clustering",  null);

		sampleClustering = Boolean.parseBoolean(sampleClusteringParam == null ? "false" : sampleClusteringParam);
		geneClustering = Boolean.parseBoolean(geneClusteringParam == null ? "false" : geneClusteringParam);

		String method = commandLine.getOptionValue("method");
		String distance = commandLine.getOptionValue("distance", "euclidean");

		if ( "kmeans".equalsIgnoreCase(method) ) {
			try {
				sampleKvalue = Integer.parseInt(commandLine.getOptionValue("sample-kvalue", "5"));
			} catch (NumberFormatException e ) {
				abort("invalidkvalue_execute_clustering", "Invalid k-value", "Invalid value (" + commandLine.getOptionValue("sample-kvalue") + ") for k-value", "Invalid value (" + commandLine.getOptionValue("sample-kvalue") + ") for k-value");
			}
			try {
				geneKvalue = Integer.parseInt(commandLine.getOptionValue("gene-kvalue", "15"));
			} catch (NumberFormatException e ) {
				abort("invalidkvalue_execute_clustering", "Invalid k-value", "Invalid value (" + commandLine.getOptionValue("gene-kvalue") + ") for k-value", "Invalid value (" + commandLine.getOptionValue("gene-kvalue") + ") for k-value");
			}
		}

		// input parameters
		//
		List<String> clusterings = new ArrayList<String>();
		if ( sampleClustering ) { clusterings.add(" samples"); }
		if ( geneClustering ) { clusterings.add(" genes"); }

		if ( clusterings.size() == 0 ) {
			abort("missingclusteringtype_execute_clustering", "Missing type of clustering: samples or/and genes", "Missing type of clustering: samples or/and genes", "Missing type of clustering: samples or/and genes");
		}

		result.addOutputItem(new Item("sampleclustering_input_param", ListUtils.toString(clusterings, ","), "Clustering of", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
		clusterings.clear();
		clusterings.add(method);
		if ( "kmeans".equalsIgnoreCase(method) ) {
			if ( sampleClustering ) { clusterings.add(" k-value = " + sampleKvalue + " for clustering of samples"); }
			if ( geneClustering ) { clusterings.add(" k-value = " + geneKvalue + " for clustering of genes"); }
		}
		result.addOutputItem(new Item("method_input_param", ListUtils.toString(clusterings, ","), "Method", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));
		result.addOutputItem(new Item("distance_input_param", distance, "Distance", Item.TYPE.MESSAGE, Arrays.asList("INPUT_PARAM"), new HashMap<String,String>(), "Input parameters"));

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
			printError("filenotfoundexception_execute_clustering", "Error", "Internal error updating job status file (file not found)");
		}

		File datasetFile = new File(commandLine.getOptionValue("dataset"));
		try {
			dataset = new Dataset(datasetFile);
			dataset.load();
			dataset.validate();
		} catch (Exception e) {
			abort("exception_execute_clustering", "Error", "Error reading dataset '" + datasetFile.getName() + "'", "");
		}		

		//		if(commandLine.hasOption("sample-filter") || commandLine.hasOption("feature-filter")) {
		//			dataset = dataset.getSubDataset(commandLine.getOptionValue("sample-filter"), "4", commandLine.getOptionValue("feature-filter"), ""); 
		//		}


		if ( !"kmeans".equalsIgnoreCase(method) && !"upgma".equalsIgnoreCase(method)  &&
				!"sota".equalsIgnoreCase(method)   && !"som".equalsIgnoreCase(method) ) {
			abort("unknownclusteringmethod_execute_clustering", "Unknown clustering method '" + method + "'", "Unknown clustering method '" + method + "'", "Unknown clustering method '" + method + "'");
		}

		List<String> tags = null;
		MultipleTree nwGenes = null, nwSamples = null;


		if ( geneClustering ) {
			try {
				jobStatus.addStatusMessage("40", "generating genes clusters");
			} catch (FileNotFoundException e) {
				printError("filenotfoundexception_execute_clustering", "Error", "Internal error updating job status file (file not found)");
			}

			try {
				nwGenes = runClustering(dataset.getDoubleMatrix(), dataset.getFeatureNames(), dataset.getSampleNames(), method, distance, geneKvalue, true);
			} catch (Exception e) {
				printError("exception_executesota_clustering", "Error", "Error running " + method + " algorithm for genes");
			}

			tags = StringUtils.toList("data,newick", ",");
			if ( nwGenes != null ) {
				try {
					String clusterFolder = outdir + "/clusters/";
					new  File(clusterFolder).mkdir();
					MultipleTreeUtils.saveClusters(nwGenes, "", clusterFolder);

					File redirectionFile = null;
					File topListFile = null;
					File bottomListFile = new File(clusterFolder + "cluster_0.txt");
					IOUtils.write(bottomListFile, nwGenes.getLabels());
					if ( bottomListFile.exists() ) {
						File[] clusterFiles = FileUtils.listFiles(new File(clusterFolder), "cluster_.+.txt", true);
						for(File clusterFile: clusterFiles) {					
							if ( !clusterFile.getName().equalsIgnoreCase(bottomListFile.getName()) ) {
								redirectionFile = new File(clusterFolder + clusterFile.getName().replace(".txt", "") + ".fatigo.redirection");
								createFatiGORedirectionFile(redirectionFile, clusterFile, bottomListFile);
							}
							redirectionFile = new File(clusterFolder + clusterFile.getName().replace(".txt", "") + ".fatigo.genome.redirection");
							createFatiGORedirectionFile(redirectionFile, clusterFile, null);							
						}					
					}

					IOUtils.write(new File(this.getOutdir() + "/genes.nw"), nwGenes.toString());		
					result.addOutputItem(new Item("gene_newick_file", "genes.nw", "Clusters of genes", TYPE.FILE, tags, new HashMap<String, String>(2), "Clusters in newick format"));
				} catch (IOException e) {
					printError("ioexception_executesota_clustering", "Error", "Error saving genes newick");
					nwGenes = null;
				}			
			}
		}

		if ( sampleClustering ) {
			try {
				jobStatus.addStatusMessage("60", "generating samples clusters");
			} catch (FileNotFoundException e) {
				printError("filenotfoundexception_execute_clustering", "Error", "Internal error updating job status file (file not found)");
			}

			try {
				nwSamples = runClustering(new DoubleMatrix(dataset.getDoubleMatrix().transpose().getData()), dataset.getSampleNames(), dataset.getFeatureNames(), method, distance, sampleKvalue);
			} catch (Exception e) {
				printError("exception_execute" + method + "_clustering", "Error", "Error running " + method + " algorithm for samples");
			}

			tags = StringUtils.toList("data,newick", ",");
			if ( nwSamples != null ) {
				try {
					IOUtils.write(new File(this.getOutdir() + "/samples.nw"), nwSamples.toString());
					result.addOutputItem(new Item("sample_newick_file", "samples.nw", "Clusters of samples", TYPE.FILE, tags, new HashMap<String, String>(2), "Clusters in newick format"));
				} catch (IOException e) {
					printError("ioexception_execute" + method + "_clustering", "Error", "Error saving samples newick");
					nwSamples = null;
				}			
			}
		}

		try {
			jobStatus.addStatusMessage("80", "generating clustering image");
		} catch (FileNotFoundException e) {
			printError("filenotfoundexception_execute_clustering", "Error", "Internal error updating job status file (file not found)");
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

		if ( dataset.getRowDimension() <= 1000 ) {
			int rowOrder[] = null;
			int columnOrder[] = null;

			try {

				imgFilename = this.getOutdir() + "/" + method + ".png";

				if ( nwGenes != null ) {
					rowOrder = getOrder(nwGenes.getLabels(), dataset.getFeatureNames());
				}

				if ( nwSamples != null ) {
					columnOrder = getOrder(nwSamples.getLabels(), dataset.getSampleNames());
				}

				//				System.out.println("row order = \n" + ArrayUtils.toString(rowOrder));
				//				System.out.println("column order = \n" + ArrayUtils.toString(columnOrder));
				//				System.out.println("nw samples labels = " + ListUtils.toString(nwSamples.getLabels(), ",") + "\ndataset sample names = " + ListUtils.toString(dataset.getSampleNames(), ","));
				DoubleMatrix matrix = orderMatrix(dataset.getDoubleMatrix(), rowOrder, columnOrder);

				ClusteringUtils.saveImageTree(matrix, nwGenes, nwSamples, (nwSamples == null ? dataset.getSampleNames() : nwSamples.getLabels()), (nwGenes == null ? dataset.getFeatureNames() : nwGenes.getLabels()), imgFilename, true);
				imgFile = new File(imgFilename);
				if ( imgFile.exists() ) {
					result.addOutputItem(new Item(method + "_clustering_image", imgFile.getName(), method.toUpperCase() + " heatmap image (png format)", TYPE.IMAGE, Arrays.asList("CLUSTER"), new HashMap<String, String>(2), "Cluster images"));					
				} else {
					printError("execute" + method + "_clustering", "Error", "Error saving clustering image");					
				}
			} catch (IOException e) {
				printError("ioexception_execute" + method + "_clustering", "error saving clustering image", e.toString(), StringUtils.getStackTrace(e));
			}
		} else {

			if ( sampleClustering &&  nwSamples != null) {
				try {

					imgFilename = this.getOutdir() + "/samples." + method + ".png";
					ClusteringUtils.saveImageTree(nwSamples, "Clusters of samples",  imgFilename, false, false);
					imgFile = new File(imgFilename);
					if ( imgFile.exists() ) {
						result.addOutputItem(new Item(method + "_clustering_image", imgFile.getName(), method.toUpperCase() + " sample clustering image (png format)", TYPE.IMAGE, new ArrayList<String>(2), new HashMap<String, String>(2), "Cluster images"));					
					} else {
						printError("execute" + method + "_clustering", "Error", "error saving sample clustering image");										
					}
				} catch (IOException e) {
					printError("ioexception_execute" + method + "_clustering", "Error", "error saving sample clustering image");
				}
			}

			printWarning("ioexception_execute" + method + "_clustering", "Warning", "This release limits the tree images to 1000 genes");
		}


		//		if ( nwGenes != null ) {
		//
		//			if ( "kmeans".equalsIgnoreCase(method) ) {
		//				String redirectionTags = null; 
		//				List<String> redirectionInputs = null;
		//				File redirectionFile = null, list1 = null;
		//				File list2 = new File(outdir + "/rownames.txt");
		//				for(int i=1 ; i<=kvalue ; i++) {
		//					// preparing significant list (top and bottom)
		//					//
		//					list1 = new File(outdir + "/cluster_" + i + ".txt");
		//					if ( list1.exists() || list2.exists() ) {
		//						redirectionFile = new File(outdir + "/cluster_" + i + "_to_fatigo.redirection");
		//
		//						redirectionInputs = new ArrayList<String>();
		//
		//						redirectionInputs.add("comparison=list2list");
		//
		//						redirectionInputs.add("list1_wum_data=true");
		//						redirectionInputs.add("list1_databox=" + list1.getName() + " (cluster " + i + " from job $JOB_NAME)");
		//						redirectionInputs.add("list1=$JOB_FOLDER/" + list1.getName());
		//
		//						redirectionInputs.add("list2_wum_data=true");
		//						redirectionInputs.add("list2_databox=" + list2.getName() + " (all genes from job $JOB_NAME)");
		//						redirectionInputs.add("list2=$JOB_FOLDER/" + list2.getName());
		//
		//						redirectionInputs.add("duplicates=ref");
		//
		//						redirectionInputs.add("tool=fatigo");
		//						redirectionInputs.add("jobname=fatigo from kmeans cluster " + 1);						
		//						redirectionInputs.add("jobdescription=redirected from job $JOB_NAME");
		//
		//						try {
		//							IOUtils.write(redirectionFile.getAbsolutePath(), redirectionInputs);
		//						} catch (IOException e) {
		//							// TODO Auto-generated catch block
		//							e.printStackTrace();
		//						}								
		//
		//						if ( redirectionFile.exists() ) {
		//							redirectionTags = "REDIRECTION(" + redirectionFile.getName() + ":Send to FatiGO tool...)";
		//							result.addOutputItem(new Item("cluster_" + i + "_to_fatigo", "", "Send cluster " + i + " to FatiGO tool", TYPE.TEXT, StringUtils.toList(redirectionTags, ","), new HashMap<String, String>(2), "Continue processing"));
		//						}
		//					}
		//				}
		//			}
		//
		//		}


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
		if ( rowOrder == null && columnOrder == null ) {
			return inputMatrix;
		}

		DoubleMatrix matrix = new DoubleMatrix(inputMatrix.getRowDimension(), inputMatrix.getColumnDimension());
		int r, c;
		for(int row=0 ; row<matrix.getRowDimension() ; row++) {
			r = ( rowOrder == null) ? row : rowOrder[row];
			for(int col=0 ; col<matrix.getColumnDimension() ; col++) {
				c = ( columnOrder == null) ? col : columnOrder[col];
				//matrix.set(row, col, inputMatrix.get(rowOrder[row], columnOrder[col]));
				matrix.set(row, col, inputMatrix.get(r, c));
			}
		}				
		return matrix;
	}


	public static void createFatiGORedirectionFile(File redirectionFile, File topListFile, File bottomListFile) {

		List<String> redirectionInputs = new ArrayList<String>();

		if ( topListFile != null && bottomListFile != null && topListFile.exists() && bottomListFile.exists() ) {
			redirectionInputs.add("comparison=list2list");

			redirectionInputs.add("list1_wum_data=true");
			redirectionInputs.add("list1_databox=" + topListFile.getName() + " (top list from job $JOB_NAME)");
			redirectionInputs.add("list1=$JOB_FOLDER/clusters/" + topListFile.getName());

			redirectionInputs.add("list2_wum_data=true");
			redirectionInputs.add("list2_databox=" + bottomListFile.getName() + " (bottom list from job $JOB_NAME)");
			redirectionInputs.add("list2=$JOB_FOLDER/clusters/" + bottomListFile.getName());

			redirectionInputs.add("duplicates=ref");
		} else if ( topListFile != null && topListFile.exists() ) {
			redirectionInputs.add("comparison=list2genome");

			redirectionInputs.add("list1_wum_data=true");
			redirectionInputs.add("list1_databox=" + topListFile.getName() + " (top list from job $JOB_NAME)");
			redirectionInputs.add("list1=$JOB_FOLDER/clusters/" + topListFile.getName());
		} else if ( bottomListFile != null && bottomListFile.exists() ) {
			redirectionInputs.add("comparison=list2genome");

			redirectionInputs.add("list1_wum_data=true");
			redirectionInputs.add("list1_databox=" + bottomListFile.getName() + " (top list from job $JOB_NAME)");
			redirectionInputs.add("list1=$JOB_FOLDER/clusters/" + bottomListFile.getName());
		}		

		redirectionInputs.add("tool=fatigo");
		redirectionInputs.add("jobname=fatigo");
		redirectionInputs.add("jobdescription=redirected from job $JOB_NAME");
		try {
			IOUtils.write(redirectionFile.getAbsolutePath(), redirectionInputs);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}		
	}

}
