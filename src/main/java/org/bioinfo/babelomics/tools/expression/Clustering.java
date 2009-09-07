package org.bioinfo.babelomics.tools.expression;

import java.awt.Color;
import java.io.File;
import java.io.FileNotFoundException; 
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.commons.exec.Command;
import org.bioinfo.commons.exec.SingleProcess;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.data.format.core.newick.NewickTree;
import org.bioinfo.data.format.io.NewickParser;
import org.bioinfo.data.format.io.exception.InvalidFormatException;
import org.bioinfo.graphics.canvas.Canvas;
import org.bioinfo.graphics.canvas.feature.ScoreFeature;
import org.bioinfo.graphics.canvas.panel.GridPanel;
import org.bioinfo.graphics.canvas.panel.NewickPanel;
import org.bioinfo.graphics.canvas.track.GridTrack;
import org.bioinfo.graphics.canvas.track.NewickTrack;
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
				IOUtils.write(new File(this.getOutdir() + "/samples.nw"), nwGenes.toString());
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
				
				saveImageTree(matrix, nwGenes, nwSamples, imgFilename);
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


	//--------------------------------------------------------------------------------------
	// The following functions must be located in another library for re-usability purposes
	//--------------------------------------------------------------------------------------
	

	private NewickTree runClustering(DoubleMatrix matrix, List<String> rowNames, List<String> colNames, String method, String distance, int kvalue) throws IOException, InvalidFormatException {
		NewickTree tree = null;

		File inputFile = File.createTempFile("input", null);
		File outputFile = File.createTempFile("output", null);

		System.out.println("(infile, outfile) = (" + inputFile.getAbsolutePath() + ", " + outputFile.getAbsolutePath() + ")");
		
		List<String> lines = new ArrayList<String>(rowNames.size() + 1);
		lines.add("#NAMES\t" + ListUtils.toString(colNames, "\t"));
		for(int i=0 ; i<rowNames.size() ; i++) {
			lines.add(rowNames.get(i) + "\t" + ListUtils.toString(ListUtils.toList(matrix.getRow(i)), "\t"));
		}
		IOUtils.write(inputFile, lines);
		
		String cmdStr = null;
		if ( "sota".equalsIgnoreCase(method) ) {
			cmdStr = System.getenv("BABELOMICS_HOME") + "/bin/sota " + inputFile.getAbsolutePath() + " " + outputFile.getAbsolutePath() + " " + distance + " -newick";
		} else if ( "som".equalsIgnoreCase(method) ) {
			
		} else if ( "upgma".equalsIgnoreCase(method) ) {
			cmdStr = System.getenv("BABELOMICS_HOME") + "/bin/cluster " + inputFile.getAbsolutePath() + " " + outputFile.getAbsolutePath() + " UPGMA " + distance;			
		} else if ( "kmeans".equalsIgnoreCase(method) ) {
			
		}
		
		Command cmd = new Command(cmdStr); 
		SingleProcess sp = new SingleProcess(cmd);
		sp.runSync();
		
		if ( outputFile.exists() && outputFile.getTotalSpace() > 0 ) {
			tree = new NewickParser().parse(IOUtils.toString(outputFile));
		}
		
		//inputFile.delete();
		//outputFile.delete();
		
		return tree;
	}
		
	
	private void saveImageTree(DoubleMatrix matrix, NewickTree vTree, NewickTree hTree, String imgFilename) throws IOException {
		
		int cellSide = 20;
		int rowLabelsWidth = getMaxStringLengh(vTree.getLabels())*9;
		int colLabelsWidth = getMaxStringLengh(hTree.getLabels())*9;
		int infoWidth = 0;
		
		
		int rowDimension = vTree.getNumberOfLeaves();
		int columnDimension = hTree.getNumberOfLeaves();
				
		System.out.println("sizes from trees: rows = " + rowDimension + ", cols = " + columnDimension);
		System.out.println("sizes from matrix: rows = " + matrix.getRowDimension() + ", cols = " + matrix.getColumnDimension());
		
		NewickPanel newickHPanel = new NewickPanel("", hTree.getNumberOfLeaves() * cellSide, 
													   hTree.getNumberOfLevels() * cellSide, 
													   rowLabelsWidth + (vTree.getNumberOfLevels() * cellSide), 
													   0);
		NewickTrack nwTrack = new NewickTrack("", "", 0, Color.RED);
		nwTrack.setNewick(hTree);
		nwTrack.setLevelSeparation(cellSide);
		nwTrack.setLeafSeparation(cellSide);
		nwTrack.setShowLabels(false);
		nwTrack.setVertical(false);
		newickHPanel.setBorder(2);
		newickHPanel.add(nwTrack);
					
		NewickPanel newickVPanel = new NewickPanel("", vTree.getNumberOfLevels() * cellSide, 
													   vTree.getNumberOfLeaves() * cellSide, 
													   0, 
													   colLabelsWidth + (hTree.getNumberOfLevels() * cellSide));
		nwTrack = new NewickTrack("", "", 0, Color.WHITE);
		nwTrack.setNewick(vTree);
		nwTrack.setLevelSeparation(cellSide);
		nwTrack.setLeafSeparation(cellSide);
		nwTrack.setShowLabels(false);
		nwTrack.setVertical(true);
		newickVPanel.add(nwTrack);

		GridPanel gridPanel = new GridPanel("", (rowDimension * cellSide) + rowLabelsWidth + infoWidth, 
												(columnDimension * cellSide) + colLabelsWidth, 
												vTree.getNumberOfLevels() * cellSide, 
												newickHPanel.getHeight());
//		GridPanel gridPanel = new GridPanel("", (rowDimension * cellSide) + rowLabelsWidth + 2 + infoWidth, 
//												(columnDimension * cellSide) + colLabelsWidth + 2, 
//												4 + x + (vTree.getNumberOfLevels()*cellSide), 
//												y + (hTree.getNumberOfLevels() * cellSide));
		GridTrack gridTrack = new GridTrack(rowDimension, columnDimension, cellSide, cellSide);
		gridTrack.setName("g r i d     t r a c k     n a m e");
		gridTrack.setColumnLabels(hTree.getLabels());
		gridTrack.setRowLabels(vTree.getLabels());
		gridTrack.setTopRegion(colLabelsWidth);
		gridTrack.setLeftRegion(rowLabelsWidth);
		gridTrack.setRightRegion(infoWidth);
		ScoreFeature feature;
		
		double mean, deviation, min, max, offset, standard;
		double[] values = new double[gridTrack.getColumnDimension()];
		for(int row=0 ; row<gridTrack.getRowDimension() ; row++) {
			mean = matrix.getRowMean(row);
			deviation = matrix.getRowStdDeviation(row);
			min = Double.MAX_VALUE;
			max = Double.MIN_NORMAL;
			for(int column=0 ; column<gridTrack.getColumnDimension() ; column++) {
				values[column] = (deviation == 0) ? Double.NaN : (matrix.get(row, column)-mean)/(deviation);
				if ( min > values[column] ) min = values[column];
				if ( max < values[column] ) max = values[column];
			}
			offset = ( min <= 0 ) ? Math.abs(min) : (-1 * min);
			for(int column=0 ; column<gridTrack.getColumnDimension() ; column++) {
				standard = (values[column] + offset) / ( max + offset);
				
				System.out.print(matrix.get(row, column) + "\t");
//				feature = new ScoreFeature("name (" + column + ", " + row + ")", "description bla, bla, bla", 0, 0, matrix.get(row, column));
				feature = new ScoreFeature("name (" + column + ", " + row + ")", "description bla, bla, bla", 0, 0, standard);
				gridTrack.setFeature(row, column, feature);
			}
			System.out.println("");
		}
		gridPanel.add(gridTrack);

		Canvas canvas = new Canvas("");
		canvas.setBorderWidth(0);
		canvas.setBorderPadding(4);
		canvas.setSpaceSeparator(0);
		canvas.setBorderColor(Color.BLACK);
		canvas.setBackGroundColor(Color.WHITE);

		int canvasHeight = gridPanel.getWidth() + newickHPanel.getHeight() + cellSide;
		int canvasWidth = gridPanel.getHeight() + newickVPanel.getWidth() + cellSide;
		System.out.println("canvas height = " + canvasHeight);
		System.out.println("canvas width  = " + canvasWidth);

		canvas.setHeight(canvasHeight);
		canvas.setWidth(canvasWidth);
		
		
		canvas.addPanel(gridPanel);
		canvas.addPanel(newickHPanel);
		canvas.addPanel(newickVPanel);
		
		canvas.render();
		canvas.save(imgFilename);		
	}

	private int getMaxStringLengh(List<String> names) {
		int max = 0;
		for(String name: names) {
			if ( name.length() > max ) {
				max = name.length();
			}
		}
		return max;
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
