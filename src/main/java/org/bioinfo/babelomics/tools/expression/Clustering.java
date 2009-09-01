package org.bioinfo.babelomics.tools.expression;

import java.awt.Color;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.math.linear.RealMatrix;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.collections.newicktree.NewickTree;
import org.bioinfo.commons.exec.Command;
import org.bioinfo.commons.exec.SingleProcess;
import org.bioinfo.commons.io.TextFileWriter;
import org.bioinfo.commons.io.utils.FileSystemUtils;
import org.bioinfo.commons.io.utils.FileUtils;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ArrayUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.graphics.canvas.Canvas;
import org.bioinfo.graphics.canvas.feature.ScoreFeature;
import org.bioinfo.graphics.canvas.panel.GridPanel;
import org.bioinfo.graphics.canvas.panel.NewickPanel;
import org.bioinfo.graphics.canvas.track.GridTrack;
import org.bioinfo.graphics.canvas.track.NewickTrack;
import org.bioinfo.io.parser.NewickParser;
import org.bioinfo.io.parser.exception.InvalidFormatException;
import org.bioinfo.math.data.DoubleMatrix;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;

public class Clustering extends BabelomicsTool {


	public Clustering() {
	}

	@Override
	public void initOptions() {
		options.addOption(OptionFactory.createOption("dataset", "the data"));
		options.addOption(OptionFactory.createOption("method", "the method, possible values: upgma, sota, som, kmeans"));
		options.addOption(OptionFactory.createOption("distance", "the distance, possible values: euclidean, spearman, pearson. Default value: euclidean, false"));
		options.addOption(OptionFactory.createOption("kvalue", "k-value for kmeans clustering. Default value: 15", false));
		options.addOption(OptionFactory.createOption("sample-filter", "class variable", false));
		options.addOption(OptionFactory.createOption("feature-filter", "class variable", false));

	}

	@Override
	public void execute() {
		Dataset dataset = null;
		String method = commandLine.getOptionValue("method");

		String distance = commandLine.getOptionValue("distance", "euclidean");
		String kvalue = commandLine.getOptionValue("time-class", "15");

		String datasetPath = commandLine.getOptionValue("dataset");
		if ( datasetPath == null ) {
			printError("missingdataset_execute_clustering", "Missing dataset", "Missing dataset");
			abort("Missing dataset", "Missing dataset");
		}

		if ( method == null ) {
			printError("missingclusteringmethod_execute_clustering", "Missing clustering method", "Missing clustering method");
			abort("Missing clustering method", "Missing clustering method");
		}

		if ( !"upgma".equalsIgnoreCase(method) && !"sota".equalsIgnoreCase(method) && !"som".equalsIgnoreCase(method) && !"kmeans".equalsIgnoreCase(method) ) {
			printError("unknowncluteringmethod_execute_clustering", "Unknown clustering method", "Unknown clustering method");
			abort("Unknown clustering method", "Unknown clustering method");			
		}
		
		try {
			jobStatus.addStatusMessage("20", "reading dataset");
		} catch (FileNotFoundException e) {
			printError("filenotfoundexception_execute_clustering", "job status file not found", e.toString(), e);
			abort("Job status file not found", e.toString());
		}

		File datasetFile = new File(commandLine.getOptionValue("dataset"));
		try {
			dataset = new Dataset(datasetFile);
		} catch (Exception e) {
			printError("exception_execute_clustering", "error reading dataset '" + datasetFile.getName() + "'", e.toString(), e);
			abort("Error reading dataset '" + datasetFile.getName() + "'", e.toString());
		}		

		if ( dataset.getDoubleMatrix() == null ) { 
			try {
				dataset.load();
			} catch (Exception e) {
				printError("exception_execute_clustering", "Error loading dataset '" + datasetFile.getName() + "'", e.toString(), e);
				abort("Error loading dataset '" + datasetFile.getName() + "'", e.toString());
			}
			dataset.validate();
		}

//		if(commandLine.hasOption("sample-filter") || commandLine.hasOption("feature-filter")) {
//			dataset = dataset.getSubDataset(commandLine.getOptionValue("sample-filter"), "4", commandLine.getOptionValue("feature-filter"), ""); 
//		}

		if ( "upgma".equalsIgnoreCase(method) ) {
			
			executeUpgma(dataset, distance);
			
		} else if ( "sota".equalsIgnoreCase(method) ) {
			
			executeSota(dataset, distance);
			
		} else if ( "som".equalsIgnoreCase(method) ) {
			
			executeSom(dataset, distance);
			
		} else if ( "kmeans".equalsIgnoreCase(method) ) {
			
			int k;
			try {
				k = Integer.parseInt(kvalue);
				executeKmeans(dataset, distance, k);
			} catch (NumberFormatException e ) {
				logger.error("Invalid k-value: " + kvalue);
			}
		}
		
		try {
			jobStatus.addStatusMessage("100", "done");
		} catch (FileNotFoundException e) {
			printError("filenotfoundexception_execute_clustering", "job status file not found", e.toString(), e);
			abort("Job status file not found", e.toString());
		}

	}


	private void executeUpgma(Dataset dataset, String distance) {
		logger.info("executing upgma, not implemented yet");
	}

	private void executeSom(Dataset dataset, String distance) {
		logger.info("executing som, not implemented yet");
	}

	private void executeKmeans(Dataset dataset, String distance, int kvalue) {
		logger.info("executing kmeans, not implemented yet");
	}

	private void executeSota(Dataset dataset, String distance) {

		NewickTree nwGenes = null, nwSamples = null;
		
		try {
			jobStatus.addStatusMessage("40", "generating genes clusters");
		} catch (FileNotFoundException e) {
			printError("filenotfoundexception_execute_clustering", "job status file not found", e.toString(), e);
			abort("Job status file not found", e.toString());
		}

		try {
			nwGenes = runSota(dataset.getDoubleMatrix(), dataset.getFeatureNames(), dataset.getSampleNames(), distance);
		} catch (Exception e) {
			printError("exception_executesota_clustering", "error running sota algorithm for genes", e.toString(), e);
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
			printError("filenotfoundexception_execute_clustering", "job status file not found", e.toString(), e);
			abort("Job status file not found", e.toString());
		}

		try {
			nwSamples = runSota(new DoubleMatrix(dataset.getDoubleMatrix().transpose().getData()), dataset.getSampleNames(), dataset.getFeatureNames(), distance);
		} catch (Exception e) {
			printError("exception_executesota_clustering", "error running sota algorithm for samples", e.toString(), e);
		}
		
		if ( nwSamples != null ) {
			try {
				IOUtils.write(new File(this.getOutdir() + "/samples.nw"), nwGenes.toString());
				result.addOutputItem(new Item("sample_newick_file", "samples.nw", "Clusters of samples (newick format file)", TYPE.FILE));
			} catch (IOException e) {
				printError("ioexception_executesota_clustering", "error saving samples newick", e.toString(), e);
				nwSamples = null;
			}			
		}
		
		try {
			jobStatus.addStatusMessage("80", "generating clustering image");
		} catch (FileNotFoundException e) {
			printError("filenotfoundexception_execute_clustering", "job status file not found", e.toString(), e);
			abort("Job status file not found", e.toString());
		}

		if ( nwGenes != null && nwSamples != null ) {
			try {
				String imgFilename = this.getOutdir() + "/sota";
				saveImageTree(dataset.getDoubleMatrix(), nwGenes, nwSamples, imgFilename);
				File imgFile = new File(imgFilename + ".png");
				if ( imgFile.exists() ) {
					result.addOutputItem(new Item("sota_clustering_image", "sota.png", "SOTA clustering image (png format)", TYPE.IMAGE));					
				} else {
					printError("executesota_clustering", "error saving clustering image", "error saving clustering image");					
				}
			} catch (IOException e) {
				printError("ioexception_executesota_clustering", "error saving clustering image", e.toString(), e);
			}
		}
	}


	
	
	
	
	//--------------------------------------------------------------------------------------
	// The following functions must be located in another library for re-usability purposes
	//--------------------------------------------------------------------------------------
	

	private NewickTree runSota(DoubleMatrix matrix, List<String> rowNames, List<String> colNames, String distance) throws IOException, InvalidFormatException {
		NewickTree tree = null;

		File inputFile = File.createTempFile("input", null);
		File outputFile = File.createTempFile("output", null);

		List<String> lines = new ArrayList<String>(rowNames.size() + 1);
		lines.add("#NAMES\t" + ListUtils.toString(colNames, "\t"));
		for(int i=0 ; i<rowNames.size() ; i++) {
			lines.add(rowNames.get(i) + "\t" + ListUtils.toString(ListUtils.toList(matrix.getRow(i)), "\t"));
		}
		IOUtils.write(inputFile, lines);
		
		String cmdStr = System.getenv("BABELOMICS_HOME") + "/bin/sota " + inputFile.getAbsolutePath() + " " + outputFile.getAbsolutePath() + " " + distance + " -newick";
		Command cmd = new Command(cmdStr); 
		SingleProcess sp = new SingleProcess(cmd);
		sp.runSync();
		
		if ( outputFile.exists() && outputFile.getTotalSpace() > 0 ) {
			tree = new NewickParser().parse(IOUtils.toString(outputFile));
		}
		
		inputFile.delete();
		outputFile.delete();
		
		return tree;
	}
		
	
	private void saveImageTree(DoubleMatrix matrix, NewickTree vTree, NewickTree hTree, String imgFilename) throws IOException {
		
		int x = 2;				
		int y = 2;
		int cellSide = 20;
		int rowLabelsWidth = 50;
		int colLabelsWidth = 50;
		int infoWidth = 0;
		
		int rowDimension = vTree.getNumberOfLeaves();
		int columnDimension = hTree.getNumberOfLeaves();
		
		Canvas canvas = new Canvas("");
		canvas.setBorderWidth(0);
		canvas.setBorderPadding(0);
		canvas.setBorderColor(Color.BLACK);
		canvas.setBackGroundColor(Color.WHITE);
		canvas.setHeight((rowDimension + hTree.getNumberOfLevels()) * cellSide + colLabelsWidth + 100);
		canvas.setWidth((columnDimension + vTree.getNumberOfLevels()) * cellSide + rowLabelsWidth + 100);
		
		System.out.println("sizes from trees: rows = " + rowDimension + ", cols = " + columnDimension);
		System.out.println("sizes from matrix: rows = " + matrix.getRowDimension() + ", cols = " + matrix.getColumnDimension());
		
		NewickPanel newickHPanel = new NewickPanel("", rowDimension*cellSide+2, (hTree.getNumberOfLevels() * cellSide)+2, 2 + x + (vTree.getNumberOfLevels()*cellSide) + rowLabelsWidth, y);
		NewickTrack nwTrack = new NewickTrack("", "", 0, Color.WHITE);
		nwTrack.setNewick(hTree);
		nwTrack.setLevelSeparation(cellSide);
		nwTrack.setLeafSeparation(cellSide);
		nwTrack.setShowLabels(false);
		nwTrack.setVertical(false);
		newickHPanel.add(nwTrack);
					
		NewickPanel newickVPanel = new NewickPanel("", (vTree.getNumberOfLevels() * cellSide)+2, (rowDimension * cellSide)+2, x, y + (hTree.getNumberOfLevels() * cellSide) + colLabelsWidth);
		nwTrack = new NewickTrack("", "", 0, Color.WHITE);
		nwTrack.setNewick(vTree);
		nwTrack.setLevelSeparation(cellSide);
		nwTrack.setLeafSeparation(cellSide);
		nwTrack.setShowLabels(false);
		nwTrack.setVertical(true);
		newickVPanel.add(nwTrack);

		GridPanel gridPanel = new GridPanel("", (rowDimension * cellSide) + rowLabelsWidth + 2 + infoWidth, 
												(columnDimension * cellSide) + colLabelsWidth + 2, 
												4 + x + (vTree.getNumberOfLevels()*cellSide), 
												y + (hTree.getNumberOfLevels() * cellSide));
		GridTrack gridTrack = new GridTrack(rowDimension, columnDimension, cellSide, cellSide);
		gridTrack.setName("g r i d     t r a c k     n a m e");
		gridTrack.setColumnLabels(hTree.getLabels());
		gridTrack.setRowLabels(vTree.getLabels());
		gridTrack.setTopRegion(colLabelsWidth);
		gridTrack.setLeftRegion(rowLabelsWidth);
		gridTrack.setRightRegion(infoWidth);
		ScoreFeature feature;
		
		
		
		for(int row=0 ; row<gridTrack.getRowDimension() ; row++) {
			for(int column=0 ; column<gridTrack.getColumnDimension() ; column++) {
				
				System.out.println("(" + row + ", " + column + ") = " + matrix.get(row, column));
//				feature = new ScoreFeature("name (" + column + ", " + row + ")", "description bla, bla, bla", 0, 0, matrix.get(row, column));
				feature = new ScoreFeature("name (" + column + ", " + row + ")", "description bla, bla, bla", 0, 0, Math.random());
				gridTrack.setFeature(row, column, feature);
			}
		}
		gridPanel.add(gridTrack);
					
		canvas.addPanel(gridPanel);
		canvas.addPanel(newickHPanel);
		canvas.addPanel(newickVPanel);
		
		canvas.render();
		canvas.save(imgFilename);		
	}	
}
