package org.bioinfo.babelomics.methods.expression.clustering;

import java.awt.Color;
import java.io.IOException;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.data.list.NamedArrayList;
import org.bioinfo.data.tree.multiway.MultipleTree;
import org.bioinfo.graphics.canvas.Canvas;
import org.bioinfo.graphics.canvas.feature.ScoreFeature;
import org.bioinfo.graphics.canvas.panel.GridPanel;
import org.bioinfo.graphics.canvas.panel.NewickPanel;
import org.bioinfo.graphics.canvas.track.GridTrack;

public class ClusteringUtils {
	

	public static void saveImageTree(MultipleTree tree, String title, String imgFilename, boolean vertical) throws IOException {
				
		int width = vertical ? (tree.getNumberOfLevels() * 10 + 30) : (tree.getNumberOfLeaves() * 10);
		int height = vertical ? (tree.getNumberOfLeaves() * 10) : (tree.getNumberOfLevels() * 10 + 30);
		
		// inferior bounds, not too small
		//
		if ( width < 256 ) width = 256;
		if (height < 256) height = 256;
		
		// superior bounds, not too big
		//
		if (width > 8000) width = 8000;
		if (height > 8000) height = 8000;
		
		NewickPanel newickPanel = new NewickPanel(title, 0, 0, width, height);
		newickPanel.setNewick(tree);
		newickPanel.setShowLabels(true);
		newickPanel.setVertical(vertical);
				
		Canvas canvas = new Canvas("");
		canvas.setBorderWidth(0);
		canvas.setBorderPadding(0);
		canvas.setSpaceSeparator(0);
		canvas.setBorderColor(Color.BLACK);
		canvas.setBackGroundColor(Color.WHITE);
		canvas.setHeight(newickPanel.getHeight());
		canvas.setWidth(newickPanel.getWidth());
		
		canvas.addPanel(newickPanel);
		
		canvas.render(false);
		
		try {
			System.out.println("saving tree in : " + imgFilename);
			System.out.println("tree labels : " + ListUtils.toString(tree.getLabels(), ","));
			canvas.save(imgFilename);
		} catch (IOException e) {
			e.printStackTrace();
		}		
	}

//	public static void saveImageTree(DoubleMatrix matrix, MultipleTree vTree, MultipleTree hTree, String imgFilename) throws IOException {
//	}

	public static Canvas generateHeatmap(Dataset dataset, String className, int[] columnOrder, int[] rowOrder, String infoName1, double[] infoList1, String infoName2, double[] infoList2) {
		int xHeatMap = 2;				
		int yHeatMap = 2;
		int rowDimension = dataset.getRowDimension();
		int columnDimension = dataset.getColumnDimension();
		int cellSide = 16;
		int rowLabelsWidth = 70;
		int colLabelsWidth = 70;
		int infoWidth = 140;
		//double min = dataset.getDoubleMatrix().getMinValue();
		//double max = dataset.getDoubleMatrix().getMaxValue();
		double offset, min, max, mean, deviation, standard;
		//System.out.println("heatmap dimensions: (rowDimension, columnDimension) = (" + rowDimension + ", " + columnDimension + ")(min, max) = (" + min + ", " + max + ")");

		Canvas canvas = new Canvas("");
		canvas.setBorderWidth(0);
		canvas.setBorderColor(Color.BLACK);
		canvas.setBackGroundColor(Color.WHITE);

		GridPanel gridPanel = new GridPanel("", xHeatMap, yHeatMap, (columnDimension * cellSide) + rowLabelsWidth + infoWidth + canvas.getBorderPadding() + canvas.getSpaceSeparator(), ((rowDimension+1) * cellSide) + colLabelsWidth + canvas.getBorderPadding() + canvas.getSpaceSeparator());
		GridTrack gridTrack = new GridTrack(rowDimension, columnDimension, cellSide, cellSide);
		gridTrack.setRowLabels(ListUtils.ordered(dataset.getFeatureNames(), rowOrder));
		List<String> columnLabels = new ArrayList<String>(dataset.getSampleNames().size());
		for(int i=0 ; i <columnOrder.length ; i++) {
			columnLabels.add(dataset.getSampleNames().get(columnOrder[i]) + " (" + dataset.getVariables().getByName(className).getValues().get(columnOrder[i]) + ")");
		}
		gridTrack.setColumnLabels(columnLabels);
		//gridTrack.setName();
		gridTrack.setTopRegion(colLabelsWidth);
		gridTrack.setLeftRegion(rowLabelsWidth);
		gridTrack.setRightRegion(infoWidth);

		double[] values = new double[columnOrder.length]; 
		int row, column;
		
		ScoreFeature feature;
		for(int i=0 ; i<rowOrder.length ; i++) {
			row = rowOrder[i];
//			System.out.println("row = " + Arrays.toString(dataset.getDoubleMatrix().getRow(row)));
//			System.out.println("row mean = " + MathUtils.mean(dataset.getDoubleMatrix().getRow(row)));
//			System.out.println("row deviation = " + MathUtils.standardDeviation(dataset.getDoubleMatrix().getRow(row)));
//			System.exit(-1);
			
			mean = dataset.getDoubleMatrix().getRowMean(row);
			deviation = dataset.getDoubleMatrix().getRowStdDeviation(row);
			min = Double.MAX_VALUE;
			max = Double.MIN_NORMAL;
			for(int j=0 ; j<columnOrder.length ; j++) {
				column = columnOrder[j];
				values[column] = (deviation == 0) ? Double.NaN : (dataset.getDoubleMatrix().get(row, column)-mean)/(deviation);
				if ( min > values[column] ) min = values[column];
				if ( max < values[column] ) max = values[column];
			}
			
			offset = ( min <= 0 ) ? Math.abs(min) : (-1 * min);
			//System.out.println("mean = " + mean + ", deviation = " + deviation + ", min = " + min + ", max = " + max + ", offset = " + offset);
			for(int j=0 ; j<columnOrder.length ; j++) {
				column = columnOrder[j];
				//System.out.print("row, column = " + row + ", " + column + ": value = "); System.out.println(dataset.getDoubleMatrix().get(row, column));
//				feature = new ScoreFeature("name (" + column + ", " + row + ")", "", 0, 0, (dataset.getDoubleMatrix().get(row, column)-min)/(max-min));
				//standard = (deviation == 0) ? Double.NaN : (dataset.getDoubleMatrix().get(row, column)-mean)/(deviation);
				standard = (values[column] + offset) / ( max + offset);
				//System.out.println("(value, standard) = (" + dataset.getDoubleMatrix().get(row, column) + ", " + standard + ")");
				feature = new ScoreFeature("name, " + dataset.getDoubleMatrix().get(row, column), "", 0, 0, standard);
				//feature = new ScoreFeature("name (" + row + ", " + column + ")", "", 0, 0, dataset.getDoubleMatrix().get(row, column));
				//feature.setJsFunction("http://www.cipf.es");
				//				gridTrack.setFeature(row, column, feature);
				gridTrack.setFeature(i, j, feature);
			}
		}

		DecimalFormat df = new DecimalFormat("##0.0000");
		List<String> aux = new ArrayList<String>(infoList1.length);
		for(int i=0 ; i<infoList1.length ; i++) {
			aux.add(df.format(infoList1[i]));
		}
		gridTrack.addInfo(new NamedArrayList(infoName1, ListUtils.ordered(aux, rowOrder)));
		aux.clear();
		for(int i=0 ; i<infoList2.length ; i++) {
			aux.add(df.format(infoList2[i]));
		}
		gridTrack.addInfo(new NamedArrayList(infoName2, ListUtils.ordered(aux, rowOrder)));

		gridPanel.add(gridTrack);
		canvas.addPanel(gridPanel);		
		canvas.render();
		
		return canvas;
	}

	
//	public static void saveImageTree(DoubleMatrix matrix, BinaryTree vTree, BinaryTree hTree, String imgFilename) throws IOException {
//		
//		int cellSide = 20;
//		int rowLabelsWidth = getMaxStringLengh(vTree.getLabels()) * 9;
//		int colLabelsWidth = getMaxStringLengh(hTree.getLabels()) * 9;
//		int infoWidth = 0;
//		
//		
//		int rowDimension = vTree.getNumberOfLeaves();
//		int columnDimension = hTree.getNumberOfLeaves();
//				
//		System.out.println("sizes from trees: rows = " + rowDimension + ", cols = " + columnDimension);
//		System.out.println("sizes from matrix: rows = " + matrix.getRowDimension() + ", cols = " + matrix.getColumnDimension());
//		
//		NewickPanel newickHPanel = new NewickPanel("", hTree.getNumberOfLeaves() * cellSide, 
//													   hTree.getNumberOfLevels() * cellSide, 
//													   rowLabelsWidth + (vTree.getNumberOfLevels() * cellSide), 
//													   0);
//		NewickTrack nwTrack = new NewickTrack("", "", 0, Color.RED);
//		nwTrack.setNewick(hTree);
//		nwTrack.setLevelSeparation(cellSide);
//		nwTrack.setLeafSeparation(cellSide);
//		nwTrack.setShowLabels(false);
//		nwTrack.setVertical(false);
////		newickHPanel.setBorder(2);
//		newickHPanel.add(nwTrack);
//					
//		NewickPanel newickVPanel = new NewickPanel("", vTree.getNumberOfLevels() * cellSide, 
//													   vTree.getNumberOfLeaves() * cellSide, 
//													   0, 
//													   colLabelsWidth + (hTree.getNumberOfLevels() * cellSide));
//		nwTrack = new NewickTrack("", "", 0, Color.WHITE);
//		nwTrack.setNewick(vTree);
//		nwTrack.setLevelSeparation(cellSide);
//		nwTrack.setLeafSeparation(cellSide);
//		nwTrack.setShowLabels(false);
//		nwTrack.setVertical(true);
//		newickVPanel.add(nwTrack);
//
//		GridPanel gridPanel = new GridPanel("", (rowDimension * cellSide) + rowLabelsWidth + infoWidth, 
//												(columnDimension * cellSide) + colLabelsWidth, 
//												vTree.getNumberOfLevels() * cellSide, 
//												newickHPanel.getHeight());
////		GridPanel gridPanel = new GridPanel("", (rowDimension * cellSide) + rowLabelsWidth + 2 + infoWidth, 
////												(columnDimension * cellSide) + colLabelsWidth + 2, 
////												4 + x + (vTree.getNumberOfLevels()*cellSide), 
////												y + (hTree.getNumberOfLevels() * cellSide));
//		GridTrack gridTrack = new GridTrack(rowDimension, columnDimension, cellSide, cellSide);
//		gridTrack.setName("g r i d     t r a c k     n a m e");
//		gridTrack.setColumnLabels(hTree.getLabels());
//		gridTrack.setRowLabels(vTree.getLabels());
//		gridTrack.setTopRegion(colLabelsWidth);
//		gridTrack.setLeftRegion(rowLabelsWidth);
//		gridTrack.setRightRegion(infoWidth);
//		ScoreFeature feature;
//		
//		double mean, deviation, min, max, offset, standard;
//		double[] values = new double[gridTrack.getColumnDimension()];
//		for(int row=0 ; row<gridTrack.getRowDimension() ; row++) {
//			mean = matrix.getRowMean(row);
//			deviation = matrix.getRowStdDeviation(row);
//			min = Double.MAX_VALUE;
//			max = Double.MIN_NORMAL;
//			for(int column=0 ; column<gridTrack.getColumnDimension() ; column++) {
//				values[column] = (deviation == 0) ? Double.NaN : (matrix.get(row, column)-mean)/(deviation);
//				if ( min > values[column] ) min = values[column];
//				if ( max < values[column] ) max = values[column];
//			}
//			offset = ( min <= 0 ) ? Math.abs(min) : (-1 * min);
//			for(int column=0 ; column<gridTrack.getColumnDimension() ; column++) {
//				standard = (values[column] + offset) / ( max + offset);
//				
////				System.out.print(matrix.get(row, column) + "\t");
////				feature = new ScoreFeature("name (" + column + ", " + row + ")", "description bla, bla, bla", 0, 0, matrix.get(row, column));
//				feature = new ScoreFeature("name (" + column + ", " + row + ")", "description bla, bla, bla", 0, 0, standard);
//				gridTrack.setFeature(row, column, feature);
//			}
////			System.out.println("");
//		}
//		gridPanel.add(gridTrack);
//
//		Canvas canvas = new Canvas("");
//		canvas.setBorderWidth(0);
//		canvas.setBorderPadding(4);
//		canvas.setSpaceSeparator(0);
//		canvas.setBorderColor(Color.BLACK);
//		canvas.setBackGroundColor(Color.WHITE);
//
//		int canvasHeight = gridPanel.getWidth() + newickHPanel.getHeight() + cellSide;
//		int canvasWidth = gridPanel.getHeight() + newickVPanel.getWidth() + cellSide;
//		System.out.println("canvas height = " + canvasHeight);
//		System.out.println("canvas width  = " + canvasWidth);
//
//		canvas.setHeight(canvasHeight);
//		canvas.setWidth(canvasWidth);
//		
//		
//		canvas.addPanel(gridPanel);
//		canvas.addPanel(newickHPanel);
//		canvas.addPanel(newickVPanel);
//		
//		canvas.render();
//		canvas.save(imgFilename);		
//	}
//
//	private static int getMaxStringLengh(List<String> names) {
//		int max = 0;
//		for(String name: names) {
//			if ( name.length() > max ) {
//				max = name.length();
//			}
//		}
//		return max;
//	}

}
