package org.bioinfo.babelomics.tools.expression.clustering;

import java.awt.Color;
import java.io.IOException;
import java.util.List;

import org.bioinfo.data.format.core.newick.NewickTree;
import org.bioinfo.graphics.canvas.Canvas;
import org.bioinfo.graphics.canvas.feature.ScoreFeature;
import org.bioinfo.graphics.canvas.panel.GridPanel;
import org.bioinfo.graphics.canvas.panel.NewickPanel;
import org.bioinfo.graphics.canvas.track.GridTrack;
import org.bioinfo.graphics.canvas.track.NewickTrack;
import org.bioinfo.math.data.DoubleMatrix;

public class ClusteringUtils {
	
	public static void saveImageTree(DoubleMatrix matrix, NewickTree vTree, NewickTree hTree, String imgFilename) throws IOException {
		
		int cellSide = 20;
		int rowLabelsWidth = getMaxStringLengh(vTree.getLabels()) * 9;
		int colLabelsWidth = getMaxStringLengh(hTree.getLabels()) * 9;
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
//		newickHPanel.setBorder(2);
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

	private static int getMaxStringLengh(List<String> names) {
		int max = 0;
		for(String name: names) {
			if ( name.length() > max ) {
				max = name.length();
			}
		}
		return max;
	}

}
