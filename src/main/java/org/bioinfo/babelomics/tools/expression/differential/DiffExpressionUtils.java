package org.bioinfo.babelomics.tools.expression.differential;

import java.awt.Color;
import java.text.DecimalFormat;
import java.util.ArrayList;
import java.util.List;

import org.bioinfo.collections.list.NamedArrayList;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.graphics.canvas.Canvas;
import org.bioinfo.graphics.canvas.feature.ScoreFeature;
import org.bioinfo.graphics.canvas.panel.GridPanel;
import org.bioinfo.graphics.canvas.track.GridTrack;

public class DiffExpressionUtils {

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

		GridPanel gridPanel = new GridPanel("", (columnDimension * cellSide) + rowLabelsWidth + infoWidth + canvas.getBorderPadding() + canvas.getSpaceSeparator(), ((rowDimension+1) * cellSide) + colLabelsWidth + canvas.getBorderPadding() + canvas.getSpaceSeparator(), xHeatMap, yHeatMap);
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

}
