package org.bioinfo.babelomics.tools;

import java.awt.Color;
import java.awt.Panel;
import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import org.apache.commons.cli.ParseException;
import org.bioinfo.babelomics.methods.functional.FatiScan;
import org.bioinfo.collections.exceptions.InvalidColumnIndexException;
import org.bioinfo.collections.matrix.DataFrame;
import org.bioinfo.commons.exec.Command;
import org.bioinfo.commons.exec.SingleProcess;
import org.bioinfo.commons.io.utils.FileUtils;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ArrayUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.data.dataset.FeatureData;
import org.bioinfo.math.data.DoubleMatrix;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;
import org.bioinfo.graphics.canvas.Canvas;
import org.bioinfo.graphics.canvas.feature.AnnotationFeature;
import org.bioinfo.graphics.canvas.feature.ScoreFeature;
import org.bioinfo.graphics.canvas.panel.AnnotationPanel;
import org.bioinfo.graphics.canvas.panel.XYPanel;
import org.bioinfo.graphics.canvas.track.AnnotationTrack;
import org.bioinfo.graphics.canvas.track.XYTrack;
import org.jfree.chart.renderer.xy.XYBarRenderer;
import org.jfree.data.xy.XYBarDataset;
import org.jfree.data.xy.XYDataset;

public class Histogram extends BabelomicsTool {

	@Override
	public void initOptions() {
		options.addOption(OptionFactory.createOption("ranked-list", "the feature data containig the ranked list"));
	}

	
		
	@Override
	protected void execute() {
		Dataset dataset = null;
		int progress = 1;
		int finalProgress = 3;
		
		try {
			jobStatus.addStatusMessage("" + (progress*100/finalProgress), "reading ranked list");
		} catch (FileNotFoundException e) {
			abort("filenotfoundexception_execute_preprocessing", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
		}

		// ranked list
		FeatureData rankedList = new FeatureData(new File(commandLine.getOptionValue("ranked-list")), true);
		
		// save id lists				
		IOUtils.write(outdir + "/ranked_list.txt", rankedList.toString());
		result.addOutputItem(new Item("ranked_list","ranked_list.txt","Ranked list",Item.TYPE.FILE,Arrays.asList("RANKED_LIST","CLEAN"),new HashMap<String,String>(),"Input data"));
		jobStatus.addStatusMessage("" + (progress*100/finalProgress), "reading ok");
		
		if ( rankedList == null ) {
			abort("rankedisnull_execute_historgam", "ranked is null", "ranked is null after reading file '", "ranked is null after reading file ");
		}
		
		double barWidth = 15;
		XYBarDataset xyBarDataset = new XYBarDataset((XYDataset) rankedList , barWidth);
		Canvas canvas = new Canvas("Histrogram");
		canvas.setBorderWidth(1);
		canvas.setBorderColor(Color.BLACK);
		canvas.setBackGroundColor(Color.WHITE);
		
		
		XYPanel xyPanel;
		XYTrack xtTrack;
		DataFrame dataFrame = rankedList.getDataFrame();
//		for (int i = 0;i < xyBarDataset.getSeriesCount(); i++){
//			 
//			//ScoreFeature sf = xyBarDataset.getYValue(i, 0);
//			ScoreFeature sf;
//			sf.setScore(xyBarDataset.getYValue(i, 0));
//			sf.setName(rankedList.getDataFrame());
//			
//			xtTrack.add(sf);
//			
//		}
		FeatureData fd= new FeatureData(dataFrame);
		
//		xtTrack.add(scoreFeature);
//		xyPanel.addXYTrack(fd);
		
		
		progress++;

		jobStatus.addStatusMessage("" + (progress*100/finalProgress), "making graph");
		
		// apply logarithm
		//
//		logger.debug("executing logarithm base " + logBase + "...\n");
//		logger.debug("end of executing logarithm base " + logBase + "\n");
			progress++;
	}
}