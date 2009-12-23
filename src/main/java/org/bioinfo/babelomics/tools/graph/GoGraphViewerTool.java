package org.bioinfo.babelomics.tools.graph;

import java.io.File;
import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.security.InvalidParameterException;
import java.util.Properties;

import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.tool.OptionFactory;

import es.blast2go.prog.graph.GetGraphApi;

public class GoGraphViewerTool  extends BabelomicsTool{

	@Override
	public void initOptions() {				
		options.addOption(OptionFactory.createOption("association-file", "the association file containing db ids"));		
		// optional
		options.addOption(OptionFactory.createOption("go-domain", "[b,m,c] biological process, molecular function, cellular component (b by default)",false,true));
		options.addOption(OptionFactory.createOption("ids-per-node-filter", "[0,1000] Number of IDs per node filter (0 by default)",false,true));
		options.addOption(OptionFactory.createOption("annot-score-node-filter", "[0,1000] Thin out graphs by the annot-score (0 by default)",false,true));
		options.addOption(OptionFactory.createOption("annot-score-parameter", "[0,1] Alpha of the annot-score formula (0.6 by default)",false,true));
		options.addOption(OptionFactory.createOption("graph-coloring", "[annot-score, id-count, input-score] Mode of graph coloring",false,true));		
	}

	@Override
	protected void execute() {		
		
		try {
			// association file		
			String associationFile = commandLine.getOptionValue("association-file");
			if(!new File(associationFile).exists()) {
				throw new FileNotFoundException(associationFile + " not found");
			}
			
			// go domain
			String goDomain = "b";		
			if(commandLine.hasOption("go-domain")) goDomain = commandLine.getOptionValue("go-domain");
			if(!goDomain.equals("b") && !goDomain.equals("m") && !goDomain.equals("c")) throw new InvalidParameterException("value " + goDomain + " for go-domain parameter is not valid");
	
			// ids per node filter
			int idsPerNodeFilter = 0;		
			if(commandLine.hasOption("ids-per-node-filter")) idsPerNodeFilter = Integer.parseInt(commandLine.getOptionValue("ids-per-node-filter"));
			if(idsPerNodeFilter<0 | idsPerNodeFilter>1000) throw new InvalidParameterException("value " + idsPerNodeFilter + " for ids-per-node-filter parameter is out of range");
			
			// annots-score node filter
			int annotScoreNodeFilter = 0;		
			if(commandLine.hasOption("annot-score-node-filter")) annotScoreNodeFilter = Integer.parseInt(commandLine.getOptionValue("annot-score-node-filter"));
			if(annotScoreNodeFilter<0 | annotScoreNodeFilter>1000) throw new InvalidParameterException("value " + annotScoreNodeFilter + " for annot-score-node-filter parameter is out of range");
			
			// annots-score node filter
			double annotScoreParameter = 0.6;
			if(commandLine.hasOption("annot-score-parameter")) annotScoreParameter = Double.parseDouble(commandLine.getOptionValue("annot-score-parameter"));
			if(annotScoreParameter<0 | annotScoreParameter>1) throw new InvalidParameterException("value " + annotScoreParameter + " for annot-score-parameter parameter is out of range");
			
			String graphColoring = "annot-score";
			if(commandLine.hasOption("graph-coloring")) graphColoring = commandLine.getOptionValue("graph-coloring");
			if(!graphColoring.equals("annot-score") && !graphColoring.equals("id-count") &&  !graphColoring.equals("input-score")) throw new InvalidParameterException("value " + graphColoring + " for graph-coloring parameter is invalid");
			
			Properties config = getConfig();
			if(config!=null){				
				GetGraphApi graph = new GetGraphApi(outdir + "/","chaval",associationFile,goDomain,idsPerNodeFilter,graphColoring,annotScoreParameter,annotScoreNodeFilter, "orange");
				graph.setDownloader(config.getProperty("JNLP_DOWNLOADER_HOST_NAME"));
				graph.setDataBase(config.getProperty("BLAST2GO_HOST_NAME"),config.getProperty("BLAST2GO_DB_NAME"),config.getProperty("BLAST2GO_DB_USER"), config.getProperty("BLAST2GO_DB_PASSWORD"));
				graph.run();
			} else {
				throw new Exception("Properties file $BABELOMICS_HOME/conf/blast2go.properties cannot be read");
			}
					
		} catch (Exception e){
			e.printStackTrace();
		}
	
		
//		GetGraphApi graph = new GetGraphApi(WORKDIR, validpref, associationFile, ontology, IdNodeFilter, modeGraphColoring, annotScoreParam, annotScoreNodeFilter, "orange");
//		graph.setDownloader(config.get("JNLP_DOWNLOADER_HOST_NAME"));										
//		graph.setDataBase(config.get("BLAST2GO_HOST_NAME"),config.get("BLAST2GO_DB_NAME"),config.get("BLAST2GO_DB_USER"),config.get("BLAST2GO_DB_PASSWORD"));
		
		
	}

//	private String makeResultXML() {		
//		//Meta nodes
//		//result.addMetaItem(new Date().toString(), "date", "date", "Date", "Metadata", "");       
//		//Result Tables group	
//		t.addOutputItem(jnlpdownloader+jnlpLink, "jnlplink", "url", "Link to interactive GoGraphViz applicacion (open with Java-Web-Start)", "Images", "");
//		t.addOutputItem(imageJpg, "imageJpg", "file", "Graph as JPG image (high resolution)", "Images", "");
//		t.addOutputItem(imagePng, "imagePng", "file", "Graph as PNG", "Images", "");
//		t.addOutputItem(imageSVG, "imageSVG", "file", "Graph as SVG", "Images", "");
//		t.addOutputItem(imageJpgsmall, "imageJpgsmall", "image", "Graph as JPG image (low resolution)", "Images", "");
//		t.addOutputItem(graphTxt, "graphTxt", "file", "Graph as textual representation", "Graph as text", "");
//		if(!error.equalsIgnoreCase(""))t.addOutputItem(error , "warning", "warning", "An error occured during execution", "Warning", "");
//		return t.generateXML();
//	}

	private Properties getConfig(){
		Properties config = new Properties();
		try {
			config.load(new FileInputStream(new File(babelomicsHomePath + "/conf/blast2go.properties")));
		} catch (Exception e) {
			e.printStackTrace();
			return null;
		} 
		return config;
	}
	
}
