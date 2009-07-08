package org.bioinfo.babelomics.tools;

import org.bioinfo.babelomics.tools.expression.Clustering;
import org.bioinfo.babelomics.tools.expression.DifferentialAnalysis;
import org.bioinfo.babelomics.tools.expression.OutlierCopa;
import org.bioinfo.babelomics.tools.expression.OutlierLrs;
import org.bioinfo.babelomics.tools.expression.Predictor;
import org.bioinfo.babelomics.tools.expression.TimeDosageSeries;
import org.bioinfo.babelomics.tools.functional.Blast2Go;
import org.bioinfo.babelomics.tools.functional.FatiGO;
import org.bioinfo.babelomics.tools.functional.FatiScan;
import org.bioinfo.babelomics.tools.functional.MarmiteScan;
import org.bioinfo.babelomics.tools.preprocessing.IdConverter;
import org.bioinfo.babelomics.tools.preprocessing.Preprocessing;



public class BabelomicsFactory {

	public static BabelomicsTool createTool(String toolName, String[] args) {
		BabelomicsTool babelomicsTool = null;

		if(toolName.equalsIgnoreCase("differential-expression")) {
			return new DifferentialAnalysis(args);
		}

		if(toolName.equalsIgnoreCase("clustering")) {
			return new Clustering(args);
		}

		if(toolName.equalsIgnoreCase("detds")) {
			return new TimeDosageSeries(args);
		}

		if(toolName.equalsIgnoreCase("copa")) {
			return new OutlierCopa(args);
		}

		if(toolName.equalsIgnoreCase("lrs")) {
			return new OutlierLrs(args);
		}

		if(toolName.equalsIgnoreCase("predictor")) {
			return new Predictor(args);
		}

		if(toolName.equalsIgnoreCase("preprocessing")) {
			return new Preprocessing(args);
		}

		if(toolName.equalsIgnoreCase("idconverter")) {
			return new IdConverter(args);
		}

		if(toolName.equalsIgnoreCase("fatigo")) {
			return new FatiGO(args);
		}

		if(toolName.equalsIgnoreCase("fatiscan")) {
			return new FatiScan(args);
		}

		if(toolName.equalsIgnoreCase("marmitescan")) {
			return new MarmiteScan(args);
		}

		if(toolName.equalsIgnoreCase("blast2go")) {
			return new Blast2Go(args);
		}

		return babelomicsTool;
	}

}
