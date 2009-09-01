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
import org.bioinfo.babelomics.tools.functional.textmining.Marmite;
import org.bioinfo.babelomics.tools.functional.textmining.MarmiteScan;
import org.bioinfo.babelomics.tools.functional.tissues.TMT;
import org.bioinfo.babelomics.tools.preprocessing.IdConverter;
import org.bioinfo.babelomics.tools.preprocessing.Preprocessing;


public class BabelomicsFactory {

	public static BabelomicsTool createTool(String toolName) {
		BabelomicsTool babelomicsTool = null;

		if(toolName.equalsIgnoreCase("differential-expression")) {
			return new DifferentialAnalysis();
		}
		
		if(toolName.equalsIgnoreCase("regression")) {
			return new DifferentialAnalysis();
		}

		if(toolName.equalsIgnoreCase("clustering")) {
			return new Clustering();
		}

		if(toolName.equalsIgnoreCase("detds")) {
			return new TimeDosageSeries();
		}

		if(toolName.equalsIgnoreCase("copa")) {
			return new OutlierCopa();
		}

		if(toolName.equalsIgnoreCase("lrs")) {
			return new OutlierLrs();
		}

		if(toolName.equalsIgnoreCase("predictor")) {
			return new Predictor();
		}

		if(toolName.equalsIgnoreCase("preprocessing")) {
			return new Preprocessing();
		}

		if(toolName.equalsIgnoreCase("idconverter")) {
			return new IdConverter();
		}

		if(toolName.equalsIgnoreCase("fatigo")) {
			return new FatiGO();
		}

		if(toolName.equalsIgnoreCase("fatiscan")) {
			return new FatiScan();
		}

		if(toolName.equalsIgnoreCase("marmite")) {
			return new Marmite();
		}

		if(toolName.equalsIgnoreCase("marmitescan")) {
			return new MarmiteScan();
		}

		if(toolName.equalsIgnoreCase("blast2go")) {
			return new Blast2Go();
		}

		if(toolName.equalsIgnoreCase("tmt")) {
			return new TMT();
		}

		return babelomicsTool;
	}

}
