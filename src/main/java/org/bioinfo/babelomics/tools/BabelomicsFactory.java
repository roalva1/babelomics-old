package org.bioinfo.babelomics.tools;

import org.bioinfo.babelomics.tools.expression.Clustering;
import org.bioinfo.babelomics.tools.expression.DifferentialAnalysis;
import org.bioinfo.babelomics.tools.expression.OutlierLrs;
import org.bioinfo.babelomics.tools.expression.Predictor;
import org.bioinfo.babelomics.tools.expression.differential.Correlation;
import org.bioinfo.babelomics.tools.expression.differential.MaSigPro;
import org.bioinfo.babelomics.tools.expression.differential.Survival;
import org.bioinfo.babelomics.tools.expression.normalization.AffyNormalization;
import org.bioinfo.babelomics.tools.functional.Blast2Go;
import org.bioinfo.babelomics.tools.functional.FatiGOTool;
import org.bioinfo.babelomics.tools.functional.FatiScanTool;
import org.bioinfo.babelomics.tools.functional.textmining.Marmite;
import org.bioinfo.babelomics.tools.functional.textmining.MarmiteScan;
import org.bioinfo.babelomics.tools.functional.tissues.AffyTmt;
import org.bioinfo.babelomics.tools.genomic.genotype.SnpAffymetrixNormalization;
import org.bioinfo.babelomics.tools.interactome.Snow;
import org.bioinfo.babelomics.tools.preprocessing.IdConverter;
import org.bioinfo.babelomics.tools.preprocessing.Preprocessing;


public class BabelomicsFactory {

	public static BabelomicsTool createTool(String toolName) {
		BabelomicsTool babelomicsTool = null;

		if(toolName.equalsIgnoreCase("affymetrix-normalization")) {
			return new AffyNormalization();
		}
		
		if(toolName.equalsIgnoreCase("affy-snp-calls")) {
			return new SnpAffymetrixNormalization();
		}

		if(toolName.equalsIgnoreCase("preprocessing")) {
			return new Preprocessing();
		}
		
		if(toolName.equalsIgnoreCase("id-converter")) {
			return new IdConverter();
		}
		
		if(toolName.equalsIgnoreCase("differential-expression")) {
			return new DifferentialAnalysis();
		}

		if(toolName.equalsIgnoreCase("correlation")) {
			return new Correlation();
		}

		if(toolName.equalsIgnoreCase("survival")) {
			return new Survival();
		}

		if(toolName.equalsIgnoreCase("masigpro")) {
			return new MaSigPro();
		}
		
		if(toolName.equalsIgnoreCase("predictor")) {
			return new Predictor();
		}

		if(toolName.equalsIgnoreCase("regression")) {
			return new DifferentialAnalysis();
		}
		
		if(toolName.equalsIgnoreCase("clustering")) {
			return new Clustering();
		}

		if(toolName.equalsIgnoreCase("fatigo")) {
			return new FatiGOTool();
		}

		if(toolName.equalsIgnoreCase("fatiscan")) {
			return new FatiScanTool();
		}

		if(toolName.equalsIgnoreCase("marmite")) {
			return new Marmite();
		}

		if(toolName.equalsIgnoreCase("marmite-scan")) {
			return new MarmiteScan();
		}
		
		if(toolName.equalsIgnoreCase("tmt-affy")) {
			return new AffyTmt();
		}
		
		if(toolName.equalsIgnoreCase("association")) {
			return new OutlierLrs();
		}
		
		if(toolName.equalsIgnoreCase("tdt")) {
			return new OutlierLrs();
		}
		
		if(toolName.equalsIgnoreCase("linkage")) {
			return new OutlierLrs();
		}
		
		if(toolName.equalsIgnoreCase("gesbap")) {
			return new OutlierLrs();
		}
		
		if(toolName.equalsIgnoreCase("snow")) {
			return new Snow();
		}

		if(toolName.equalsIgnoreCase("blast2go")) {
			return new Blast2Go();
		}

		return babelomicsTool;
	}

}
