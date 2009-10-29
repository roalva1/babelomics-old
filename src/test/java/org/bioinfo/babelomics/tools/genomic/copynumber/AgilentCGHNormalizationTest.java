package org.bioinfo.babelomics.tools.genomic.copynumber;

import java.io.File;

import org.bioinfo.babelomics.BabelomicsMain;
import org.junit.Test;


public class AgilentCGHNormalizationTest {
	
	@Test
	public void test() {
		String outdirname = "/tmp/agilent-cgh-normalization";
		new File(outdirname).mkdir();
		
		System.out.println("----- copynumberanalysis dnacopy from directory ------");
		String []args = {"-tool", "copy-number-normalization", "--raw-dir", "/mnt/commons/test/biodata/example/cgh/agilent/normalization/dataset1", "-o", outdirname, "--bg-correction", "normexp", "--wa-normalization", "loess", "--ba-normalization", "quantile", "--design", "1"};
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}

}
