package org.bioinfo.babelomics.tools.genomic.copynumber;

import org.bioinfo.babelomics.BabelomicsMain;
import org.junit.Test;


public class CopyNumberAnalysisTest {
	
	@Test
	public void test() {
		System.out.println("----- copynumberanalysis dnacopy ------");
		String []args = {"-tool", "copy-number", "-normalized-file", "/mnt/commons/test/biodata/example/cgh/agilent/segmentation/data1.txt", "-o", "/tmp/copynumber-dnacopy", "-segmentation-method", "dnacopy", "-cgh-mcr", "-gap-allowed", "400"};
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}

	public void test1() {
		System.out.println("----- copynumberanalysis glad ------");
		String []args = {"-tool", "copy-number", "-normalized-file", "/mnt/commons/test/biodata/example/cgh/agilent/segmentation/data1.txt", "-o", "/tmp/copynumber-glad", "-segmentation-method", "glad"};
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}
}
