package org.bioinfo.babelomics.tools.expression.normalization;

import static org.junit.Assert.fail;

import java.io.File;
import java.util.Arrays;

import org.bioinfo.babelomics.BabelomicsMain;
import org.junit.Test;


public class AgilentExpression2CNormalizationTest {
	public void Test() {
		String outDirName = "/tmp/AgilentExpression2CNormalizationTest";
		new File(outDirName).mkdir();
		String rawDirName = "/mnt/commons/test/biodata/example/GSE11968_RAW/";
		
		String []args = { "--tool", "agilent-expression-two-colors-normalization","--log-level", "2", "--raw-dir", rawDirName, "-o", outDirName};

		System.out.println("command line parameters --> " + Arrays.toString(args));
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.toString());
		}
	}	

	public void Test1() {
		String outDirName = "/tmp/AgilentExpression2CNormalizationTest1";
		new File(outDirName).mkdir();
		String rawDirName = "/mnt/commons/test/biodata/example/GSE11968_RAW/agilent.tar.gz";
		//String rawDirName = "/mnt/commons/test/biodata/example/GSE11968_RAW/agilent.zip";
		
		String []args = { "--tool", "agilent-expression-two-colors-normalization","--log-level", "2", "--compressed-file", rawDirName, "-o", outDirName};

		System.out.println("command line parameters --> " + Arrays.toString(args));
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.toString());
		}
	}	

	@Test
	public void Test2() {
		String outDirName = "/tmp/AgilentExpression2CNormalizationTest2";
		new File(outDirName).mkdir();
		String dataset = "/mnt/commons/test/biodata/example/GSE11968_RAW/agilent.zip";
		
		String []args = { "--tool", "agilent-expression-two-colors-normalization","--log-level", "2", "--compressed-file", dataset, "-o", outDirName, "--bg-correction", "normexp", "--wa-normalization", "loess", "--ba-normalization", "quantile"};

		System.out.println("command line parameters --> " + Arrays.toString(args));
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.toString());
		}
	}	
	
}
