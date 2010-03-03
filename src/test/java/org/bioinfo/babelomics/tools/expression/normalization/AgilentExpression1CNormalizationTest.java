package org.bioinfo.babelomics.tools.expression.normalization;


import static org.junit.Assert.*;

import java.io.File;
import java.util.Arrays;

import org.bioinfo.babelomics.BabelomicsMain;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;
import org.junit.runner.notification.Failure;

public class AgilentExpression1CNormalizationTest {

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}
	
	@Test
	public void Test() {
		String outDirName = "/tmp/AgilentExpression1CNormalizationTest";
		new File(outDirName).mkdir();
		String dataset = "/mnt/commons/test/biodata/example/agilent_one_color_expression.zip";
		
		String []args = { "--tool", "agilent-expression-one-color-normalization","--log-level", "2", "--compressed-file", dataset, "-o", outDirName, "--bg-correction", "none", "--ba-normalization", "quantiles"};

		System.out.println("command line parameters --> " + Arrays.toString(args));
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.toString());
		}
	}	

}
