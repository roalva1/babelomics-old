package org.bioinfo.babelomics.tools.expression.normalization;


import static org.junit.Assert.fail;

import java.io.File;
import java.util.Arrays;

import org.bioinfo.babelomics.BabelomicsMain;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class GenePixExpression2CNormalizationTest {

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void Test() {
		String outDirName = "/tmp/GenePixExpression2CNormalizationTest";
		new File(outDirName).mkdir();
		String dataset = "/mnt/commons/test/biodata/example/twocolor_genepix_expression.zip";
		
		String []args = { "--tool", "genepix-expression-two-colors-normalization","--log-level", "2", "--compressed-file", dataset, "-o", outDirName, "--bg-correction", "normexp", "--ba-normalization", "scale", "--wa-normalization", "loess"};

		System.out.println("command line parameters --> " + Arrays.toString(args));
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.toString());
		}
	}		
}
