package org.bioinfo.babelomics.tools.expression;


import java.io.File;
import java.util.Arrays;

import org.bioinfo.babelomics.BabelomicsMain;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class RawExpressionViewerTest {

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}

	public void Test() {
		String dataset = "/mnt/commons/test/biodata/example/onecolor_agilent_expression.zip";
		String outdir = "/tmp/RawExpressionViewerTest";
		new File(outdir).mkdir();

		int channels = 1;
		String technology = "agilent";
		
		String []args = { "--tool", "raw-expression-viewer","--log-level", "2", "--compressed-file", dataset, "-o", outdir, "--technology", technology, "--channels", ""+channels};

		System.out.println("----- " + technology + " " + channels + " color ----------------> " + Arrays.toString(args));
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}
	}	

	public void Test1() {
		String dataset = "/mnt/commons/test/biodata/example/twocolor_agilent_expression.zip";
		String outdir = "/tmp/RawExpressionViewerTest1";
		new File(outdir).mkdir();

		int channels = 2;
		String technology = "agilent";
		
		String []args = { "--tool", "raw-expression-viewer","--log-level", "2", "--compressed-file", dataset, "-o", outdir, "--technology", technology, "--channels", ""+channels};

		System.out.println("----- " + technology + " " + channels + " color ----------------> " + Arrays.toString(args));
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}
	}	

	public void Test2() {
		String dataset = "/mnt/commons/test/biodata/example/onecolor_genepix_expression.zip";
		String outdir = "/tmp/RawExpressionViewerTest2";
		new File(outdir).mkdir();

		int channels = 1;
		String technology = "genepix";
		
		String []args = { "--tool", "raw-expression-viewer","--log-level", "2", "--compressed-file", dataset, "-o", outdir, "--technology", technology, "--channels", ""+channels};

		System.out.println("----- " + technology + " " + channels + " color ----------------> " + Arrays.toString(args));
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}
	}	

	@Test
	public void Test3() {
		String dataset = "/mnt/commons/test/biodata/example/twocolor_genepix_expression.zip";
		String outdir = "/tmp/RawExpressionViewerTest3";
		new File(outdir).mkdir();

		int channels = 2;
		String technology = "genepix";
		
		String []args = { "--tool", "raw-expression-viewer","--log-level", "2", "--compressed-file", dataset, "-o", outdir, "--technology", technology, "--channels", ""+channels};

		System.out.println("----- " + technology + " " + channels + " color ----------------> " + Arrays.toString(args));
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}
	}	
}
