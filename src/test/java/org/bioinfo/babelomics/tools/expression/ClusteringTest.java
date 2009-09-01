package org.bioinfo.babelomics.tools.expression;


import static org.junit.Assert.fail;

import java.io.File;
import java.util.Arrays;

import org.bioinfo.babelomics.BabelomicsMain;
import org.bioinfo.commons.io.utils.IOUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class ClusteringTest {

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void notest() {		
	}

	public void Test6() {
		String dataset = "/home/joaquin/tests/preprocessed.txt";
		String outdir = "/home/joaquin/tests/clustering";
		String []args = { "--tool", "clustering","--log-level", "2", "--dataset", dataset, "-o", outdir, "--method", "sota", "--distance", "euclidean"};

		System.out.println("----------------> " + Arrays.toString(args));
		
		try {
			BabelomicsMain.main(args); 
			System.out.println("input dataset:\n" + IOUtils.toString(new File(dataset)));
			System.out.println("cluster of genes:\n" + IOUtils.toString(new File(outdir + "/genes.nw")));
			System.out.println("cluster of samples:\n" + IOUtils.toString(new File(outdir + "/samples.nw")));
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.toString());
			//System.out.println(e.toString());
		}
	}	

}
