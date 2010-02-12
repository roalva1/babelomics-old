package org.bioinfo.babelomics.tools.preprocessing;


import java.io.File;
import java.util.Arrays;

import org.bioinfo.babelomics.BabelomicsMain;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class CreateAnnotationTest {

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void Test1() {
		String outdir = "/tmp/CreateAnnotationTest1";
		new File(outdir).mkdir();

		String []args = { "--tool", "create-annotation","--log-level", "2", "--species", "hsa", "--all-genome", "true", "--go-cc", "--go-bp", "--go-mf", "-o", outdir};

		System.out.println("CreateAnnotationTest1, args : " + Arrays.toString(args));
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}
	}	

}
