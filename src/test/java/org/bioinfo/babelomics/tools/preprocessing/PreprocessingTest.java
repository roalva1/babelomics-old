package org.bioinfo.babelomics.tools.preprocessing;


import static org.junit.Assert.fail;

import java.io.File;
import java.util.Arrays;

import org.bioinfo.babelomics.BabelomicsMain;
import org.bioinfo.commons.io.utils.IOUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class PreprocessingTest {

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}
	
	@Test
	public void notest(){
		System.out.println("----------------> ");
	}
	
<<<<<<< HEAD:src/test/java/org/bioinfo/babelomics/tools/preprocessing/PreprocessingTest.java
	public void Test3() {
		String dataset = "/mnt/commons/test/biodata/example/dataset1.txt";
		String outdir = "/tmp/preprocessing";
		String []args = { "--tool", "preprocessing","--log-level", "2", "--dataset", dataset, "-o", outdir,"--merge-replicates", "mean", "--logarithm-base", "2", "--impute-missing", "zero", "--filter-missing", "90", "--report", "pdf"};
=======
	@Test
	public void Test() {
		String dataset = "/mnt/commons/test/biodata/example/dataset_example.txt";
		String outdir = "/tmp/PreprocessingTest";
		new File(outdir).mkdir();
		String []args = { "--tool", "preprocessing","--log-level", "2", "--dataset", dataset, "-o", outdir,"--merge-replicates", "mean", "--logarithm-base", "2", "--impute-missing", "zero", "--filter-missing", "90"};
>>>>>>> a78add3d31a7884921cb78eab15152ae82b8b5b8:src/test/java/org/bioinfo/babelomics/tools/preprocessing/PreprocessingTest.java

		System.out.println("----------------> " + Arrays.toString(args));
		
		try {
			BabelomicsMain.main(args); 
			System.out.println("input dataset:\n" + IOUtils.toString(new File(dataset)));
			System.out.println("output dataset:\n" + IOUtils.toString(new File(outdir + "/preprocessed.txt")));
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.toString());
			//System.out.println(e.toString());
		}
	}	

	@Test
	public void Test1() {
		String dataset = "/mnt/commons/test/biodata/example/dataset_example.txt";
		String outdir = "/tmp/PreprocessingTest1";
		new File(outdir).mkdir();
		String []args = { "--tool", "preprocessing","--log-level", "2", "--dataset", dataset, "-o", outdir, "--impute-missing", "zero"};

		System.out.println("----------------> " + Arrays.toString(args));
		
		try {
			BabelomicsMain.main(args);
			System.out.println("input dataset:\n" + IOUtils.toString(new File(dataset)));
			System.out.println("output dataset:\n" + IOUtils.toString(new File(outdir + "/preprocessed.txt")));
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.toString());
			//System.out.println(e.toString());
		}
	}	


	@Test
	public void Test2() {
		String dataset = "/mnt/commons/test/biodata/example/dataset_example.txt";
		String outdir = "/tmp/PreprocessingTest2";
		new File(outdir).mkdir();
		String filename = "/mnt/commons/test/biodata/example/known_genes.txt";
		String []args = { "--tool", "preprocessing","--log-level", "2", "--dataset", dataset, "-o", outdir, "--gene-file-filter", filename};

		System.out.println("----------------> " + Arrays.toString(args));
		
		try {
			BabelomicsMain.main(args);
			System.out.println("input dataset:\n" + IOUtils.toString(new File(dataset)));
			System.out.println("output dataset:\n" + IOUtils.toString(new File(outdir + "/preprocessed.txt")));
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.toString());
			//System.out.println(e.toString());
		}
	}	

	@Test
	public void Test3() {
		String dataset = "/mnt/commons/test/biodata/example/dataset_example.txt";
		String outdir = "/tmp/PreprocessingTest3";
		new File(outdir).mkdir();
		String []args = { "--tool", "preprocessing","--log-level", "2", "--dataset", dataset, "-o", outdir, "--impute-missing", "knn", "--kvalue", "3"};

		System.out.println("----------------> " + Arrays.toString(args));
		
		try {
			BabelomicsMain.main(args); 
			System.out.println("input dataset:\n" + IOUtils.toString(new File(dataset)));
			System.out.println("output dataset:\n" + IOUtils.toString(new File(outdir + "/preprocessed.txt")));
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.toString());
			//System.out.println(e.toString());
		}
	}	
}
