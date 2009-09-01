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
	public void Test0() {
	}

//	public void Test1() {
//		System.out.println("-----     ------");
//		//String []args = {"-dataset", "/mnt/commons/test/biodata/example/dataset_example1.txt", "-o", "/tmp", "--logarithm-base", "10"};
//		//String []args = {"-dataset", "/mnt/commons/test/biodata/example/dataset_example1.txt", "-o", "/tmp", "--impute-missing", "zero"};
//		//String []args = {"-dataset", "/mnt/commons/test/biodata/example/dataset_example1.txt", "-o", "/tmp", "--merge-replicates", "mean"};
//		String []args = {"-dataset", "/mnt/commons/test/biodata/example/dataset_example1.txt", "-o", "/tmp", "--logarithm-base", "10", "--impute-missing", "zero", "--filter-missing", "90", "--merge-replicates", "mean"};
//		
//		try {
//			Preprocessing prep = new Preprocessing(args);
//			prep.execute();
//		} catch (Exception e) {
//			e.printStackTrace();
//			//System.out.println(e.toString());
//		}
//	}
	@Test
	public void notest(){
		
	}
	
	public void Test3() {
		String dataset = "/home/joaquin/tests/dataset1.txt";
		String outdir = "/home/joaquin/tests/preprocessing";
		String []args = { "--tool", "preprocessing","--log-level", "2", "--dataset", dataset, "-o", outdir,"--merge-replicates", "mean", "--logarithm-base", "2", "--impute-missing", "zero", "--filter-missing", "90"};

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

	public void Test4() {
		String dataset = "/home/joaquin/tests/dataset1.txt";
		String outdir = "/home/joaquin/tests/preprocessing";
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


	public void Test5() {
		String dataset = "/home/joaquin/tests/dataset1.txt";
		String outdir = "/home/joaquin/tests/preprocessing";
		String filename = "/home/joaquin/tests/known.txt";
		String []args = { "--tool", "preprocessing","--log-level", "2", "--dataset", dataset, "-o", outdir, "--gene-list-filter", filename};

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

	public void Test6() {
		String dataset = "/home/joaquin/tests/dataset1.txt";
		String outdir = "/home/joaquin/tests/preprocessing";
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
