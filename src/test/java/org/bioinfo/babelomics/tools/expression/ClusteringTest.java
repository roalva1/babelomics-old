package org.bioinfo.babelomics.tools.expression;


import static org.junit.Assert.fail;

import java.io.File;
import java.util.Arrays;

import org.bioinfo.babelomics.BabelomicsMain;
import org.bioinfo.commons.io.utils.FileUtils;
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

	@Test
	public void Test1() {
//		String dataset = "/mnt/commons/test/biodata/example/preprocessed.txt";
		String dataset = "/mnt/commons/test/biodata/example/cyano.txt";
		String outdir = "/tmp/upgma";
		String []args = { "--tool", "clustering","--log-level", "2", "--dataset", dataset, "-o", outdir, "--method", "upgma", "--distance", "euclidean"};

		System.out.println("UPGMA ----------------> " + Arrays.toString(args));
		try {
			FileUtils.createDirectory(outdir);
			BabelomicsMain.main(args); 
			System.out.println("input dataset:\n" + IOUtils.toString(new File(dataset)));
			System.out.println("cluster of genes:\n" + IOUtils.toString(new File(outdir + "/genes.nw")));
			System.out.println("cluster of samples:\n" + IOUtils.toString(new File(outdir + "/samples.nw")));
		} catch (Exception e) {
			e.printStackTrace();
			//fail(e.toString());
			//System.out.println(e.toString());
		}
	}	

	public void Test2() {
//		String dataset = "/mnt/commons/test/biodata/example/preprocessed.txt";
		String dataset = "/mnt/commons/test/biodata/example/cyano.txt";
		String outdir = "/tmp/sota";
		String []args = { "--tool", "clustering","--log-level", "2", "--dataset", dataset, "-o", outdir, "--method", "sota", "--distance", "euclidean"};

		System.out.println("SOTA ----------------> " + Arrays.toString(args));
		
		try {
			FileUtils.createDirectory(outdir);
			BabelomicsMain.main(args); 
			System.out.println("input dataset:\n" + IOUtils.toString(new File(dataset)));
//			System.out.println("cluster of genes:\n" + IOUtils.toString(new File(outdir + "/genes.nw")));
//			System.out.println("cluster of samples:\n" + IOUtils.toString(new File(outdir + "/samples.nw")));
		} catch (Exception e) {
			e.printStackTrace();
			//fail(e.toString());
			//System.out.println(e.toString());
		}
	}	

	public void Test3() {
		String dataset = "/mnt/commons/test/biodata/example/preprocessed.txt";
		String outdir = "/tmp/som";
		String []args = { "--tool", "clustering","--log-level", "2", "--dataset", dataset, "-o", outdir, "--method", "som", "--distance", "euclidean"};

		System.out.println("SOM ----------------> " + Arrays.toString(args));
		
		try {
			FileUtils.createDirectory(outdir);
			BabelomicsMain.main(args); 
			System.out.println("input dataset:\n" + IOUtils.toString(new File(dataset)));
			System.out.println("cluster of genes:\n" + IOUtils.toString(new File(outdir + "/genes.nw")));
			System.out.println("cluster of samples:\n" + IOUtils.toString(new File(outdir + "/samples.nw")));
		} catch (Exception e) {
			e.printStackTrace();
			//fail(e.toString());
			//System.out.println(e.toString());
		}
	}	

	public void Test4() {
		String dataset = "/mnt/commons/test/biodata/example/preprocessed.txt";
//		String dataset = "/mnt/commons/test/biodata/example/cyano.txt";
		String outdir = "/tmp/kmeans";
		String []args = { "--tool", "clustering","--log-level", "2", "--dataset", dataset, "-o", outdir, "--method", "kmeans", "--distance", "euclidean", "--kvalue", "4"};

		System.out.println("KMEANS ----------------> " + Arrays.toString(args));
		
		try {
			FileUtils.createDirectory(outdir);
			BabelomicsMain.main(args); 
//			System.out.println("input dataset:\n" + IOUtils.toString(new File(dataset)));
			System.out.println("cluster of genes:\n" + IOUtils.toString(new File(outdir + "/genes.nw")));
			System.out.println("cluster of samples:\n" + IOUtils.toString(new File(outdir + "/samples.nw")));
		} catch (Exception e) {
			e.printStackTrace();
			//fail(e.toString());
			//System.out.println(e.toString());
		}
	}	
}
