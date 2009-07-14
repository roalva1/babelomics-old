package org.bioinfo.babelomics.tool.preprocessing;

import static org.junit.Assert.fail;

import org.bioinfo.babelomics.BabelomicsMain;
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
	
//	@Test
//	public void Test0() {
//	}
//
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
	
	public void Test3() {
		System.out.println("-----     ------");
		//String []args = {"--tool", "preprocessing", "-dataset", "/mnt/commons/test/biodata/example/dataset_example1.txt", "-o", "/tmp", "--logarithm-base", "10"};
		//String []args = {"-dataset", "/mnt/commons/test/biodata/example/dataset_example1.txt", "-o", "/tmp", "--impute-missing", "zero"};
		//String []args = {"-dataset", "/mnt/commons/test/biodata/example/dataset_example1.txt", "-o", "/tmp", "--merge-replicates", "mean"};
		String []args = { "--tool", "preprocessing","--log-level", "2", "--dataset", "/mnt/commons/test/biodata/example/dataset_example1.txt", "-o", "/tmp","--merge-replicates", "mean", "--logarithm-base", "2", "--impute-missing", "zero", "--filter-missing", "90", "--jobname", "sample name", "--submit_button Run", "--session-id", "QIVn13S0VRgTvluw1dwbOxcYXfVusjoyhINlcuouYIZc7jPT5IWc8nRDVPfAYZrP"};
		
		try {
			BabelomicsMain.main(args); 
			System.out.println("holaaa");
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.toString());
			//System.out.println(e.toString());
		}
	}


}
