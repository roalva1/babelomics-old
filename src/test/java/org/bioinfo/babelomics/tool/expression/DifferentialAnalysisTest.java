package org.bioinfo.babelomics.tool.expression;


import org.bioinfo.babelomics.tools.expression.DifferentialAnalysis;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;


public class DifferentialAnalysisTest {

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void test1() {
		System.out.println("-----     ------");
		//String []args = {"-dataset", "/mnt/commons/test/biodata/example/dataset_example1.txt", "-o", "/tmp", "--logarithm-base", "10"};
		//String []args = {"-dataset", "/mnt/commons/test/biodata/example/dataset_example1.txt", "-o", "/tmp", "--impute-missing", "zero"};
		//String []args = {"-dataset", "/mnt/commons/test/biodata/example/dataset_example1.txt", "-o", "/tmp", "--merge-replicates", "mean"};
		String []args = {"-dataset", "/mnt/commons/test/biodata/example/dataset_example2.txt", "-o", "/mnt/commons/test/biodata/example/out", "-test", "t-test", "-class", "sex"};
		
		try {
			DifferentialAnalysis diffexpr = new DifferentialAnalysis(args);
			diffexpr.execute();
		} catch (Exception e) {
			e.printStackTrace();
			//System.out.println(e.toString());
		}		
	}
}
