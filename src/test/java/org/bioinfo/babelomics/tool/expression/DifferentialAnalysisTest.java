package org.bioinfo.babelomics.tool.expression;


import java.util.ArrayList;
import java.util.Collection;
import java.util.LinkedHashMap;
import java.util.List;
import java.util.Map;
import java.util.TreeMap;

import org.bioinfo.babelomics.BabelomicsMain;
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

	public void test() {
		System.out.println("-----     ------");
		//String []args = {"-dataset", "/mnt/commons/test/biodata/example/dataset_example1.txt", "-o", "/tmp", "--logarithm-base", "10"};
		//String []args = {"-dataset", "/mnt/commons/test/biodata/example/dataset_example1.txt", "-o", "/tmp", "--impute-missing", "zero"};
		//String []args = {"-dataset", "/mnt/commons/test/biodata/example/dataset_example1.txt", "-o", "/tmp", "--merge-replicates", "mean"};
		String []args = {"-dataset", "/mnt/commons/test/biodata/example/dataset_example2.txt", "-o", "/mnt/commons/test/biodata/example/out", "-test", "t-test", "-class", "sex"};
		
		try {
			DifferentialAnalysis diffexpr = new DifferentialAnalysis();
			diffexpr.execute();
		} catch (Exception e) {
			e.printStackTrace();
			//System.out.println(e.toString());
		}		
	}
	
	
	@Test
	public void test1() {
				
	    
		System.out.println("-----   pearson correlation  ------");
		//String []args = {"-dataset", "/mnt/commons/test/biodata/example/dataset_example1.txt", "-o", "/tmp", "--logarithm-base", "10"};
		//String []args = {"-dataset", "/mnt/commons/test/biodata/example/dataset_example1.txt", "-o", "/tmp", "--impute-missing", "zero"};
		//String []args = {"-dataset", "/mnt/commons/test/biodata/example/dataset_example1.txt", "-o", "/tmp", "--merge-replicates", "mean"};
		String []args = {"-tool", "differential-expression", "-dataset", "/mnt/commons/test/biodata/example/dataset_example.txt", "-o", "/tmp/correlationPearsonFinal", "-test", "pearson", "-class", "only_males"};
		
		try {
			BabelomicsMain.main(args); 
//			DifferentialAnalysis diffexpr = new DifferentialAnalysis(args);
//			diffexpr.execute();
		} catch (Exception e) {
			e.printStackTrace();
			//System.out.println(e.toString());
		}		
	}
	
	public void test2() {
		System.out.println("-----     ------");
		//String []args = {"-dataset", "/mnt/commons/test/biodata/example/dataset_example1.txt", "-o", "/tmp", "--logarithm-base", "10"};
		//String []args = {"-dataset", "/mnt/commons/test/biodata/example/dataset_example1.txt", "-o", "/tmp", "--impute-missing", "zero"};
		//String []args = {"-dataset", "/mnt/commons/test/biodata/example/dataset_example1.txt", "-o", "/tmp", "--merge-replicates", "mean"};
		String []args = {"-tool", "differential-expression", "-dataset", "/mnt/commons/test/biodata/example/dataset_example.txt", "-o", "/tmp/pearson", "-test", "pearson", "-class", "only_males"};
		
		try {
			BabelomicsMain.main(args); 
//			DifferentialAnalysis diffexpr = new DifferentialAnalysis(args);
//			diffexpr.execute();
		} catch (Exception e) {
			e.printStackTrace();
			//System.out.println(e.toString());
		}		
	}
	
	public void test3() {
		System.out.println("-----     ------");
		//String []args = {"-dataset", "/mnt/commons/test/biodata/example/dataset_example1.txt", "-o", "/tmp", "--logarithm-base", "10"};
		//String []args = {"-dataset", "/mnt/commons/test/biodata/example/dataset_example1.txt", "-o", "/tmp", "--impute-missing", "zero"};
		//String []args = {"-dataset", "/mnt/commons/test/biodata/example/dataset_example1.txt", "-o", "/tmp", "--merge-replicates", "mean"};
		String []args = {"-tool", "differential-expression", "-dataset", "/mnt/commons/test/biodata/example/dataset_example.txt", "-o", "/tmp", "-test", "regression", "-class", "O2_conc"};
		
		try {
			BabelomicsMain.main(args); 
//			DifferentialAnalysis diffexpr = new DifferentialAnalysis(args);
//			diffexpr.execute();
		} catch (Exception e) {
			e.printStackTrace();
			//System.out.println(e.toString());
		}		
	}

}
