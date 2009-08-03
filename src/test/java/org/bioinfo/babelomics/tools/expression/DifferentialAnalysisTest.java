package org.bioinfo.babelomics.tools.expression;


import org.bioinfo.babelomics.BabelomicsMain;
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
	public void test() {
		System.out.println("-----     ------");
		String []args = {"-tool", "differential-expression","-dataset", "/mnt/commons/test/biodata/example/twoclasses.txt", "-o", "/tmp/ttest", "-test", "t-test", "-class", "class"};
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}
	
	public void test1() {
					    
		System.out.println("-----   pearson ------");
		//String []args = {"-tool", "differential-expression", "-dataset", "/mnt/commons/test/biodata/example/dataset_example.txt", "-o", "/tmp/pearson", "-test", "pearson", "-class", "O2_conc"};
		String []args = {"-tool", "differential-expression", "-dataset", "/mnt/commons/test/biodata/example/correlation.txt", "-o", "/tmp/pearson", "-test", "pearson", "-class", "indep"};
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}
	
	@Test
	public void test2() {
		System.out.println("-----  spearman ------");
		//String []args = {"-tool", "differential-expression", "-dataset", "/mnt/commons/test/biodata/example/dataset_example.txt", "-o", "/tmp/spearman", "-test", "spearman", "-class", "O2_conc"};
		String []args = {"-tool", "differential-expression", "-dataset", "/mnt/commons/test/biodata/example/correlation.txt", "-o", "/tmp/spearman", "-test", "spearman", "-class", "indep"};
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}
	
	public void test3() {
		System.out.println("-----  regression  ------");
		//String []args = {"-tool", "differential-expression", "-dataset", "/mnt/commons/test/biodata/example/dataset_example.txt", "-o", "/tmp/regression", "-test", "regression", "-class", "O2_conc"};
		String []args = {"-tool", "differential-expression", "-dataset", "/mnt/commons/test/biodata/example/correlation.txt", "-o", "/tmp/regression", "-test", "regression", "-class", "indep"};
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}
	
	
	public void test4() {
		System.out.println("-----     ------");
		String []args = {"-tool", "differential-expression", "-dataset", "/mnt/commons/test/biodata/example/survival.txt", "-o", "/tmp/cox", "-test", "cox", "-time-class", "time", "-censored-class", "censored"};
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}

	public void test5() {
		System.out.println("-----     ------");
		String []args = {"-tool", "differential-expression","-dataset", "/mnt/commons/test/biodata/example/multiclasses.txt", "-o", "/tmp/anova", "-test", "anova", "-class", "class"};
		
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
