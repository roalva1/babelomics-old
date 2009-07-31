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

	public void test() {
		System.out.println("-----     ------");
		//String []args = {"-dataset", "/mnt/commons/test/biodata/example/dataset_example1.txt", "-o", "/tmp", "--logarithm-base", "10"};
		//String []args = {"-dataset", "/mnt/commons/test/biodata/example/dataset_example1.txt", "-o", "/tmp", "--impute-missing", "zero"};
		//String []args = {"-dataset", "/mnt/commons/test/biodata/example/dataset_example1.txt", "-o", "/tmp", "--merge-replicates", "mean"};
		String []args = {"-tool", "differential-expression","-dataset", "/mnt/commons/test/biodata/example/dataset_example.txt", "-o", "/tmp/ttest", "-test", "t-test", "-class", "sex"};
		
		try {
			BabelomicsMain.main(args); 
//			DifferentialAnalysis diffexpr = new DifferentialAnalysis(args);
//			diffexpr.execute();
		} catch (Exception e) {
			e.printStackTrace();
			//System.out.println(e.toString());
		}		
	}
	
	@Test
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
		//String []args = {"-dataset", "/mnt/commons/test/biodata/example/dataset_example1.txt", "-o", "/tmp", "--logarithm-base", "10"};
		//String []args = {"-dataset", "/mnt/commons/test/biodata/example/dataset_example1.txt", "-o", "/tmp", "--impute-missing", "zero"};
		//String []args = {"-dataset", "/mnt/commons/test/biodata/example/dataset_example1.txt", "-o", "/tmp", "--merge-replicates", "mean"};
		String []args = {"-tool", "differential-expression", "-dataset", "/mnt/commons/test/biodata/example/dataset_example.txt", "-o", "/tmp/coxAnalisys", "-test", "cox", "-class", "only_males"};
		
		try {
			BabelomicsMain.main(args); 
//			DifferentialAnalysis diffexpr = new DifferentialAnalysis(args);
//			diffexpr.execute();
		} catch (Exception e) {
			e.printStackTrace();
			//System.out.println(e.toString());
		}		
	}
	
	public void test5() {
		System.out.println("-----     ------");
		//String []args = {"-dataset", "/mnt/commons/test/biodata/example/dataset_example1.txt", "-o", "/tmp", "--logarithm-base", "10"};
		//String []args = {"-dataset", "/mnt/commons/test/biodata/example/dataset_example1.txt", "-o", "/tmp", "--impute-missing", "zero"};
		//String []args = {"-dataset", "/mnt/commons/test/biodata/example/dataset_example1.txt", "-o", "/tmp", "--merge-replicates", "mean"};
		String []args = {"-tool", "differential-expression","-dataset", "/mnt/commons/test/biodata/example/dataset_example.txt", "-o", "/tmp/anova", "-test", "anova", "-class", "only_males"};
		
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
