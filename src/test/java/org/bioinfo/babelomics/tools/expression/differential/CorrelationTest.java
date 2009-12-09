package org.bioinfo.babelomics.tools.expression.differential;

import org.bioinfo.babelomics.BabelomicsMain;
import org.junit.Test;


public class CorrelationTest {

	public void test1() {
	    
		System.out.println("-----   pearson ------");
		//String []args = {"-tool", "differential-expression", "-dataset", "/mnt/commons/test/biodata/example/dataset_example.txt", "-o", "/tmp/pearson", "-test", "pearson", "-class", "O2_conc"};
		String []args = {"-tool", "correlation", "-dataset", "/mnt/commons/test/biodata/example/correlation.txt", "-o", "/tmp/pearson", "-test", "pearson", "-class", "indep"};
		
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
		String []args = {"-tool", "correlation", "-dataset", "/mnt/commons/test/biodata/example/correlation.txt", "-o", "/tmp/spearman", "-test", "spearman", "-class", "indep"};
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}
	
	public void test3() {
		System.out.println("-----  regression  ------");
		//String []args = {"-tool", "differential-expression", "-dataset", "/mnt/commons/test/biodata/example/dataset_example.txt", "-o", "/tmp/regression", "-test", "regression", "-class", "O2_conc"};
		String []args = {"-tool", "correlation", "-dataset", "/mnt/commons/test/biodata/example/correlation.txt", "-o", "/tmp/regression", "-test", "regression", "-class", "indep"};
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}

}
