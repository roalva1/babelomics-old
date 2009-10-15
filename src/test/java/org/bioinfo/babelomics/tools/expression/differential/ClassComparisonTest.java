package org.bioinfo.babelomics.tools.expression.differential;

import org.bioinfo.babelomics.BabelomicsMain;
import org.junit.Test;


public class ClassComparisonTest {
	
	public void test1() {
	    
		System.out.println("----- two classes - limma ------");
		//String []args = {"-tool", "differential-expression", "-dataset", "/mnt/commons/test/biodata/example/dataset_example.txt", "-o", "/tmp/pearson", "-test", "pearson", "-class", "O2_conc"};
		String []args = {"-tool", "class-comparison", "-dataset", "/mnt/commons/test/biodata/example/twoclasses100.txt", "-o", "/tmp/limmatwo", "-test", "limma", "-class-name", "class", "-class-values", "basal,luminal"};
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}

	public void test2() {
	    
		System.out.println("----- one class - limma ------");
		String []args = {"-tool", "class-comparison", "-dataset", "/mnt/commons/test/biodata/example/twoclasses100.txt", "-o", "/tmp/limmaone", "-test", "limma", "-class-name", "class", "-class-values", "luminal"};
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}

	public void test3() {
	    
		System.out.println("----- fold change ------");
		//String []args = {"-tool", "differential-expression", "-dataset", "/mnt/commons/test/biodata/example/dataset_example.txt", "-o", "/tmp/pearson", "-test", "pearson", "-class", "O2_conc"};
		String []args = {"-tool", "class-comparison", "-dataset", "/mnt/commons/test/biodata/example/twoclasses100.txt", "-o", "/tmp/foldchange", "-test", "fold-change", "-class-name", "class", "-class-values", "basal,luminal"};
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}
	
	@Test
	public void test4() {
	    
		System.out.println("----- multi class - limma ------");
		//String []args = {"-tool", "differential-expression", "-dataset", "/mnt/commons/test/biodata/example/dataset_example.txt", "-o", "/tmp/pearson", "-test", "pearson", "-class", "O2_conc"};
		String []args = {"-tool", "class-comparison", "-dataset", "/mnt/commons/test/biodata/example/multiclasses.txt", "-o", "/tmp/limmamulti", "-test", "limma", "-class-name", "class", "-class-values", "A,B,C"};
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}
}
