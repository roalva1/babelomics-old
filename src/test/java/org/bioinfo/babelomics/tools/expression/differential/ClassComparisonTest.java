package org.bioinfo.babelomics.tools.expression.differential;

import java.io.File;

import org.bioinfo.babelomics.BabelomicsMain;
import org.junit.Test;


public class ClassComparisonTest {
	
	@Test
	public void notest() {
	}

	@Test
	public void Test() {
		String dataset = "/mnt/commons/test/biodata/example/twoclasses100.txt";
		String outdir = "/tmp/ClassComparisonTest";
		new File(outdir).mkdir();

		System.out.println("----- one class - limma ------");
		String []args = {"--tool", "class-comparison", "--dataset", dataset, "-o", outdir, "--test", "limma", "--class-name", "class", "--class-values", "luminal", "--correction", "fdr"};
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}

	@Test
	public void Test1() {
	    
		String dataset = "/mnt/commons/test/biodata/example/twoclasses100.txt";
		String outdir = "/tmp/ClassComparisonTest1";
		new File(outdir).mkdir();

		System.out.println("----- two classes - limma ------");
		String []args = {"--tool", "class-comparison", "--dataset", dataset, "-o", outdir, "--test", "limma", "--class-name", "class", "--class-values", "basal,luminal", "--correction", "bonferroni"};
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}

	@Test
	public void Test2() {
	    
		String dataset = "/mnt/commons/test/biodata/example/twoclasses100.txt";
		String outdir = "/tmp/ClassComparisonTest2";
		new File(outdir).mkdir();

		System.out.println("----- two classes - ttest ------");
		String []args = {"--tool", "class-comparison", "--dataset", dataset, "-o", outdir, "--test", "t", "--class-name", "class", "--class-values", "basal,luminal", "--correction", "hochberg"};
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}

	@Test
	public void Test3() {
	    
		String dataset = "/mnt/commons/test/biodata/example/twoclasses100.txt";
		String outdir = "/tmp/ClassComparisonTest3";
		new File(outdir).mkdir();

		System.out.println("----- two classes - fold change ------");
		String []args = {"--tool", "class-comparison", "--dataset", dataset, "-o", outdir, "--test", "fold-change", "--class-name", "class", "--class-values", "basal,luminal", "--correction", "holm"};
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}
		
	@Test
	public void Test4() {
	    
		String dataset = "/mnt/commons/test/biodata/example/multiclasses.txt";
		String outdir = "/tmp/ClassComparisonTest4";
		new File(outdir).mkdir();

		System.out.println("----- multi classes - limma ------");
		String []args = {"--tool", "class-comparison", "--dataset", dataset, "-o", outdir, "--test", "limma", "--class-name", "indep", "--class-values", "1,3,5,7,9", "--correction", "by"};
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}

	@Test
	public void Test5() {
	    
		String dataset = "/mnt/commons/test/biodata/example/multiclasses.txt";
		String outdir = "/tmp/ClassComparisonTest5";
		new File(outdir).mkdir();

		System.out.println("----- multi classes - anova ------");
		String []args = {"--tool", "class-comparison", "--dataset", dataset, "-o", outdir, "--test", "anova", "--class-name", "indep", "--class-values", "1,3,5,7,9", "--correction", "and now what?"};
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}
}
