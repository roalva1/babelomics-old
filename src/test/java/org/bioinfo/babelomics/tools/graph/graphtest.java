package org.bioinfo.babelomics.tools.graph;


import static org.junit.Assert.fail;

import java.io.File;
import java.util.Arrays;

import org.bioinfo.babelomics.BabelomicsMain;
import org.bioinfo.commons.io.utils.IOUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class graphtest {

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}
	
	@Test
	public void notest(){
		System.out.println("----------------> ");
	}
	
	
	@Test
	public void Test1() {
		//String dataset = "/mnt/commons/test/fatiscanmini2.txt";
		
		//String dataset = "/mnt/commons/test/biodata/example/dataset_example4.txt";
		String dataset = "/mnt/commons/test/biodata/newick1.nw";
		String outdir = "/tmp/histogram";
		new File(outdir).mkdir();
		//String []args = { "--tool", "histogram","--log-level", "2", "--datalist", dataset,"--class", "sex","--boxplot","true","--histogram","true", "-o", outdir};
		String []args = { "--tool", "descriptive-statistics","--log-level", "2", "--datalist", dataset,"-tree","true", "-o", outdir};

		System.out.println("----------------> " + args.toString());
		
		try {
			BabelomicsMain.main(args); 
//			System.out.println("input dataset:\n" + IOUtils.toString(new File(dataset)));
//			System.out.println("output dataset:\n" + IOUtils.toString(new File(outdir + "/histogram.txt")));
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.toString());
			//System.out.println(e.toString());
		}
	}	
	
	
	
	
	
	
	
	
	
	
}
