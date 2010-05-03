package org.bioinfo.babelomics.tools.functional;


import java.io.File;
import java.io.IOException;

import org.apache.commons.cli.ParseException;
import org.bioinfo.babelomics.tools.BabelomicsFactory;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class FatiGOTest {

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}
	
	@Test
	public void Test0() {
		String outdir = "/tmp/fatigo";
		new File(outdir).mkdir();
		
		String list1 = "/mnt/commons/test/example.motor";
		String list2 = "/mnt/commons/test/example.apoptosis";
//		String []args = {"--tool", "fatigo" ,"--list1", list1, "--list2", list2, "-o", "/tmp/fatigo"};
		String []args = {"--tool", "fatigo" ,"--list1", list1, "--list2", list2, "--go-bp", "--kegg", "-o", outdir, "--species", "hsa", "--home", System.getenv("BABELOMICS_HOME")};
		try {
			FatiGOTool fatigo = (FatiGOTool)BabelomicsFactory.createTool("fatigo");
			fatigo.parse(args);
			fatigo.run();
		} catch (ParseException e) {
			e.printStackTrace();
		} catch (IOException e) {			
			e.printStackTrace();
		}		
	}

	public void Test1() {
		////
		String []args = {"-list1", "/mnt/commons/test/biodata/example/list3.txt", "-o", "/tmp", "-rest-of-genome", "1", "--go-bp-db", "true", "--go-bp-min-level", "5", "--go-bp-max-level", "12", "--remove-duplicates", "ref", "--home", System.getenv("BABELOMICS_HOME")};
		//String []args = {"-list1", "/mnt/commons/test/biodata/example/list3.txt", "-o", "/tmp", "-list2", "/mnt/commons/test/biodata/example/list4.txt", "--go-bp-db", "true", "--go-bp-min-level", "5", "--go-bp-max-level", "12", "--remove-duplicates", "ref"};
		//String []args = {"-list1", "/mnt/commons/test/biodata/example/list3.txt", "-o", "/tmp", "-list2", "/mnt/commons/test/biodata/example/list4.txt", "--go-bp-db", "true", "--go-bp-min-level", "5", "--go-bp-max-level", "5"};
//		try {
//			FatiGOTool fatigo = new FatiGOTool();
//			fatigo.execute();			
//		} catch (Exception e) {
//			e.printStackTrace();
//			//System.out.println(e.toString());
//		}
	}

}
