package org.bioinfo.babelomics.tools.preprocessing;


import java.io.File;
import java.util.Arrays;

import org.bioinfo.babelomics.BabelomicsMain;
import org.bioinfo.commons.io.utils.FileUtils;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class IDConverterTest {

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}

	@Test
	public void Test() {
		String outdir = "/tmp/IDConverterTest";
		new File(outdir).mkdir();

		String []args = { "--tool", "id-converter","--log-level", "2", "--species", "hsa", "--list", "AATK,BIRC7,PARM1,kkkkkk", "--go", "--entrezgene", "--interpro", "-o", outdir};

		System.out.println("ID Converter Test, args : " + Arrays.toString(args));
		try {
			FileUtils.createDirectory(outdir);
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}
	}	

	@Test
	public void Test1() {
		//String dataset = "/mnt/commons/test/biodata/example/genes.txt";
		String dataset = "/mnt/commons/test/tools/tmt/list2_liver.txt";
		String outdir = "/tmp/IDConverterTest1";
		new File(outdir).mkdir();
		
		String []args = { "--tool", "id-converter","--log-level", "2", "--species", "hsa", "--listfile", dataset, "--go", "--kegg", "-o", outdir};
		//String []args = { "--tool", "id-converter","--log-level", "2", "--listfile", dataset, "--go", "--entrezgene", "--interpro", "--havana_gene", "--ensembl_gene", "--codelink", "--kegg", "-o", outdir};

		System.out.println("ID Converter Test1, args : " + Arrays.toString(args));
		try {
			FileUtils.createDirectory(outdir);
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}
	}	

}
