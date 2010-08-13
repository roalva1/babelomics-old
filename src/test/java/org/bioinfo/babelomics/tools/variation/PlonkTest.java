package org.bioinfo.babelomics.tools.variation;

import org.bioinfo.babelomics.BabelomicsMain;
import org.junit.Test;

public class PlonkTest {
	@Test
	public void Test() {
		
//		String []args = {"--tool", "ped", "-o", "/home/ralonso/ped/", "--input-ped-file", "/home/ralonso/ped/job.ped", "--input-map-file", "/home/ralonso/ped/job.map", "--format", "tped", "--output-file","/home/ralonso/ped/outJob", "--home", System.getenv("BABELOMICS_HOME")};
//		String []args = {"--tool", "ped", "-o", "/home/ralonso/ped/", "--input-ped-file", "/home/ralonso/ped/test.ped", "--input-map-file", "/home/ralonso/ped/test.map", 
//				"--input-snp-file", "/home/ralonso/ped/filterFiles/snpList", "--input-base-file","/home/ralonso/ped/filterFiles/bpList","--home", System.getenv("BABELOMICS_HOME")};
//		String []args = {"--tool", "ped", "-o", "/home/ralonso/ped/", "--input-ped-file", "/home/ralonso/ped/test.ped", "--input-map-file", "/home/ralonso/ped/test.map", 
//				"--input-ch-file", "/home/ralonso/ped/filterFiles/chList", "--input-gd-file","/home/ralonso/ped/filterFiles/gdList",
//				"--ch","17","--from-base","41419601","--to-base","41419605","--home", System.getenv("BABELOMICS_HOME")};
//		String []args = {"--tool", "plonk", "-o", "/home/ralonso/ped/", "--input-ped-file", "/home/ralonso/ped/test01.ped", "--input-map-file", "/home/ralonso/ped/test.map", 
//				"--family", "2", "--home", System.getenv("BABELOMICS_HOME")};
//		String []args = {"--tool", "plonk", "-o", "/home/ralonso/ped/", "--input-ped-file", "/home/ralonso/ped/test.ped", "--input-map-file", "/home/ralonso/ped/test.map", 
//				"--snp-format", "0,1", "--output-ped-file","/home/ralonso/ped/test01.ped","--home", System.getenv("BABELOMICS_HOME")};
		String []args = {"--tool", "plonk", "-o", "/home/ralonso/ped/", "--input-ped-file", "/home/ralonso/ped/test.ped", "--input-map-file", "/home/ralonso/ped/test.map", 
				"parser", "tped", "--output-file","/home/ralonso/ped/test","--home", System.getenv("BABELOMICS_HOME")};

		try {
//			for(int i=0; i < args.length; i++)
//				System.out.println(args[i]);
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
