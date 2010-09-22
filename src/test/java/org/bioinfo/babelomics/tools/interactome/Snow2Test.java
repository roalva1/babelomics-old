package org.bioinfo.babelomics.tools.interactome;

import org.bioinfo.babelomics.BabelomicsMain;
import org.junit.Test;

public class Snow2Test {

	@Test
	public void Test0() {
	}
	
	@Test
	public void Test1() {

		String []args = {
				"--tool", "snow2", 
				"-o", "/tmp/snow2/test1", 
//				"--log-file","/tmp/snow2/test1/result.log",
//				"--log-level","1",
//				"-s", "/mnt/commons/babelomics/tests/snow2/ej1/sce_alldb_proteins_interactome_nr.sif", 
//				"-t", "/home/ralonso/appl/babelomics/sce/proteins/sce_alldb_proteins_interactome_nr_topo.txt",
//				"--randoms", "10",
//				"--randoms-size","2", 
				"--o-name","result",
				"--list1","/mnt/commons/babelomics/tests/snow2/ej1/UPYDOWN_HET_list_uniq",
				"--list2","/mnt/commons/babelomics/tests/snow2/ej1/UPYDOWN_HET_list_uniq_mitad",
				"--interactome","sce",
				"--type","proteins",
				"--side", "less",
//				"--intermediate",
//				"--components",
//				"--bicomponents",
				"--images",
				"--json",
				"--home", System.getenv("BABELOMICS_HOME")};
//		String comando="./babelomics.sh --tool snow2 -o /tmp/snow2/test1 --o-name result --list1 /mnt/commons/babelomics/tests/snow2/ej1/UPYDOWN_HET_list_uniq --list2 /mnt/commons/babelomics/tests/snow2/ej1/list2 --interactome sce --side less";
		main(args);
	}
	
//	@Test
	public void Test2() {

		String []args = {
				"--tool", "snow2", 
				"--sif-file", "/mnt/commons/babelomics/tests/snow2/ej8/ej8.sif",
//				"-t", "/mnt/commons/babelomics/tests/snow2/ej8/ej8.topo", 
				"-o", "/tmp/snow2/test2",
//				"--o-sif-topo-file",
//				"-topo-file", "/home/ralonso/appl/babelomics/sce/proteins/sce_alldb_proteins_interactome_nr_topo.txt",
				"--randoms", "10",
				"--randoms-size","2", 
				"--o-name","resultSmall",
				"--interactome","own",
				"--type", "proteins",
				"--list1","/mnt/commons/babelomics/tests/snow2/ej8/list3.txt",
//				"--list2","/mnt/commons/babelomics/tests/snow2/ej8/list2",
				"--side", "less",
//				"--intermediate",
				"--json",
				"--images",
				"--home", System.getenv("BABELOMICS_HOME")};

		main(args);
	}
	
	//@Test
	public void DotSvgJsonTest() {

		String []args = {
				"--tool", "snow2", 
				"-o", "/tmp/ej8/",  
				"-s", "/mnt/commons/babelomics/tests/snow2/ej8/ej8.sif",
				"-t", "/mnt/commons/babelomics/tests/snow2/ej8/ej8.topo",
				"--node-file1","/mnt/commons/babelomics/tests/snow2/ej8/list1",
				"--node-file2","/mnt/commons/babelomics/tests/snow2/ej8/list2",
//				"--dot",
				"--o-name", "result",
				"--intermediate",
//				"--bicomponents",
				"--json",
				"--home", System.getenv("BABELOMICS_HOME")};

		main(args);
	}
	
	@Test
	public void Test3() {

		String []args = {
				"--tool", "snow2", 
				"-o", "/tmp/snow2/test3", 
//				"--log-file","/tmp/snow2/test1/result.log",
//				"--log-level","1",
//				"-s", "/mnt/commons/babelomics/tests/snow2/ej1/sce_alldb_proteins_interactome_nr.sif", 
//				"-t", "/home/ralonso/appl/babelomics/sce/proteins/sce_alldb_proteins_interactome_nr_topo.txt",
//				"--randoms", "10",
//				"--randoms-size","2", 
				"--o-name","result",
				"--list1","/mnt/commons/babelomics/tests/snow2/ej1/UPYDOWN_HET_list_uniq",
				"--list2","/mnt/commons/babelomics/tests/snow2/ej1/UPYDOWN_HET_list_uniq_mitad",
				"--interactome","sce",
				"--type","proteins",
				"--side", "less",
//				"--intermediate",
				"--components",
//				"--bicomponents",
				"--images",
				"--json",
				"--home", System.getenv("BABELOMICS_HOME")};
//		String comando="./babelomics.sh --tool snow2 -o /tmp/snow2/test1 --o-name result --list1 /mnt/commons/babelomics/tests/snow2/ej1/UPYDOWN_HET_list_uniq --list2 /mnt/commons/babelomics/tests/snow2/ej1/list2 --interactome sce --side less";
		main(args);
	}

	public void main(String []args){
		try {
//			for(String arg : args)
//				System.out.println(arg);
			BabelomicsMain.main(args);
				
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
