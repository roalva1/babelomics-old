package org.bioinfo.babelomics.tools.interactome;

import org.bioinfo.babelomics.BabelomicsMain;
import org.junit.Test;

public class Snow2Test {

	
	@Test
	public void Test1() {
		//Comando de Alicia
		String []args = {
				"--tool", "snow2", 
				"-o", "/tmp/snow2/test1/",  
//				"-s", "/mnt/commons/babelomics/tests/snow2/ej1/sce_alldb_proteins_interactome_nr.sif", 
//				"-t", "/home/ralonso/appl/babelomics/sce/proteins/sce_alldb_proteins_interactome_nr_topo.txt",
				"--randoms", "10",
				"--randoms-size","2", 
				"--o-name","result",
//				"--list1","/home/ralonso/appl/babelomics/sce/proteins/Ejemplo1/UPYDOWN_HET_list_uniq",
				"--list1","/mnt/commons/babelomics/tests/snow2/ej1/list2",
//				"--species","sce",
				"--interactome","sce",
				"--side", "less",
				"--intermediate",
				"--images",
				"--json", "1",
				//"--svg", "1",
				"--home", System.getenv("BABELOMICS_HOME")};
		try {
//			for(String arg : args)
//				System.out.println(arg);
			BabelomicsMain.main(args);
				
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
	
	@Test
	public void Test2() {

		//Comando de Alicia
		String []args2 = {
				"--tool", "snow2", 
				"-o", "/tmp/snow2/test2/",
//				"-sif-file", "/mnt/commons/babelomics/tests/snow2/ej1/sce_alldb_proteins_interactome_nr.sif", 
//				"-topo-file", "/home/ralonso/appl/babelomics/sce/proteins/sce_alldb_proteins_interactome_nr_topo.txt",
				"--randoms", "10",
				"--randoms-size","5", 
				"--o-name","result",
				"--interactome","sce",
//				"--o-topo-file","/tmp/snow2/result",
//				"--o-means-file","/tmp/snow2/result",
				"--list1","/mnt/commons/babelomics/tests/snow2/ej1/UPYDOWN_HET_list_uniq",
				"--list2","/mnt/commons/babelomics/tests/snow2/ej1/list2",
//				"--o-components-file","/tmp/snow2/result",
//				"--o-stats-file","/tmp/snow2/result",
//				"--o-images-file","/tmp/snow2/result",
				"--side", "less",
				"--intermediate",
				"--json", "1",
				//"--svg", "1",
				"--home", System.getenv("BABELOMICS_HOME")};
		//String comando = "./babelomics.sh --tool snow2 -o /tmp/ -s /home/ralonso/appl/babelomics/sce/proteins/sce_alldb_proteins_interactome_nr.sif -t /home/ralonso/appl/babelomics/sce/proteins/sce_alldb_proteins_interactome_nr_topo.txt --randoms 1000 --randoms-size 50 --o-topo-file /home/ralonso/appl/babelomics/sce/proteins/Ejemplo1/result --o-means-file /home/ralonso/appl/babelomics/sce/proteins/Ejemplo1/result --node-file1 /home/ralonso/appl/babelomics/sce/proteins/Ejemplo1/UPYDOWN_HET_list_uniq --o-components-file /home/ralonso/appl/babelomics/sce/proteins/Ejemplo1/result --o-stats-file /home/ralonso/appl/babelomics/sce/proteins/Ejemplo1/result --side less";

		try {
			BabelomicsMain.main(args2);
				
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

	@Test
	public void Test3() {

		//Comando de Alicia
		String []args2 = {
				"--tool", "snow2", 
				"-o", "/tmp/snow2/test3/",
				"-sif-file", "/mnt/commons/babelomics/tests/snow2/ej1/sce_alldb_proteins_interactome_nr.sif", 
//				"-topo-file", "/home/ralonso/appl/babelomics/sce/proteins/sce_alldb_proteins_interactome_nr_topo.txt",
				"--randoms", "10",
				"--randoms-size","5", 
				"--o-name","result",
				"--interactome","own",
//				"--o-topo-file","/tmp/snow2/result",
//				"--o-means-file","/tmp/snow2/result",
				"--list1","/mnt/commons/babelomics/tests/snow2/ej1/UPYDOWN_HET_list_uniq",
				"--list2","/mnt/commons/babelomics/tests/snow2/ej1/list2",
//				"--o-components-file","/tmp/snow2/result",
//				"--o-stats-file","/tmp/snow2/result",
//				"--o-images-file","/tmp/snow2/result",
				"--side", "less",
				"--intermediate",
				"--json", "1",
				//"--svg", "1",
				"--home", System.getenv("BABELOMICS_HOME")};
		//String comando = "./babelomics.sh --tool snow2 -o /tmp/ -s /home/ralonso/appl/babelomics/sce/proteins/sce_alldb_proteins_interactome_nr.sif -t /home/ralonso/appl/babelomics/sce/proteins/sce_alldb_proteins_interactome_nr_topo.txt --randoms 1000 --randoms-size 50 --o-topo-file /home/ralonso/appl/babelomics/sce/proteins/Ejemplo1/result --o-means-file /home/ralonso/appl/babelomics/sce/proteins/Ejemplo1/result --node-file1 /home/ralonso/appl/babelomics/sce/proteins/Ejemplo1/UPYDOWN_HET_list_uniq --o-components-file /home/ralonso/appl/babelomics/sce/proteins/Ejemplo1/result --o-stats-file /home/ralonso/appl/babelomics/sce/proteins/Ejemplo1/result --side less";

		try {
			BabelomicsMain.main(args2);
				
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
