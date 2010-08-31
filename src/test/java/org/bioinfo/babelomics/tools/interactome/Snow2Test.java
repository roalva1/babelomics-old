package org.bioinfo.babelomics.tools.interactome;

import org.bioinfo.babelomics.BabelomicsMain;
import org.junit.Test;

public class Snow2Test {

	@Test
	public void Test() {

		//Comando de Alicia
		String []args2 = {
				"--tool", "snow2", 
				"-o", "/tmp/",  
				"-s", "/mnt/commons/babelomics/tests/snow2/ej1/sce_alldb_proteins_interactome_nr.sif", 
//				"-t", "/home/ralonso/appl/babelomics/sce/proteins/sce_alldb_proteins_interactome_nr_topo.txt",
				"--randoms", "1",
				"--randoms-size","5", 
				"--o-topo-file","/mnt/commons/babelomics/tests/snow2/ej1/result",
				"--o-means-file","/mnt/commons/babelomics/tests/snow2/ej1/result",
//				"--node-file1","/home/ralonso/appl/babelomics/sce/proteins/Ejemplo1/UPYDOWN_HET_list_uniq",
				"--node-file1","/mnt/commons/babelomics/tests/snow2/ej1/list2",
				"--o-components-file","/mnt/commons/babelomics/tests/snow2/ej1/result",
				"--o-stats-file","/mnt/commons/babelomics/tests/snow2/ej1/result",
				"--o-images-file","/mnt/commons/babelomics/tests/snow2/ej1/result",
				"--side", "less",
				"--intermediate",
				"--home", System.getenv("BABELOMICS_HOME")};
		String comando = "./babelomics.sh --tool snow2 -o /tmp/ -s /home/ralonso/appl/babelomics/sce/proteins/sce_alldb_proteins_interactome_nr.sif -t /home/ralonso/appl/babelomics/sce/proteins/sce_alldb_proteins_interactome_nr_topo.txt --randoms 1000 --randoms-size 50 --o-topo-file /home/ralonso/appl/babelomics/sce/proteins/Ejemplo1/result --o-means-file /home/ralonso/appl/babelomics/sce/proteins/Ejemplo1/result --node-file1 /home/ralonso/appl/babelomics/sce/proteins/Ejemplo1/UPYDOWN_HET_list_uniq --o-components-file /home/ralonso/appl/babelomics/sce/proteins/Ejemplo1/result --o-stats-file /home/ralonso/appl/babelomics/sce/proteins/Ejemplo1/result --side less";

		try {
			BabelomicsMain.main(args2);
				
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
