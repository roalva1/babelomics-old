package org.bioinfo.babelomics.tools.interactome;

import org.bioinfo.babelomics.BabelomicsMain;
import org.junit.Test;

public class Snow2Test {

	@Test
	public void Test() {
//		#1:If you have sif-file and you want to create a .topo file from the full sif-file (obviously no randoms)
//		String []args = {"--tool", "snow2r", "-o", "/tmp/",  "-s", "/home/ralonso/appl/babelomics/sce/proteins/sce_alldb_proteins_interactome_nr.sif", "--o-sif-topo-file", "/home/ralonso/appl/babelomics/sce/proteins/sce_alldb_proteins_interactome_nr_topo.txt", "--home", System.getenv("BABELOMICS_HOME")};

		// Creating randoms
//		String []args = {"--tool", "snow2r", "-o", "/tmp/",  "-s", "/home/ralonso/appl/babelomics/ej8/ej8.sif", "--randoms", "4","--randoms-size","3", 
//				"--o-topo-file","/home/ralonso/appl/babelomics/ej8/ej8","--o-median-file","/home/ralonso/appl/babelomics/ej8/ej8","--home", System.getenv("BABELOMICS_HOME")};

//		 Creating randoms + 2 subnets with
		String []args = {
				"--tool", "snow2", 
				"-o", "/tmp/",  
				"-s", "/home/ralonso/appl/babelomics/ej8/ej8.sif", 
				"-t","/home/ralonso/appl/babelomics/ej8/ej8.topo",
//				"--randoms", "4",
//				"--randoms-size","3", 
				"--o-topo-file","/home/ralonso/appl/babelomics/ej8/ej8",
				"--o-means-file","/home/ralonso/appl/babelomics/ej8/ej8",
				"--node-file1","/home/ralonso/appl/babelomics/ej8/list1",
//				"--node-file2","/home/ralonso/appl/babelomics/ej8/list2",
				"--o-components-file","/home/ralonso/appl/babelomics/ej8/ej8",
				"--o-stats-file","/home/ralonso/appl/babelomics/ej8/ej8",
				"--o-images-file","/home/ralonso/appl/babelomics/ej8/ej8",
				"--side", "less",
//				"--o-json-file","/home/ralonso/appl/babelomics/ej8/ej8",
//				"--o-dot-file","/home/ralonso/appl/babelomics/ej8/ej8dot",
//				"--o-svg-file","/home/ralonso/appl/babelomics/ej8/ej8svg",
//				"--intermediate",
				"--home", System.getenv("BABELOMICS_HOME")};

		//Comando de Alicia
		String []args2 = {
				"--tool", "snow2", 
				"-o", "/tmp/",  
				"-s", "/home/ralonso/appl/babelomics/sce/proteins/sce_alldb_proteins_interactome_nr.sif", 
//				"-t", "/home/ralonso/appl/babelomics/sce/proteins/sce_alldb_proteins_interactome_nr_topo.txt",
				"--randoms", "1",
				"--randoms-size","5", 
				"--o-topo-file","/home/ralonso/appl/babelomics/sce/proteins/Ejemplo1/result",
				"--o-means-file","/home/ralonso/appl/babelomics/sce/proteins/Ejemplo1/result",
//				"--node-file1","/home/ralonso/appl/babelomics/sce/proteins/Ejemplo1/UPYDOWN_HET_list_uniq",
				"--node-file1","/home/ralonso/appl/babelomics/sce/proteins/Ejemplo1/list2",
				"--o-components-file","/home/ralonso/appl/babelomics/sce/proteins/Ejemplo1/result",
				"--o-stats-file","/home/ralonso/appl/babelomics/sce/proteins/Ejemplo1/result",
				"--o-images-file","/home/ralonso/appl/babelomics/sce/proteins/Ejemplo1/result",
				"--side", "less",
				"--intermediate",
				"--home", System.getenv("BABELOMICS_HOME")};
		String comando = "./babelomics.sh --tool snow2 -o /tmp/ -s /home/ralonso/appl/babelomics/sce/proteins/sce_alldb_proteins_interactome_nr.sif -t /home/ralonso/appl/babelomics/sce/proteins/sce_alldb_proteins_interactome_nr_topo.txt --randoms 1000 --randoms-size 50 --o-topo-file /home/ralonso/appl/babelomics/sce/proteins/Ejemplo1/result --o-means-file /home/ralonso/appl/babelomics/sce/proteins/Ejemplo1/result --node-file1 /home/ralonso/appl/babelomics/sce/proteins/Ejemplo1/UPYDOWN_HET_list_uniq --o-components-file /home/ralonso/appl/babelomics/sce/proteins/Ejemplo1/result --o-stats-file /home/ralonso/appl/babelomics/sce/proteins/Ejemplo1/result --side less";

		try {
//			for(int i=0; i < args.length; i++)
//				System.out.println(args[i]);
			//for (int i=0; i < 100; i++){
				//System.out.println("Numero:"+i);
				BabelomicsMain.main(args2);
//			}
				
		} catch (Exception e) {
			e.printStackTrace();
		}
	}

}
