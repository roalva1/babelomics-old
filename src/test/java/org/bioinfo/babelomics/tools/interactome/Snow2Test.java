package org.bioinfo.babelomics.tools.interactome;
import org.bioinfo.babelomics.BabelomicsMain;
import org.junit.Test;

public class Snow2Test {
	
	@Test
	public void Test() {
		
//		#1:If you have sif-file and you want to create a .topo file from the full sif-file (obviously no randoms)
//		String []args = {"--tool", "snow2", "-o", "/tmp/",  "--sif-file", "/home/ralonso/appl/babelomics/ej8/ej8.sif", "--o-file-topo", "/home/ralonso/appl/babelomics/ej8/ej8", "--home", System.getenv("BABELOMICS_HOME")};

// 		#2:If you have sif-file and you want to create a .topo and/or .mcn file random list without a given size
//		String []args = {"--tool", "snow2", "-o", "/tmp/",  "--sif-file", "/home/ralonso/appl/bioinfo-networks/mine.sif", "--file-topo", "/home/ralonso/appl/bioinfo-networks/mineTopo", "--file-mcn", "/home/ralonso/appl/bioinfo-networks/mineMcn", "--randoms", "2","--home", System.getenv("BABELOMICS_HOME")};

//		#3:If you have sif-file and you want to create a .topo and/or .mcn file random list with a given size
//		String []args = {"--tool", "snow2", "-o", "/tmp/",  "--sif-file", "/home/ralonso/appl/babelomics/sce/proteins/sce_alldb_proteins_interactome_nr.sif", "--o-file-mcn", "/home/ralonso/appl/bioinfo-networks/mineMcn", "--randoms", "10", "--random-size", "30","--home", System.getenv("BABELOMICS_HOME")};

// 		#4:If you have sif-file and you want to create a .topo and/or .mcn file from a given node-list
//		String []args = {"--tool", "snow2", "-o", "/tmp/", "--sif-file", "/home/ralonso/appl/bioinfo-networks/mine.sif", "--node-list", "v1,v2,v0",  "--file-topo", "/home/ralonso/appl/bioinfo-networks/mineTopo", "--file-mcn", "/home/ralonso/appl/bioinfo-networks/mineMcn", "--home", System.getenv("BABELOMICS_HOME")};

//		#5:If you have file-topo-values and you want to create a .mcn file from a given node-list but with the values of the file-topo-values
//		String []args = {"--tool", "snow2", "-o", "/tmp/file-topo-values", "--file-topo-values", "/home/ralonso/appl/bioinfo-networks/mineTopo.topo", "--node-list", "v1,v2,v0", "--file-mcn", "/home/ralonso/appl/bioinfo-networks/mineMcn", "--home", System.getenv("BABELOMICS_HOME")};

//		Case: You get the values for a network, but with the values of the subnet
//		String []args = {"--tool", "snow2", "-o",".", "--sif-file", "sce_alldb_proteins_interactome_nr.sif", "--node-file", "UPYDOWN_HET_list_uniq", "--file-mcn", "fileMcn", "--file-topo", "fileTopo","--home", System.getenv("BABELOMICS_HOME")};
		
//		String []args = {"--tool", "snow2", "-o","/opt/babelomics/", "--sif-file", "/opt/babelomics/sce_alldb_proteins_interactome_nr.sif", "--randoms", "2", "--random-size", "50", "--home", System.getenv("BABELOMICS_HOME")};

//		Case: We want to do an statistic Kolmogorov test 
//		String []args = {"--tool", "snow2", "-o","/home/ralonso/appl/babelomics/", "--sif-file", "/home/ralonso/appl/babelomics/ej8/ej8.sif", "--file-topo-values", "/home/ralonso/appl/babelomics/ej8/ej8.topo", "--node-file1", "/home/ralonso/appl/babelomics/ej8/list1", "--side", "greater", "--randoms", "10", "--random-size", "3", "--o-file-components", "/home/ralonso/appl/babelomics/ej8/ej8","--images","--home", System.getenv("BABELOMICS_HOME")};
//		String []args = {"--tool", "snow2", "-o","/home/ralonso/appl/babelomics/", "--sif-file", "/home/ralonso/appl/babelomics/sce/proteins/sce_alldb_proteins_interactome_nr.sif", "--file-topo-values", "/home/ralonso/appl/babelomics/sce/proteins/sce_alldb_proteins_interactome_nr.topo", "--node-file1", "/home/ralonso/appl/babelomics/sce/proteins/Ejemplo1/UPYDOWN_HET_list_uniq", "--side", "less", "--randoms", "1000", "--random-size", "100", "--images","--home", System.getenv("BABELOMICS_HOME")};

//		Case: We want to do an statistic Kolmogorov and Wilcoxon test with 2 lists 
//		String []args = {"--tool", "snow2", "-o","/home/ralonso/appl/babelomics/", "--sif-file", "/home/ralonso/appl/babelomics/sce/proteins/sce_alldb_proteins_interactome_nr.sif", "--file-topo-values", "/home/ralonso/appl/babelomics/sce/proteins/sce_alldb_proteins_interactome_nr.topo", "--node-file1", "/home/ralonso/appl/babelomics/sce/proteins/Ejemplo1/UPYDOWN_HET_list_uniq", "--node-file2", "/home/ralonso/appl/babelomics/sce/proteins/Ejemplo2/UPYDOWN_HOM_list_uniq", "--side", "greater", "--images","--home", System.getenv("BABELOMICS_HOME")};
		String []args = {"--tool", "snow2", "-o","/home/ralonso/appl/babelomics/", "--sif-file", "/opt/babelomics/hsa_alldb_genes_interactome_nr.sif", "--node-file1", "/opt/babelomics/snpbreastcancer_down_114.txt","--o-file-mcn","/opt/babelomics/snpbreastcancer_down_114", "--o-file-components", "/opt/babelomics/snpbreastcancer_down_114", "--no-bicomponents", "--home", System.getenv("BABELOMICS_HOME")};

		try {
			for(int i=0; i < args.length; i++)
				System.out.println(args[i]);
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
