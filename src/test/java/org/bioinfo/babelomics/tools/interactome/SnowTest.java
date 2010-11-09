package org.bioinfo.babelomics.tools.interactome;

import java.io.File;

import org.bioinfo.babelomics.BabelomicsMain;
import org.junit.Test;

public class SnowTest {

	@Test
	public void newStatsTest(){
		//./babelomics.sh --tool snow2 --outdir /tmp/ --sif-file /mnt/commons/babelomics/tests/snow2/ej8/ej8.sif  --o-sif-topo-file --o-name prueba --interactome own
		String outdir = "/tmp/snow2/newStatsTest";
		new File(outdir).mkdirs();
		String []args = {
				"--tool", "snow2", 
				"-o", outdir, 
				"--sif-file", "/mnt/commons/babelomics/tests/snow2/ej8/ej9.sif",
				"--topo-file", "/mnt/commons/babelomics/tests/snow2/ej8/ej9_topo.txt",
				"--list1","/mnt/commons/babelomics/tests/snow2/ej8/list2",
				"--interactome","own",
				"--type","proteins",
				"--randoms", "500",
				"--components", "1",
				"--bicomponents", "1",
				"--intermediate", "1",
				"--images", "1",
				"--xml",
				"--sif",
				"--side", "less",
				"--home", System.getenv("BABELOMICS_HOME")};
		main(args);
	}
//	@Test
	public void createTopoFile(){
		//./babelomics.sh --tool snow2 --outdir /tmp/ --sif-file /mnt/commons/babelomics/tests/snow2/ej8/ej8.sif  --o-sif-topo-file --o-name prueba --interactome own
		String outdir = "/mnt/commons/babelomics/tests/snow2/ej8/";
		new File(outdir).mkdirs();
		String []args = {
				"--tool", "snow2", 
				"-o", outdir, 
				"--sif-file", "/mnt/commons/babelomics/tests/snow2/ej8/ej9.sif",
				"--o-sif-topo-file",
				"--interactome","own",
				"--home", System.getenv("BABELOMICS_HOME")};
		main(args);
	}

	//@Test
	public void TestGenes() {

		String outdir = "/tmp/snow2/testGenes";
		new File(outdir).mkdirs();

		String []args = {
				"--tool", "snow2", 
				"-o", outdir, 
				"--randoms", "1",
				"--o-name","result",
				"--list1","/mnt/commons/babelomics/tests/snow2/listas/hsa/genes/list1.txt",
				"--interactome","hsa",
				"--type","genes",
				"--side", "less",
				"--components", "1",
				"--intermediate", "0",
				"--images", "1",
				"--xml",
				"--home", System.getenv("BABELOMICS_HOME")};
		main(args);
	}
	
	
	
	//@Test
	public void Test1() {

		String outdir = "/tmp/snow2/test1";
		new File(outdir).mkdirs();

		String []args = {
				"--tool", "snow2", 
				"-o", outdir, 
				"--o-name","result",
				"--list1","/mnt/commons/babelomics/tests/snow2/ej1/UPYDOWN_HET_list_uniq",
				"--list2","/mnt/commons/babelomics/tests/snow2/ej1/UPYDOWN_HET_list_uniq_mitad",
				"--interactome","sce",
				"--type","proteins",
				"--side", "less",
				"--components", "1",
				"--intermediate", "0",
				"--images", "1",
				"--xml",
				"--home", System.getenv("BABELOMICS_HOME")};
		main(args);
	}


	//@Test
	public void TestMouse() {

		String outdir = "/tmp/snow2/testMouse";
		new File(outdir).mkdirs();

		String []args = {
				"--tool", "snow2", 
				"-o", outdir, 
				"--sif-file", "/mnt/commons/babelomics/tests/snow2/listas/own/proteins/mmu_alldb_proteins_interactome_nr.sif",
				//				"-t", "/mnt/commons/babelomics/tests/snow2/ej8/ej8.topo", 
				//				"--o-sif-topo-file",
				//				"-topo-file", "/home/ralonso/appl/babelomics/sce/proteins/sce_alldb_proteins_interactome_nr_topo.txt",
				"--randoms", "500",
				"--o-name","resultSmall",
				"--interactome","own",
				"--type", "proteins",
				"--list1","/mnt/commons/babelomics/tests/snow2/listas/own/proteins/mmu_example.txt",
				//				"--list2","/mnt/commons/babelomics/tests/snow2/ej8/list2",
				"--side", "less",
				//				"--intermediate",
				"--json",
				"--images",
				"--home", System.getenv("BABELOMICS_HOME")};
		//		String comando="./babelomics.sh --tool snow2 -o /tmp/snow2/test1 --o-name result --list1 /mnt/commons/babelomics/tests/snow2/ej1/UPYDOWN_HET_list_uniq --list2 /mnt/commons/babelomics/tests/snow2/ej1/list2 --interactome sce --side less";
		main(args);
	}

	public void Test2() {

		String outdir = "/tmp/snow2/test2";
		new File(outdir).mkdirs();

		String []args = {
				"--tool", "snow2", 
				"-o", outdir, 
				"--sif-file", "/mnt/commons/babelomics/tests/snow2/ej8/ej8.sif",
				//				"-t", "/mnt/commons/babelomics/tests/snow2/ej8/ej8.topo", 
				//				"--o-sif-topo-file",
				//				"-topo-file", "/home/ralonso/appl/babelomics/sce/proteins/sce_alldb_proteins_interactome_nr_topo.txt",
				"--randoms", "10",
				"--o-name","resultSmall",
				"--interactome","own",
				"--type", "proteins",
				"--list1","/mnt/commons/babelomics/tests/snow2/ej8/list1",
				//				"--list2","/mnt/commons/babelomics/tests/snow2/ej8/list2",
				"--side", "less",
				"--intermediate", "1",
				"--json",
				"--images",
				"--home", System.getenv("BABELOMICS_HOME")};

		main(args);
	}

	//	@Test
	public void SnowExampleOneList(){
		///httpd/bioinfo/babelomics/babelomics.sh --tool snow2 --outdir /httpd/bioinfo/wum_sessions_v0.7/4164/jobs/2964 --log-file /httpd/bioinfo/wum_sessions_v0.7/4164/jobs/2964/job.log --list2 none --randoms-size 2 --json 1 --list1 /httpd/bioinfo/wum_sessions_v0.7/4164/data/27302/chr_9_block6.txt --side less --images  --randoms 10 --interactome hsa --o-name result

		String outdir = "/tmp/snow2/test2";
		new File(outdir).mkdirs();

		String []args = {
				"--tool", "snow2", 
				"-o", outdir, 
				//					"-t", "/mnt/commons/babelomics/tests/snow2/ej8/ej8.topo", 
				//					"--o-sif-topo-file",
				//					"-topo-file", "/home/ralonso/appl/babelomics/sce/proteins/sce_alldb_proteins_interactome_nr_topo.txt",
				"--randoms", "500",
				//					"--randoms-size","2", 
				"--o-name","result",
				"--interactome","hsa",
				"--type", "proteins",
				"--list1","/mnt/commons/babelomics/tests/snow2/listas/hsa/proteins/chr_10_block225.txt",
				"--side", "less",
				"--intermediate",
				"--json",
				"--images",
				"--home", System.getenv("BABELOMICS_HOME")};
		main(args);
	}

	//@Test
	public void SnowExampleTwoLists(){
		///httpd/bioinfo/babelomics/babelomics.sh --tool snow2 --outdir /httpd/bioinfo/wum_sessions_v0.7/4164/jobs/2964 --log-file /httpd/bioinfo/wum_sessions_v0.7/4164/jobs/2964/job.log --list2 none --randoms-size 2 --json 1 --list1 /httpd/bioinfo/wum_sessions_v0.7/4164/data/27302/chr_9_block6.txt --side less --images  --randoms 10 --interactome hsa --o-name result

		String outdir = "/tmp/snow2/testTwoLists";
		new File(outdir).mkdirs();

		String []args = {
				"--tool", "snow2", 
				"-o", outdir, 
				"--o-name","result",
				"--components", "1",
				"--interactome","hsa",
				"--type", "proteins",
				"--list1","/mnt/commons/babelomics/tests/snow2/listas/hsa/proteins/chr_10_block225.txt",
				"--list2","/mnt/commons/babelomics/tests/snow2/listas/hsa/proteins/chr_11_block80.txt",
				"--side", "less",
				"--intermediate", "0",
				"--xml",
				"--images",
				"--home", System.getenv("BABELOMICS_HOME")};
		main(args);
	}

	public void Test3() {

		String outdir = "/tmp/snow2/test3";
		new File(outdir).mkdirs();

		String []args = {
				"--tool", "snow2", 
				"-o", outdir, 
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
				"--intermediate", "0",
				"--components", "1",
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
