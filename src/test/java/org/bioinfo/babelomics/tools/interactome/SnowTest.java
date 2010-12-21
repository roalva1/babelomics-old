package org.bioinfo.babelomics.tools.interactome;

import java.io.File;

import org.bioinfo.babelomics.BabelomicsMain;
import org.junit.Test;

public class SnowTest {

	
	//@Test
	public void SnowExampleOneList(){
		///httpd/bioinfo/babelomics/babelomics.sh --tool snow2 --outdir /httpd/bioinfo/wum_sessions_v0.7/4164/jobs/2964 --log-file /httpd/bioinfo/wum_sessions_v0.7/4164/jobs/2964/job.log --list2 none --randoms-size 2 --json 1 --list1 /httpd/bioinfo/wum_sessions_v0.7/4164/data/27302/chr_9_block6.txt --side less --images  --randoms 10 --interactome hsa --o-name result

		String outdir = "/tmp/snow2/test1";
		new File(outdir).mkdirs();

		String []args = {
				"--tool", "snow2", 
				"-o", outdir, 
				"--randoms", "1",
				"--o-name","result",
				"--interactome","hsa",
				"--group","all",
				"--type", "proteins",
				"--list1","/mnt/commons/babelomics/tests/snow2/listas/hsa/proteins/chr_10_block225.txt",
				"--side", "less",
				"--intermediate","0",
				"--components","1",
				"--xml",
				"--images",
				"--home", System.getenv("BABELOMICS_HOME")};
		main(args);
	}
	
	//@Test
	public void genesTest(){

		String outdir = "/tmp/snow2/test2";
		new File(outdir).mkdirs();

		String []args = {
				"--tool", "snow2", 
				"-o", outdir, 
				"--randoms", "50",
				"--o-name","result",
				"--interactome","hsa",
				"--group","all",
				"--type", "genes",
				"--list1","/mnt/commons/babelomics/tests/snow2/listas/hsa/genes/list1.txt",
				"--side", "less",
				"--intermediate","1",
				"--components","1",
				"--xml",
				"--images",
				"--home", System.getenv("BABELOMICS_HOME")};
		main(args);
	}
	@Test
	public void transcriptsTest(){

		String outdir = "/tmp/snow2/test3";
		new File(outdir).mkdirs();

		String []args = {
				"--tool", "snow2", 
				"-o", outdir, 
				"--randoms", "50",
				"--o-name","result",
				"--interactome","hsa",
				"--group","all",
				"--type", "transcripts",
				"--list1","/mnt/commons/babelomics/tests/snow2/listas/hsa/transcripts/list1.txt",
				"--side", "less",
				"--intermediate","1",
//				"--components","1",
				"--xml",
				"--images",
				"--home", System.getenv("BABELOMICS_HOME")};
		main(args);
	}

//	@Test
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
				"--json",
				//"--sif",
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
	public void ownInteractomeTest(){
		///httpd/bioinfo/babelomics/babelomics.sh --tool snow2 --outdir /httpd/bioinfo/wum_sessions_v0.7/4164/jobs/2964 --log-file /httpd/bioinfo/wum_sessions_v0.7/4164/jobs/2964/job.log --list2 none --randoms-size 2 --json 1 --list1 /httpd/bioinfo/wum_sessions_v0.7/4164/data/27302/chr_9_block6.txt --side less --images  --randoms 10 --interactome hsa --o-name result

		String outdir = "/tmp/snow2/ownTest";
		new File(outdir).mkdirs();

		String []args = {
				"--tool", "snow2", 
				"-o", outdir, 
				"--randoms", "2",
				"--o-name","result",
				"--interactome","own",
				"--group","all",
				"--sif-file","/mnt/commons/babelomics/tests/snow2/ej8/ej8.sif",
				"--list1","/mnt/commons/babelomics/tests/snow2/ej8/list1",
				"--side", "less",
				"--intermediate","1",
				"--components","1",
				"--xml",
				"--sif",
				"--images",
				"--home", System.getenv("BABELOMICS_HOME")};
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
