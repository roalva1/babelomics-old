package org.bioinfo.babelomics.tools.interactome;

import java.io.File;

import org.bioinfo.babelomics.BabelomicsMain;
import org.junit.Test;

public class GSnowTest {
	@Test
	public void example1(){
		// Essential_genes_in_cancer_cell_line_K562
		// /opt/babelomics/babelomics.sh --tool network-miner --outdir /httpd/bioinfo/wum_sessions_v0.7/4164/jobs/5880 --log-file /httpd/bioinfo/wum_sessions_v0.7/4164/jobs/5880/job.log --order ascendant --significant-value 0.05 --list /opt/babelomics/example/K562_symbol.txt --randoms 1000 --components true --interactome hsa --intermediate 0 --group curated --type genes --o-name result
		String outdir = "/tmp/gsnow/example1";
		new File(outdir).mkdirs();

		String []args = {
				"--tool", "network-miner", 
				"-o", outdir, 
				"--o-name","result",
				"--interactome","hsa",
				"--type", "genes",
				"--group", "curated",
				"--intermediate","0",
				"--randoms","1000",
				"--components","true",
				"--order","ascendant",
				"--significant-value","0.05",
				"--list","/opt/babelomics/example/K562_symbol.txt",
				"--home", System.getenv("BABELOMICS_HOME")};
		main(args);
	}
	//@Test
	public void example2(){
		// Genes_up_in_control_Vs_case_Hirschsprung_disease
		// /opt/babelomics/babelomics.sh --tool network-miner --outdir /httpd/bioinfo/wum_sessions_v0.7/4164/jobs/5879 --log-file /httpd/bioinfo/wum_sessions_v0.7/4164/jobs/5879/job.log --order descendant --significant-value 0.05 --list /opt/babelomics/example/Genes_up_in_control_Vs_case_Hirschsprung_disease.txt --randoms 1000 --interactome hsa --components true --intermediate 1 --group curated --type genes --o-name result
		String outdir = "/tmp/gsnow/example2";
		new File(outdir).mkdirs();

		String []args = {
				"--tool", "network-miner", 
				"-o", outdir, 
				"--o-name","result",
				"--interactome","hsa",
				"--type", "genes",
				"--group", "curated",
				"--intermediate","1",
				"--randoms","10",
				"--components","true",
				"--order","descendant",
				"--significant-value","0.05",
				"--list","/opt/babelomics/example/Genes_up_in_control_Vs_case_Hirschsprung_disease.txt",
				"--home", System.getenv("BABELOMICS_HOME")};
		main(args);
		
	}
	//@Test
	public void test2(){

		String outdir = "/tmp/gsnow/test2";
		new File(outdir).mkdirs();

		String []args = {
				"--tool", "network-miner", 
				"-o", outdir, 
				"--o-name","output",
				"--interactome","hsa",
				"--type", "genes",
				"--group", "curated",
//				"--select-mcn","abs-min",
				"--intermediate","1",
//				"--side","less",
				"--randoms","2",
				"--components","true",
//				"--number-items","50",
				"--order","ascendant",
				"--significant-value","1000000",
//				"--cut-off","5",
//				"--list","/mnt/commons/babelomics/tests/snow2/listas/hsa/transcripts/list1.txt",
//				"--list","/mnt/commons/babelomics/tests/snow2/listas/eco/genes/list1.txt",
				"--list","/mnt/commons/babelomics/tests/snow2/listas/hsa/proteins/snow_viz.txt",/*chr_11_block80.txt,chr_13_block95.txt, chr_20_block115.txt */
				"--home", System.getenv("BABELOMICS_HOME")};
		
		//ensg_translated_JURKAT.txt, K562.txt, ensg_translated_REH.txt, ensg_K562.txt
		//sce_prots_mit_uniq.txt,UPYDOWN_HET_list_uniq
		main(args);
		///opt/babelomics/babelomics.sh --tool gsnow --outdir /httpd/bioinfo/wum_sessions_v0.7/4164/jobs/3865 --log-file /httpd/bioinfo/wum_sessions_v0.7/4164/jobs/3865/job.log --significant-value 0.05 --list /httpd/bioinfo/wum_sessions_v0.7/4164/data/27734/hsa.txt --randoms 1000 --components true --interactome hsa --intermediate 1 --group all --number-items 200 --type proteins --o-name result

			
	}

	//	@Test
	public void test1(){

		String outdir = "/tmp/gsnow/test";
		new File(outdir).mkdirs();

		String []args = {
				"--tool", "network-miner", 
				"-o", outdir, 
				"--o-name","output",
				"--interactome","hsa",
				"--type", "proteins",
				"--group", "all",
				"--intermediate","1",
//				"--side","less",
				"--randoms","1000",
				"--components","true",
				"--number-items","200",
//				"--order","descendant",
				"--significant-value","0.05",
//				"--cut-off","5",
				"--list","/httpd/bioinfo/wum_sessions_v0.7/4164/data/27734/hsa.txt",
				//"--list","/mnt/commons/babelomics/tests/snow2/listas/hsa/proteins/aprocesar4.txt",
				"--home", System.getenv("BABELOMICS_HOME")};
		main(args);
		///opt/babelomics/babelomics.sh --tool gsnow --outdir /httpd/bioinfo/wum_sessions_v0.7/4164/jobs/3865 --log-file /httpd/bioinfo/wum_sessions_v0.7/4164/jobs/3865/job.log --significant-value 0.05 --list /httpd/bioinfo/wum_sessions_v0.7/4164/data/27734/hsa.txt --randoms 1000 --components true --interactome hsa --intermediate 1 --group all --number-items 200 --type proteins --o-name result

			
	}
	
//	@Test
	public void testGsnowRandomsGenerator(){

		String outdir = "/tmp/gsnow/testGenerator";
		new File(outdir).mkdirs();

		String []args = {
				"--tool", "network-miner", 
				"-o", outdir, 
				"--o-name","valores.txt",
				"--interactome","hsa",
				"--type", "proteins",
				"--group", "all",
				"--intermediate","0",
				"--size-min","30",
				"--size-max","40",
				"--randoms","5",
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

