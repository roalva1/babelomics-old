package org.bioinfo.babelomics.tools.interactome;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.bioinfo.babelomics.BabelomicsMain;
import org.bioinfo.commons.io.TextFileReader;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.math.stats.inference.KolmogorovSmirnovTest;

public class Snow2Test {
	
//	@Test
//	public void testR() throws REXPMismatchException, NumberFormatException, IOException, REngineException{
//		RConnection c = new RConnection();
//		// Create test data (Java array, R Vector)
//		TextFileReader tfrX = new TextFileReader("x.txt");
//		String line = null;
//		String[] fields;
//		List<Double> xList = new ArrayList<Double>();
//		while((line = tfrX.readLine()) != null) {
//			if(!line.startsWith("#")){
//				fields = line.split("\t");
//				xList.add(Double.parseDouble(fields[0]));
//			}
//		}
//		tfrX.close();	
//		TextFileReader tfrY = new TextFileReader("y.txt");
//		line = null;
//		List<Double> yList = new ArrayList<Double>();
//		while((line = tfrY.readLine()) != null) {
//			if(!line.startsWith("#")){
//				fields = line.split("\t");
//				yList.add(Double.parseDouble(fields[0]));
//			}
//		}
//		tfrY.close();
//		double[] x = ListUtils.toArray(xList);
//		double[] y = ListUtils.toArray(yList);
//		c.assign("x",x);
//		c.assign("y",y);
//		RList l=c.eval("ks.test(x,y, alternative = \"greater\")").asList();
//		REXPDouble result = (REXPDouble) l.get(1);
//	    System.out.println("exp: "+result.toString());
//	}
//	@Test
	public void testReadXY() throws IOException {
		TextFileReader tfrX = new TextFileReader("x.txt");
		String line = null;
		String[] fields;
		List<Double> x = new ArrayList<Double>();
		while((line = tfrX.readLine()) != null) {
			if(!line.startsWith("#")){
				fields = line.split("\t");
				x.add(Double.parseDouble(fields[0]));
			}
		}
		tfrX.close();	
		TextFileReader tfrY = new TextFileReader("y.txt");
		line = null;
		List<Double> y = new ArrayList<Double>();
		while((line = tfrY.readLine()) != null) {
			if(!line.startsWith("#")){
				fields = line.split("\t");
				y.add(Double.parseDouble(fields[0]));
			}
		}
		tfrY.close();
		
		System.out.println("result :\n" +  new KolmogorovSmirnovTest().compute(ListUtils.toArray(x), ListUtils.toArray(y)).toString()+"\n");	
		System.out.println("result :\n" +  KolmogorovSmirnovTest.kstest(ListUtils.toArray(x), ListUtils.toArray(y),"greater"));	
		System.out.println("result :\n" +  KolmogorovSmirnovTest.kstest(ListUtils.toArray(x), ListUtils.toArray(y),"less"));
		System.out.println("result :\n" +  KolmogorovSmirnovTest.kstest(ListUtils.toArray(x), ListUtils.toArray(y),"two-sided"));	
	}
	//	@Test
	public void TestKolmogorov() throws IOException {
		
		/*double[] d1 = {75.95259857177734, 10.248600006103516, 13.988200187683105, 1.3875499963760376, 0.9268500208854675, 1020.8099975585938, 3.442150115966797, 11.628800392150879, 24.642000198364258, 14.097299575805664, 17.189599990844727, 3.263200044631958, 0.24850499629974365, -0.34060800075531006, 0.3484489917755127, 2.0719199180603027, 7.16609001159668, 8.231969833374023, 36.24879837036133, 3.2494900226593018, 39.06919860839844, 19.201799392700195, 7.988969802856445, -0.0522489994764328, 4.554649829864502, 2.5156800746917725, 1.8990099430084229, 16.530899047851562, 9.296680450439453, 2.2606000900268555, 16.273399353027344, 6.935699939727783, 0.868369996547699, 2.863409996032715, 6.438650131225586, 6.1732001304626465};
		double[] d2 = {17.322099685668945, 4.429460048675537, 2.3720200061798096, 6.197969913482666, 7.515659809112549, 79.73970031738281, 23.520999908447266, 3.796639919281006, 31.923900604248047, 3.196429967880249, 2.2524900436401367, -0.5383650064468384, -1.014140009880066, 3.161329984664917, 5.246250152587891, 12.421500205993652, 40.577301025390625, 119.0530014038086, 14.748200416564941, 2.790839910507202, 1.9048099517822266, 7.025420188903809, 91.17900085449219, 9.022159576416016, 13.212699890136719, -0.9780979752540588, 5.8169097900390625, 1.3014099597930908, 77.50810241699219, 9.051329612731934, 4.987810134887695, 0.46449100971221924, 1.644760012626648, -2.1945300102233887, 15.038100242614746, 11.851099967956543, 319.0639953613281, 19.84709930419922, 8.98723030090332, 22.35810089111328, 75.95259857177734, -3.5373899936676025, 10.966500282287598, 6.023519992828369, 7.452789783477783, 37.92369842529297, 4.633339881896973, 11.421799659729004, 25.256099700927734};
		*/
		BufferedWriter bwd3 = new BufferedWriter(new FileWriter("d3.txt"));
		BufferedWriter bwd2 = new BufferedWriter(new FileWriter("d2.txt"));
		
		StringBuilder d2String = new StringBuilder();
		double[] d2 = new double[5000];
		for (int i=0; i < d2.length; i++){
			double x = (double) (Math.random()*100.0);
			d2[i] = x;
			d2String.append(x+"\n");
		}
		StringBuilder d3String = new StringBuilder();
		double[] d3 = new double[5000];
		for(int i=0; i<5000; i++){
			double x = (double) (Math.random()*100.0);
			d3[i] = x;
			d3String.append(x+"\n");
		}
		bwd3.write(d3String.toString());
		bwd2.write(d2String.toString());
		bwd3.close();
		bwd2.close();
//		System.out.println("result :\n" +  new KolmogorovSmirnovTest().compute(d1, d2).toString()+"\n");	
//		System.out.println("result :\n" +  KolmogorovSmirnovTest.kstest(d1, d2,"greater"));	
//		System.out.println("result :\n" +  KolmogorovSmirnovTest.kstest(d1, d2,"less"));
//		System.out.println("result :\n" +  KolmogorovSmirnovTest.kstest(d1, d2,"two-sided"));	
		
		System.out.println("result :\n" +  new KolmogorovSmirnovTest().compute(d3, d2).toString()+"\n");	
		System.out.println("result :\n" +  KolmogorovSmirnovTest.kstest(d3, d2,"greater"));	
		System.out.println("result :\n" +  KolmogorovSmirnovTest.kstest(d3, d2,"less"));
		System.out.println("result :\n" +  KolmogorovSmirnovTest.kstest(d3, d2,"two-sided"));	
		
	
	}
//	@Test
	public void Test() {
		
//		#1:If you have sif-file and you want to create a .topo file from the full sif-file (obviously no randoms)
//		String []args = {"--tool", "snow2", "-o", "/tmp/",  "--sif-file", "/home/ralonso/appl/babelomics/hsa/proteins/hsa_alldb_proteins_interactome_nr.sif", "--o-file-topo", "/home/ralonso/appl/babelomics/hsa/proteins/hsa_alldb_proteins_interactome_nr.topo", "--home", System.getenv("BABELOMICS_HOME")};

// 		#2:If you have sif-file and you want to create a .topo and/or .mcn file random list without a given size
//		String []args = {"--tool", "snow2", "-o", "/tmp/",  "--sif-file", "/home/ralonso/appl/bioinfo-networks/mine.sif", "--file-topo", "/home/ralonso/appl/bioinfo-networks/mineTopo", "--file-mcn", "/home/ralonso/appl/bioinfo-networks/mineMcn", "--randoms", "2","--home", System.getenv("BABELOMICS_HOME")};

//		#3:If you have sif-file and you want to create a .topo and/or .mcn file random list with a given size
//		String []args = {"--tool", "snow2", "-o", "/tmp/",  "--sif-file", "/home/ralonso/appl/bioinfo-networks/hsa_alldb_proteins_interactome_nr.sif", "--o-file-mcn", "/home/ralonso/appl/bioinfo-networks/mineMcn", "--randoms", "1000", "--random-size", "3","--home", System.getenv("BABELOMICS_HOME")};

// 		#4:If you have sif-file and you want to create a .topo and/or .mcn file from a given node-list
//		String []args = {"--tool", "snow2", "-o", "/tmp/", "--sif-file", "/home/ralonso/appl/bioinfo-networks/mine.sif", "--node-list", "v1,v2,v0",  "--file-topo", "/home/ralonso/appl/bioinfo-networks/mineTopo", "--file-mcn", "/home/ralonso/appl/bioinfo-networks/mineMcn", "--home", System.getenv("BABELOMICS_HOME")};

//		#5:If you have file-topo-values and you want to create a .mcn file from a given node-list but with the values of the file-topo-values
//		String []args = {"--tool", "snow2", "-o", "/tmp/file-topo-values", "--file-topo-values", "/home/ralonso/appl/bioinfo-networks/mineTopo.topo", "--node-list", "v1,v2,v0", "--file-mcn", "/home/ralonso/appl/bioinfo-networks/mineMcn", "--home", System.getenv("BABELOMICS_HOME")};

//		Case: You get the values for a network, but with the values of the subnet
//		String []args = {"--tool", "snow2", "-o",".", "--sif-file", "sce_alldb_proteins_interactome_nr.sif", "--node-file", "UPYDOWN_HET_list_uniq", "--file-mcn", "fileMcn", "--file-topo", "fileTopo","--home", System.getenv("BABELOMICS_HOME")};
		
//		String []args = {"--tool", "snow2", "-o","/opt/babelomics/", "--sif-file", "/opt/babelomics/sce_alldb_proteins_interactome_nr.sif", "--randoms", "2", "--random-size", "50", "--home", System.getenv("BABELOMICS_HOME")};

//		Case: We want to do an statistic Kolmogorov test 
		String []args = {"--tool", "snow2", "-o","/home/ralonso/appl/babelomics/", "--sif-file", "/home/ralonso/appl/babelomics/hsa/proteins/hsa_alldb_proteins_interactome_nr.sif", "--file-topo-values", "/home/ralonso/appl/babelomics/hsa/proteins/protein_interactome_params.hsa_alldb_proteins_interactome_nr", "--node-file1", "/home/ralonso/appl/babelomics/hsa/proteins/Ejemplo5/crs1_block559_uniprot.txt", "--randoms", "1", "--random-size", "1","--home", System.getenv("BABELOMICS_HOME")};

		try {
			for(int i=0; i < args.length; i++)
				System.out.println(args[i]);
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}
	}
}
