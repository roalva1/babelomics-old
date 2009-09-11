package org.bioinfo.babelomics.tools.functional;


import org.bioinfo.babelomics.BabelomicsMain;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class MarmiteTest {

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}
	
	@Test
	public void test0() {
		
	}
	
	@Test
	public void test1() {
		System.out.println("----- marmite test  ------");
//		String []args = {"-tool", "marmite", "-list1", "/mnt/commons/test/tools/marmite/marmite1.txt",  "-list2", "/mnt/commons/test/tools/marmite/marmite2.txt", "-o", "/mnt/commons/test/tools/marmite/out/", "-bioentity-name", "diseases", "-bioentity-score-filter", "5"};
		String []args = {"-tool", "marmite", "-list1", "/mnt/commons/test/tools/marmite/marmite1a.txt",  "-list2", "/mnt/commons/test/tools/marmite/marmite2a.txt", "-o", "/mnt/commons/test/tools/marmite/out/", "-bioentity-name", "diseases", "-bioentity-score-filter", "5"};
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}

}
