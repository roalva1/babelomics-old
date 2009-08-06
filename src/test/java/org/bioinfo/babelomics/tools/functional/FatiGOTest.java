package org.bioinfo.babelomics.tools.functional;


import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class FatiGOTest {

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}
	
	@Test
	public void Test0() {
		
	}

	public void Test1() {
		////
		String []args = {"-list1", "/mnt/commons/test/biodata/example/list3.txt", "-o", "/tmp", "-rest-of-genome", "1", "--go-bp-db", "true", "--go-bp-min-level", "5", "--go-bp-max-level", "12", "--remove-duplicates", "ref"};
		//String []args = {"-list1", "/mnt/commons/test/biodata/example/list3.txt", "-o", "/tmp", "-list2", "/mnt/commons/test/biodata/example/list4.txt", "--go-bp-db", "true", "--go-bp-min-level", "5", "--go-bp-max-level", "12", "--remove-duplicates", "ref"};
		//String []args = {"-list1", "/mnt/commons/test/biodata/example/list3.txt", "-o", "/tmp", "-list2", "/mnt/commons/test/biodata/example/list4.txt", "--go-bp-db", "true", "--go-bp-min-level", "5", "--go-bp-max-level", "5"};
		try {
			FatiGO fatigo = new FatiGO(args);
			fatigo.execute();
		} catch (Exception e) {
			e.printStackTrace();
			//System.out.println(e.toString());
		}
	}

}
