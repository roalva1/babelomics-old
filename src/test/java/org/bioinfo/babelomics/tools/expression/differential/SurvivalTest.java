package org.bioinfo.babelomics.tools.expression.differential;

import org.bioinfo.babelomics.BabelomicsMain;
import org.junit.Test;


public class SurvivalTest {
	
	@Test
	public void test() {
		System.out.println("----- cox ------");
		String []args = {"-tool", "survival", "-dataset", "/mnt/commons/test/biodata/example/survival.txt", "-o", "/tmp/cox", "-test", "cox", "-time-class", "time", "-censored-class", "censored"};
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}

}
