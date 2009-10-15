package org.bioinfo.babelomics.tools.expression.differential;

import org.bioinfo.babelomics.BabelomicsMain;
import org.junit.Test;


public class TimeSeriesTest {
	
	@Test
	public void test() {
		System.out.println("----- masigpro ------");
		String []args = {"--tool", "time-series", "--dataset", "/mnt/commons/test/biodata/example/masigpro.dataset", "-o", "/tmp/masigpro", "--test", "masigpro", "--contin-class", "contin", "--series-class", "series"};
		
		try {
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
		}		
	}
}
