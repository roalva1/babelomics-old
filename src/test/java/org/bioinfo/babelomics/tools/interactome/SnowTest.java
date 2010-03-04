package org.bioinfo.babelomics.tools.interactome;


import static org.junit.Assert.fail;

import java.io.File;
import java.io.IOException;
import java.util.Arrays;

import org.bioinfo.babelomics.BabelomicsMain;
import org.bioinfo.commons.Config;
import org.bioinfo.commons.io.utils.FileUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.tool.OptionFactory;
import org.junit.After;
import org.junit.Before;
import org.junit.Test;

public class SnowTest {

	@Before
	public void setUp() throws Exception {
	}

	@After
	public void tearDown() throws Exception {
	}

	//	options.addOption(OptionFactory.createOption("list1", "the list1"));		
	//	options.addOption(OptionFactory.createOption("list2", "the list2", false));		
	//	options.addOption(OptionFactory.createOption("interactome", "Select interactome: join (for human interactome with all ppis), intersect (for human interacionts with only ppis detected by two methods) and own (for your own interactions)", false));
	//	options.addOption(OptionFactory.createOption("own-interactions", "submit a file with interactome ",false));
	//	options.addOption(OptionFactory.createOption("check-interactions", "Set this option if proteins in interactions and list are in same id",false,false));
	//	options.addOption(OptionFactory.createOption("id-nature", "Nature of your lists: genes or proteins", false));
	//	options.addOption(OptionFactory.createOption("interactions-number", "Maximum number of external proteins introduced: 0, 1, 2 or 3", false));


	public void Test1() {
		String list1 = "/mnt/commons/test/tools/snow/brca1_overexp_dn.txt";
		String outdir = "/tmp/SnowTest1";
		new File(outdir).mkdir();

		String []args = { "--tool", "snow","--log-level", "2", "--list1", list1, "-o", outdir};

		System.out.println("executing ----------------> " + Arrays.toString(args));
		try {
			FileUtils.createDirectory(outdir);
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.toString());
		}

	}

	@Test
	public void Test2() {
		String title;
		Config config;
		try {
			config = new Config("/mnt/commons/test/tools/snow/output.properties");

			String key;
			String keyStr = "INTERACTOME.BETWEENNESS,INTERACTOME.COEFFICIENT,INTERACTOME.COEFFICIENT,NETWORK.BETWEENNESS,NETWORK.COEFFICIENT,NETWORK.COEFFICIENT";
			for(String item: StringUtils.toList(keyStr, ",")) {
				key = item + ".PVALUE";
				System.out.println(key + " = " + config.getProperty(key));

				key = item + ".SIDE";
				System.out.println(key + " = " + config.getProperty(key));

				System.out.println("");
			}

		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}


}
