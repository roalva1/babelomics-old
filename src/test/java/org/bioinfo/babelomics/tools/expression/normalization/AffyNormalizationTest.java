package org.bioinfo.babelomics.tools.expression.normalization;

import static org.junit.Assert.fail;

import java.io.File;
import java.util.Arrays;

import org.bioinfo.babelomics.BabelomicsMain;
import org.bioinfo.commons.io.utils.FileUtils;
import org.bioinfo.commons.io.utils.IOUtils;
import org.junit.Test;


public class AffyNormalizationTest {

	@Test
	public void Test() {
		String dataset = "/mnt/commons/test/biodata/example/CEL.tar.gz";
		String outdir = "/tmp/affy-normalization";
		String []args = { "--tool", "affy-normalization","--log-level", "2", "--compressed-file", dataset, "-o", outdir, "--rma"};

		System.out.println("----------------> " + Arrays.toString(args));
		try {
			FileUtils.createDirectory(outdir);
			BabelomicsMain.main(args); 
		} catch (Exception e) {
			e.printStackTrace();
			fail(e.toString());
		}
	}	

}
