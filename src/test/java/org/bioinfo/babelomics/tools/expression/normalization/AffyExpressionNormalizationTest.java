package org.bioinfo.babelomics.tools.expression.normalization;

import static org.junit.Assert.fail;

import java.util.Arrays;

import org.bioinfo.babelomics.BabelomicsMain;
import org.bioinfo.commons.io.utils.FileUtils;
import org.junit.Test;


public class AffyExpressionNormalizationTest {
	@Test
	public void Test() {
		String dataset = "/mnt/commons/test/biodata/example/CEL.tar.gz";
		//String dataset = "/mnt/commons/test/biodata/example/cel.zip";
		String outdir = "/tmp/affy-normalization";
		String []args = { "--tool", "affy-expression-normalization","--log-level", "2", "--compressed-file", dataset, "-o", outdir, "--rma", "--plier", "--calls"};

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
