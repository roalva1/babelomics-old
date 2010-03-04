package org.bioinfo.babelomics.tools.expression.normalization;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.bioinfo.commons.io.utils.IOUtils;

public class NormalizationUtils {

	public static void createPreprocessingRedirectionFile(File redirectionFile, File fileToRedirect) {
		List<String> redirectionInputs = new ArrayList<String>();
		redirectionInputs.add("tool=preprocessing");
		redirectionInputs.add("jobname=preprocessing");
		redirectionInputs.add("job_description=redirected from job $JOB_NAME");
		redirectionInputs.add("dataset=" + fileToRedirect.getName() + "(from job $JOB_NAME)::::$JOB_FOLDER/" + fileToRedirect.getName());
		redirectionInputs.add("merge-replicates=mean");
		redirectionInputs.add("impute-missing=mean");
		try {
			IOUtils.write(redirectionFile.getAbsolutePath(), redirectionInputs);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}

	}
}
