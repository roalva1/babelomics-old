package org.bioinfo.babelomics.methods.genomic.copynumber;

import java.util.ArrayList;
import java.util.List;

import org.bioinfo.commons.exec.Command;
import org.bioinfo.commons.exec.SingleProcess;
import org.bioinfo.commons.utils.ListUtils;

public class CopyNumberUtils {

	public static void normalization(String binPath, String bgCorrect, String waNormalization, String baNormalization, String design, String inDirname, String outDirname) {
	
		List<String> env = new ArrayList<String>();
		env.add("datadir=" + inDirname);
		env.add("source=agilent");
		env.add("methodbackcorrect=" + bgCorrect);
		env.add("methodnormWA=" +  waNormalization);
		env.add("methodnormBA=" + baNormalization);
		env.add("design=" + design);
		env.add("outfile=cgh_normalization.txt");
		env.add("outpositions=cgh_positions.txt");
		env.add("outdir=" + outDirname);

		Command cmd = new Command("R CMD BATCH --no-save --no-restore " + binPath + " " + outDirname + "/cgh.normalization.r.log", env);
		SingleProcess sp = new SingleProcess(cmd);
		sp.runSync();
		
		System.out.println("cmd line = " + cmd.getCommandLine());
		System.out.println("cmd env = " + ListUtils.toString(env, " "));		
	}
}
