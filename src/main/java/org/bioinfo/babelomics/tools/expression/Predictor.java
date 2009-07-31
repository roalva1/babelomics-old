package org.bioinfo.babelomics.tools.expression;

import java.io.File;
import java.io.IOException;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.ParseException;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.tool.OptionFactory;

public class Predictor extends BabelomicsTool {


	public Predictor(String[] args) {
		initOptions();
	}

	@Override
	public void initOptions() {
		
		// data 
		options.addOption(OptionFactory.createOption("dataset", "the data"));
		options.addOption(OptionFactory.createOption("class", "class variable"));
		options.addOption(OptionFactory.createOption("sample-filter", "class variable", false));
		options.addOption(OptionFactory.createOption("feature-filter", "class variable", false));
		
		// feature selection (gene selection)
		options.addOption(OptionFactory.createOption("gene-selection", "the gene selection, valid values: f-ratio, wilcoxon", false));
		//options.addOption(OptionFactory.createOption("trainning-size", "number of genes to use in trainning separated by commas, default:2,5,10,20,35,50", false));		

		// classifiers
		  // knn
		options.addOption(OptionFactory.createOption("knn", "Knn classifier"));
		options.addOption(OptionFactory.createOption("knn-neighbors", " number of neighbors (5 by default) or auto for automated tunning", true));
		  // svn
		options.addOption(OptionFactory.createOption("svn", "Support vector machine"));

		

	}

	@Override
	public void execute() {
		logger.info("hasta aqui llegó mi ejecución");
	}

	private void executeSvm(Dataset dataset, String className, String nbOfGenes, String geneSelection) {
		logger.info("executing svm, not implemented yet");
	}

	private void executeKnn(Dataset dataset, String className, String nbOfGenes, String geneSelection) {
		logger.info("executing knn, not implemented yet");
	}
	
	private void executeDlda(Dataset dataset, String className, String nbOfGenes, String geneSelection) {
		logger.info("executing dlda, not implemented yet");
	}

	private void executePam(Dataset dataset, String className, String nbOfGenes, String geneSelection) {
		logger.info("executing pam, not implemented yet");
	}

	private void executeSom(Dataset dataset, String className, String nbOfGenes, String geneSelection) {
		logger.info("executing som, not implemented yet");
	}
}
