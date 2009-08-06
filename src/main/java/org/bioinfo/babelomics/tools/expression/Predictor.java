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
		getOptions().addOption(OptionFactory.createOption("dataset", "the data"));
		getOptions().addOption(OptionFactory.createOption("algorithms", "the algorithms separated by commas, possible values: svn, knn, dlda, pam, som"));
		getOptions().addOption(OptionFactory.createOption("class", "class variable"));
		getOptions().addOption(OptionFactory.createOption("trainning-size", "number of genes to use in trainning separated by commas, default:2,5,10,20,35,50", false));
		getOptions().addOption(OptionFactory.createOption("gene-selection", "the gene selection, valid values: f-ratio, wilcoxon", false));
		getOptions().addOption(OptionFactory.createOption("sample-filter", "class variable", false));
		getOptions().addOption(OptionFactory.createOption("feature-filter", "class variable", false));

	}

	@Override
	public void execute() {
		try {
			Dataset dataset = new Dataset(new File(commandLine.getOptionValue("dataset")));
			String className = commandLine.getOptionValue("class");
			String algorithms = commandLine.getOptionValue("algorithms");

			String nbOfGenes = commandLine.getOptionValue("trainning-size", "2,5,10,20,35,50");
			String geneSelection = commandLine.getOptionValue("gene-selection", "none");
				
			if(commandLine.hasOption("sample-filter") || commandLine.hasOption("feature-filter")) {
				dataset = dataset.getSubDataset(commandLine.getOptionValue("sample-filter"), "4", commandLine.getOptionValue("feature-filter"), ""); 
			}
			System.out.println(dataset.toString()+"\n");
			

			String[] algs = algorithms.split(",");
			for(String algorithm: algs) {
				if (algorithm.equals("svm")) {
					executeSvm(dataset, className, nbOfGenes, geneSelection);
				}
				if (algorithm.equals("knn")) {
					executeKnn(dataset, className, nbOfGenes, geneSelection);
				}
				if (algorithm.equals("dlda")) {
					executeDlda(dataset, className, nbOfGenes, geneSelection);
				}
				if (algorithm.equals("pam")) {
					executePam(dataset, className, nbOfGenes, geneSelection);
				}				
				if (algorithm.equals("som")) {
					executeSom(dataset, className, nbOfGenes, geneSelection);
				}				
			}
			// create global graph
			

			logger.warn("que raroo....");
		} catch (Exception e) {
			logger.error("Error opening the dataset", e.toString());
		} 
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
