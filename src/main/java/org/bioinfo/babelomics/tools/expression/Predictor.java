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
		super(args);
		initOptions();
	}

	@Override
	public void initOptions() {
		options.addOption(OptionFactory.createOption("dataset", "the data"));
		options.addOption(OptionFactory.createOption("algorithms", "the algorithms separated by commas, possible values: svn, knn, dlda, pam, som"));
		options.addOption(OptionFactory.createOption("class", "class variable"));
		options.addOption(OptionFactory.createOption("trainning-size", "number of genes to use in trainning separated by commas, default:2,5,10,20,35,50", false));
		options.addOption(OptionFactory.createOption("gene-selection", "the gene selection, valid values: f-ratio, wilcoxon", false));
		options.addOption(OptionFactory.createOption("sample-filter", "class variable", false));
		options.addOption(OptionFactory.createOption("feature-filter", "class variable", false));

	}

	@Override
	public void execute() {
		try {
			CommandLine cmd = parse(args);

			Dataset dataset = new Dataset(new File(cmd.getOptionValue("dataset")));
			String className = cmd.getOptionValue("class");
			String algorithms = cmd.getOptionValue("algorithms");

			String nbOfGenes = cmd.getOptionValue("trainning-size", "2,5,10,20,35,50");
			String geneSelection = cmd.getOptionValue("gene-selection", "none");
				
			if(cmd.hasOption("sample-filter") || cmd.hasOption("feature-filter")) {
				dataset = dataset.getSubDataset(cmd.getOptionValue("sample-filter"), "4", cmd.getOptionValue("feature-filter"), ""); 
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
		} catch (ParseException e) {
			logger.error("Error parsing command line", e.toString());
			System.out.println("\n");
			printUsage();
		} catch (IOException e) {
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
