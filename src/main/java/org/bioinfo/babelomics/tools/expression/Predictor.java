package org.bioinfo.babelomics.tools.expression;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.apache.commons.cli.OptionGroup;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.mlpr.classifier.Knn;
import org.bioinfo.mlpr.classifier.RForest;
import org.bioinfo.mlpr.classifier.Svm;
import org.bioinfo.mlpr.evaluation.result.EvaluationResult;
import org.bioinfo.mlpr.utils.InstancesBuilder;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;

import weka.core.Attribute;
import weka.core.Instances;

public class Predictor extends BabelomicsTool {

	private double progressStep;
	private double progressCurrent = 0;
	
	public Predictor() {
		initOptions();
	}

	@Override
	public void initOptions() {
		// data
		//
		options.addOption(OptionFactory.createOption("dataset", "the data"));
		options.addOption(OptionFactory.createOption("dataset-arff", "dataset is in arff format",false,false));
		
		// class attribute
		//
		options.addOption(OptionFactory.createOption("class", "corresponding class attribute in the dataset"));
		options.addOption(OptionFactory.createOption("class-file", "class variable file",false));
		
		// sample and feature filters
		//
		options.addOption(OptionFactory.createOption("sample-filter", "Sample filter", false));
		options.addOption(OptionFactory.createOption("feature-filter", "Feature filter", false));

		// classifiers algorithms
		//
		OptionGroup classifiers = new OptionGroup();
		classifiers.setRequired(false);

		// KNN
		classifiers.addOption(OptionFactory.createOption("knn", "Classify dataset with a KNN classifier",false,false));
		options.addOption(OptionFactory.createOption("knn-tune", "Perform automated number of neighbors tunning",false,false));
		options.addOption(OptionFactory.createOption("knn-neighbors", "Knn number of neighbors (5 by default)", false, true));
		
		// SVM
		classifiers.addOption(OptionFactory.createOption("svm", "Classify dataset with a SVM classifier", false,false));
		options.addOption(OptionFactory.createOption("svm-tune", "Perform automated parameter tunning",false,false));
		options.addOption(OptionFactory.createOption("svm-cost", "----", false));
		
		// Random forest
		classifiers.addOption(OptionFactory.createOption("random-forest", "Classify dataset with a Random Forest classifier", false,false));
		options.addOption(OptionFactory.createOption("random-forest-tune", "Perform automated number of tree tunning",false,false));
		options.addOption(OptionFactory.createOption("random-forest-trees", "Number of trees", false));
		
		// DLDA
		classifiers.addOption(OptionFactory.createOption("dlda", "Classify dataset with a DLDA classifier", false,false));
		options.addOption(OptionFactory.createOption("dlda-tune", "Perform automated parameter tunning",false,false));
		
		// SOM
		classifiers.addOption(OptionFactory.createOption("som", "Classify dataset with a SOM classifier", false,false));
		options.addOption(OptionFactory.createOption("som-tune", "Perform automated parameter tunning",false,false));
		
		// PAM
		classifiers.addOption(OptionFactory.createOption("pam", "Classify dataset with a PAM classifier", false,false));
		options.addOption(OptionFactory.createOption("pam-tune", "Perform automated parameter tunning",false,false));		
		options.addOptionGroup(classifiers);

		// feature selection (gene selection), and other options
		//
		options.addOption(OptionFactory.createOption("gene-selection", "the gene selection, valid values: f-ratio, wilcoxon", false));
		//options.addOption(OptionFactory.createOption("trainning-size", "number of genes to use in trainning separated by commas, default:2,5,10,20,35,50", false));		

	}

	@Override
	public void execute() {
		initProgress();
		logger.info("Welcome to prophet...");

		try {

			Instances instances = null;
			File datasetFile = new File(commandLine.getOptionValue("dataset"));
			String className = commandLine.getOptionValue("class");

			// data file checking
			//
			if(datasetFile.exists()) {

				try {
					jobStatus.addStatusMessage(StringUtils.decimalFormat(progressCurrent), "reading dataset");
					progressCurrent += progressStep;
				} catch (FileNotFoundException e) {
					abort("filenotfoundexception_execute_predictor", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
				}

				// data set loading 
				if(commandLine.hasOption("dataset-arff")){
					instances = InstancesBuilder.getInstancesFromArrfFile(datasetFile,"sample_name");
				} else {
					Dataset dataset = new Dataset(datasetFile, true);
					List<List<String>> data = new ArrayList<List<String>>(dataset.getColumnDimension());
					for(int i = 0 ; i<dataset.getColumnDimension() ; i++) {
						//						System.out.println("row -> " + Arrays.toString(dataset.getDoubleMatrix().getRow(i)));
						//						System.out.println("list -> " + ListUtils.toString(ListUtils.toStringList(dataset.getDoubleMatrix().getRow(i)), "\t"));
						data.add(ListUtils.toStringList(dataset.getDoubleMatrix().getColumn(i)));
					}
					System.out.println("0,0 = " + data.get(0).get(0));
					System.out.println("variables value = " + ListUtils.toString(dataset.getVariables().getByName(className).getValues()));

					instances = InstancesBuilder.getInstancesFromDataList(data, dataset.getFeatureNames(), dataset.getSampleNames(), dataset.getVariables().getByName(className).getValues());
				}				
			} else {
				abort("execute_predictor", "dataset " + datasetFile.getName() + "not found", "dataset " + datasetFile.getName() + "not found", "dataset " + datasetFile.getName() + "not found");
			}

			// class definition 
			//
			Attribute classAttr = instances.attribute(className);
			instances.setClassIndex(classAttr.index());					

			// classifiers
			//
			if(commandLine.hasOption("svm")) executeSvm(instances);
			if(commandLine.hasOption("knn")) executeKnn(instances);
			if(commandLine.hasOption("random-forest")) executeRandomForest(instances);
			if(commandLine.hasOption("dlda")) executeDlda(instances);
			if(commandLine.hasOption("som")) executeSom(instances);
			if(commandLine.hasOption("pam")) executePam(instances);

		} catch (Exception e) {
			logger.error(e.getMessage());
			e.printStackTrace(logger.getLogWriter());
		}

	}


	private void executeKnn(Instances instances) {
		try {
			jobStatus.addStatusMessage(StringUtils.decimalFormat(progressCurrent), "executing KNN classifier");
			progressCurrent += progressStep;
		} catch (FileNotFoundException e) {
			abort("filenotfoundexception_executeknn_predictor", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
		}
		
		// init classifier
		Knn knn = new Knn();
		// init params
		if(commandLine.hasOption("knn-tune")) knn.setTuneParameters(true);
		else {
			if(commandLine.hasOption("knn-neighbors")) knn.setKnn(Integer.parseInt(commandLine.getOptionValue("knn-neighbors")));						
		}
		// train
		knn.train(instances);
		// results
		int best = knn.getEvaluationResultList().getBestRootMeanSquaredErrorIndex();
		EvaluationResult bestRMSE = knn.getEvaluationResultList().get(best);
		logger.println("Best RMSE classification");
		logger.println("Knn: " + knn.getKnnValues()[best] + " neighbors");
		logger.println(bestRMSE.toString());

		try {
			IOUtils.write(new File(outdir + "/knn.txt"), "Best RMSE classification\n\nKnn: " + knn.getKnnValues()[best] + " neighbors\n\n" + bestRMSE.toString());
			result.addOutputItem(new Item("knn_result_file", "knn.txt", "KNN result file", TYPE.FILE));
		} catch (IOException e) {
			printError("ioexception_executeknn_predictor", "Error saving KNN results", "Error saving KNN results");
		}
		
		
		// ??????
	}

	private void executeSvm(Instances instances) {
		try {
			jobStatus.addStatusMessage(StringUtils.decimalFormat(progressCurrent), "executing SVM classifier");
			progressCurrent += progressStep;
		} catch (FileNotFoundException e) {
			abort("filenotfoundexception_executesvm_predictor", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
		}
		
		// init classifier
		Svm svm = new Svm();
		// init params
		if(commandLine.hasOption("svm-tune")) svm.setTuneParameters(true);
		else {
			if(commandLine.hasOption("svm-cost")) svm.setCost(Integer.parseInt(commandLine.getOptionValue("svm-cost")));						
		}
		// train
		svm.train(instances);
		// results
		int best = svm.getEvaluationResultList().getBestRootMeanSquaredErrorIndex();
		EvaluationResult bestRMSE = svm.getEvaluationResultList().get(best);
		logger.println("Beset RMSE classification");
		logger.println("Cost: " + svm.getNumCostArray()[best]);
		logger.println(bestRMSE.toString());

		try {
			IOUtils.write(new File(outdir + "/svm.txt"), "Best RMSE classification\n\nCost: " + svm.getNumCostArray()[best] + "\n\n" + bestRMSE.toString());
			result.addOutputItem(new Item("svm_result_file", "svm.txt", "SVM result file", TYPE.FILE));
		} catch (IOException e) {
			printError("ioexception_executesvm_predictor", "Error saving SVM results", "Error saving SVM results");
		}
		// ??????
	}

	private void executeRandomForest(Instances instances) {
		try {
			jobStatus.addStatusMessage(StringUtils.decimalFormat(progressCurrent), "executing Random Forest classifier");
			progressCurrent += progressStep;
		} catch (FileNotFoundException e) {
			abort("filenotfoundexception_executerandomforest_predictor", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
		}

		// init classifier
		RForest randomForest = new RForest();
		// init params
		if(commandLine.hasOption("random-forest-tune")) randomForest.setTuneParameters(true);
		else {
			if(commandLine.hasOption("random-forest-trees")) randomForest.setNumTrees(Integer.parseInt(commandLine.getOptionValue("random-forest-trees")));						
		}
		// train
		randomForest.train(instances);
		// results
		int best = randomForest.getEvaluationResultList().getBestRootMeanSquaredErrorIndex();
		EvaluationResult bestRMSE = randomForest.getEvaluationResultList().get(best);
		logger.println("Best RMSE classification");
		logger.println("Number of trees: " + randomForest.getNumTreesArray()[best]);
		logger.println(bestRMSE.toString());

		try {
			IOUtils.write(new File(outdir + "/random_forest.txt"), "Best RMSE classification\n\nNumber of trees: " + randomForest.getNumTreesArray()[best] + "\n\n" + bestRMSE.toString());
			result.addOutputItem(new Item("random_forest_result_file", "random_forest.txt", "Random Forest result file", TYPE.FILE));
		} catch (IOException e) {
			printError("ioexception_executerandomforest_predictor", "Error saving Random Forest results", "Error saving Random Forest results");
		}
		// ??????
	}

	private void executeDlda(Instances instances) {		
		printError("executeDlda_predictor", "DLDA is not implemented yet", "DLDA is not implemented yet");
	}

	private void executeSom(Instances instances) {
		printError("executeDlda_predictor", "SOM is not implemented yet", "DLDA is not implemented yet");
	}

	private void executePam(Instances instances) {
		printError("executeDlda_predictor", "PAM is not implemented yet", "DLDA is not implemented yet");
	}	
	
	private void initProgress() {
		int steps = 2; // it includes reading dataset and done
		
		if(commandLine.hasOption("svm")) steps++;
		if(commandLine.hasOption("knn")) steps++;
		if(commandLine.hasOption("random-forest")) steps++;
		if(commandLine.hasOption("dlda")) steps++;
		if(commandLine.hasOption("som")) steps++;
		if(commandLine.hasOption("pam")) steps++;
		
		progressStep = 100.0 / steps;
		progressCurrent = progressStep;
	}
}
