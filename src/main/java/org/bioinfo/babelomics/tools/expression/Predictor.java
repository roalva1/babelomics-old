package org.bioinfo.babelomics.tools.expression;

import java.io.File;

import org.apache.commons.cli.OptionGroup;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.mlpr.classifier.Knn;
import org.bioinfo.mlpr.classifier.RForest;
import org.bioinfo.mlpr.classifier.Svm;
import org.bioinfo.mlpr.evaluation.result.EvaluationResult;
import org.bioinfo.mlpr.utils.InstancesBuilder;
import org.bioinfo.tool.OptionFactory;

import weka.core.Attribute;
import weka.core.Instances;

public class Predictor extends BabelomicsTool {


	public Predictor() {
		initOptions();
	}

	@Override
	public void initOptions() {
		// data 
		  // dataset
		options.addOption(OptionFactory.createOption("dataset", "the data"));
		options.addOption(OptionFactory.createOption("dataset-arff", "dataset is in arff format",false,false));
		  // class attribute
		options.addOption(OptionFactory.createOption("class", "corresponding class attribute in the dataset",false));
		options.addOption(OptionFactory.createOption("class-file", "class variable file",false));
		  // sample filter
		options.addOption(OptionFactory.createOption("sample-filter", "Sample filter", false));
		  // feature filter
		options.addOption(OptionFactory.createOption("feature-filter", "Feature filter", false));
		
		// feature selection (gene selection)
		options.addOption(OptionFactory.createOption("gene-selection", "the gene selection, valid values: f-ratio, wilcoxon", false));
		//options.addOption(OptionFactory.createOption("trainning-size", "number of genes to use in trainning separated by commas, default:2,5,10,20,35,50", false));		

		// classifiers
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
		

	}

	@Override
	public void execute() {
				
		logger.info("Welcome to prophet...");
				
		try {
			
			// data file checking
			File arff = new File(commandLine.getOptionValue("dataset"));
			if(arff.exists()) {
				
				// data set loading 
				Instances instances;
				if(commandLine.hasOption("dataset-arff")){
					instances = InstancesBuilder.getInstancesFromArrfFile(arff,"sample_name");
				} else {
					throw new Exception("Babelomics dataset reading not yet implemented");
				}
				
				// class definition 
				if(commandLine.hasOption("class")){
					Attribute classAttr = instances.attribute(commandLine.getOptionValue("class"));
					if(classAttr==null) throw new Exception("class attribute not found in dataset");
					else instances.setClassIndex(classAttr.index());					
				} else {
					throw new Exception("Class file reading not yet implemented");					
				}
									
				// classifiers
				if(commandLine.hasOption("svm")) executeSvm(instances);
				if(commandLine.hasOption("knn")) executeKnn(instances);
				if(commandLine.hasOption("random-forest")) executeRandomForest(instances);
				if(commandLine.hasOption("dlda")) executeDlda(instances);
				if(commandLine.hasOption("som")) executeSom(instances);
				if(commandLine.hasOption("pam")) executePam(instances);
				
				
			} else {
				logger.error("dataset " + commandLine.getOptionValue("dataset") + "not found");
			}
			
			
		} catch (Exception e) {
			logger.error(e.getMessage());
			e.printStackTrace(logger.getLogWriter());
		}
		
	}
		
	
	private void executeKnn(Instances instances) {
		logger.println("");
		logger.println("KNN classifier..................................................................");
		logger.println("");
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
		logger.println("Beset RMSE classification");
		logger.println("Knn: " + knn.getKnnValues()[best] + " neighbors");
		logger.println(bestRMSE.toString());
		// ??????
	}
	
	private void executeSvm(Instances instances) {
		logger.println("");
		logger.println("SVM classifier...................................................................");
		logger.println("");
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
		// ??????
	}
	
	private void executeRandomForest(Instances instances) {
		logger.println("");
		logger.println("Random Forest classifier.........................................................");
		logger.println("");
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
		logger.println("Beset RMSE classification");
		logger.println("Number of trees: " + randomForest.getNumTreesArray()[best]);
		logger.println(bestRMSE.toString());
		// ??????
	}
	
	private void executeDlda(Instances instances) {		
		logger.error("DLDA is not implemented yet");
	}
	
	private void executeSom(Instances instances) {
		logger.error("SOM is not implemented yet (and perhaps never!!)");
	}
	
	private void executePam(Instances instances) {
		logger.error("PAM is not implemented yet (and perhaps never!!)");
	}
	
}
