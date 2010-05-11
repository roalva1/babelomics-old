package org.bioinfo.babelomics.tools.expression;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.HashMap;
import java.util.List;

import org.apache.commons.cli.OptionGroup;
import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.chart.XYLineChart;
import org.bioinfo.commons.io.utils.FileUtils;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ArrayUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.mlpr.classifier.AbstractFeatureSelector;
import org.bioinfo.mlpr.classifier.CfsFeatureSelector;
import org.bioinfo.mlpr.classifier.GenericClassifier;
import org.bioinfo.mlpr.classifier.Knn;
import org.bioinfo.mlpr.classifier.RForest;
import org.bioinfo.mlpr.classifier.Svm;
import org.bioinfo.mlpr.classifier.result.ClassificationResult;
import org.bioinfo.mlpr.evaluation.KFoldCrossValidation;
import org.bioinfo.mlpr.evaluation.result.EvaluationResult;
import org.bioinfo.mlpr.utils.InstancesBuilder;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.exception.InvalidParemeterException;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;

import weka.core.Attribute;
import weka.core.Instances;

public class Predictor extends BabelomicsTool {

	private double progressStep;
	private double progressCurrent = 0;
	private StringBuilder combinedTable;
	private int numberOfBestSelectedClassifications = 5;
	
	public Predictor() {		
		initOptions();
	}


	@Override
	public void initOptions() {
		// data
		options.addOption(OptionFactory.createOption("dataset", "the data"));
		options.addOption(OptionFactory.createOption("dataset-arff", "dataset is in arff format",false,false));
		
		// data test
		options.addOption(OptionFactory.createOption("test", "the data test",false));
		options.addOption(OptionFactory.createOption("test-arff", "dataset test is in arff format",false,false));
		
		// class attribute
		options.addOption(OptionFactory.createOption("class", "corresponding class attribute in the dataset",false,true));
		options.addOption(OptionFactory.createOption("class-file", "class variable file",false,true));
		
		// sample and feature filters
		options.addOption(OptionFactory.createOption("sample-filter", "Sample filter", false));
		options.addOption(OptionFactory.createOption("feature-filter", "Feature filter", false));

		// classifiers algorithms
		OptionGroup classifiers = new OptionGroup();
		classifiers.setRequired(false);
		
			// Evaluation
			options.addOption(OptionFactory.createOption("cross-validation", "Perform cross validation analysis",false,false));			
			options.addOption(OptionFactory.createOption("cross-validation-folds", "Number of folds in cross validation evaluation",false,true));
			options.addOption(OptionFactory.createOption("validation-repeats", "Number of repeat each randomization",false,true));
			options.addOption(OptionFactory.createOption("feature-selection", "Feature selection",false,true));			

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
			
//			// DLDA
//			classifiers.addOption(OptionFactory.createOption("dlda", "Classify dataset with a DLDA classifier", false,false));
//			options.addOption(OptionFactory.createOption("dlda-tune", "Perform automated parameter tunning",false,false));
//			
//			// SOM
//			classifiers.addOption(OptionFactory.createOption("som", "Classify dataset with a SOM classifier", false,false));
//			options.addOption(OptionFactory.createOption("som-tune", "Perform automated parameter tunning",false,false));
//			
//			// PAM
//			classifiers.addOption(OptionFactory.createOption("pam", "Classify dataset with a PAM classifier", false,false));
//			options.addOption(OptionFactory.createOption("pam-tune", "Perform automated parameter tunning",false,false));
			
		options.addOptionGroup(classifiers);

		// feature selection (gene selection), and other options
//		options.addOption(OptionFactory.createOption("gene-selection", "the gene selection, valid values: f-ratio, wilcoxon", false));
		//options.addOption(OptionFactory.createOption("trainning-size", "number of genes to use in trainning separated by commas, default:2,5,10,20,35,50", false));		

	}

	@Override
	public void execute() {
		
		logger.info("Welcome to prophet...");
		
		// init status
		combinedTable = new StringBuilder();
		initStatus();
		
		try {

			Instances instances = null;
			File datasetFile = new File(commandLine.getOptionValue("dataset"));

			// data file checking
			if(datasetFile.exists()) {

				// update status
				updateStatus(progressCurrent,"Reading dataset");
				
				// data set loading 
				if(commandLine.hasOption("dataset-arff")){
					instances = InstancesBuilder.getInstancesFromArrfFile(datasetFile,"sample_name");
					instances.setClassIndex(instances.numAttributes()-1);
				} else {
					// convert Dataset to Instances format (data is trasposed!!)
					Dataset dataset = new Dataset(datasetFile, true);
					List<List<String>> data = new ArrayList<List<String>>(dataset.getColumnDimension());
					for(int i = 0 ; i<dataset.getColumnDimension() ; i++) {
						data.add(ArrayUtils.toStringList(dataset.getDoubleMatrix().getColumn(i)));
					}
					// class values
					List<String> classValues = null;
					if(commandLine.hasOption("class-file")){
						String classFileName = commandLine.getOptionValue("class-file");
						FileUtils.checkFile(classFileName);
						classValues = StringUtils.toList(IOUtils.toString(classFileName));
					} else {
						if(commandLine.hasOption("class")){
							String className = commandLine.getOptionValue("class");
							if(dataset.getVariables()!=null && dataset.getVariables().getByName(className)!=null && dataset.getVariables().getByName(className).getValues()!=null){
								classValues = dataset.getVariables().getByName(className).getValues();	
							} else {
								throw new Exception("class not found in dataset");
							}
						} else {
							throw new Exception("class or class-file is not defined"); 
						}
					}
					System.err.println("class values: " + classValues);
					// create instances
					instances = InstancesBuilder.getInstancesFromDataList(data, dataset.getFeatureNames(), dataset.getSampleNames(), classValues);
					// define class attribute
					Attribute classAttr = instances.attribute("class");
					instances.setClassIndex(classAttr.index());		
				}
			} else {
				abort("execute_predictor", "dataset " + datasetFile.getName() + "not found", "dataset " + datasetFile.getName() + "not found", "dataset " + datasetFile.getName() + "not found");
			}
			 

			Instances test = null;
			
			if(commandLine.hasOption("test")){
				File testFile = new File(commandLine.getOptionValue("test"));				
				// data set loading 
				if(commandLine.hasOption("test-arff")){
					instances = InstancesBuilder.getInstancesFromArrfFile(datasetFile,"sample_name");				
				} else {
					Dataset dataset = new Dataset(testFile, true);
					List<List<String>> data = new ArrayList<List<String>>(dataset.getColumnDimension());
					for(int i = 0 ; i<dataset.getColumnDimension() ; i++) {
						data.add(ArrayUtils.toStringList(dataset.getDoubleMatrix().getColumn(i)));
					}				
					// create instances
					test = InstancesBuilder.getTestInstancesFromDataList(data, dataset.getFeatureNames(), dataset.getSampleNames());				
				}
			}
			
			
			// classifiers
			if(commandLine.hasOption("svm")) executeSvm(instances,test);
			if(commandLine.hasOption("knn")) executeKnn(instances,test);
			if(commandLine.hasOption("random-forest")) executeRandomForest(instances,test);
//			if(commandLine.hasOption("dlda")) executeDlda(instances);
//			if(commandLine.hasOption("som")) executeSom(instances);
//			if(commandLine.hasOption("pam")) executePam(instances);

			
		} catch (Exception e) {
			logger.error(e.getMessage());
			e.printStackTrace(logger.getLogWriter());
		}
		
		// classification table
		try {
			IOUtils.write(new File(outdir + "/combined_table.txt"), GenericClassifier.getResultsTableHeader() + "\n" +  combinedTable.toString());	
			result.getOutputItems().add(0, new Item("combined_table", "combined_table.txt", " Combined results (best " + numberOfBestSelectedClassifications + " per classifier)", TYPE.FILE,Arrays.asList("TABLE","PREDICTOR_COMBINED_TABLE"),new HashMap<String,String>(),"Summary",""));
		} catch (IOException e) {
			printError("ioexception_predictor", "Error saving combined table", "Error saving combined table");
		}
		
	}


	private void executeKnn(Instances instances, Instances test) throws Exception {
		
		// update status
		updateStatus(progressCurrent,"executing KNN classifier");
				
		// init classifier
		Knn knn = new Knn();
		
		// params
		if(commandLine.hasOption("knn-tune")) {
			knn.setTuneParameters(true);
		} else {
			if(commandLine.hasOption("knn-neighbors")) knn.setKnn(Integer.parseInt(commandLine.getOptionValue("knn-neighbors")));						
		}
		
		// validation
		setValidation(knn);
		
		// train
		knn.train(instances);
		
		// output results
		int best = knn.getEvaluationResultList().getBestRootMeanSquaredErrorIndex();
		EvaluationResult bestRMSE = knn.getEvaluationResultList().get(best);
		
		Knn testKnn = new Knn(best);
		List<ClassificationResult> classificationResult = knn.test(instances, test);
//		logger.println("Best RMSE classification");
//		logger.println("Knn: " + knn.getKnnValues()[best] + " neighbors");
//		logger.println(bestRMSE.toString());
		
		// save results
		saveClassifierResults(knn);
		
		
	}
	
	
	private void executeSvm(Instances instances, Instances test) throws Exception {
		
		// update status
		updateStatus(progressCurrent,"executing SVM classifier");
		
		// init classifier
		Svm svm = new Svm();
		// init params
		if(commandLine.hasOption("svm-tune")) svm.setTuneParameters(true);
		else {
			if(commandLine.hasOption("svm-cost")) svm.setCost(Integer.parseInt(commandLine.getOptionValue("svm-cost")));						
		}
		
		// validation
		setValidation(svm);
				
		// train
		svm.train(instances);
		// results
		int best = svm.getEvaluationResultList().getBestRootMeanSquaredErrorIndex();
		EvaluationResult bestRMSE = svm.getEvaluationResultList().get(best);
//		logger.println("Beset RMSE classification");
//		logger.println("Cost: " + svm.getCostValues()[best]);
//		logger.println(bestRMSE.toString());

		try {
			IOUtils.write(new File(outdir + "/svm.txt"), "Best RMSE classification\n\nCost: " + svm.getCostValues()[best] + "\n\n" + bestRMSE.toString());
			result.addOutputItem(new Item("svm_result_file", "svm.txt", "SVM result file", TYPE.FILE));
		} catch (IOException e) {
			printError("ioexception_executesvm_predictor", "Error saving SVM results", "Error saving SVM results");
		}
		// ??????
		
		// save results
		saveClassifierResults(svm);
	}

	private void executeRandomForest(Instances instances, Instances test) throws Exception {
		
		// update status
		updateStatus(progressCurrent,"executing Random Forest classifier");

		// init classifier
		RForest randomForest = new RForest();
		
		// init params
		if(commandLine.hasOption("random-forest-tune")) randomForest.setTuneParameters(true);
		else {
			if(commandLine.hasOption("random-forest-trees")) randomForest.setNumTrees(Integer.parseInt(commandLine.getOptionValue("random-forest-trees")));						
		}
		
		// validation
		setValidation(randomForest);
		
		// train
		randomForest.train(instances);
		
		// results
		int best = randomForest.getEvaluationResultList().getBestRootMeanSquaredErrorIndex();
		EvaluationResult bestRMSE = randomForest.getEvaluationResultList().get(best);
//		logger.println("Best RMSE classification");
//		logger.println("Number of trees: " + randomForest.getNumTreesArray()[best]);
//		logger.println(bestRMSE.toString());

		saveClassifierResults(randomForest);
		
	}

private void saveClassifierResults(GenericClassifier classifier){
		
		int best = classifier.getEvaluationResultList().getBestAreaUnderRocIndex();
		String name = classifier.getClassifierName();
				
		
		try {
			IOUtils.write(new File(outdir + "/" + name + ".txt"), "Best AUC classification\n\nParameters: " + classifier.getParamList().get(best) + "\n\n" + classifier.getEvaluationResultList().getBestAreaUnderRoc().toString());			
			result.addOutputItem(new Item(name + "_result_file", name + ".txt", name + " result file", TYPE.FILE,name + " results",""));
		} catch (IOException e) {
			printError("ioexception_execute_" + name + "_predictor", "Error saving " + name + " results", "Error saving " + name + " results");
		}
		
		// classification table
		try {
			IOUtils.write(new File(outdir + "/" + name + "_table.txt"), classifier.getResultsTable(true));
			result.addOutputItem(new Item(name + "_table", name + "_table.txt", name + " classifications", TYPE.FILE,Arrays.asList("TABLE","PREDICTOR_TABLE"),new HashMap<String,String>(),name + " results",""));
		} catch (IOException e) {
			printError("ioexception_" + name + "_predictor", "Error saving " + name + " classifications table", "Error saving " + name + " classifications table");
		}
				
		// Comparative plot
		try {
			XYLineChart comparativeplot = classifier.getComparativeMetricsPlot();
			comparativeplot.save(outdir + "/" + name + "_comparative_plot.png",500,350,"png");
			result.addOutputItem(new Item(name + "_comparative_plot", name + "_comparative_plot.png", name + " comparative plot", TYPE.IMAGE, Arrays.asList(""),new HashMap<String,String>(),name + " results",""));
		} catch (IOException e) {
			printError("ioexception_" + name + "_predictor", "Error saving " + name + " comparative plot", "Error saving " + name + " comparative plot");
			e.printStackTrace();
		}
		
//		// Sample classification rate plot
//		try {
//			Canvas sampleClassficationPlot = classifier.getSampleClassificationRatePlot();			 
//			sampleClassficationPlot.save(outdir + "/" + name + "_classification_rate_plot.png");
//			result.addOutputItem(new Item(name + "_classification_rate_plot", name + "_classification_rate_plot.png", name + " classification rate plot", TYPE.IMAGE, Arrays.asList(""),new HashMap<String,String>(),name + " results",""));
//		} catch (IOException e) {			
//			printError("ioexception_" + name + "_predictor", "Error saving " + name + " classification rate plot", "Error saving " + name + " classification rate plot");
//			e.printStackTrace();
//		} catch (InvalidColumnIndexException e) {
//			// TODO Auto-generated catch block
//			e.printStackTrace();
//		}
		
		combinedTable.append(classifier.getSortedResultsTable(numberOfBestSelectedClassifications));
	}
	

	private void setValidation(GenericClassifier classifier){
		logger.println("SETTING VALIDATION!!!!!!!!!!!!!!!!!!!!!!!");
		
		// feature selection
		AbstractFeatureSelector featureSelector = null;
		if(commandLine.hasOption("feature-selection")){	
			if("cfs".equalsIgnoreCase(commandLine.getOptionValue("feature-selection"))){
				featureSelector = new CfsFeatureSelector();
			} else if("pca".equalsIgnoreCase(commandLine.getOptionValue("feature-selection"))){
				featureSelector = new CfsFeatureSelector(); // change by PCA
			} else if("ga".equalsIgnoreCase(commandLine.getOptionValue("feature-selection"))){
				featureSelector = new CfsFeatureSelector(); // change by Genetic algorithm
			}
		}
		// validation
//		if(commandLine.hasOption("cross-validation")) {
			int repeats = Integer.parseInt(commandLine.getOptionValue("validation-repeats", "10"));
			int folds = Integer.parseInt(commandLine.getOptionValue("cross-validation-folds", "5"));
			classifier.setClassifierEvaluation(new KFoldCrossValidation(repeats, folds, featureSelector));
//			if(commandLine.hasOption("feature-selection")){
//				classifier.setClassifierEvaluation(new KFoldCrossValidation(repeats, folds, featureSelector));
//				logger.println("setting cross-validation evaluation (repeats = " + repeats + ", folds = " + folds + ", feature-selection = " + commandLine.getOptionValue("feature-selection") + ")");	
//			} else {
//				classifier.setClassifierEvaluation(new KFoldCrossValidation(repeats, folds));
//				logger.println("setting cross-validation evaluation (repeats = " + repeats + ", folds = " + folds + ", no feature selection");		
//			}					
//		}
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
	
	
	/*
	 * 
	 * STATUS MANAGEMENT (TO REMOVE!!!!!!!!)
	 * 
	 */
	
	private void initStatus() {
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
	
	private void updateStatus(double progressCurrent, String message){
		try {
			jobStatus.addStatusMessage(StringUtils.decimalFormat(progressCurrent),message);
			progressCurrent += progressStep;
		} catch (FileNotFoundException e) {
			abort("filenotfoundexception_executeknn_predictor", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
		}
	}
}
