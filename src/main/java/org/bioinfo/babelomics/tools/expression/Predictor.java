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
	}
	
}
