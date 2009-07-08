package org.bioinfo.babelomics.tools.preprocessing;

import java.io.File;
import java.io.IOException;

import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;

public class Preprocessing extends BabelomicsTool {


	public Preprocessing(String[] args) {
		super(args);
		initOptions();
	}

	@Override
	public void initOptions() {
		//super.initOptions();
		options.addOption(OptionFactory.createOption("dataset", "the data"));
		options.addOption(OptionFactory.createOption("logarithm-base", "the logarithm base to apply transformation, possible values: e, 2 10 for log base 2, log base 2 and log base 10"));
		options.addOption(OptionFactory.createOption("merge-replicates", "method to merge replicates, valid values are: mean, median", false));
		options.addOption(OptionFactory.createOption("filter-missing", "minimum percentage of existing values, from 0 to 100", false));
		options.addOption(OptionFactory.createOption("impute-missing", "method to impute missing values, valid values are: zero, mean, median, knn", false));
		options.addOption(OptionFactory.createOption("kvalue", "kvalue for knn impute method, default 15", false));
		options.addOption(OptionFactory.createOption("sample-filter", "class variable", false));
		options.addOption(OptionFactory.createOption("feature-filter", "class variable", false));
	}


	@Override
	public void execute() {

		Dataset dataset;
		try {
			dataset = new Dataset(new File(commandLine.getOptionValue("dataset")));

			String logBase = commandLine.getOptionValue("logarithm-base", null);
			String mergeMethod = commandLine.getOptionValue("merge-replicates", null);
			String filterPercentage = commandLine.getOptionValue("filter-missing", null);
			String imputeMethod = commandLine.getOptionValue("impute-missing", null);
			String kvalue = commandLine.getOptionValue("impute-missing", "15");

			if(commandLine.hasOption("sample-filter") || commandLine.hasOption("feature-filter")) {
				dataset = dataset.getSubDataset(commandLine.getOptionValue("sample-filter"), "4", commandLine.getOptionValue("feature-filter"), ""); 
			}

			//System.out.println("before executing preprocessing...\n" + dataset.toString() + "\n");

			// apply logarithm
			//
			if ( logBase != null ) {
				try {
					if ( dataset.getDoubleMatrix() == null ) { 
						dataset.load();
						dataset.validate();
					}
					
					jobStatus.addStatusMessage("25", "reading data");
					logger.debug("executing logarithm...\n");
					dataset.setDoubleMatrix(dataset.getDoubleMatrix().applyLogarithm(logBase));
					
					logger.debug("validating...\n");
					jobStatus.addStatusMessage("40", "reading data");
					dataset.validate();
					
					logger.debug("saving...\n");
					jobStatus.addStatusMessage("80", "saving data");
					dataset.save(outdir+"/preprocessed.txt");
					result.addOutputItem(new Item("prepocessed_file",outdir+"/preprocessed.txt", "The preprocessed file is: ", TYPE.FILE));
					jobStatus.addStatusMessage("100", "done");
				} 
				//				catch (InvalidParameterException e) {				
				//					logger.error("Invalid logarithm base '" + logBase + "', valid values are e, 2, 10");
				//					System.out.println("Invalid logarithm base '" + logBase + "', valid values are e, 2, 10\n");
				//					printUsage();
				//					return;
				//				}
				catch (CloneNotSupportedException e) {
					e.printStackTrace();
				}
			}

			// merge replicated rows
			//
			if ( mergeMethod != null ) {
				if ( dataset.getDoubleMatrix() == null ) { 
					dataset.load();
					dataset.validate();
				}

				System.out.println("merge replicated (" + mergeMethod + "), input :\n" + dataset.toString());
				dataset = dataset.mergeReplicatedFeatures(mergeMethod);
				System.out.println("merge replicated (" + mergeMethod + "), output :\n" + dataset.toString());
			}


			// filter missing values according to the given percentage
			//
			if ( filterPercentage != null ) {
				try {
					double perc = Double.parseDouble(filterPercentage);					

					System.out.println("input :\n" + dataset.toString());
					System.out.println("executing filter missing values " + filterPercentage + "% ...\n");
					dataset = dataset.filterRowsByPercOfMissingValues(perc);
					dataset.validate();
					System.out.println("output :\n" + dataset.toString());
					System.out.println("end of filter missing values " + filterPercentage + "% ...\n");

				} catch (Exception e) {
					e.printStackTrace();
					logger.error("Invalid percentage value '" + filterPercentage + "' to filter missing values, it must range from 0 to 100");
					System.out.println("Invalid percentage value '" + filterPercentage + "' to filter missing values, it must range from 0 to 100\n");
					printUsage();
					return;
				}
			}


			// imputing missing values
			//
			if ( imputeMethod != null ) {
				try {
					if ( dataset.getDoubleMatrix() == null ) { 
						dataset.load();
						dataset.validate();
					}

					System.out.println("input :\n" + dataset.toString());
					System.out.println("executing impute missing values...\n");
					dataset.setDoubleMatrix(dataset.getDoubleMatrix().imputeMissingValuesInRows(imputeMethod));
					dataset.validate();
					System.out.println("output :\n" + dataset.toString());
					System.out.println("end of impute missing values...\n");

				} 
				//				catch (InvalidParameterException e) {				
				//					logger.error("Invalid impute missing method '" + imputeMethod + "', valid values are zero, average, median");
				//					System.out.println("Invalid logarithm base '" + logBase + "', valid values are e, 2, 10\n");
				//					printUsage();
				//					return;
				//				}
				catch (CloneNotSupportedException e) {
					e.printStackTrace();
				}	
				/*				
				int k;
				try {
					k = Integer.parseInt(kvalue);
					dMatrix = dMatrix.imputeMissingValuesInRows(imputeMethod);
				} catch (NumberFormatException e) {
					logger.error("Invalid k-value '" + kvalue + "'");
					System.out.println("Invalid k-value '" + kvalue + "'\n");
					printUsage();
					return;
				}
				 */
			}
		} catch (IOException e1) {
			e1.printStackTrace();
		}

	}
}
