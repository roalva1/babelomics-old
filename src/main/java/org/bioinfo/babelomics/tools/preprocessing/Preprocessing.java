package org.bioinfo.babelomics.tools.preprocessing;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.commons.io.utils.FileUtils;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ArrayUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.math.data.DoubleMatrix;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;

public class Preprocessing extends BabelomicsTool {


	public Preprocessing() {
		initOptions();
	}

	@Override
	public void initOptions() {
		options.addOption(OptionFactory.createOption("dataset", "the data"));
		options.addOption(OptionFactory.createOption("logarithm-base", "the logarithm base to apply transformation, possible values: e, 2 10 for log base 2, log base 2 and log base 10", false));
		options.addOption(OptionFactory.createOption("merge-replicates", "method to merge replicates, valid values are: mean, median", false));
		options.addOption(OptionFactory.createOption("filter-missing", "minimum percentage of existing values, from 0 to 100", false));
		options.addOption(OptionFactory.createOption("impute-missing", "method to impute missing values, valid values are: zero, mean, median, knn", false));
		options.addOption(OptionFactory.createOption("kvalue", "kvalue for knn impute method, default 15", false));
		options.addOption(OptionFactory.createOption("gene-list-filter", "This option will remove all the patterns of the genes that are not present in this gene list", false));
//		options.addOption(OptionFactory.createOption("sample-filter", "class variable", false));
//		options.addOption(OptionFactory.createOption("feature-filter", "class variable", false));
	}


	@Override
	public void execute() {

		Dataset dataset;
		try {

			String logBase = commandLine.getOptionValue("logarithm-base", null);
			String mergeMethod = commandLine.getOptionValue("merge-replicates", null);
			String filterPercentage = commandLine.getOptionValue("filter-missing", null);
			String imputeMethod = commandLine.getOptionValue("impute-missing", null);
			String kvalue = commandLine.getOptionValue("kvalue", "15");
			String filterFile = commandLine.getOptionValue("gene-list-filter", null);
			
			int progress = 1, finalProgress = 3;
			if ( imputeMethod != null && !("none".equalsIgnoreCase(imputeMethod)) ) finalProgress++; 
			if ( filterPercentage != null && !("none".equalsIgnoreCase(filterPercentage)) ) finalProgress++;
			if ( mergeMethod != null && !("none".equalsIgnoreCase(mergeMethod)) ) finalProgress++;
			if ( logBase != null && !("none".equalsIgnoreCase(logBase)) ) finalProgress++;
			
			jobStatus.addStatusMessage("" + (progress*100/finalProgress), "reading dataset");
			dataset = new Dataset(new File(commandLine.getOptionValue("dataset")));			
			if(commandLine.hasOption("sample-filter") || commandLine.hasOption("feature-filter")) {
				dataset = dataset.getSubDataset(commandLine.getOptionValue("sample-filter"), "4", commandLine.getOptionValue("feature-filter"), ""); 
			}
			progress++;

			//System.out.println("before executing preprocessing...\n" + dataset.toString() + "\n");

			// apply logarithm
			//
			if ( logBase != null && !("none".equalsIgnoreCase(logBase)) ) {
				
				jobStatus.addStatusMessage("" + (progress*100/finalProgress), "applying logarithm base " + logBase);
				logger.debug("executing logarithm base " + logBase + "...\n");
				
				try {
					if ( dataset.getDoubleMatrix() == null ) { 
						dataset.load();
						dataset.validate();
					}
					
					dataset.setDoubleMatrix(dataset.getDoubleMatrix().applyLogarithm(logBase));
					
//					logger.debug("validating...\n");
//					jobStatus.addStatusMessage("40", "reading data");
					dataset.validate();
					
//					logger.debug("saving...\n");
//					jobStatus.addStatusMessage("80", "saving data");
//					dataset.write(new File(outdir+"/preprocessed.txt"));
//					result.addOutputItem(new Item("prepocessed_file",outdir+"/preprocessed.txt", "The preprocessed file is: ", TYPE.FILE));
//					jobStatus.addStatusMessage("100", "done");
					
					//System.out.println("end of logarithm (base " + logBase + ")");
					
				} 
				catch (CloneNotSupportedException e) {
					printError("cloneexecption_logbase_execute_preprocessing", "logarithm base preprocessing error", e.toString(), e);
				}
				logger.debug("end of executing logarithm base " + logBase + "\n");
				progress++;
			}
			
			// merge replicated rows
			//
			if ( mergeMethod != null && !("none".equalsIgnoreCase(mergeMethod)) ) {
				jobStatus.addStatusMessage("" + (progress*100/finalProgress), "merging replicated rows (" + mergeMethod + ")");
				logger.debug("merging replicated rows (" + mergeMethod + ")...\n");

				if ( dataset.getDoubleMatrix() == null ) { 
					dataset.load();
					dataset.validate();
				}

				//System.out.println("merge replicated (" + mergeMethod + "), input :\n" + dataset.toString());
				dataset = dataset.mergeReplicatedFeatures(mergeMethod);
				//System.out.println("merge replicated (" + mergeMethod + "), output :\n" + dataset.toString());
				
				logger.debug("end of merging replicated rows (" + mergeMethod + ")\n");
				progress++;
			}


			// filter missing values according to the given percentage
			//
			if ( filterPercentage != null && !("none".equalsIgnoreCase(filterPercentage)) ) {
				jobStatus.addStatusMessage("" + (progress*100/finalProgress), "filtering missing values by percentage " + filterPercentage);
				logger.debug("filtering missing values by percentage " + filterPercentage + "...\n");
				try {
					double perc = Double.parseDouble(filterPercentage);					

					//System.out.println("input :\n" + dataset.toString());
					//System.out.println("executing filter missing values " + filterPercentage + "% ...\n");
					dataset = dataset.filterRowsByPercOfMissingValues(perc);
					dataset.validate();
					//System.out.println("output :\n" + dataset.toString());
					//System.out.println("end of filter missing values " + filterPercentage + "% ...\n");

				} catch (Exception e) {
					printError("execption_filterpercentage_execute_preprocessing", "filter percentage preprocessing error", e.toString(), e);
					e.printStackTrace();
				}
				logger.debug("end of filtering missing values by percentage " + filterPercentage + "\n");
				progress++;
			}


			// imputing missing values
			//
			if ( imputeMethod != null && !("none".equalsIgnoreCase(imputeMethod)) ) {
				jobStatus.addStatusMessage("" + (progress*100/finalProgress), "imputing missing values (" + imputeMethod + ")");
				logger.debug("imputing missing values (" + imputeMethod + ")...\n");
				try {
					if ( dataset.getDoubleMatrix() == null ) { 
						dataset.load();
						dataset.validate();
					}

					//System.out.println("input :\n" + dataset.toString());
					//System.out.println("executing impute missing values...\n");
					if ( "knn".equalsIgnoreCase(imputeMethod) ) {
						Dataset newDataset = knnImpute(dataset, Integer.parseInt(kvalue));
						if ( newDataset != null ) {
							dataset = newDataset;
						}
					} else {
						dataset.setDoubleMatrix(dataset.getDoubleMatrix().imputeMissingValuesInRows(imputeMethod));
					}
					dataset.validate();
					//System.out.println("output :\n" + dataset.toString());
					//System.out.println("end of impute missing values...\n");

				} 
				//				catch (InvalidParameterException e) {				
				//					logger.error("Invalid impute missing method '" + imputeMethod + "', valid values are zero, average, median");
				//					System.out.println("Invalid logarithm base '" + logBase + "', valid values are e, 2, 10\n");
				//					printUsage();
				//					return;
				//				}
				catch (CloneNotSupportedException e) {
					printError("cloneexecption_impute_execute_preprocessing", "impute missing values preprocessing error", e.toString(), e);
					e.printStackTrace();
				}	
				logger.debug("imputing missing values (" + imputeMethod + ")\n");
				progress++;
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
			
			if ( filterFile != null ) {
				if ( dataset.getDoubleMatrix() == null ) { 
					dataset.load();
					dataset.validate();
				}

				List<String> validGenes = IOUtils.grep(filterFile, "[^#]");
								
				if ( validGenes != null && validGenes.size() > 0 ) {
					List<String> genes = dataset.getFeatureNames();					
					List<Integer> rows = new ArrayList<Integer>();
					for (int row=0 ; row<genes.size() ; row++) {
						if ( !validGenes.contains(genes.get(row))) {
							rows.add(row);
						}
					}
					
					dataset = dataset.filterRows(rows);
					dataset.validate();					
				}
			}
			
			
			
			logger.debug("saving dataset...\n");
			jobStatus.addStatusMessage("" + (progress*100/finalProgress), "saving data");
			dataset.write(new File(this.getOutdir()+"/preprocessed.txt"));
			
			result.addOutputItem(new Item("prepocessed_file","preprocessed.txt", "The preprocessed file is: ", TYPE.FILE));
			jobStatus.addStatusMessage("100", "done");
		} catch (Exception e1) {
			try {
				jobStatus.addStatusMessage("100", "error");
			} catch (FileNotFoundException e) {
				e.printStackTrace();
			}
			e1.printStackTrace();
		}

	}
	
	
	private Dataset knnImpute(Dataset dataset, int kvalue) {
		Dataset newDataset = null;
		try {
			String wd = outdir + "/knn";
			File inputFile = new File(wd + "/in.txt");
			File outputFile = new File(wd + "/out.txt");
			if ( new File(wd).isDirectory() || FileUtils.createDirectory(wd) ) {
				dataset.write(inputFile);
				String cmdStr = System.getenv("BABELOMICS_HOME") + "/bin/KNNimpute -K=" + kvalue + " " + inputFile.getAbsolutePath() + " " + outputFile.getAbsolutePath();
				Runtime.getRuntime().exec(cmdStr);
				System.out.println("cmd = " + cmdStr);				
				if ( outputFile.exists() ) {
					List<String> values;
					List<String> lines = IOUtils.grep(outputFile, "[^#]");
					DoubleMatrix matrix = new DoubleMatrix(dataset.getRowDimension(), dataset.getColumnDimension());
					for(int row=0 ; row<lines.size() ; row++) {
						values = StringUtils.stringToList("\t", lines.get(row));
						values.remove(0);
						matrix.setRow(row, ArrayUtils.toDoubleArray(ListUtils.toStringArray(values)));
					}
					newDataset = dataset;
					newDataset.setDoubleMatrix(matrix);
				}
			}
		} catch (IOException e) {
			newDataset = null;
		}
		return newDataset;
	}
}
