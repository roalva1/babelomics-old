package org.bioinfo.babelomics.tools.preprocessing;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;

import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.commons.exec.Command;
import org.bioinfo.commons.exec.SingleProcess;
import org.bioinfo.commons.io.utils.FileUtils;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ArrayUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.data.dataset.Dataset;
import org.bioinfo.data.dataset.SampleVariable;
import org.bioinfo.math.data.DoubleMatrix;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;
import org.bioinfo.tool.result.Item.TYPE;

public class Preprocessing extends BabelomicsTool {

	public Preprocessing() {
	}

	@Override
	public void initOptions() {
		options.addOption(OptionFactory.createOption("dataset", "the data"));
		options.addOption(OptionFactory.createOption("logarithm-base", "the logarithm base to apply transformation, possible values: e, 2 10 for log base 2, log base 2 and log base 10", false));
		options.addOption(OptionFactory.createOption("merge-replicates", "method to merge replicates, valid values are: mean or median", false));
		options.addOption(OptionFactory.createOption("filter-missing", "minimum percentage of existing values, from 0 to 100", false));
		options.addOption(OptionFactory.createOption("impute-missing", "method to impute missing values, valid values are: zero, mean, median, knn", false));
		options.addOption(OptionFactory.createOption("kvalue", "kvalue for knn impute method, default 15", false));
		options.addOption(OptionFactory.createOption("gene-file-filter", "This option will remove all the patterns of the genes that are not present in this gene file", false));
		//		options.addOption(OptionFactory.createOption("sample-filter", "class variable", false));
		//		options.addOption(OptionFactory.createOption("feature-filter", "class variable", false));
	}

	public void execute2() {
		
	}
	
	@Override
	public void execute() {
		Dataset dataset = null;

		int kvalue = 15;
		String imputeMethodMsg = "";
		String logBase = commandLine.getOptionValue("logarithm-base", null);
		String mergeMethod = commandLine.getOptionValue("merge-replicates", null);
		String filterPercentage = commandLine.getOptionValue("filter-missing", null);
		String imputeMethod = commandLine.getOptionValue("impute-missing", null);
		String filterFilename = commandLine.getOptionValue("gene-file-filter", null);

		int progress = 1;
		int finalProgress = 3;
		if ( logBase != null && !("none".equalsIgnoreCase(logBase)) ) finalProgress++;
		if ( mergeMethod != null && !("none".equalsIgnoreCase(mergeMethod)) ) finalProgress++;
		if ( filterPercentage != null && !("none".equalsIgnoreCase(filterPercentage)) ) finalProgress++;
		if ( imputeMethod != null && !("none".equalsIgnoreCase(imputeMethod)) ) finalProgress++; 
		if ( filterFilename != null && new File(filterFilename).exists() ) finalProgress++; 

		imputeMethodMsg = imputeMethod;
		if ( "knn".equalsIgnoreCase(imputeMethod) ) {
			try {
				kvalue = Integer.parseInt(commandLine.getOptionValue("k-value", "15"));
			} catch (NumberFormatException e) {
				abort("numberformatexception_execute_preprocessing", "Invalid k-value '" + commandLine.getOptionValue("kvalue", "15") + "' for knn imputation", e.toString(), StringUtils.getStackTrace(e));
			}
			imputeMethodMsg = imputeMethodMsg + ", k-value = " + kvalue;
		}
		
		try {
			jobStatus.addStatusMessage("" + (progress*100/finalProgress), "reading dataset");
		} catch (FileNotFoundException e) {
			abort("filenotfoundexception_execute_preprocessing", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
		}

		String datasetPath = commandLine.getOptionValue("dataset");
		if ( datasetPath == null ) {
			abort("missingdataset_execute_preprocessing", "Missing dataset", "Missing dataset", "Missing dataset");
		}
		File datasetFile = new File(commandLine.getOptionValue("dataset"));
		try {
			dataset = new Dataset(datasetFile);
		} catch (Exception e) {
			abort("exception_execute_preprocessing", "error reading dataset '" + datasetFile.getName() + "'", e.toString(), StringUtils.getStackTrace(e));
		}		

		//		if(commandLine.hasOption("sample-filter") || commandLine.hasOption("feature-filter")) {
		//			dataset = dataset.getSubDataset(commandLine.getOptionValue("sample-filter"), "4", commandLine.getOptionValue("feature-filter"), "");
		//		}
		
		if ( dataset == null ) {
			abort("datasetisnull_execute_preprocessing", "dataset is null", "Dataset is null after reading file '" + datasetFile.getName() + "'", "Dataset is null after reading file " + datasetFile.getAbsolutePath());
		}
		
		progress++;


		
		// apply logarithm
		//
		if ( logBase != null && !("none".equalsIgnoreCase(logBase)) ) {

			try {
				jobStatus.addStatusMessage(StringUtils.decimalFormat((double)progress*100/finalProgress, "##.00"), "applying logarithm base " + logBase);
			} catch (FileNotFoundException e) {
				abort("filenotfoundexception_execute_preprocessing", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
			}
			logger.debug("executing logarithm base " + logBase + "...\n");

			try {
				if ( dataset.getDoubleMatrix() == null ) { 
					try {
						dataset.load();
					} catch (Exception e) {
						abort("exception_logbase_execute_preprocessing", "Error loading dataset '" + datasetFile.getName() + "' before applying logarithm base " + logBase, e.toString(), StringUtils.getStackTrace(e));
					}
					dataset.validate();
				}

				dataset.setDoubleMatrix(dataset.getDoubleMatrix().applyLogarithm(logBase));				
				dataset.validate();					
			} 
			catch (CloneNotSupportedException e) {
				abort("cloneexecption_logbase_execute_preprocessing", "logarithm base preprocessing error", e.toString(), StringUtils.getStackTrace(e));
			}
			logger.debug("end of executing logarithm base " + logBase + "\n");
			progress++;
		}

		// merge replicated rows
		//
		if ( mergeMethod != null && !("none".equalsIgnoreCase(mergeMethod)) ) {
			try {
				jobStatus.addStatusMessage("" + (progress*100/finalProgress), "merging replicated rows (" + mergeMethod + ")");
			} catch (FileNotFoundException e) {
				abort("filenotfoundexception_execute_preprocessing", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
			}
			logger.debug("merging replicated rows (" + mergeMethod + ")...\n");

			if ( dataset.getDoubleMatrix() == null ) { 
				try {
					dataset.load();
				} catch (Exception e) {
					abort("exception_mergereplicated_execute_preprocessing", "Error loading dataset '" + datasetFile.getName() + "' before merging replicated rows with method " + mergeMethod, e.toString(), StringUtils.getStackTrace(e));
				}
				dataset.validate();
			}

			try {
				dataset = dataset.mergeReplicatedFeatures(mergeMethod);
			} catch (Exception e) {
				abort("exception_mergereplicated_execute_preprocessing", "Error merging replicated rows with method " + mergeMethod, e.toString(), StringUtils.getStackTrace(e));
			}

			logger.debug("end of merging replicated rows (" + mergeMethod + ")\n");
			progress++;
		}


		// filter missing values according to the given percentage
		//
		if ( filterPercentage != null && !("none".equalsIgnoreCase(filterPercentage)) ) {
			try {
				jobStatus.addStatusMessage("" + (progress*100/finalProgress), "filtering missing values by percentage " + filterPercentage);
			} catch (FileNotFoundException e) {
				abort("filenotfoundexception_execute_preprocessing", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
			}
			logger.debug("filtering missing values by percentage " + filterPercentage + "...\n");

			if ( dataset.getDoubleMatrix() == null ) { 
				try {
					dataset.load();
				} catch (Exception e) {
					abort("exception_filterpermissingvalues_execute_preprocessing", "Error loading dataset '" + datasetFile.getName() + "' before filtering by missing values, percentage " + filterPercentage + "%", e.toString(), StringUtils.getStackTrace(e));
				}
				dataset.validate();
			}
						
			try {
				double perc = Double.parseDouble(filterPercentage);					

				dataset = dataset.filterRowsByPercOfMissingValues(perc);
				dataset.validate();

			} catch (Exception e) {
				abort("exception_filterbymissingvalues_execute_preprocessing", "Error filtering by missing values, percentage value of " + filterPercentage + "%", e.toString(), StringUtils.getStackTrace(e));
			}
			logger.debug("end of filtering missing values by percentage " + filterPercentage + "\n");
			progress++;
		}


		// imputing missing values
		//
		if ( imputeMethod != null && !("none".equalsIgnoreCase(imputeMethod)) ) {
			try {
				jobStatus.addStatusMessage("" + (progress*100/finalProgress), "imputing missing values (" + imputeMethodMsg + ")");
			} catch (FileNotFoundException e) {
				abort("filenotfoundexception_execute_preprocessing", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
			}
			logger.debug("imputing missing values (" + imputeMethodMsg + ")...\n");
			try {
				if ( dataset.getDoubleMatrix() == null ) { 
					try {
						dataset.load();
					} catch (Exception e) {
						abort("exception_imputingmissingvalues_execute_preprocessing", "Error loading dataset '" + datasetFile.getName() + "' before imputing missing values with method " + imputeMethodMsg, e.toString(), StringUtils.getStackTrace(e));
					}
					dataset.validate();
				}

				if ( "knn".equalsIgnoreCase(imputeMethod) ) {
					Dataset newDataset = knnImpute(dataset, kvalue);
					if ( newDataset != null ) {
						dataset = newDataset;
					}
				} else {
					dataset.setDoubleMatrix(dataset.getDoubleMatrix().imputeMissingValuesInRows(imputeMethod));
				}
				dataset.validate();

			} catch (CloneNotSupportedException e) {
				abort("cloneexeption_imputemissingvalues_execute_preprocessing", "Error imputing missing values with method " + imputeMethodMsg, e.toString(), StringUtils.getStackTrace(e));
			}	
			logger.debug("end of imputing missing values (" + imputeMethodMsg + ")\n");
			progress++;
		}

		// filter by gene names
		//
		if ( filterFilename != null && !("none".equalsIgnoreCase(filterFilename)) ) {
			try {
				jobStatus.addStatusMessage("" + (progress*100/finalProgress), "filtering by names");
			} catch (FileNotFoundException e) {
				abort("filenotfoundexception_execute_preprocessing", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
			}
			logger.debug("filtering by names...\n");

			if ( dataset.getDoubleMatrix() == null ) { 
				try {
					dataset.load();
				} catch (Exception e) {
					abort("exception_filterbynames_execute_preprocessing", "Error loading dataset '" + datasetFile.getName() + "' before filtering by names", e.toString(), StringUtils.getStackTrace(e));
				}
				dataset.validate();
			}

			List<String> validGenes = null;
			try {
				System.out.println("---> " + filterFilename + ", content = " + IOUtils.toString(filterFilename));
				
				validGenes = IOUtils.grep(filterFilename, "[^#].+", false);
//				validGenes = IOUtils.grep(filterFilename, "f");
				
				System.out.println("size of valid genes = " + validGenes.size());
				if ( validGenes.size() > 0 ) {
					System.out.println("list = " + ListUtils.toString(validGenes));
				}
				
				
			} catch (IOException e) {
				abort("ioexecption_filterbynames_execute_preprocessing", "Error reading names from file '" + new File(filterFilename).getName() + "' when filtering by names", e.toString(), StringUtils.getStackTrace(e));
			}
			
			//System.out.println("valid genes = " + ListUtils.toString(validGenes));

			if ( validGenes != null && validGenes.size() > 0 ) {
				List<String> genes = dataset.getFeatureNames();					
				List<Integer> rows = new ArrayList<Integer>();
				for (int row=0 ; row<genes.size() ; row++) {
					if ( !validGenes.contains(genes.get(row))) {
						rows.add(row);
					}
				}

				try {
					dataset = dataset.filterRows(rows);
				} catch (Exception e) {
					abort("exception_filterbynames_execute_preprocessing", "Error filtering rows by names", e.toString(), StringUtils.getStackTrace(e));
				}
				dataset.validate();	
				
				
				System.out.println("------------> feature data is null ? " + (dataset.getFeatureData() == null));
			}
			
			logger.debug("end of filtering by names\n");
			progress++;
		}



		logger.debug("saving dataset...\n");
		try {
			jobStatus.addStatusMessage("" + (progress*100/finalProgress), "saving results");
		} catch (FileNotFoundException e) {
			abort("filenotfoundexception_execute_preprocessing", "job status file not found", e.toString(), StringUtils.getStackTrace(e));
		}
		
		try {
			File file = new File(this.getOutdir() + "/preprocessed.txt");
			dataset.save(file);
			
			
			if ( file.exists() ) {				
				String tags = "DATA,DATAMATRIX,EXPRESSION";
				
				
				File redirectionFile = new File(outdir + "/clustering.redirection");
				createClusteringRedirectionFile(redirectionFile, file);
				if ( redirectionFile.exists() ) {
					tags = tags + ",REDIRECTION(" + redirectionFile.getName() + ":Send to Clustering tool...)";
				}
				
				if (dataset.getVariables() != null && dataset.getVariables().size() > 0 ) {
					redirectionFile = new File(outdir + "/classcomparison.redirection");
					createClassComparisonRedirectionFile(redirectionFile, file);
					if ( redirectionFile.exists() ) {
						tags = tags + ",REDIRECTION(" + redirectionFile.getName() + ":Send to Class-comparison tool...)";
					}
					
					redirectionFile = new File(outdir + "/correlation.redirection");
					createCorrelationRedirectionFile(redirectionFile, file);
					if ( redirectionFile.exists() ) {
						tags = tags + ",REDIRECTION(" + redirectionFile.getName() + ":Send to Correlation tool...)";
					}
					
					redirectionFile = new File(outdir + "/classprediction.redirection");
					createClassPredictionRedirectionFile(redirectionFile, file);
					if ( redirectionFile.exists() ) {
						tags = tags + ",REDIRECTION(" + redirectionFile.getName() + ":Send to Class-prediction tool...)";
					}
				}
				
				result.addOutputItem(new Item("prepocessed_file", file.getName(), "Preprocessed file", TYPE.FILE, StringUtils.toList(tags, ","), new HashMap<String, String>(2), "Preprocessed data"));
			}
		} catch (IOException e) {
			abort("ioexception_savingresults_execute_preprocessing", "error saving output file", e.toString(), StringUtils.getStackTrace(e));
		}		
	}


	private void createClusteringRedirectionFile(File redirectionFile, File fileToRedirect) {
		List<String> redirectionInputs = new ArrayList<String>();
		redirectionInputs.add("tool=clustering");
		redirectionInputs.add("jobname=clustering");
		redirectionInputs.add("jobdescription=redirected from job $JOB_NAME");
		redirectionInputs.add("dataset_databox=" + fileToRedirect.getName() + " (from job $JOB_NAME)");
		redirectionInputs.add("dataset=$JOB_FOLDER/" + fileToRedirect.getName());
		redirectionInputs.add("dataset_wum_data=true");
		redirectionInputs.add("method=upgma");
		redirectionInputs.add("distance=euclidean");
		try {
			IOUtils.write(redirectionFile.getAbsolutePath(), redirectionInputs);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private void createClassComparisonRedirectionFile(File redirectionFile, File fileToRedirect) {
		List<String> redirectionInputs = new ArrayList<String>();
		redirectionInputs.add("tool=class-comparison");
		redirectionInputs.add("jobname=class-comparison");
		redirectionInputs.add("jobdescription=redirected from job $JOB_NAME");
		redirectionInputs.add("dataset_databox=" + fileToRedirect.getName() + " (from job $JOB_NAME)");
		redirectionInputs.add("dataset=$JOB_FOLDER/" + fileToRedirect.getName());
		redirectionInputs.add("dataset_wum_data=true");
		redirectionInputs.add("correction=fdr");
		try {
			IOUtils.write(redirectionFile.getAbsolutePath(), redirectionInputs);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private void createCorrelationRedirectionFile(File redirectionFile, File fileToRedirect) {
		List<String> redirectionInputs = new ArrayList<String>();
		redirectionInputs.add("tool=correlation");
		redirectionInputs.add("jobname=correlation");
		redirectionInputs.add("jobdescription=redirected from job $JOB_NAME");
		redirectionInputs.add("dataset_databox=" + fileToRedirect.getName() + " (from job $JOB_NAME)");
		redirectionInputs.add("dataset=$JOB_FOLDER/" + fileToRedirect.getName());
		redirectionInputs.add("dataset_wum_data=true");
		redirectionInputs.add("test=pearson");
		redirectionInputs.add("correction=fdr");
		try {
			IOUtils.write(redirectionFile.getAbsolutePath(), redirectionInputs);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}

	private void createClassPredictionRedirectionFile(File redirectionFile, File fileToRedirect) {
		List<String> redirectionInputs = new ArrayList<String>();
		redirectionInputs.add("tool=class-prediction");
		redirectionInputs.add("jobname=class-prediction");
		redirectionInputs.add("jobdescription=redirected from job $JOB_NAME");
		redirectionInputs.add("dataset_databox=" + fileToRedirect.getName() + " (from job $JOB_NAME)");
		redirectionInputs.add("dataset=$JOB_FOLDER/" + fileToRedirect.getName());
		redirectionInputs.add("dataset_wum_data=true");
		redirectionInputs.add("svm=svm");
		try {
			IOUtils.write(redirectionFile.getAbsolutePath(), redirectionInputs);
		} catch (IOException e) {
			// TODO Auto-generated catch block
			e.printStackTrace();
		}
	}
	
	/**
	 * 
	 * @param dataset
	 * @param kvalue
	 * @return
	 */
	private Dataset knnImpute(Dataset dataset, int kvalue) {
		Dataset newDataset = null;
		String wd = outdir + "/knn";
		try {
			File inputFile = new File(wd + "/in.txt");
			File outputFile = new File(wd + "/out.txt");
			if ( new File(wd).isDirectory() || FileUtils.createDirectory(wd) ) {
				dataset.save(inputFile);
				
				String cmdStr = System.getenv("BABELOMICS_HOME") + "/bin/KNNimpute -K=" + kvalue + " " + inputFile.getAbsolutePath() + " " + outputFile.getAbsolutePath();
				Command cmd = new Command(cmdStr); 
				SingleProcess sp = new SingleProcess(cmd);
				sp.runSync();
				
				//System.out.println("cmd = " + cmdStr);				
				if ( outputFile.exists() ) {
					List<String> values;
					List<String> lines = IOUtils.grep(outputFile, "[^#].+");
					DoubleMatrix matrix = new DoubleMatrix(dataset.getRowDimension(), dataset.getColumnDimension());
					for(int row=0 ; row<lines.size() ; row++) {
						values = StringUtils.toList(lines.get(row), "\t");
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
		FileUtils.deleteDirectory(new File(wd));
		return newDataset;
	}
}
