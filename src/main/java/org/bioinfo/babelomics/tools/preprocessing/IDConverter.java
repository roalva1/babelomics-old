package org.bioinfo.babelomics.tools.preprocessing;

import java.io.File;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;

import org.bioinfo.babelomics.tools.BabelomicsTool;
import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.commons.utils.ListUtils;
import org.bioinfo.commons.utils.MapUtils;
import org.bioinfo.commons.utils.StringUtils;
import org.bioinfo.infrared.common.dbsql.DBConnector;
import org.bioinfo.infrared.common.feature.FeatureList;
import org.bioinfo.infrared.core.XRef;
import org.bioinfo.infrared.core.dbsql.XRefDBManager;
import org.bioinfo.tool.OptionFactory;
import org.bioinfo.tool.result.Item;

public class IDConverter extends BabelomicsTool {

	public IDConverter() {
		initOptions();
	}

	@Override
	public void initOptions() {
		getOptions().addOption(OptionFactory.createOption("listfile", "File containning the IDs to convert", false));
		getOptions().addOption(OptionFactory.createOption("list", "IDs to convert separated by commas", false));
		
		getOptions().addOption(OptionFactory.createOption("affy_focus", "Affymx Microarray Focus", false, false)); 
		getOptions().addOption(OptionFactory.createOption("affy_hcg110", "Affymx Microarray HCG110", false, false)); 
		getOptions().addOption(OptionFactory.createOption("affy_hugenefl", "Affymx Microarray HuGeneFL", false, false)); 
		getOptions().addOption(OptionFactory.createOption("affy_human_exon_1.0_st_v2", "Affymx Microarray Human Exon 1.0 ST v2", false, false)); 
		getOptions().addOption(OptionFactory.createOption("affy_human_gene_1.0_st", "Affymx Microarray Human Gene 1.0 ST", false, false));
		getOptions().addOption(OptionFactory.createOption("affy_u133", "Affymx Microarray U133", false, false));
		getOptions().addOption(OptionFactory.createOption("affy_u95", "Affymx Microarray U95", false, false));
		getOptions().addOption(OptionFactory.createOption("agilent_cgh", "Agilent CGH", false, false));
		getOptions().addOption(OptionFactory.createOption("agilent_probe", "Agilent Probe", false, false));
		getOptions().addOption(OptionFactory.createOption("biocarta", "Biocarta", false, false));
		getOptions().addOption(OptionFactory.createOption("ccds", "CCDS", false, false));
		getOptions().addOption(OptionFactory.createOption("clone-based_(ensembl)", "Clone-based (Ensembl)", false, false));
		getOptions().addOption(OptionFactory.createOption("clone-based_(vega)", "Clone-based (Vega)", false, false));
		getOptions().addOption(OptionFactory.createOption("codelink", "GE Healthcare/Amersham Codelink WGA", false, false));
		getOptions().addOption(OptionFactory.createOption("embl", "EMBL", false, false));
		getOptions().addOption(OptionFactory.createOption("ensembl_gene", "Ensembl gene", false, false));
		getOptions().addOption(OptionFactory.createOption("ensembl_protein", "Ensembl protein", false, false));
		getOptions().addOption(OptionFactory.createOption("ensembl_transcript", "Ensembl transcript", false, false));
		getOptions().addOption(OptionFactory.createOption("ensembl_transcript_having_same_cds", "Ensembl transcript having same CDS", false, false));
		getOptions().addOption(OptionFactory.createOption("entrezgene", "EntrezGene", false, false));
		getOptions().addOption(OptionFactory.createOption("go", "GO", false, false));
		getOptions().addOption(OptionFactory.createOption("havana_gene", "Havana gene", false, false));
		getOptions().addOption(OptionFactory.createOption("havana_transcript", "Havana transcript", false, false));
		getOptions().addOption(OptionFactory.createOption("havana_transcript_having_same_cds", "Havana transcript having same CDS", false, false));
		getOptions().addOption(OptionFactory.createOption("havana_translation", "Havana translation", false, false));
		getOptions().addOption(OptionFactory.createOption("hgnc_(automatic)", "HGNC (automatic)", false, false));
		getOptions().addOption(OptionFactory.createOption("hgnc_(curated)", "HGNC (curated)", false, false));
		getOptions().addOption(OptionFactory.createOption("hgnc_symbol", "HGNC symbol", false, false));
		getOptions().addOption(OptionFactory.createOption("illumina_v1", "Illumina V1", false, false));
		getOptions().addOption(OptionFactory.createOption("illumina_v2", "Illumina V2", false, false));
		getOptions().addOption(OptionFactory.createOption("imgt/gene-db", "IMGT/GENE-DB", false, false));
		getOptions().addOption(OptionFactory.createOption("imgt/ligm-db", "IMGT/LIGM-DB", false, false));
		getOptions().addOption(OptionFactory.createOption("interpro", "Interpro", false, false));
		getOptions().addOption(OptionFactory.createOption("ipi", "IPI", false, false));
		getOptions().addOption(OptionFactory.createOption("kegg", "KEGG", false, false));
		getOptions().addOption(OptionFactory.createOption("mirbase", "miRBase", false, false));
		getOptions().addOption(OptionFactory.createOption("mirna_gene", "miRNA gene", false, false));
		getOptions().addOption(OptionFactory.createOption("omim_disease", "MIM disease", false, false));
		getOptions().addOption(OptionFactory.createOption("omim_gene", "MIM gene", false, false));
		getOptions().addOption(OptionFactory.createOption("pdb", "PDB", false, false));
		getOptions().addOption(OptionFactory.createOption("projected_go", "Projected GO", false, false));
		getOptions().addOption(OptionFactory.createOption("protein_id", "Protein ID", false, false));
		getOptions().addOption(OptionFactory.createOption("refseq_dna", "RefSeq DNA", false, false));
		getOptions().addOption(OptionFactory.createOption("refseq_dna_predicted", "RefSeq DNA predicted", false, false));
		getOptions().addOption(OptionFactory.createOption("refseq_peptide", "RefSeq peptide", false, false));
		getOptions().addOption(OptionFactory.createOption("refseq_peptide_predicted", "RefSeq peptide predicted", false, false));
		getOptions().addOption(OptionFactory.createOption("rfam", "RFAM", false, false));
		getOptions().addOption(OptionFactory.createOption("transcript_having_exact_match_between_ensembl_and_havana", "Transcript having exact match between ENSEMBL and HAVANA", false, false));
		getOptions().addOption(OptionFactory.createOption("ucsc_stable_id", "UCSC stable ID", false, false));
		getOptions().addOption(OptionFactory.createOption("unigene", "UniGene", false, false));
		getOptions().addOption(OptionFactory.createOption("uniprotkb/splicevariant", "UniprotKB/SpliceVariant", false, false));
		getOptions().addOption(OptionFactory.createOption("uniprotkb/swissprot", "UniProtKB/Swiss-Prot", false, false));
		getOptions().addOption(OptionFactory.createOption("uniprotkb/trembl", "UniProtKB/TrEMBL", false, false));
	}

	@Override
	protected void execute() {
		//try {	
		List<String> ids = null;
		List<String> dbNames = new ArrayList<String> ();

		// check input IDs
		//
		String inputIds = commandLine.getOptionValue("list", null);			
		if ( inputIds == null ) {
			String fileName = commandLine.getOptionValue("listfile", null);
			if ( fileName != null ) {
				try {
					ids = IOUtils.column(new File(fileName), 0);
				} catch (IOException e) {
					abort("ioexception_execute_idconverter", e.getMessage(), e.getMessage(), e.getMessage());
				}
			} else {
				abort("error_execute_idconverter", "Missing input file", "Missing input file", "Missing input file");					
			}
		} else {
			ids = StringUtils.toList(inputIds, ",");
		}

		if ( ids == null || ids.size() == 0) {
			abort("error_execute_idconverter", "No IDs found", "No IDs found", "No IDs found");
		}

		// check output refs
		//
//		if ( commandLine.hasOption("go") ) dbNames.add("go");
//		if ( commandLine.hasOption("entrezgene") ) dbNames.add("entrezgene");
//		if ( commandLine.hasOption("interpro") ) dbNames.add("interpro");

		if ( commandLine.hasOption("affy_focus") ) dbNames.add("affy_focus"); 
		if ( commandLine.hasOption("affy_hcg110") ) dbNames.add("affy_hcg110"); 
		if ( commandLine.hasOption("affy_hugenefl") ) dbNames.add("affy_hugenefl"); 
		if ( commandLine.hasOption("affy_human_exon_1.0_st_v2") ) dbNames.add("affy_human_exon_1.0_st_v2"); 
		if ( commandLine.hasOption("affy_human_gene_1.0_st") ) dbNames.add("affy_human_gene_1.0_st");
		if ( commandLine.hasOption("affy_u133") ) dbNames.add("affy_u133");
		if ( commandLine.hasOption("affy_u95") ) dbNames.add("affy_u95");
		if ( commandLine.hasOption("agilent_cgh") ) dbNames.add("agilent_cgh");
		if ( commandLine.hasOption("agilent_probe") ) dbNames.add("agilent_probe");
		if ( commandLine.hasOption("biocarta") ) dbNames.add("biocarta");
		if ( commandLine.hasOption("ccds") ) dbNames.add("ccds");
		if ( commandLine.hasOption("clone-based_(ensembl)") ) dbNames.add("clone-based_(ensembl)");
		if ( commandLine.hasOption("clone-based_(vega)") ) dbNames.add("clone-based_(vega)");
		if ( commandLine.hasOption("codelink") ) dbNames.add("codelink");
		if ( commandLine.hasOption("embl") ) dbNames.add("embl");
		if ( commandLine.hasOption("ensembl_gene") ) dbNames.add("ensembl_gene");
		if ( commandLine.hasOption("ensembl_protein") ) dbNames.add("ensembl_protein");
		if ( commandLine.hasOption("ensembl_transcript") ) dbNames.add("ensembl_transcript");
		if ( commandLine.hasOption("ensembl_transcript_having_same_cds") ) dbNames.add("ensembl_transcript_having_same_cds");
		if ( commandLine.hasOption("entrezgene") ) dbNames.add("entrezgene");
		if ( commandLine.hasOption("go") ) dbNames.add("go");
		if ( commandLine.hasOption("havana_gene") ) dbNames.add("havana_gene");
		if ( commandLine.hasOption("havana_transcript") ) dbNames.add("havana_transcript");
		if ( commandLine.hasOption("havana_transcript_having_same_cds") ) dbNames.add("havana_transcript_having_same_cds");
		if ( commandLine.hasOption("havana_translation") ) dbNames.add("havana_translation");
		if ( commandLine.hasOption("hgnc_(automatic)") ) dbNames.add("hgnc_(automatic)");
		if ( commandLine.hasOption("hgnc_(curated)") ) dbNames.add("hgnc_(curated)");
		if ( commandLine.hasOption("hgnc_symbol") ) dbNames.add("hgnc_symbol");
		if ( commandLine.hasOption("illumina_v1") ) dbNames.add("illumina_v1");
		if ( commandLine.hasOption("illumina_v2") ) dbNames.add("illumina_v2");
		if ( commandLine.hasOption("imgt/gene-db") ) dbNames.add("imgt/gene-db");
		if ( commandLine.hasOption("imgt/ligm-db") ) dbNames.add("imgt/ligm-db");
		if ( commandLine.hasOption("interpro") ) dbNames.add("interpro");
		if ( commandLine.hasOption("ipi") ) dbNames.add("ipi");
		if ( commandLine.hasOption("kegg") ) dbNames.add("kegg");
		if ( commandLine.hasOption("mirbase") ) dbNames.add("mirbase");
		if ( commandLine.hasOption("mirna_gene") ) dbNames.add("mirna_gene");
		if ( commandLine.hasOption("omim_disease") ) dbNames.add("omim_disease");
		if ( commandLine.hasOption("omim_gene") ) dbNames.add("omim_gene");
		if ( commandLine.hasOption("pdb") ) dbNames.add("pdb");
		if ( commandLine.hasOption("projected_go") ) dbNames.add("projected_go");
		if ( commandLine.hasOption("protein_id") ) dbNames.add("protein_id");
		if ( commandLine.hasOption("refseq_dna") ) dbNames.add("refseq_dna");
		if ( commandLine.hasOption("refseq_dna_predicted") ) dbNames.add("refseq_dna_predicted");
		if ( commandLine.hasOption("refseq_peptide") ) dbNames.add("refseq_peptide");
		if ( commandLine.hasOption("refseq_peptide_predicted") ) dbNames.add("refseq_peptide_predicted");
		if ( commandLine.hasOption("rfam") ) dbNames.add("rfam");
		if ( commandLine.hasOption("transcript_having_exact_match_between_ensembl_and_havana") ) dbNames.add("transcript_having_exact_match_between_ensembl_and_havana");
		if ( commandLine.hasOption("ucsc_stable_id") ) dbNames.add("ucsc_stable_id");
		if ( commandLine.hasOption("unigene") ) dbNames.add("unigene");
		if ( commandLine.hasOption("uniprotkb/splicevariant") ) dbNames.add("uniprotkb/splicevariant");
		if ( commandLine.hasOption("uniprotkb/swissprot") ) dbNames.add("uniprotkb/swissprot");
		if ( commandLine.hasOption("uniprotkb/trembl") ) dbNames.add("uniprotkb/trembl");
		
		if ( dbNames == null || dbNames.size() == 0 ) {
			abort("error_execute_idconverter", "Missing DB names", "Missing DB names", "Missing DB names");
		}


		try {
			DBConnector dbConnector = new DBConnector(species, new File(babelomicsHomePath + "/conf/infrared.properties"));
			
			StringBuilder line = new StringBuilder();
			List<String> alls = new ArrayList<String>();
			List<String> outIds = new ArrayList<String>(); 
			Map<String, List<String>> idsMap = new HashMap<String, List<String>>(); 
			
			System.out.println("-> infrared config = " + System.getenv("BABELOMICS_HOME") + "/conf/infrared.conf");
			
			System.out.println("-> db connector = " + dbConnector.toString());
			System.out.println("-> db names = " + ListUtils.toString(dbNames, ","));
			List<Map<String, FeatureList<XRef>>> list = new XRefDBManager(dbConnector).getListByDBNames(ids, dbNames);

			if ( list != null && list.size() > 0 ) {
				
				line.append("#NAMES\t").append(ListUtils.toString(MapUtils.getKeys(list.get(0)), "\t"));
				alls.add(line.toString());
				
				for(int i=0 ; i<list.size() ; i++) {
					line.delete(0, line.length());
					line.append(ids.get(i));
					System.out.println("Converting " + ids.get(i));
					for(String key: MapUtils.getKeys(list.get(i))) {
						System.out.print("to " + key + " ---> ");
						outIds.clear();
						for(XRef xref: list.get(i).get(key)) {
							outIds.add(xref.getId());
						}
						if ( ! idsMap.containsKey(key) ) {
							idsMap.put(key, new ArrayList<String>());
							idsMap.get(key).add("#NAMES\t" + key);
						}
						idsMap.get(key).add(ids.get(i) + "\t" + ListUtils.toString(outIds,","));
						line.append("\t").append(ListUtils.toString(outIds,","));
					}
					alls.add(line.toString());
				}
				
				// save results
				//
				String fileName = "alls.txt";
				IOUtils.write(new File(outdir + "/" + fileName), alls);
				result.addOutputItem(new Item("All IDs", fileName, "ID conversion", Item.TYPE.FILE, new ArrayList<String>(), new HashMap<String,String>(), "Results in a single file"));
				
				for(String key: MapUtils.getKeys(idsMap)) {
					fileName = key + ".txt";
					IOUtils.write(new File(outdir + "/" + fileName), idsMap.get(key));					
					result.addOutputItem(new Item(key + " IDs", fileName, "ID conversion", Item.TYPE.FILE, new ArrayList<String>(), new HashMap<String,String>(), "Results per file"));
				}
			}
		} catch (Exception e) {
			this.printError("exception_execute_idconverter", e.toString(), e.toString(), e);
		}
	}

}
