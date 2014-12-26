package org.bioinfo.babelomics.utils;


import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.infrared.core.common.FeatureList;
import org.bioinfo.infrared.core.funcannot.AnnotationItem;
import org.bioinfo.infrared.funcannot.filter.Filter;
import org.bioinfo.infrared.funcannot.filter.GOFilter;

import org.bioinfo.babelomics.utils.AnnotationManager;
import org.bioinfo.babelomics.utils.XrefManager;

import java.io.IOException;
import java.util.*;

/**
 * Created by ralonso on 12/26/14.
 */
public class Prueba {

    private List<String> list;
    private String species;
    private String db;
    private Filter filter;

    private String annotationFile;

    private List<String> rawAnnots;

    private String babelomicsHome = System.getenv("BABELOMICS_HOME");

    //    public Prueba(){}
    public Prueba(List<String> list, String species, String db, Filter filter) {
        this.list = list;
        this.species = species;
        this.db = db;
        this.filter = filter;
        String annotationFile = this.babelomicsHome + "/conf/annotations/" + this.species + "/" + db + ".txt";
        try {
            this.rawAnnots = IOUtils.readLines(annotationFile);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public FeatureList<AnnotationItem> getAnnotations() {

        /** Parse Xref **/
        //    public XrefManager(List<String> list, String species) {

        XrefManager xrefManager = new XrefManager(list, species);
        Map<String, Set<String>> list1Xref = xrefManager.getXrefs();

        /** Get map with this structure: [featureId in ENSTRANSCRIPT: anotationsList] **/
        Map<String, List<String>> featAnnot = parseFeatureAnnotations();

        /** Get map with the annotations related to the original id with this structure: [annotationId: originalIdList ]**/
        Map<String, Set<String>> annotFeat = this.parseAnnotationsOriginalId(list1Xref, featAnnot);

        /** Filter annotations **/
        FeatureList<AnnotationItem> annotations = this.filter(annotFeat);
        return annotations;
    }

    public Map<String, Set<String>> parseAnnotationsOriginalId(Map<String, Set<String>> list1Xref, Map<String, List<String>> featAnnot) {
        Map<String, Set<String>> annotFeat = new HashMap<String, Set<String>>();

        for (String id : list1Xref.keySet()) {
            Set<String> xrefs = list1Xref.get(id);
            for (String xref : xrefs) {
                if (featAnnot.containsKey(xref)) {
                    for (String annot : featAnnot.get(xref)) {
                        if (!annotFeat.containsKey(annot)) {
                            annotFeat.put(annot, new HashSet<String>());
                        }
                        Set<String> aux = annotFeat.get(annot);
                        aux.add(id);
                        annotFeat.put(annot, aux);
                    }
                }

            }
        }
        return annotFeat;
    }

    public FeatureList<AnnotationItem> filter(Map<String, Set<String>> annotFeature) {
        FeatureList<AnnotationItem> annotations = new FeatureList<AnnotationItem>();
        if (filter instanceof GOFilter) {
            int maxNumberGenes = ((GOFilter) filter).getMaxNumberGenes();
            int minNumberGenes = ((GOFilter) filter).getMinNumberGenes();
            for (String annot : annotFeature.keySet()) {
                Set<String> features = annotFeature.get(annot);
                /** Size filter **/
                if (minNumberGenes <= features.size() && features.size() <= maxNumberGenes) {
                    for (String feature : features) {
                        annotations.add(new AnnotationItem(feature, annot));
                    }
                }
            }
        }
        return annotations;
    }

    /**
     * Get map with this structure: [featureId:anotationsList] *
     */
    public Map<String, List<String>> parseFeatureAnnotations() {
        Map<String, List<String>> featureAnnot = new HashMap<String, List<String>>();
        for (String raw : rawAnnots) {
            if (raw.startsWith("#") || !raw.contains("\t")) {
                continue;
            }
            String fields[] = raw.split("\t");
            String feature = fields[0];
            String annot = fields[1];
            List<String> annotations = new ArrayList<String>();
            if (!featureAnnot.containsKey(feature)) {
                featureAnnot.put(feature, annotations);
            } else {
                annotations = featureAnnot.get(feature);

            }
            annotations.add(annot);
            featureAnnot.put(feature, annotations);

        }
        return featureAnnot;
    }
}
