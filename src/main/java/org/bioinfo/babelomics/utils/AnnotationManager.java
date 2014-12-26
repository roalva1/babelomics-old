package org.bioinfo.babelomics.utils;

import org.bioinfo.commons.io.utils.IOUtils;
import org.bioinfo.infrared.core.common.FeatureList;
import org.bioinfo.infrared.core.funcannot.AnnotationItem;
import org.bioinfo.infrared.funcannot.filter.Filter;
import org.bioinfo.infrared.funcannot.filter.GOFilter;

import java.io.IOException;
import java.util.*;

/**
 * Created by ralonso on 12/22/14.
 */
public class AnnotationManager {

    private String fileName;
    List<String> rawAnnots;
    Set<String> uniqAnnots;
    Map<String, List<String>> annotFeature;

    public AnnotationManager(String fileName) {
        this.fileName = fileName;
        try {
            this.rawAnnots = IOUtils.readLines(fileName);
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    public void process() {
        uniqAnnots = removeDuplicates();
        annotFeature = getAnnotationFeatures();
    }

    public Set<String> removeDuplicates() {
        /** Remove duplicates annotation**/
        Set<String> uniqAnnots = new HashSet<String>();

        for (String annot : rawAnnots) {
            if (annot.startsWith("#") || !annot.contains("\t")) {
                continue;
            }
            if (!uniqAnnots.contains(annot))
                uniqAnnots.add(annot);
        }
        return uniqAnnots;
    }

    public Map<String, List<String>> getAnnotationFeatures() {
        Map<String, List<String>> annotFeature = new HashMap<String, List<String>>();
        for (String uniqs : uniqAnnots) {
            String uniqsArr[] = uniqs.split("\t");
            String feature = uniqsArr[0];
            String annot = uniqsArr[1];
            List<String> features = new ArrayList<String>();
            if (!annotFeature.containsKey(annot)) {
                annotFeature.put(annot, features);
            } else {
                features = annotFeature.get(annot);

            }
            features.add(feature);
            annotFeature.put(annot, features);
        }
        return annotFeature;
    }

    public FeatureList<AnnotationItem> filter(Filter filter) {
        FeatureList<AnnotationItem> annotations = new FeatureList<AnnotationItem>();
        if (filter instanceof GOFilter) {
            int maxNumberGenes = ((GOFilter) filter).getMaxNumberGenes();
            int minNumberGenes = ((GOFilter) filter).getMinNumberGenes();
            for (String annot : annotFeature.keySet()) {
                List<String> features = annotFeature.get(annot);
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
    public List<String> getRawAnnots(){
        return this.rawAnnots;
    }

}
