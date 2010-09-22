package org.bioinfo.babelomics.methods.functional;

import java.util.List;

import org.bioinfo.infrared.common.feature.FeatureList;
import org.bioinfo.infrared.funcannot.AnnotationItem;

public abstract class FunctionalTest {

	
	public abstract void test(List<String> list1, List<String> list2, FeatureList<AnnotationItem> items, int testMode);
	
	
}
