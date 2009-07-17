package org.bioinfo.babelomics.tools.functional;

public class Score {
	private String name;
	private float value;

	public Score(String name, Float value) {
		this.setName(name);
		this.value = value;
	}

	public void setName(String name) {
		this.name = name;
	}

	public String getName() {
		return name;
	}
	
	public void setValue(float value) {
		this.value = value;
	}

	public float getValue() {
		return value;
	}
	
	public String toString() {
		return "[" + name + ", " + value + "]";
	}
	
}
