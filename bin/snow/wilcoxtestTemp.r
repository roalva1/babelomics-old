## kstest.r
## 2010-07-23 aamadoz@cipf.es
## kolmogorov test with input data from java test library (testing snow versions)

setwd ("/home/ralonso/appl/babelomics/")

## load functions
#source (file.path (.job$funcdir, "functions.r"))

## load ks input data
listfile <- "connections_list1"
listdata <- scan(listfile, what = double(0), sep = '\n', quote = NULL)
intfile <- "connections_list2"
intdata <- scan(intfile, what = double(0), sep = '\n', quote = NULL)

# ks test
#ks.test(listdata, intdata, alternative = "two.sided")
#ks.test(listdata, intdata, alternative = "less")
wilcox.test(listdata, intdata, alternative = "greater")
