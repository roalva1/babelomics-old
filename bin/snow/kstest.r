## ks.test.r
## 2010-07-23 ralonso@cipf.es
## ks test with input data from java test library (testing snow versions)

setwd ("/home/ralonso/appl/babelomics/")

outfile <- Sys.getenv("outfile")
list1degree <- Sys.getenv("list1degree")
list2degree <- Sys.getenv("list2degree")
list1relBetweenness <- Sys.getenv("list1relBetweenness")
list2relBetweenness <- Sys.getenv("list2relBetweenness")
list1clustering <- Sys.getenv("list1clustering")
list2clustering <- Sys.getenv("list2clustering")
side <- Sys.getenv("side")

list1degree <- scan(list1degree, what = double(0), sep = '\n', quote = NULL)
list2degree <- scan(list2degree, what = double(0), sep = '\n', quote = NULL)
list1relBetweenness <- scan(list1relBetweenness, what = double(0), sep = '\n', quote = NULL)
list2relBetweenness <- scan(list2relBetweenness, what = double(0), sep = '\n', quote = NULL)
list1clustering <- scan(list1clustering, what = double(0), sep = '\n', quote = NULL)
list2clustering <- scan(list2clustering, what = double(0), sep = '\n', quote = NULL)

# ks test
resultDegree <- ks.test(list1degree,list2degree, alternative = side)
resultRelBetweenness <- ks.test(list1relBetweenness,list2relBetweenness, alternative = side)
resultClustering <- ks.test(list1clustering,list2clustering, alternative = side)

result <- paste("#parameter\tpval\tside\nconnections\t",resultDegree['p.value'], "\t",resultDegree['alternative'], "\nbetweenness\t", resultRelBetweenness['p.value'], "\t",resultRelBetweenness['alternative'],"\nclustering\t",resultClustering['p.value'], "\t",resultClustering['alternative'])

write.table(file=outfile, result, quote = FALSE, row.names = FALSE, col.names = FALSE)
#write.table(file=outfile,#parameter\tpval\tside)
#write("#parameter\tpval\tside",file=outfile)




