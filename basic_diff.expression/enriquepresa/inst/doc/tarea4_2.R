## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE)

## ---- echo=FALSE--------------------------------------------------------------
doAll=FALSE

## ---- eval=doAll, results="hide"----------------------------------------------
#  pacman::p_load(BiocGenerics, dplyr, EnrichmentBrowser, enriquepresa, ReportingTools, tami)

## ---- eval=doAll, warning=FALSE, message=FALSE--------------------------------
#  data(PRJNA517180, package = "enriquepresa")
#  
#  cel_go_bp = getGenesets(org = "cel", db = "go", go.onto = "BP")
#  cel_go_cc = getGenesets(org = "cel", db = "go", go.onto = "CC")
#  cel_go_mf = getGenesets(org = "cel", db = "go", go.onto = "MF")

## ---- eval=doAll, message=FALSE-----------------------------------------------
#  PRJNA517180_bp_sc =
#    GeneSetTest(x = PRJNA517180,y="type",
#                test = edgercommon, association="pvalue",
#                correction="BH",
#                GeneNullDistr = "randomization",
#                GeneSetNullDistr = "self-contained",
#                alternative="two-sided",nmax = 100,
#                id = "ENTREZID",gsc=cel_go_bp,descriptive=maxmean,
#                foutput = "PRJNA517180_bp_sc")
#  
#  PRJNA517180_bp_sc = tidy(PRJNA517180_bp_sc)
#  PRJNA517180_bp_sc = (PRJNA517180_bp_sc[with(PRJNA517180_bp_sc, order(adjp)), ])
#  
#  PRJNA517180_bp_c =
#    GeneSetTest(x = PRJNA517180,y="type",
#                test = edgercommon,association="pvalue",
#                correction="BH",
#                GeneNullDistr = "randomization",
#                GeneSetNullDistr = "competitive",
#                alternative="two-sided",nmax = 100,
#                id = "ENTREZID",gsc=cel_go_bp,descriptive=maxmean,
#                foutput = "PRJNA517180_bp_c")
#  
#  PRJNA517180_bp_c = tidy(PRJNA517180_bp_c)
#  PRJNA517180_bp_c = PRJNA517180_bp_c[with(PRJNA517180_bp_c, order(adjp)), ]
#  
#  PRJNA517180_cc_sc =
#    GeneSetTest(x = PRJNA517180,y="type",
#                test = edgercommon,association="pvalue",
#                correction="BH",
#                GeneNullDistr = "randomization",
#                GeneSetNullDistr = "self-contained",
#                alternative="two-sided",nmax = 100,
#                id = "ENTREZID",gsc=cel_go_cc,descriptive=maxmean,
#                foutput = "PRJNA517180_cc_sc")
#  
#  PRJNA517180_cc_sc = tidy(PRJNA517180_cc_sc)
#  PRJNA517180_cc_sc = PRJNA517180_cc_sc[with(PRJNA517180_cc_sc, order(adjp)), ]
#  
#  PRJNA517180_cc_c =
#    GeneSetTest(x = PRJNA517180,y="type",
#                test = edgercommon,association="pvalue",
#                correction="BH",
#                GeneNullDistr = "randomization",
#                GeneSetNullDistr = "competitive",
#                alternative="two-sided",nmax = 100,
#                id = "ENTREZID",gsc=cel_go_cc,descriptive=maxmean,
#                foutput = "PRJNA517180_cc_c")
#  
#  PRJNA517180_cc_c = tidy(PRJNA517180_cc_c)
#  PRJNA517180_cc_c = PRJNA517180_cc_c[with(PRJNA517180_cc_c, order(adjp)), ]
#  
#  PRJNA517180_mf_sc =
#    GeneSetTest(x = PRJNA517180,y="type",
#                test = edgercommon,association="pvalue",
#                correction="BH",
#                GeneNullDistr = "randomization",
#                GeneSetNullDistr = "self-contained",
#                alternative="two-sided",nmax = 100,
#                id = "ENTREZID",gsc=cel_go_mf, descriptive=maxmean,
#                foutput = "PRJNA517180_mf_sc")
#  
#  PRJNA517180_mf_sc = tidy(PRJNA517180_mf_sc)
#  PRJNA517180_mf_sc = PRJNA517180_mf_sc[with(PRJNA517180_mf_sc, order(adjp)), ]
#  
#  PRJNA517180_mf_c =
#    GeneSetTest(x = PRJNA517180,y="type",
#                test = edgercommon,association="pvalue",
#                correction="BH",
#                GeneNullDistr = "randomization",
#                GeneSetNullDistr = "competitive",
#                alternative="two-sided",nmax = 100,
#                id = "ENTREZID",gsc=cel_go_mf,descriptive=maxmean,
#                foutput = "PRJNA517180_mf_c")
#  
#  PRJNA517180_mf_c = tidy(PRJNA517180_mf_c)
#  PRJNA517180_mf_c = PRJNA517180_mf_c[with(PRJNA517180_mf_c, order(adjp)), ]

## ---- eval=doAll, warning=FALSE, message=FALSE--------------------------------
#  PRJNA517180_bp_sc = enriquepresa::dfGO(PRJNA517180_bp_sc)
#  PRJNA517180_bp_c = enriquepresa::dfGO(PRJNA517180_bp_c)
#  PRJNA517180_mf_sc = enriquepresa::dfGO(PRJNA517180_mf_sc)
#  PRJNA517180_mf_c = enriquepresa::dfGO(PRJNA517180_mf_c)
#  PRJNA517180_cc_sc = enriquepresa::dfGO(PRJNA517180_cc_sc)
#  PRJNA517180_cc_c = enriquepresa::dfGO(PRJNA517180_cc_c)
#  
#  head(PRJNA517180_bp_sc)
#  
#  foutput = "PRJNA517180_bp_sc"
#  htmlRep1 = HTMLReport(shortName = foutput,title = foutput,
#                        reportDirectory = "./reports/PRJNA517180_gsa")
#  
#  publish(PRJNA517180_bp_sc, htmlRep1)
#  finish(htmlRep1)
#  
#  foutput = "PRJNA517180_bp_c"
#  htmlRep1 = HTMLReport(shortName = foutput,title = foutput,
#                        reportDirectory = "./reports/PRJNA517180_gsa")
#  
#  publish(PRJNA517180_bp_c, htmlRep1)
#  finish(htmlRep1)
#  
#  foutput = "PRJNA517180_cc_sc"
#  htmlRep1 = HTMLReport(shortName = foutput,title = foutput,
#                        reportDirectory = "./reports/PRJNA517180_gsa")
#  
#  publish(PRJNA517180_cc_sc, htmlRep1)
#  finish(htmlRep1)
#  
#  foutput = "PRJNA517180_cc_c"
#  htmlRep1 = HTMLReport(shortName = foutput,title = foutput,
#                        reportDirectory = "./reports/PRJNA517180_gsa")
#  
#  publish(PRJNA517180_cc_c, htmlRep1)
#  finish(htmlRep1)
#  
#  foutput = "PRJNA517180_mf_sc"
#  htmlRep1 = HTMLReport(shortName = foutput,title = foutput,
#                        reportDirectory = "./reports/PRJNA517180_gsa")
#  
#  publish(PRJNA517180_mf_sc, htmlRep1)
#  finish(htmlRep1)
#  
#  foutput = "PRJNA517180_mf_c"
#  htmlRep1 = HTMLReport(shortName = foutput,title = foutput,
#                        reportDirectory = "./reports/PRJNA517180_gsa")
#  
#  publish(PRJNA517180_mf_c, htmlRep1)
#  finish(htmlRep1)

## ---- eval=doAll, warning=FALSE, message=FALSE--------------------------------
#  data(Eset50467, package = "enriquepresa")
#  
#  hsa_go_bp = getGenesets(org = "hsa", db = "go", go.onto = "BP")
#  hsa_go_cc = getGenesets(org = "hsa", db = "go", go.onto = "CC")
#  hsa_go_mf = getGenesets(org = "hsa", db = "go", go.onto = "MF")

## ---- eval=doAll, message=FALSE-----------------------------------------------
#  Eset50467_bp_sc =
#    GeneSetTest(x = Eset50467,y="FactorValue..compound.",
#                test = rowtmod, association="pvalue",
#                correction="BH",
#                GeneNullDistr = "randomization",
#                GeneSetNullDistr = "self-contained",
#                alternative="two-sided",nmax = 100,
#                id = "ENTREZID",gsc=hsa_go_bp,descriptive=maxmean,
#                foutput = "Eset50467_bp_sc")
#  
#  Eset50467_bp_sc = tidy(Eset50467_bp_sc)
#  Eset50467_bp_sc = Eset50467_bp_sc[with(Eset50467_bp_sc, order(adjp)), ]
#  
#  Eset50467_bp_c =
#    GeneSetTest(x = Eset50467,y="FactorValue..compound.",
#                test = rowtmod,association="pvalue",
#                correction="BH",
#                GeneNullDistr = "randomization",
#                GeneSetNullDistr = "competitive",
#                alternative="two-sided",nmax = 100,
#                id = "ENTREZID",gsc=hsa_go_bp,descriptive=maxmean,
#                foutput = "Eset50467_bp_c")
#  
#  Eset50467_bp_c = tidy(Eset50467_bp_c)
#  Eset50467_bp_c = Eset50467_bp_c[with(Eset50467_bp_c, order(adjp)), ]
#  
#  Eset50467_cc_sc =
#    GeneSetTest(x = Eset50467,y="FactorValue..compound.",
#                test = rowtmod,association="pvalue",
#                correction="BH",
#                GeneNullDistr = "randomization",
#                GeneSetNullDistr = "self-contained",
#                alternative="two-sided",nmax = 100,
#                id = "ENTREZID",gsc=hsa_go_cc,descriptive=maxmean,
#                foutput = "Eset50467_cc_sc")
#  
#  Eset50467_cc_sc = tidy(Eset50467_cc_sc)
#  Eset50467_cc_sc = Eset50467_cc_sc[with(Eset50467_cc_sc, order(adjp)), ]
#  
#  Eset50467_cc_c =
#    GeneSetTest(x = Eset50467,y="FactorValue..compound.",
#                test = rowtmod,association="pvalue",
#                correction="BH",
#                GeneNullDistr = "randomization",
#                GeneSetNullDistr = "competitive",
#                alternative="two-sided",nmax = 100,
#                id = "ENTREZID",gsc=hsa_go_cc,descriptive=maxmean,
#                foutput = "Eset50467_cc_c")
#  
#  Eset50467_cc_c = tidy(Eset50467_cc_c)
#  Eset50467_cc_c = Eset50467_cc_c[with(Eset50467_cc_c, order(adjp)), ]
#  
#  Eset50467_mf_sc =
#    GeneSetTest(x = Eset50467,y="FactorValue..compound.",
#                test = rowtmod,association="pvalue",
#                correction="BH",
#                GeneNullDistr = "randomization",
#                GeneSetNullDistr = "self-contained",
#                alternative="two-sided",nmax = 100,
#                id = "ENTREZID",gsc=hsa_go_mf, descriptive=maxmean,
#                foutput = "Eset50467_mf_sc")
#  
#  Eset50467_mf_sc = tidy(Eset50467_mf_sc)
#  Eset50467_mf_sc = Eset50467_mf_sc[with(Eset50467_mf_sc, order(adjp)), ]
#  
#  Eset50467_mf_c =
#    GeneSetTest(x = Eset50467,y="FactorValue..compound.",
#                test = rowtmod,association="pvalue",
#                correction="BH",
#                GeneNullDistr = "randomization",
#                GeneSetNullDistr = "competitive",
#                alternative="two-sided",nmax = 100,
#                id = "ENTREZID",gsc=hsa_go_mf,descriptive=maxmean,
#                foutput = "Eset50467_mf_c")
#  
#  Eset50467_mf_c = tidy(Eset50467_mf_c)
#  Eset50467_mf_c = Eset50467_mf_c[with(Eset50467_mf_c, order(adjp)), ]

## ---- eval=doAll, warning=FALSE, message=FALSE--------------------------------
#  Eset50467_bp_sc = enriquepresa::dfGO(Eset50467_bp_sc)
#  Eset50467_bp_c = enriquepresa::dfGO(Eset50467_bp_c)
#  Eset50467_mf_sc = enriquepresa::dfGO(Eset50467_mf_sc)
#  Eset50467_mf_c = enriquepresa::dfGO(Eset50467_mf_c)
#  Eset50467_cc_sc = enriquepresa::dfGO(Eset50467_cc_sc)
#  Eset50467_cc_c = enriquepresa::dfGO(Eset50467_cc_c)
#  
#  head(Eset50467_bp_sc)
#  
#  foutput = "Eset50467_bp_sc"
#  htmlRep1 = HTMLReport(shortName = foutput,title = foutput,
#                        reportDirectory = "./reports/Eset50467_gsa")
#  
#  publish(Eset50467_bp_sc, htmlRep1)
#  finish(htmlRep1)
#  
#  foutput = "Eset50467_bp_c"
#  htmlRep1 = HTMLReport(shortName = foutput,title = foutput,
#                        reportDirectory = "./reports/Eset50467_gsa")
#  
#  publish(Eset50467_bp_c, htmlRep1)
#  finish(htmlRep1)
#  
#  foutput = "Eset50467_cc_sc"
#  htmlRep1 = HTMLReport(shortName = foutput,title = foutput,
#                        reportDirectory = "./reports/Eset50467_gsa")
#  
#  publish(Eset50467_cc_sc, htmlRep1)
#  finish(htmlRep1)
#  
#  foutput = "Eset50467_cc_c"
#  htmlRep1 = HTMLReport(shortName = foutput,title = foutput,
#                        reportDirectory = "./reports/Eset50467_gsa")
#  
#  publish(Eset50467_cc_c, htmlRep1)
#  finish(htmlRep1)
#  
#  foutput = "Eset50467_mf_sc"
#  htmlRep1 = HTMLReport(shortName = foutput,title = foutput,
#                        reportDirectory = "./reports/Eset50467_gsa")
#  
#  publish(Eset50467_mf_sc, htmlRep1)
#  finish(htmlRep1)
#  
#  foutput = "Eset50467_mf_c"
#  htmlRep1 = HTMLReport(shortName = foutput,title = foutput,
#                        reportDirectory = "./reports/Eset50467_gsa")
#  
#  publish(Eset50467_mf_c, htmlRep1)
#  finish(htmlRep1)

## ---- eval=doAll--------------------------------------------------------------
#  lista = list(PRJNA517180_bp_c, PRJNA517180_bp_sc, PRJNA517180_cc_c, PRJNA517180_cc_sc, PRJNA517180_mf_c, PRJNA517180_mf_sc)
#  names(lista) = c("PRJNA517180_bp_c", "PRJNA517180_bp_sc", "PRJNA517180_cc_c", "PRJNA517180_cc_sc", "PRJNA517180_mf_c", "PRJNA517180_mf_sc")
#  
#  lista2 = list(Eset50467_bp_c, Eset50467_bp_sc, Eset50467_cc_c, Eset50467_cc_sc, Eset50467_mf_c, Eset50467_mf_sc)
#  names(lista2) = c("Eset50467_bp_c", "Eset50467_bp_sc", "Eset50467_cc_c", "Eset50467_cc_sc", "Eset50467_mf_c", "Eset50467_mf_sc")
#  
#  sapply(lista, function(x) length(which(x$adjp<0.05)))
#  sapply(lista2, function(x) length(which(x$adjp<0.05)))

