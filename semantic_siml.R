#!/usr/bin/env Rscript

# TODO: add onassis to env specification 
suppressPackageStartupMessages(require(Onassis))
suppressPackageStartupMessages(require(optparse))
suppressPackageStartupMessages(require(workflowscriptscommon))
suppressPackageStartupMessages(require(org.Hs.eg.db))









# create ConceptMapper dictionary 
obo = "./cl-basic.obo"
