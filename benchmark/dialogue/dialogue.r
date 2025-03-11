library('DIALOGUE')

rA<-readRDS("test.example.rds")

# save the PCA components
write.csv(rA$A@X, file = paste0("./cca_output/", "X_pca_A.csv"))
write.csv(rA$B@X, file = paste0("./cca_output/", "X_pca_B.csv"))
write.csv(rA$C@X, file = paste0("./cca_output/", "X_pca_C.csv"))

param <- DLG.get.param(k = 3,
                       results.dir = "results/",
                       conf = c("gender","sample.quality","cellQ"), # Confounding factors
                       pheno = "pathology") # Phenotype (optional)

res<-DIALOGUE.run(rA = rA, # list of cell.type objects
                  main = "toyExample",
                  param = param)

# save the MCP scores for comparison in python
for (ct in names(res$cca.scores)) {
   write.csv(res$cca.scores[[ct]], file = paste0("./cca_output/", ct, "_cca_scores.csv"))
}
