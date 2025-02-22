library('DIALOGUE')

rA<-readRDS("/ictstr01/home/icb/yuge.ji/projects/DIALOGUE/Data/test.example.rds")

# save the PCA components
write.csv(rA$A@X, file = paste0("./cca_output/", "X_pca_A.csv"))
write.csv(rA$B@X, file = paste0("./cca_output/", "X_pca_B.csv"))
write.csv(rA$C@X, file = paste0("./cca_output/", "X_pca_C.csv"))

res<-DIALOGUE.run(rA = rA,main = "toyExample", conf=c("cellQ"))

# save the MCP scores for comparison in python
for (ct in names(res$cca.scores)) {
   write.csv(res$cca.scores[[ct]], file = paste0("./cca_output/", ct, "_cca_scores.csv"))
}