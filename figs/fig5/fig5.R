library("Seurat")
library("mitch")

load("seurat_data.RData")

#unique(immune.combined@meta.data$celltype.stim)

cd14.mono<-FindMarkers(immune.combined, ident.1 = "CD14 Mono_STIM", ident.2 = "CD14 Mono_CTRL", print.bar = FALSE)
pdc<-FindMarkers(immune.combined, ident.1 = "pDC_STIM", ident.2 = "pDC_CTRL", print.bar = FALSE)
cd4.mem <- FindMarkers(immune.combined, ident.1 = "CD4 Memory T_STIM", ident.2 = "CD4 Memory T_CTRL", print.bar = FALSE)
t.act <- FindMarkers(immune.combined, ident.1 = "T activated_STIM", ident.2 = "T activated_CTRL", print.bar = FALSE)
cd4.nat <- FindMarkers(immune.combined, ident.1 = "CD4 Naive T_STIM", ident.2 = "CD4 Naive T_CTRL", print.bar = FALSE)
cd8.t <- FindMarkers(immune.combined, ident.1 = "CD8 T_STIM", ident.2 = "CD8 T_CTRL", print.bar = FALSE)
mk <- FindMarkers(immune.combined, ident.1 = "Mk_STIM", ident.2 = "Mk_CTRL", print.bar = FALSE)
b.act <- FindMarkers(immune.combined, ident.1 = "B activated_STIM", ident.2 = "B activated_CTRL", print.bar = FALSE)
b <- FindMarkers(immune.combined, ident.1 = "B_STIM", ident.2 = "B_CTRL", print.bar = FALSE)
dc <- FindMarkers(immune.combined, ident.1 = "DC_STIM", ident.2 = "DC_CTRL", print.bar = FALSE)
cd16.mono <- FindMarkers(immune.combined, ident.1 = "CD16 Mono_STIM", ident.2 = "CD16 Mono_CTRL", print.bar = FALSE)
nk <- FindMarkers(immune.combined, ident.1 = "NK_STIM", ident.2 = "NK_CTRL", print.bar = FALSE)
eryth <- FindMarkers(immune.combined, ident.1 = "Eryth_STIM", ident.2 = "Eryth_CTRL", print.bar = FALSE)


library("mitch")
cell_type_de<-list("cd14.mono"=cd14.mono,"pdc"=pdc,"cd4.mem"=cd4.mem,"t.act"=t.act,
  "cd4.nat"=cd4.nat,"cd8.t"=cd8.t,"mk.b"=mk,"b.act"=b.act,"b"=b,"dc"=dc,
  "cd16.mono"=cd16.mono,"nk"=nk,"eryth"=eryth)

x<-mitch_import(cell_type_de,DEtype="seurat")

# get gene sets from Reactome
download.file("https://reactome.org/download/current/ReactomePathways.gmt.zip","ReactomePathways.gmt.zip")
unzip("ReactomePathways.gmt.zip")
genesets<-gmt_import("ReactomePathways.gmt")

res<-mitch_calc(x,genesets,resrows=200,bootstraps=1000,priority="confidence")
mitch_plots(res,outfile="scrnaseq.pdf")
mitch_report(res,"scrnaseq.html")
