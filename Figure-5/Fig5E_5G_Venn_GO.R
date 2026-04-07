{
  setwd("~/Library/CloudStorage/OneDrive-Personal/HKU PhD/Figure 5 new/Venn Diagram")
  library(VennDiagram)
  library(grid)
  library(gridSVG)
  library(DESeq2)
  library(tidyverse)
  library(tidyr)
  library(tibble)
  library(dplyr)
  library(stringr)
  library(forcats)
  library(ggnewscale)
  library(clusterProfiler)
  library(AnnotationDbi)
  library(org.Hs.eg.db)
  library(biomaRt)
  library(svglite)
}


# Data import --------------------------------------------------
{
  gene_list_XIST_binding <- c(
    "A3GALT2", "ACKR1", "ADORA1", "ADORA2BP1", "AIM2", "AL359753.1", "AL451006.1", "AL513327.1", "AL583832.1", "AL591704.5", "AL591704.7", "AL591704.9",
    "ALDH4A1", "ARHGEF10L", "ASCL5", "ASTN1", "BTG2", "C1orf132", "C1orf147", "C1orf68", "C1orf95", "C1QB", "CACNA1E", "CADM3", "CADM3-AS1", "CAMK1G",
    "CAMTA1", "CD1A", "CD1B", "CD1D", "CD1E", "CD34", "CDK18", "CHI3L1", "CHIT1", "CNTN2", "CSMD2", "CSMD2-AS1", "ELL2P1", "EPHB2", "FCAMR", "FMOD", "FTLP18",
    "G0S2", "GM140", "GNPAT", "GOLT1A", "GRIK3", "HHAT", "HIVEP3", "HMGB4", "HSD11B1", "HSPD1P14", "IFI16", "IGFN1", "IKBKE", "IPO9-AS1", "ITGB3BP", "IVL",
    "KCNA10", "KCND3", "KCNH1", "KCNH1-IT1", "KISS1", "KPRP", "LAD1", "LAMB3", "LARP7P1", "LCE1B", "LCE2B", "LCE2C", "LCE2D", "LCE3A", "LCE3C", "LCE3E",
    "LCE4A", "LELP1", "LINC00302", "LINC00862", "LINC01136", "LINC01138", "LINC01353", "LINC01527", "LOR", "LRRN2", "MEGF6", "MIR4255", "MIR4260", "MIR4695",
    "MIR488", "MIR551A", "MNDA", "MYBPH", "MYOG", "NAV1", "NFASC", "NPM1P40", "OPTC", "OR10AA1P", "OR10R1P", "OR10R3P", "OR2AQ1P", "OR6K1P", "OR6K2", "OR6K3",
    "OR6K4P", "OR6K5P", "OR6K6", "OR6N1", "OR6N2", "PAX7", "PGLYRP3", "PGLYRP4", "PGM1", "PHC2", "PIGR", "PKP1", "PLEKHA6", "PLXNA2", "PPFIA4", "PRELP", "PRR9",
    "PYHIN1", "RAD1P2", "RASSF5", "RBBP5", "RN7SKP16", "RN7SL130P", "RNA5SP43", "RNA5SP60", "RNA5SP70", "RNA5SP75", "RNU1-8P", "RNU5A-8P", "RNU6-160P", "RNVU1-3",
    "ROR1", "ROR1-AS1", "RP1-140J1.1", "RP1-140J1.4", "RP1-20N18.4", "RP1-272L16.1", "RP1-28O10.1", "RP1-43O17.1", "RP1-43O17.2", "RP1-52J10.9", "RP1-8B22.2",
    "RP11-144L1.4", "RP11-144L1.8", "RP11-203F10.1", "RP11-203F10.5", "RP11-203F10.6", "RP11-23I7.1", "RP11-25B7.1", "RP11-2P2.1", "RP11-328D5.1", "RP11-335O13.7",
    "RP11-335O13.8", "RP11-404O13.1", "RP11-415J8.3", "RP11-415J8.5", "RP11-430C7.4", "RP11-430C7.5", "RP11-435P24.2", "RP11-451O13.1", "RP11-494K3.2",
    "RP11-504A18.1", "RP11-532L16.1", "RP11-532L16.3", "RP11-534L20.5", "RP11-536L3.4", "RP11-567E21.3", "RP11-75I2.3", "RP11-88H9.2", "RP13-279N23.2",
    "RP4-549F15.1", "RP4-614N24.1", "RP5-1180C18.1", "RP5-879K22.1", "RP5-965F6.2", "RPS29P6", "S100A12", "S100A7", "S100A7A", "S100A7L2", "S100A8", "S100A9",
    "SMCP", "SMIM12", "SPRR1A", "SPRR1B", "SPRR2A", "SPRR2B", "SPRR2C", "SPRR2D", "SPRR2E", "SPRR2G", "SPRR3", "SPRR4", "TAS1R2", "TMEM81", "TMEM9", "TNNT2",
    "TNR", "TRAF3IP3", "TTC34", "ADARB2", "ADARB2-AS1", "ADRA2A", "AL138925.1", "AL731561.1", "AL731568.1", "C10orf105", "C10orf11", "C10orf120", "C10orf128",
    "C10orf71", "C10orf71-AS1", "CDH23", "CDH23-AS1", "COL13A1", "COX6CP15", "CPXM2", "CXCL12", "DMBT1", "GDF10", "GDF2", "KCNMA1", "KCNMA1-AS1", "KCNMA1-AS2",
    "KCNMA1-AS3", "KLF6", "LINC00702", "LINC00856", "MTND1P20", "MTND2P15", "NEUROG3", "RBP3", "RP11-184A2.2", "RP11-184A2.3", "RP11-318C4.3", "RP11-343J3.2",
    "RP11-343J3.5", "RP11-343J3.6", "RP11-343J3.7", "RP11-432J9.3", "RP11-432J9.6", "RP11-479A21.1", "RP11-482E14.1", "RP11-89K18.1", "RP11-90J7.4", "RPL5P26",
    "SLC29A3", "SORCS3", "UNC5B", "WDFY4", "ZMIZ1-AS1", "ABCC8", "AC074191.1", "AC084859.1", "ANKK1", "ANO1", "AP000770.1", "AP000797.1", "AP000797.2",
    "AP000797.3", "AP000797.4", "AP000908.1", "AP000997.2", "AP000997.3", "ARNTL", "BTBD10", "CADM1", "CD5", "CD6", "CDON", "CENPUP1", "CEP164", "COPB1",
    "CTC-497E21.3", "CTC-774J1.1", "CTD-2210P24.1", "CTD-2234N14.1", "DRD2", "DSCAML1", "FAR1", "FLI1", "FXYD2", "GRIK4", "HTR3A", "HTR3B", "IFITM9P", "INSC",
    "KCNC1", "KCNJ11", "KIRREL3", "KIRREL3-AS1", "KIRREL3-AS2", "KIRREL3-AS3", "LINC00900", "LINC00958", "MICAL2", "MICALCL", "MIR139", "MIR3167", "MIR4694",
    "MIR6124", "MYEOV", "MYOD1", "NAV2", "NAV2-AS2", "NAV2-AS4", "NAV2-AS5", "NAV2-IT1", "NCAM1", "NCAM1-AS1", "NNMT", "OPCML", "OTOG", "PDE2A", "PKNOX2",
    "PLEKHA7", "POU2AF1", "PVRL1", "RASSF10", "RN7SKP151", "RNA5SP331", "RNA5SP332", "RNU6-933P", "RNU7-49P", "RP11-115C10.1", "RP11-140L24.3", "RP11-159N11.3",
    "RP11-168K9.1", "RP11-168K9.2", "RP11-196E1.3", "RP11-21L19.1", "RP11-231N3.1", "RP11-23B7.4", "RP11-265D17.2", "RP11-27G22.1", "RP11-358H18.2",
    "RP11-358H18.3", "RP11-359E10.1", "RP11-396O20.1", "RP11-396O20.2", "RP11-413N13.1", "RP11-50B3.2", "RP11-531H8.2", "RP11-583F24.3", "RP11-583F24.4",
    "RP11-627G23.1", "RP11-643C9.2", "RP11-64D24.2", "RP11-688I9.2", "RP11-688I9.4", "RP11-688I9.5", "RP11-744N12.3", "RP11-756D7.1", "RP11-839D17.3",
    "RP11-881M11.4", "RP11-881M11.8", "RP11-98J9.1", "RP11-98J9.2", "RP11-98J9.3", "RRAS2", "SDHCP4", "SENCR", "SERGEF", "SLC15A3", "SPON1", "ST3GAL4", "TEAD1",
    "TENM4", "TMEM109", "TMEM132A", "TMPRSS4", "TMPRSS4-AS1", "TRIM29", "TTC12", "USH1C", "VPS37C", "ZBTB16", "AC007618.3", "ANO2", "BTBD11", "C12orf80",
    "CACNA1C", "CACNA1C-AS1", "CACNA1C-AS2", "CACNA1C-AS3", "CCND2-AS1", "CRACR2A", "DYNLL1P4", "FLJ12825", "GALNT8", "GLULP5", "HAUS8P1", "KCNA1", "KCNA5",
    "KCNA6", "KRT7", "KRT80", "KRT81", "KRT86", "KRT87P", "KSR2", "LINC00592", "LINC01234", "METTL7AP1", "MYO1H", "NOS1", "NTF3", "PARP11", "PDE1B", "RBM19",
    "RIMBP2", "RN7SKP216", "RNA5SP371", "RNU6-174P", "RP11-100F15.2", "RP11-1028N23.2", "RP11-1028N23.3", "RP11-1028N23.4", "RP11-1038A11.1", "RP11-116D17.2",
    "RP11-139B1.1", "RP11-143E21.2", "RP11-143E21.3", "RP11-234B24.4", "RP11-364C11.3", "RP11-429A20.2", "RP11-429A20.3", "RP11-429A20.4", "RP11-664D1.1",
    "RP11-834C11.4", "RP11-834C11.5", "RP11-845M18.7", "RP3-377H17.2", "RP3-416H24.1", "RPH3A", "SRRM4", "TBX5", "TBX5-AS1", "TMEM132B", "UBA52P7", "AC007375.1",
    "AE000661.37", "AL133167.1", "CMTM5", "CTD-2201G16.1", "EFS", "ESRRB", "GPATCH2L", "IL25", "MIR208A", "MIR208B", "MYH6", "MYH7", "NRXN3", "RP11-124D2.3",
    "RP11-185P18.2", "RP11-187O7.3", "RP11-361H10.3", "RP11-361H10.5", "RP11-463C8.4", "RP11-516J2.1", "SLC22A17", "TMEM63C", "TRAC", "TRDC", "TRDD3", "TRDJ1",
    "TRDJ2", "TRDJ3", "TRDJ4", "TUNAR", "AC007950.2", "C15orf59", "CCDC33", "CORO2B", "CTD-2071N1.1", "CTD-2071N1.2", "CTD-2632K10.1", "HCN4", "ISLR", "ITGA11",
    "LINGO1", "MED28P6", "NTRK3", "PCAT29", "RP11-208K4.1", "RP11-272D12.1", "RP11-279F6.1", "RP11-60L3.1", "RP11-60L3.2", "RP11-60L3.3", "RP11-665J16.1",
    "RP11-8P11.3", "RP11-8P11.4", "STRA6", "THSD4", "TLE3", "USP3", "AC007339.1", "AC007339.2", "ADAM3B", "ADAMTS18", "CES1P2", "CHP2", "CLEC3A", "CNGB1",
    "CTA-345G4.1", "CTC-391G2.1", "CTD-2015G9.1", "CTD-2015G9.2", "CTD-2385L22.2", "CTD-2600O9.1", "GNAO1", "GPR114", "GPR56", "GPR97", "GRIN2A", "KIFC3",
    "KRT8P22", "MIR1273H", "MIR6772", "MMP15", "PDILT", "PRKCB", "RBFOX1", "RP11-158I3.1", "RP11-26L20.3", "RP11-281J9.2", "RP11-297M9.1", "RP11-358L22.2",
    "RP11-405F3.5", "RP11-420N3.2", "RP11-420N3.3", "RP11-58A18.2", "RP11-805I24.1", "RP11-895K13.2", "SCNN1G", "SLC6A2", "UMOD", "USB1", "VAT1L", "WWOX",
    "ZNF423", "AC019349.5", "AC125421.1", "ASIC2", "BTBD17", "CCL2", "CCL7", "CTC-304I17.1", "CTC-304I17.3", "CTC-304I17.4", "CTD-2514K5.4", "CTD-2532D12.4",
    "CTD-2532D12.5", "CTD-2582D11.1", "DNAI2", "GPR142", "KIF19", "KRT13", "KRT14", "KRT15", "KRT16", "KRT17", "KRT19", "KRT32", "KRT35", "KRT36", "KRT38",
    "KRT42P", "KRT43P", "KRT9", "LINC00469", "LINC00974", "MIR6510", "PITPNM3", "RBFOX3", "RNA5SP438", "RNU2-32P", "RP11-101O21.1", "RP11-17M24.1", "RP11-17M24.2",
    "RP11-398J5.1", "RP11-449L23.2", "RP11-449L23.3", "RP11-647F2.2", "SDK2", "SLC13A5", "SOX9-AS1", "TTYH2", "ALPK2", "BCL2", "CPLX4", "CTD-2130O13.1",
    "MIR122", "MIR3591", "NEDD4L", "RNF165", "RNU4-17P", "RP11-1151B14.1", "RP11-1151B14.5", "RP11-627G18.4", "RP11-749H17.2", "RP11-795H16.2", "RP11-795H16.3",
    "RPS3AP49", "SLC14A2", "ST8SIA5", "AC011513.4", "ADAMTS10", "CACNA1A", "CEA", "CEACAM3", "CEACAM5", "CEACAM6", "CEACAM7", "CTC-525D6.1", "CTC-525D6.2",
    "CTC-525D6.3", "CTD-3187F8.11", "CTD-3187F8.12", "CTD-3187F8.14", "CTD-3187F8.2", "CTD-3187F8.7", "MUC16", "OR2Z1", "RPL23AP78", "SIGLEC17P", "SIGLEC20P",
    "SIGLECL1", "VSTM2B", "ZNF536", "AC007395.3", "AC007395.4", "AC007743.1", "AC016764.1", "AC017084.1", "AC018866.1", "AC022201.4", "AC022201.5", "AC073257.2",
    "ACSL3", "ADD2", "ALK", "BRD7P6", "CCDC85A", "CD207", "CLEC4F", "FAM136A", "FIGLA", "GLI2", "HMGN2P21", "INHBB", "LINC01101", "MIR6809", "MYT1L", "RP11-297J22.1",
    "RP11-481J13.1", "RP11-482H16.1", "RUFY4", "SH3RF3", "TGFA", "TGFA-IT1", "TNS1", "AL049647.1", "AL121588.1", "AL121756.1", "AL139429.1", "AL354984.1",
    "AL354984.3", "BPIFB2", "BPIFB3", "BPIFB4", "BPIFB6", "CDH4", "CTCFL", "DEFB129", "ELMO2", "EYA2", "HDHD1P3", "KCNB1", "LINC00659", "LINC01272", "MIR4532",
    "NTSR1", "PHACTR3", "PI3", "PMEPA1", "PTPRT", "RNU7-173P", "RP11-103J8.1", "RP11-429E11.2", "RP11-560A15.3", "RP13-379L11.2", "RP13-379L11.4",
    "RP3-453C12.14", "RP3-461P17.9", "RP4-705O1.1", "RP4-719C8.1", "RP5-827E24.1", "RPL5P2", "RPS2P7", "SDC4", "SLC12A5", "SLC24A3", "SLX4IP", "SPINT3", "SUN5",
    "TOX2", "WFDC10A", "WFDC11", "WFDC2", "WFDC8", "WFDC9", "ZBP1", "ZNF663P", "AP000261.1", "AP000265.1", "AP000266.7", "AP001615.9", "C21orf128", "HUNK",
    "HUNK-AS1", "LINC00111", "LINC00112", "LINC00159", "LINC00479", "MIR6814", "MIS18A", "MIS18A-AS1", "MRAP", "MX1", "RIPK4", "TPT1P1", "UMODL1", "URB1",
    "AC244157.1", "CACNA1I", "CTA-929C8.6", "CTA-929C8.7", "CTA-929C8.8", "ENTHD1", "GRAP2", "IGLL5", "IGLV2-11", "LINC01422", "LL22NC03-84E4.13", "MIATNB",
    "MPPED1", "MYO18B", "RP1-172B20.6", "RP1-205F14P.1", "RP1-40G4P.1", "RP11-46E17.6", "AC080008.1", "AC117401.1", "ATP2B2", "ATP2B2-IT1", "CRYGS", "DGKG",
    "EIF2B5", "EPHB1", "EPHB3", "ETV5", "ETV5-AS1", "KALRN", "KY", "MIR5002", "MIR885", "PLXNA1", "RP11-407B7.1", "RP11-443P15.2", "RP11-48F14.1", "RP11-48F14.2",
    "RP11-78H24.1", "RPL7P15", "SRGAP3", "TBCCD1", "AC080003.1", "AFAP1", "AFAP1-AS1", "C4orf50", "JAKMIP1", "MIR4274", "MIR4798", "PPP2R2C", "PSAPL1", "RBPJ",
    "RN7SKP113", "RNF150", "RP11-5K16.3", "RPL14P3", "SCOC", "SCOC-AS1", "SORCS2", "TBC1D9", "TNRC18P1", "AC005592.1", "AC005592.2", "AC005592.3", "AC011343.1",
    "ARHGAP26", "CTB-105N12.2", "CTB-178M22.1", "CTB-178M22.2", "CTB-37A13.1", "CTC-558O2.1", "DOCK2", "FAM196B", "FGF1", "HMP19", "KCNIP1", "MIR218-2", "MIR378E",
    "MIR585", "RP11-619L12.3", "RP11-619L12.4", "RPS12P10", "SLIT3", "SPRY4", "SPRY4-IT1_1", "SPRY4-IT1_2", "TENM2", "WWC1", "C6orf223", "DAAM2", "GLP1R",
    "KCNK16", "KCNK17", "KIF6", "LINC00951", "LRFN2", "MOCS1", "MTRF1L", "PHACTR1", "RGS17", "RNU6-250P", "RP1-137F1.3", "RP1-202I21.3", "RP11-121P10.1",
    "RP11-256G5.1", "RP11-552E20.1", "RP11-552E20.4", "RP11-570K4.1", "RP11-61I13.3", "RP3-462C17.1", "RP5-1120P11.1", "TDRG1", "TRERF1", "AC004540.5",
    "AC004947.2", "ADCYAP1R1", "AUTS2", "C7orf71", "GIMAP1", "GIMAP5", "HOXA10", "HOXA10-AS", "HOXA10-HOXA9", "HOXA9", "KIAA0087", "PLXNA4", "RP5-1007F24.1",
    "RP5-978I12.1", "SKAP2", "ADRA1A", "CHRNA2", "CLU", "EBAG9", "EPHX2", "GFRA2", "GULOP", "IMPA1P", "MIR6842", "NIPA2P4", "PAG1", "PTK2B", "RNU7-85P",
    "RP11-1149M10.1", "RP11-1149M10.2", "RP11-172E10.1", "RP11-26J3.1", "RP11-27N21.3", "RP11-486M23.2", "RP11-486M23.3", "RPS26P34", "SLC10A5P1", "SYBU",
    "AKAP2", "AL161784.1", "AL353141.1", "AL357936.1", "AL365274.1", "AL390240.1", "AL390240.2", "C5", "CNTRL", "COL27A1", "COL5A1", "DAB2IP", "GGTA1P", "GSN",
    "GSN-AS1", "MIR455", "PALM2", "PALM2-AKAP2", "PAPPA", "PAPPA-AS2", "RAB14", "RN7SL181P", "RN7SL187P", "RNY4P18", "RP11-162D16.2", "RP11-244O19.1",
    "RP11-402G3.3", "RP11-402G3.4", "RP11-406O23.2", "RP11-428F18.2", "RP11-45A16.4", "RP11-477J21.2", "RP11-477J21.6", "RP11-524G24.2", "STOM", "TNFSF15",
    "TTLL11", "TTLL11-IT1", "ZNF618"
  )
  genes_binding <- read.csv("~/Library/CloudStorage/OneDrive-Personal/HKU PhD/Figure 5 new/Venn Diagram/target_genes_avg_v2.csv")
}

{
  # 1) Make a clean ENSG key (only for rows that are ENSG)
  genes_binding2 <- genes_binding %>%
    mutate(
      ensembl_gene_id = ifelse(grepl("^ENSG", x), sub("\\..*$", "", x), NA_character_)
    )
  
  # 2) Old way: org.Hs.eg.db mapping
  genes_binding2 <- genes_binding2 %>%
    mutate(
      gene_symbol_org = mapIds(
        org.Hs.eg.db,
        keys = ensembl_gene_id,
        column = "SYMBOL",
        keytype = "ENSEMBL",
        multiVals = "first"
      )
    )
  
  # 3) biomaRt mapping table
  mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  
  bm <- getBM(
    attributes = c("ensembl_gene_id", "hgnc_symbol"),
    filters    = "ensembl_gene_id",
    values     = unique(na.omit(genes_binding2$ensembl_gene_id)),
    mart       = mart
  ) %>%
    filter(hgnc_symbol != "") %>%
    distinct(ensembl_gene_id, .keep_all = TRUE)
  
  # 4) Join biomaRt symbols back (as a second column)
  genes_binding2 <- genes_binding2 %>%
    left_join(bm, by = "ensembl_gene_id") %>%
    rename(gene_symbol_biomart = hgnc_symbol)

  # quick checks
  head(genes_binding2)
  table(is.na(genes_binding2$gene_symbol_org), is.na(genes_binding2$gene_symbol_biomart))
}
{
  genes_binding2 <- genes_binding2 %>%
    mutate(
      gene_symbol_org = as.character(gene_symbol_org),
      gene_symbol_final = coalesce(
        gene_symbol_biomart,
        gene_symbol_org
      )
    )
  
  genes_binding_clean <- genes_binding2 %>%
    mutate(
      gene_symbol_final = case_when(
        # case 1: already a gene symbol (not ENSG)
        !grepl("^ENSG", x) ~ x,
        
        # case 2: ENSG but mapped by either method
        grepl("^ENSG", x) & !is.na(gene_symbol_final) ~ gene_symbol_final,
        
        # case 3: ENSG but unmapped by both
        TRUE ~ NA_character_
      )
    ) %>%
    # discard ENSG that could not be annotated
    filter(!is.na(gene_symbol_final)) %>%
    # remove middle columns
    dplyr::select(x, gene_symbol_final)
  
  annotated_changed <- genes_binding_clean %>%
    filter(x != gene_symbol_final)
}


# Data cleaning -----------------------------------------------------------
{
  set_xist <- gene_list_XIST_binding %>%
    as.character() %>%
    trimws() %>%
    toupper() %>%
    unique()
  
  set_binding <- genes_binding_clean$gene_symbol_final %>%
    as.character() %>%
    trimws() %>%
    toupper() %>%
    unique()
}
{
  # get chromosome info for genes in set_binding
  mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  
  binding_chr <- getBM(
    attributes = c(
      "hgnc_symbol",
      "chromosome_name",
      "start_position",
      "end_position",
      "strand"
    ),
    filters = "hgnc_symbol",
    values = unique(set_binding),
    mart = mart
  ) %>%
    filter(hgnc_symbol != "") %>%
    distinct(hgnc_symbol, .keep_all = TRUE)
  
  # make names consistent with your set_binding
  binding_chr <- binding_chr %>%
    mutate(
      hgnc_symbol = toupper(trimws(hgnc_symbol))
    )
  
  # join back to full set_binding list to see mapped / unmapped genes
  binding_chr_full <- tibble(gene_symbol_final = unique(set_binding)) %>%
    left_join(binding_chr, by = c("gene_symbol_final" = "hgnc_symbol"))
  
  # quick checks
  head(binding_chr_full)
  table(is.na(binding_chr_full$chromosome_name))
}
{
  set_binding_autosome <- binding_chr_full %>%
    filter(
      !chromosome_name %in% c("X", "Y"),
      gene_symbol_final != "CETN2"
    ) %>%
    pull(gene_symbol_final) %>%
    unique()
  
  # quick checks
  length(set_binding)
  length(set_binding_autosome)
  
  removed_genes <- set_binding[set_binding %in% c(
    binding_chr_full %>%
      filter(chromosome_name %in% c("X", "Y")) %>%
      pull(gene_symbol_final),
    "CETN2"
  )] %>%
    unique() %>%
    sort()
  
  removed_genes
}


# Venn plot ---------------------------------------------------------------
venn_plot_xist_binding <- venn.diagram(
  x = list(
    XIST_binding = set_xist,
    genes_binding = set_binding_autosome
  ),
  category.names = c("XIST binding", "Genes binding"),
  filename = NULL,
  output = TRUE,
  fill = c("#FFCDBA", "#FFC566"), alpha = 0.5,
  lwd = 2, lty = "solid",
  cex = 1.8, fontface = "bold", fontfamily = "sans",
  cat.cex = 2.4, cat.fontface = "bold", cat.fontfamily = "sans",
  cat.pos = c(-40, 40), cat.dist = c(0.08, 0.08),
  print.mode = "raw",
  margin = 0.1
)

# Save SVG
{
  svglite("fig5_venn_XISTbinding_vs_genesbinding_edit.svg", 
          bg = "transparent", 
          width = 8, height = 6)
  grid.newpage()
  pushViewport(viewport(width = unit(1, "npc"), 
                        height = unit(0.7, "npc")))
  grid.draw(venn_plot_xist_binding)
  popViewport()
  dev.off()
}

# Save the gene lists behind the Venn
{
  XIST_only <- setdiff(set_xist, set_binding)
  binding_only <- setdiff(set_binding, set_xist)
  XIST_and_binding <- intersect(set_xist, set_binding)
  
  max_len <- max(length(XIST_only), length(binding_only), length(XIST_and_binding))
  
  venn_genes <- data.frame(
    XIST_only = c(XIST_only, rep(NA, max_len - length(XIST_only))),
    dataset_only = c(binding_only, rep(NA, max_len - length(binding_only))),
    overlap = c(XIST_and_binding, rep(NA, max_len - length(XIST_and_binding)))
  )
}

write.csv(venn_genes, "fig5_venn_XISTbinding_vs_genesbinding_edit.csv", row.names = FALSE)


# GO ----------------------------------------------------------------------
{
  # 1) pull overlap genes from the dataframe
  overlap_genes <- venn_genes$overlap %>%
    as.character() %>%
    trimws() %>%
    toupper() %>%              # keep consistent with how we made the Venn sets
    .[!is.na(.) & . != ""] %>%
    unique()
  
  # 2) enrichGO
  GO_overlap <- enrichGO(
    gene = overlap_genes,
    OrgDb = org.Hs.eg.db,
    keyType = "SYMBOL",
    ont = "ALL",
    pAdjustMethod = "BH"
  )
  
  # 3) make dataframe + your extra columns
  GO_overlap_df <- as.data.frame(GO_overlap) %>%
    mutate(
      Description_raw = Description,                  # keep original
      FoldEnrichment = as.numeric(FoldEnrichment),
      log2_OE = log2(FoldEnrichment),
      neglog10p = -log10(pvalue)
    )
  
}

# 4) GO plot
ggplot(
  GO_overlap_df %>% slice_max(order_by = neglog10p, n = 11),
  aes(
    x = neglog10p,
    y = reorder(Description, neglog10p),
    fill = log2_OE
    )
  ) +
  geom_col(color = "black") +
  labs(x = "-log10(p-value)", y = NULL, fill = "log2(O/E)") +
  scale_fill_gradient(
    low = "#E6E0DC",
    high = "salmon",
    labels = scales::number_format(accuracy = 0.1),
    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")
  ) +
  theme_classic() +
  theme(
    axis.text.y = element_text(size = 13),
    axis.title.x = element_text(size = 13),
    axis.text.x = element_text(size = 13),
    legend.title = element_text(size = 13),
    legend.text = element_text(size = 13)
  )

{
  # 5) select terms
  patterns_keep <- c(
    "cell fate specification",
    "motor neuron apoptotic",
    "skeletal muscle organ development",
    "synaptic membrane",
    "pore complex"
  )
  
  GO_overlap_df_sel <- GO_overlap_df %>%
    filter(stringr::str_detect(Description_raw, paste(patterns_keep, collapse = "|"))) %>%
    mutate(Description = stringr::str_wrap(Description_raw, width = 30))   # wrap AFTER filtering
}

ggplot(
  GO_overlap_df_sel %>% slice_max(order_by = neglog10p, n = 5),
  aes(
    x = neglog10p,
    y = reorder(Description, neglog10p),
    fill = log2_OE
    )
  ) +
  geom_col(color = "black") +
  labs(x = "-log10 (p-value)", y = NULL, fill = "log2 (O/E)") +
  scale_fill_gradient(
    low = "#EFE6D8",
    high = "#D55E00", 
    labels = scales::number_format(accuracy = 0.1),
    guide = guide_colorbar(frame.colour = "black", ticks.colour = "black")
  ) +
  theme_classic() +
  theme(
    axis.line = element_line(linewidth = 1.4, colour = "black"),
    axis.text.y = element_text(size = 24, face = "bold"),
    axis.title.x = element_text(size = 28, face = "bold"),
    axis.text.x = element_text(size = 18),
    legend.title = element_text(size = 24),
    legend.text = element_text(size = 24)
  )

ggsave(
  "fig5_BrainOrganoid_vs_Public_GO_3.tiff",
  plot = last_plot(),
  width = 10,    
  height = 5,      
  units = "in",
  dpi = 600,
  limitsize = FALSE,
  compression = "lzw"
)

