# 06-4_gene_sets.R - Builds the 14 prespecified gene sets and the KEGG/GO sets (before any gene-set result is read). See README.md.
#      KEGG_cfa04710 Circadian rhythm / KEGG_cfa04713 Circadian entrainment / GO_0007623 circadian rhythm / GO_0042752 regulation of circadian rhythm
#      GWAS_chronotype_Jones2019 (PMID 30696823; chronotype + morning person)
#      GWAS_sleepdur_Dashti2019 (PMID 30846698; sleep duration + short + long sleep)
#      GWAS_insomnia_Jansen2019 (PMID 30804565; insomnia)
#      GWAS_cognition_Davies2018 (PMID 29844566; general cognitive ability)

library(data.table)
options(width = 220)

base_dir <- "/data/02_ByProject/DogAgingProject"
qc_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023")   
out_prefix <- file.path(qc_dir, "06-4")
ref_dir <- paste0(out_prefix, "_ref")
res63_dir <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023/06-3_magma")
hs_info_63 <- file.path(base_dir, "GeneticData/analysis/04_GWAS_2023/06-3_ref/Homo_sapiens.gene_info.gz")
v2022_dir <- file.path(base_dir, "GeneticData/analysis/03_GWAS_Pure")
ref22 <- file.path(v2022_dir, "06-1_ref")
sets2022_file <- file.path(v2022_dir, "06-1_sets_explore.txt")
setinfo2022_file <- file.path(v2022_dir, "06-1_set_info.csv")
ortho_file <- file.path(ref22, "gene_orthologs_human_dog.tsv")
g2g_hum_file <- file.path(ref22, "gene2go_human.tsv")
g2g_dog_file <- file.path(ref22, "gene2go_dog.tsv")
obo_file <- file.path(ref22, "go-basic.obo")
gene_loc_file <- file.path(v2022_dir, "06-1_gene_loc.txt")
ft_file <- file.path(ref22, "GCF_000002285.3_CanFam3.1_feature_table.txt.gz")
hs_info_url <- "https://ftp.ncbi.nlm.nih.gov/gene/DATA/GENE_INFO/Mammalia/Homo_sapiens.gene_info.gz"
gwas_url <- "https://ftp.ebi.ac.uk/pub/databases/gwas/releases/latest/gwas-catalog-associations_ontology-annotated-full.zip"
palmer_url <- "https://www.ebi.ac.uk/europepmc/webservices/rest/PMC7906136/supplementaryFiles"
gse_url <- "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE304nnn/GSE304125/suppl/GSE304125_Cortex_EV_NormCounts.xlsx"
q_age <- 0.05
univ_files <- c(w10 = file.path(res63_dir, "main.genes.out"), w50 = file.path(res63_dir, "main_w50.genes.out"))
for (f in c(sets2022_file, setinfo2022_file, ortho_file, g2g_hum_file, g2g_dog_file, obo_file, gene_loc_file, ft_file, univ_files))
  if (!file.exists(f)) stop("파일 없음: ", f)
dir.create(ref_dir, showWarnings = FALSE, recursive = TRUE)

human_taxon <- 9606; dog_taxon <- 9615
p_gwas <- 5e-8
p_gwas_sens <- 1e-6
set_min <- 10; aux_max <- 500
kegg_exclude <- "^cfa0(5|7)"

circ_sets <- c("KEGG_cfa04710", "KEGG_cfa04713", "GO_0007623", "GO_0042752")
splice_sets <- c("GO_0000398", "GO_0048024", "KEGG_cfa03040")
gwas_spec <- data.table(
  set = c("GWAS_chronotype_Jones2019", "GWAS_sleepdur_Dashti2019", "GWAS_insomnia_Jansen2019",
          "GWAS_AD_Bellenguez2022", "GWAS_cognition_Davies2018"),
  family = c(rep("sleep_gwas", 3), rep("AD_cognition_gwas", 2)),
  pmid = c("30696823", "30846698", "30804565", "35379992", "29844566"),
  traits = list(c("Chronotype", "Morning person"),
                c("Sleep duration", "Sleep duration (short sleep)", "Sleep duration (long sleep)"),
                "Insomnia", "Alzheimer's disease", "General cognitive ability"),
  name = c("Chronotype GWAS genes (Jones 2019)", "Sleep duration GWAS genes (Dashti 2019)", "Insomnia GWAS genes (Jansen 2019)",
           "Alzheimer's disease GWAS genes (Bellenguez 2022)", "General cognitive ability GWAS genes (Davies 2018)"))
neuro_splice <- c("NOVA1", "NOVA2", "RBFOX1", "RBFOX2", "RBFOX3", "PTBP1", "PTBP2", "SRRM3", "SRRM4",
                  "ELAVL2", "ELAVL3", "ELAVL4", "CELF1", "CELF2", "CELF3", "CELF4", "CELF5", "CELF6",
                  "MBNL1", "MBNL2", "KHDRBS1", "KHDRBS2", "KHDRBS3", "QKI")
primary_order <- c(circ_sets, gwas_spec$set, splice_sets, "SPLICE_neuronal_regulators", "AGE_brain_Palmer2021")
family_of <- c(setNames(rep("circadian", 4), circ_sets), setNames(gwas_spec$family, gwas_spec$set),
               setNames(rep("splicing", 4), c(splice_sets, "SPLICE_neuronal_regulators")), AGE_brain_Palmer2021 = "brain_aging_expr")

gloc <- fread(gene_loc_file, header = FALSE, col.names = c("GENE", "CHR", "START", "STOP", "STRAND", "SYMBOL"))
gloc[, GENE := as.character(GENE)]
univ <- lapply(univ_files, function(f) as.character(fread(f, skip = "GENE", select = "GENE")$GENE))
cat(sprintf("개 유전자 %s개, 06-3 분석 유전자 ±10 kb %s개 · ±50 kb %s개\n",
            format(nrow(gloc), big.mark = ","), format(length(univ$w10), big.mark = ","), format(length(univ$w50), big.mark = ",")))

hs_file <- if (file.exists(hs_info_63)) hs_info_63 else file.path(ref_dir, "Homo_sapiens.gene_info.gz")
if (!file.exists(hs_file)) { message("사람 gene_info 내려받기"); download.file(hs_info_url, hs_file, mode = "wb", quiet = TRUE) }
hs <- fread(cmd = sprintf("gzip -dc '%s'", hs_file), sep = "\t", quote = "",
            select = c("GeneID", "Symbol", "Synonyms", "type_of_gene", "description", "Full_name_from_nomenclature_authority"))
hs[, GeneID := as.character(GeneID)]
syn <- hs[Synonyms != "-", .(alias = unlist(strsplit(Synonyms, "|", fixed = TRUE))), by = GeneID]
syn <- syn[!alias %in% hs$Symbol][, if (uniqueN(GeneID) == 1) .SD[1], by = alias]
sym2id <- function(s) {
  id <- hs$GeneID[match(s, hs$Symbol)]; how <- ifelse(is.na(id), NA_character_, "symbol")
  j <- is.na(id); id[j] <- syn$GeneID[match(s[j], syn$alias)]; how[j & !is.na(id)] <- "synonym"
  list(id = id, how = how)
}

orth <- fread(ortho_file, sep = "\t")
setnames(orth, c("tax_id", "GeneID", "relationship", "Other_tax_id", "Other_GeneID"))
orth <- unique(rbind(orth[tax_id == human_taxon, .(human = as.character(GeneID), dog = as.character(Other_GeneID))],
                     orth[tax_id == dog_taxon, .(human = as.character(Other_GeneID), dog = as.character(GeneID))]))
n_orth_raw <- nrow(orth)
orth <- orth[dog %in% gloc$GENE][, if (.N == 1) .SD, by = dog][, if (.N == 1) .SD, by = human]
orth[, by := "NCBI"]
same <- merge(hs[type_of_gene == "protein-coding", .(human = GeneID, SYMBOL = Symbol)], gloc[, .(dog = GENE, SYMBOL)], by = "SYMBOL")
same <- same[, if (.N == 1) .SD, by = SYMBOL][!human %in% orth$human & !dog %in% orth$dog]
fam_expand <- "^(ZNF|ZSCAN|ZKSCAN|ZFP|ZIK)[0-9]|^OR[0-9]"
ft <- fread(cmd = sprintf("gzip -dc '%s' | awk -F'\\t' '$1==\"mRNA\"' | cut -f14,16", ft_file), sep = "\t", quote = "",
            header = FALSE, col.names = c("name", "dog"), colClasses = "character")
ft <- unique(ft[, name := tolower(sub(", transcript variant .*$", "", name))][dog %in% gloc$GENE & name != ""])
hn <- hs[type_of_gene == "protein-coding" & Full_name_from_nomenclature_authority != "-" &
         !grepl(fam_expand, Symbol), .(human = GeneID, name = tolower(Full_name_from_nomenclature_authority))][name != ""]
fulln <- merge(ft[, if (uniqueN(dog) == 1) .SD, by = name], hn[, if (uniqueN(human) == 1) .SD, by = name], by = "name")
fulln <- fulln[, if (.N == 1) .SD, by = dog][, if (.N == 1) .SD, by = human]
fulln <- fulln[!human %in% c(orth$human, same$human) & !dog %in% c(orth$dog, same$dog)]
omap <- rbind(orth[, .(human, dog, by)], same[, .(human, dog, by = "same_symbol")], fulln[, .(human, dog, by = "same_fullname")])
omap[, `:=`(hs_symbol = hs$Symbol[match(human, hs$GeneID)], dog_symbol = gloc$SYMBOL[match(dog, gloc$GENE)])]
cat(sprintf("사람 → 개: NCBI 표 %s쌍 중 CanFam3.1 유전자를 가리키는 1:1 %s쌍 + 같은 이름 %s쌍 + 같은 전체 이름 %s쌍 = 개 유전자 %s개\n",
            format(n_orth_raw, big.mark = ","), format(nrow(orth), big.mark = ","), format(nrow(same), big.mark = ","),
            format(nrow(fulln), big.mark = ","), format(uniqueN(omap$dog), big.mark = ",")))
if (anyDuplicated(omap$dog) || anyDuplicated(omap$human)) stop("사람 ↔ 개 대응이 1:1이 아님")
to_dog <- function(hid) omap[human %in% hid, unique(dog)]

read_g2g <- function(f) {
  d <- fread(f, sep = "\t", quote = "")
  setnames(d, c("tax_id", "GeneID", "GO_ID", "Evidence", "Qualifier", "GO_term", "PubMed", "Category"))
  unique(d[Category == "Process" & !grepl("NOT", Qualifier), .(GeneID = as.character(GeneID), GO_ID)])
}
g2g_dog <- read_g2g(g2g_dog_file)[GeneID %in% gloc$GENE]
g2g_hum <- merge(read_g2g(g2g_hum_file), omap[, .(human, dog)], by.x = "GeneID", by.y = "human")[, .(GeneID = dog, GO_ID)]
g2g <- unique(rbind(g2g_dog, g2g_hum))

obo <- readLines(obo_file)
blk <- cumsum(obo == "[Term]") * (cumsum(obo == "[Term]") > 0)
blk[cumsum(grepl("^\\[Typedef\\]", obo)) > 0] <- 0
ob <- data.table(b = blk, line = obo)[b > 0]
terms <- ob[, .(id = sub("^id: ", "", line[grepl("^id: ", line)][1]),
                name = sub("^name: ", "", line[grepl("^name: ", line)][1]),
                namespace = sub("^namespace: ", "", line[grepl("^namespace: ", line)][1]),
                obsolete = any(line == "is_obsolete: true")), by = b]
parents <- ob[grepl("^is_a: GO:|^relationship: part_of GO:", line), .(b, parent = sub(".*(GO:[0-9]{7}).*", "\\1", line))]
parents <- merge(parents, terms[, .(b, child = id)], by = "b")[, .(child, parent)]
par_list <- split(parents$parent, parents$child)
anc_cache <- new.env(hash = TRUE)
ancestors <- function(t) {
  if (!is.null(r <- anc_cache[[t]])) return(r)
  p <- par_list[[t]]
  r <- if (is.null(p)) character(0) else unique(c(p, unlist(lapply(p, ancestors))))
  assign(t, r, envir = anc_cache); r
}
term_anc <- rbindlist(lapply(unique(g2g$GO_ID), function(t) data.table(GO_ID = t, term = c(t, ancestors(t)))))
go_prop <- unique(merge(g2g, term_anc, by = "GO_ID", allow.cartesian = TRUE)[, .(GeneID, term)])
go_prop <- merge(go_prop, terms[obsolete == FALSE & namespace == "biological_process", .(term = id, name)], by = "term")
go_sets <- go_prop[, .(genes = list(unique(GeneID))), by = .(term, name)]
go_sets[, set := paste0("GO_", sub("GO:", "", term))]
go_list <- setNames(go_sets$genes, go_sets$set)
cat(sprintf("GO: 개 자체 %s + 사람에서 옮긴 %s 연결 → 유전자 %s개, BP term %s개 (상위 term까지 채운 뒤)\n",
            format(nrow(g2g_dog), big.mark = ","), format(nrow(g2g_hum), big.mark = ","), format(uniqueN(g2g$GeneID), big.mark = ","),
            format(length(go_list), big.mark = ",")))

gw_file <- file.path(ref_dir, "gwas_catalog_5studies.tsv")
if (!file.exists(gw_file)) {
  zf <- file.path(ref_dir, "gwas_catalog_assoc.zip")
  if (!file.exists(zf) || system2("unzip", c("-tq", shQuote(zf)), stdout = FALSE, stderr = FALSE) != 0) {
    message("GWAS Catalog 연관 표 내려받기 (약 70 MB): ", gwas_url)
    st <- system2("curl", c("-sL", "--fail", "--retry", "3", "-o", shQuote(zf), shQuote(gwas_url)))
    if (st != 0 || !file.exists(zf) || file.size(zf) < 1e7) stop("GWAS Catalog 내려받기 실패")
  }
  inner <- system2("unzip", c("-Z1", shQuote(zf)), stdout = TRUE)[1]
  g <- fread(cmd = sprintf("unzip -p '%s' '%s'", zf, inner), sep = "\t", quote = "", colClasses = "character")
  g <- g[PUBMEDID %in% gwas_spec$pmid]
  fwrite(g, gw_file, sep = "\t")
  writeLines(c(paste("downloaded:", format(file.mtime(zf), "%Y-%m-%d %H:%M %Z")), paste("url:", gwas_url), paste("file in zip:", inner)),
             file.path(ref_dir, "gwas_catalog_version.txt"))
  unlink(zf)
}
gw_ver <- readLines(file.path(ref_dir, "gwas_catalog_version.txt"))[1]
g <- fread(gw_file, sep = "\t", quote = "", colClasses = "character")
need <- c("PUBMEDID", "DISEASE/TRAIT", "REPORTED GENE(S)", "MAPPED_GENE", "P-VALUE", "SNPS", "STUDY ACCESSION")
if (length(miss <- setdiff(need, names(g)))) stop("GWAS Catalog 열 없음: ", paste(miss, collapse = ", "))
g[, p := suppressWarnings(as.numeric(`P-VALUE`))]

not_gene <- c("", "NR", "nr", "intergenic", "Intergenic", "NA")
split_genes <- function(x, mapped) {
  s <- if (mapped) strsplit(x, " - |, |,|; |;| x ") else strsplit(x, ", |,|; |;")
  v <- trimws(gsub('"', "", unlist(s)))
  v[!v %in% not_gene & !grepl("^LOC|^LINC|^MIR[0-9]|^RNU|^SNOR", v)]
}
gene_sources <- list(reported    = list(mapped = FALSE, p = p_gwas),
                     mapped      = list(mapped = TRUE,  p = p_gwas),
                     mapped_p1e6 = list(mapped = TRUE,  p = p_gwas_sens))
gw_genes <- rbindlist(lapply(seq_len(nrow(gwas_spec)), function(i) {
  sp <- gwas_spec[i]
  d0 <- g[PUBMEDID == sp$pmid & `DISEASE/TRAIT` %in% sp$traits[[1]] & !is.na(p)]
  if (!nrow(d0[p <= p_gwas])) stop("GWAS 연관이 없음: ", sp$set, " — 형질 이름 확인: ", paste(unique(g[PUBMEDID == sp$pmid, `DISEASE/TRAIT`]), collapse = " / "))
  rbindlist(lapply(names(gene_sources), function(src) {
    m <- gene_sources[[src]]$mapped
    d <- d0[p <= gene_sources[[src]]$p]
    col <- if (m) "MAPPED_GENE" else "REPORTED GENE(S)"
    x <- d[, .(hs_symbol = split_genes(get(col), m)), by = .(SNPS, p)]
    if (!nrow(x)) return(NULL)
    x[, .(n_assoc = uniqueN(SNPS), min_p = min(p)), by = hs_symbol][, `:=`(set = sp$set, gene_source = src)]
  }))
}))
gwas_add <- data.table(set = "GWAS_AD_Bellenguez2022", gene_source = c("mapped", "mapped_p1e6"), hs_symbol = "APOE",
                       n_assoc = NA_integer_, min_p = NA_real_)
gw_genes <- rbind(gw_genes, gwas_add[!gw_genes, on = .(set, gene_source, hs_symbol)], use.names = TRUE)
mp <- sym2id(gw_genes$hs_symbol)
gw_genes[, `:=`(hs_gene = mp$id, id_by = mp$how)]
gw_genes[, dog_gene := omap$dog[match(hs_gene, omap$human)]]
gw_genes[, `:=`(map_by = omap$by[match(hs_gene, omap$human)], dog_symbol = gloc$SYMBOL[match(dog_gene, gloc$GENE)])]
gw_genes[, why_missing := fifelse(is.na(hs_gene), "사람 GeneID 없음", fifelse(is.na(dog_gene), "개 유전자로 못 옮김", ""))]
gw_study <- g[PUBMEDID %in% gwas_spec$pmid, .(n_assoc_all = .N, n_assoc_5e8 = sum(p <= p_gwas, na.rm = TRUE), n_assoc_1e6 = sum(p <= p_gwas_sens, na.rm = TRUE),
                                              accession = paste(unique(`STUDY ACCESSION`), collapse = ";")), by = .(PUBMEDID, trait = `DISEASE/TRAIT`)]
gw_study[, used := mapply(function(pm, tr) tr %in% gwas_spec[pmid == pm, traits][[1]], PUBMEDID, trait)]

read_xlsx_sheet <- function(f, sheet) {
  rd <- function(p) paste(system2("unzip", c("-p", shQuote(f), shQuote(p)), stdout = TRUE), collapse = "\n")
  unesc <- function(x) { x <- gsub("&lt;", "<", x, fixed = TRUE); x <- gsub("&gt;", ">", x, fixed = TRUE)
                         x <- gsub("&quot;", "\"", x, fixed = TRUE); x <- gsub("&apos;", "'", x, fixed = TRUE); gsub("&amp;", "&", x, fixed = TRUE) }
  wb <- rd("xl/workbook.xml"); rels <- rd("xl/_rels/workbook.xml.rels")
  sh <- regmatches(wb, gregexpr("<sheet [^>]*>", wb))[[1]]
  sh_name <- unesc(sub('.*name="([^"]+)".*', "\\1", sh)); sh_rid <- sub('.*r:id="([^"]+)".*', "\\1", sh)
  rl <- regmatches(rels, gregexpr("<Relationship [^>]*>", rels))[[1]]
  tgt <- setNames(sub('.*Target="([^"]+)".*', "\\1", rl), sub('.*Id="([^"]+)".*', "\\1", rl))
  if (!sheet %in% sh_name) stop("시트 없음: ", sheet, " (", f, ")")
  p <- sub("^/", "", tgt[[sh_rid[match(sheet, sh_name)]]]); if (!startsWith(p, "xl/")) p <- paste0("xl/", p)
  ss <- rd("xl/sharedStrings.xml")
  ss <- unesc(gsub("<[^>]+>", "", regmatches(ss, gregexpr("<si>.*?</si>", ss, perl = TRUE))[[1]]))
  x <- rd(p)
  cells <- regmatches(x, gregexpr("<c [^>]*?(/>|>.*?</c>)", x, perl = TRUE))[[1]]
  ref <- sub('.*\\br="([A-Z]+)([0-9]+)".*', "\\1 \\2", sub(">.*", "", cells))
  val <- ifelse(grepl("<v>", cells), sub(".*<v>(.*?)</v>.*", "\\1", cells), "")
  isstr <- grepl('t="s"', sub(">.*", "", cells))
  val[isstr & val != ""] <- ss[as.integer(val[isstr & val != ""]) + 1]
  col <- sub(" .*", "", ref); row <- as.integer(sub(".* ", "", ref))
  cn <- vapply(strsplit(col, ""), function(l) sum((match(l, LETTERS)) * 26^(rev(seq_along(l)) - 1)), 0)
  d <- data.table(row, cn, val)
  m <- dcast(d, row ~ cn, value.var = "val", fill = "")
  m[, row := NULL]; m
}

palmer_file <- file.path(ref_dir, "palmer2021_brain_signature.tsv")
if (!file.exists(palmer_file)) {
  zf <- file.path(ref_dir, "palmer2021_supp.zip")
  if (!file.exists(zf) || system2("unzip", c("-tq", shQuote(zf)), stdout = FALSE, stderr = FALSE) != 0) {
    message("Palmer 2021 보충 자료 내려받기 (Europe PMC, 약 15 MB; 몇 분 걸릴 수 있음)")
    st <- system2("curl", c("-sL", "--fail", "--retry", "3", "-m", "900", "-o", shQuote(zf), shQuote(palmer_url)))
    if (st != 0 || system2("unzip", c("-tq", shQuote(zf)), stdout = FALSE, stderr = FALSE) != 0) stop("Palmer 2021 보충 자료 내려받기 실패")
  }
  td <- file.path(ref_dir, "palmer2021_tmp"); dir.create(td, showWarnings = FALSE)
  system2("unzip", c("-oq", shQuote(zf), "aging-13-202648-s004.xlsx", "aging-13-202648-s005.xlsx", "-d", shQuote(td)))
  rd_sig <- function(f, sh, dirn) {
    m <- read_xlsx_sheet(file.path(td, f), sh)
    hdr <- which(m[[1]] == "Gene")[1]; if (is.na(hdr)) stop("Palmer 표 머리글 없음: ", sh)
    d <- m[(hdr + 1):.N]; setnames(d, as.character(unlist(m[hdr]))[seq_len(ncol(d))])
    d[Gene != "", .(direction = dirn, Gene, Entrez, NumSucc, NumTrials, PVal, QValue, sheet = sh)]
  }
  ps <- rbind(rd_sig("aging-13-202648-s004.xlsx", "Over.Brain", "up"), rd_sig("aging-13-202648-s005.xlsx", "Under.Brain", "down"))
  fwrite(ps, palmer_file, sep = "\t")
  writeLines(c(paste("downloaded:", format(Sys.time(), "%Y-%m-%d %H:%M %Z")), paste("url:", palmer_url),
               "files: aging-13-202648-s004.xlsx (Over.Brain = Table S4), aging-13-202648-s005.xlsx (Under.Brain = Table S8)",
               "Palmer D et al. (2021) Ageing transcriptome meta-analysis reveals similarities and differences between key mammalian tissues. Aging 13:3313-3341"),
             file.path(ref_dir, "palmer2021_version.txt"))
  unlink(c(zf, td), recursive = TRUE)
}
ps <- fread(palmer_file, colClasses = c(Entrez = "character"))
ps[, QValue := as.numeric(QValue)]
age_sig <- ps[QValue < q_age]
mp <- sym2id(age_sig$Gene)
age_sig[, hs_gene := fifelse(Entrez %in% hs$GeneID, Entrez, mp$id)]
age_sig[, `:=`(dog_gene = omap$dog[match(hs_gene, omap$human)], map_by = omap$by[match(hs_gene, omap$human)])]
age_sig[, dog_symbol := gloc$SYMBOL[match(dog_gene, gloc$GENE)]]
age_up_set <- list(AGE_brain_up_Palmer2021 = to_dog(age_sig[direction == "up", hs_gene]))
cat(sprintf("Palmer 2021 뇌 signature (q < %g): 증가 %d · 감소 %d개 → 개 유전자 %d · %d개\n", q_age,
            age_sig[direction == "up", .N], age_sig[direction == "down", .N],
            age_sig[direction == "up", uniqueN(na.omit(dog_gene))], age_sig[direction == "down", uniqueN(na.omit(dog_gene))]))

gse_file <- file.path(ref_dir, "gse304125_cortex_normcounts.tsv.gz")
if (!file.exists(gse_file)) {
  xf <- file.path(ref_dir, "GSE304125_Cortex_EV_NormCounts.xlsx")
  if (!file.exists(xf)) {
    message("GSE304125 정규화 count 내려받기 (약 9 MB)")
    st <- system2("curl", c("-sL", "--fail", "--retry", "3", "-m", "900", "-o", shQuote(xf), shQuote(gse_url)))
    if (st != 0 || !file.exists(xf) || file.size(xf) < 1e6) stop("GSE304125 내려받기 실패")
  }
  m <- read_xlsx_sheet(xf, "Cortex Norm Counts")
  hdr <- as.character(unlist(m[1])); d <- m[-1]
  setnames(d, c("symbol", hdr[-1]))
  if (!all(grepl("^[yoc][0-9]+$", names(d)[-1]))) stop("GSE304125 열 이름이 예상과 다름: ", paste(names(d), collapse = ","))
  fwrite(d, gse_file, sep = "\t")
  writeLines(c(paste("downloaded:", format(Sys.time(), "%Y-%m-%d %H:%M %Z")), paste("url:", gse_url),
               "sheet: Cortex Norm Counts (DESeq2 normalized; y = young, o = old, c = CCD; canFam6 annotation)"),
             file.path(ref_dir, "gse304125_version.txt"))
  unlink(xf)
}
gse <- fread(cmd = sprintf("gzip -dc '%s'", gse_file))
grp <- c(y = "young", o = "old", c = "CCD")
smp <- names(gse)[-1]
cnt <- as.matrix(gse[, ..smp]); storage.mode(cnt) <- "double"
mean_g <- function(k) rowMeans(cnt[, substr(smp, 1, 1) == k, drop = FALSE])
ex <- data.table(symbol = gse$symbol, m_y = mean_g("y"), m_o = mean_g("o"), m_c = mean_g("c"))
ex <- ex[(m_y * 3 + m_o * 3 + m_c * 3) / 9 >= 10]
ex[, `:=`(lfc_old_young = log2((m_o + 0.5) / (m_y + 0.5)), lfc_ccd_young = log2((m_c + 0.5) / (m_y + 0.5)),
          lfc_ccd_old = log2((m_c + 0.5) / (m_o + 0.5)))]
ex <- ex[, if (.N == 1) .SD, by = symbol]
conc <- rbindlist(lapply(c("lfc_old_young", "lfc_ccd_young", "lfc_ccd_old"), function(v) rbindlist(lapply(c("up", "down"), function(dirn) {
  g <- unique(age_sig[direction == dirn, Gene]); x <- ex[symbol %in% g, get(v)]; bg <- ex[!symbol %in% age_sig$Gene, get(v)]
  p_bg <- if (dirn == "up") mean(bg > 0) else mean(bg < 0)
  k <- if (dirn == "up") sum(x > 0) else sum(x < 0)
  data.table(comparison = sub("lfc_", "", v), palmer = dirn, n_expressed = length(x), same_direction = k, frac = round(k / length(x), 3),
             frac_background = round(p_bg, 3), p_binom = signif(binom.test(k, length(x), p_bg, alternative = "greater")$p.value, 3),
             median_lfc = round(median(x), 3), median_lfc_background = round(median(bg), 3),
             p_wilcox = signif(wilcox.test(x, bg, alternative = if (dirn == "up") "greater" else "less")$p.value, 3))
}))))

read_sets <- function(f) { l <- strsplit(readLines(f), "[ \t]+"); setNames(lapply(l, `[`, -1), vapply(l, `[`, "", 1)) }
s22 <- read_sets(sets2022_file)
info22 <- fread(setinfo2022_file)
kegg_list <- s22[grepl("^KEGG_", names(s22))]
if (length(miss <- setdiff(grep("^KEGG", c(circ_sets, splice_sets), value = TRUE), names(kegg_list)))) stop("2022 KEGG 세트 없음: ", paste(miss, collapse = ", "))
if (length(miss <- setdiff(grep("^GO", c(circ_sets, splice_sets), value = TRUE), names(go_list)))) stop("GO 세트 없음: ", paste(miss, collapse = ", "))
db_list <- c(kegg_list, go_list)

gwas_set <- function(src) setNames(lapply(gwas_spec$set, function(s) gw_genes[set == s & gene_source == src & !is.na(dog_gene), unique(dog_gene)]), gwas_spec$set)
ns_hs <- sym2id(neuro_splice)
if (anyNA(ns_hs$id)) stop("신경 스플라이싱 유전자 사람 GeneID 없음: ", paste(neuro_splice[is.na(ns_hs$id)], collapse = ", "))
fixed <- c(db_list[c(circ_sets, splice_sets)], list(SPLICE_neuronal_regulators = to_dog(ns_hs$id), AGE_brain_Palmer2021 = to_dog(age_sig[direction %in% c("up", "down"), hs_gene])))
variants <- list(primary = c(fixed, gwas_set("mapped")), reported = c(fixed, gwas_set("reported")), p1e6 = c(fixed, gwas_set("mapped_p1e6")))
variants <- lapply(variants, function(v) v[primary_order])

name_of <- c(setNames(info22$name, info22$set), setNames(go_sets$name, go_sets$set), setNames(gwas_spec$name, gwas_spec$set),
             SPLICE_neuronal_regulators = "Neuronal splicing regulators (curated, 24 genes)",
            AGE_brain_Palmer2021 = "Brain ageing expression signature, up + down (Palmer 2021 meta-analysis, q < 0.05)")
name_of <- name_of[!duplicated(names(name_of))]
source_of <- c(setNames(rep("KEGG (2022)", length(kegg_list)), names(kegg_list)), setNames(rep("GO_BP", length(go_list)), names(go_list)),
               setNames(rep("GWAS Catalog", 5), gwas_spec$set), SPLICE_neuronal_regulators = "curated",
               AGE_brain_Palmer2021 = "Palmer 2021 (Aging 13:3313) Table S4·S8")
n_in <- function(v, u) vapply(v, function(x) sum(x %in% u), 0L)
si <- data.table(set = primary_order, family = family_of[primary_order], source = source_of[primary_order], name = name_of[primary_order])
si[set %in% gwas_spec$set, n_hs_input := gw_genes[gene_source == "mapped", uniqueN(hs_symbol), by = set][match(si[set %in% gwas_spec$set, set], set), V1]]
si[set == "SPLICE_neuronal_regulators", n_hs_input := length(neuro_splice)]
si[set == "AGE_brain_Palmer2021", n_hs_input := nrow(age_sig)]
for (v in names(variants)) {
  si[, (paste0("n_", v)) := lengths(variants[[v]])]
  si[, (paste0("n_", v, "_w10")) := n_in(variants[[v]], univ$w10)]
  si[, (paste0("n_", v, "_w50")) := n_in(variants[[v]], univ$w50)]
}
si[, n_by_same_symbol := vapply(variants$primary, function(x) sum(x %in% omap[by == "same_symbol", dog]), 0L)]
si[, n_by_same_fullname := vapply(variants$primary, function(x) sum(x %in% omap[by == "same_fullname", dog]), 0L)]
si[grepl("^KEGG_", set), `:=`(n_by_same_symbol = NA_integer_, n_by_same_fullname = NA_integer_)]
si[, p1e6_same_as_primary := mapply(function(a, b) setequal(a, b), variants$primary, variants$p1e6)]
si[, `:=`(n_2022 = vapply(set, function(s) if (s %in% names(s22)) length(s22[[s]]) else NA_integer_, 0L),
          n_shared_2022 = vapply(set, function(s) if (s %in% names(s22)) length(intersect(s22[[s]], variants$primary[[s]])) else NA_integer_, 0L))]
si[, kept := n_primary_w10 >= set_min]
si[, dog_symbols := vapply(variants$primary, function(x) paste(sort(gloc$SYMBOL[match(x, gloc$GENE)]), collapse = ";"), "")]

aux <- db_list[lengths(db_list) >= set_min & lengths(db_list) <= aux_max & !grepl(paste0("^KEGG_", sub("\\^", "", kegg_exclude)), names(db_list))]
aux <- c(aux, db_list[setdiff(circ_sets, names(aux))])
aux_info <- data.table(set = names(aux), source = source_of[names(aux)], name = name_of[names(aux)], n_genes = lengths(aux),
                       n_w10 = n_in(aux, univ$w10))
aux_info[, `:=`(n_2022 = vapply(set, function(s) if (s %in% names(s22)) length(s22[[s]]) else NA_integer_, 0L),
                in_2022 = set %in% names(s22))]

write_sets <- function(v, f) {
  v <- v[n_in(v, univ$w10) >= set_min]
  writeLines(vapply(names(v), function(s) paste(c(s, v[[s]]), collapse = " "), ""), f)
  length(v)
}
p1e6_changed <- si[p1e6_same_as_primary == FALSE, set]
p1e6_out <- setNames(variants$p1e6[p1e6_changed], paste0(p1e6_changed, "_p1e6"))
n_written <- c(primary  = write_sets(variants$primary,  paste0(out_prefix, "_sets_primary.txt")),
               reported = write_sets(variants$reported, paste0(out_prefix, "_sets_primary_reported.txt")),
               p1e6     = write_sets(p1e6_out,          paste0(out_prefix, "_sets_primary_p1e6.txt")),
               age_up   = write_sets(age_up_set,        paste0(out_prefix, "_sets_age_up.txt")),
               aux      = write_sets(aux,               paste0(out_prefix, "_sets_aux.txt")))
fwrite(age_sig[order(direction, QValue)], paste0(out_prefix, "_age_genes.csv"))
fwrite(conc, paste0(out_prefix, "_age_dog_check.csv"))
fwrite(si, paste0(out_prefix, "_set_info.csv"))
fwrite(aux_info, paste0(out_prefix, "_aux_info.csv"))
fwrite(omap[order(by, hs_symbol)], paste0(out_prefix, "_ortholog_map.csv"))
fwrite(gw_genes[order(set, gene_source, hs_symbol), .(set, gene_source, hs_symbol, n_assoc, min_p = signif(min_p, 3), hs_gene, id_by,
                                                     dog_gene, dog_symbol, map_by, why_missing)], paste0(out_prefix, "_gwas_genes.csv"))
ov <- CJ(a = primary_order, b = primary_order)[a < b]
ov[, `:=`(n_a = lengths(variants$primary[a]), n_b = lengths(variants$primary[b]),
          n_shared = mapply(function(x, y) length(intersect(variants$primary[[x]], variants$primary[[y]])), a, b))]
ov[, jaccard := round(n_shared / (n_a + n_b - n_shared), 3)]
fwrite(ov[order(-jaccard)], paste0(out_prefix, "_set_overlap.csv"))

sink(paste0(out_prefix, "_summary.txt"), split = TRUE)
cat("06-4 MAGMA 유전자 세트 (우리 GWAS의 z·p는 읽지 않음 — 결과를 보기 전에 고정)\n")
cat(sprintf("GWAS Catalog %s; 연관 p ≤ %g (민감도 p ≤ %g)\n", gw_ver, p_gwas, p_gwas_sens))
cat(sprintf("사람 → 개: NCBI 1:1 %s쌍 + 같은 이름 %s쌍 + 같은 전체 이름 %s쌍 (NCBI 표 %s쌍 중 %s쌍은 CanFam3.1에 없는 개 ID를 가리키거나 1:1이 아님)\n",
            format(nrow(orth), big.mark = ","), format(nrow(same), big.mark = ","), format(nrow(fulln), big.mark = ","),
            format(n_orth_raw, big.mark = ","), format(n_orth_raw - nrow(orth), big.mark = ",")))
cat("  같은 전체 이름으로 옮긴 것 중 이름이 LOC가 아닌 개 유전자: ",
    paste(omap[by == "same_fullname" & !grepl("^LOC", dog_symbol), paste0(dog_symbol, "↔", hs_symbol)], collapse = ", "), "\n\n", sep = "")
cat("GWAS 논문별 형질 (used = 세트에 씀)\n"); print(gw_study[order(PUBMEDID, -used)])
cat("\nGWAS 세트: 사람 유전자 → 개 유전자 (mapped = 위치 유전자, 주 분석 / reported = 저자 유전자, 민감도)\n")
print(gw_genes[, .(n_hs = .N, hs_id = sum(!is.na(hs_gene)), dog = sum(!is.na(dog_gene)), dog_by_same_symbol = sum(map_by == "same_symbol", na.rm = TRUE),
                   dog_by_same_fullname = sum(map_by == "same_fullname", na.rm = TRUE),
                   hs_by_synonym = sum(id_by == "synonym", na.rm = TRUE), not_moved = sum(why_missing != "")),
               by = .(set, gene_source)][order(set, gene_source)])
cat("  못 옮긴 사람 유전자 예 (mapped, 단백질 코딩만, 앞 20개; 전체는 _gwas_genes.csv):\n")
pc <- hs[type_of_gene == "protein-coding", GeneID]
print(gw_genes[gene_source == "mapped" & why_missing != "" & hs_gene %in% pc,
               .(n_protein_coding = .N, examples = paste(head(hs_symbol, 20), collapse = ", ")), by = set])
cat(sprintf("\n뇌 노화 발현 세트 (Palmer 2021, q < %g; 사람 Entrez → 개): 증가 %d → 개 %d · 감소 %d → 개 %d (못 옮긴 것: %s)\n", q_age,
            age_sig[direction == "up", .N], age_sig[direction == "up", uniqueN(na.omit(dog_gene))],
            age_sig[direction == "down", .N], age_sig[direction == "down", uniqueN(na.omit(dog_gene))],
            paste(age_sig[is.na(dog_gene), Gene], collapse = ", ")))
cat(sprintf("개 피질 방향 확인 (GSE304125: 어린 개 3 · 늙은 개 3 · CCD 3, 9마리 평균 정규화 count ≥ 10인 유전자 %s개; 설명용, 세트 정의에 쓰지 않음)\n",
            format(nrow(ex), big.mark = ",")))
cat("  same_direction = Palmer와 같은 방향으로 바뀐 유전자 수, frac_background = 나머지 발현 유전자에서 같은 방향 비율 (이항 검정의 기대값)\n")
print(conc)
cat(sprintf("\n주 가설 세트 (±10 kb 분석 유전자 안 %d개 이상만 파일에 씀): primary %d개 · reported %d개 · p1e6 %d개 · age_up %d개\n",
            set_min, n_written["primary"], n_written["reported"], n_written["p1e6"], n_written["age_up"]))
print(si[, .(set, family, n_hs_input, n_primary, n_primary_w10, n_primary_w50, n_reported_w10, n_p1e6_w10, p1e6_same = p1e6_same_as_primary,
             n_by_same_symbol, n_by_same_fullname, n_2022, n_shared_2022, kept)])
cat("\n주 가설 세트끼리 겹침 (Jaccard ≥ 0.1)\n"); print(ov[jaccard >= 0.1][order(-jaccard)])
cat(sprintf("\n보조 세트 %s개 (KEGG %d · GO %s; 2022 %s개)\n", format(n_written["aux"], big.mark = ","), sum(aux_info$source != "GO_BP" & aux_info$n_w10 >= set_min),
            format(sum(aux_info$source == "GO_BP" & aux_info$n_w10 >= set_min), big.mark = ","), format(length(s22), big.mark = ",")))
cat(sprintf("  GO 세트 중 2022에도 있던 것 %s개, 유전자 수 중앙값 2022 → 이번 %g → %g\n",
            format(aux_info[source == "GO_BP", sum(in_2022)], big.mark = ","),
            aux_info[source == "GO_BP" & in_2022, median(n_2022)], aux_info[source == "GO_BP" & in_2022, median(n_genes)]))
cat("\n세트별 개 유전자 (primary)\n")
for (i in seq_len(nrow(si))) cat(sprintf("  %s (%d): %s\n", si$set[i], si$n_primary[i], strtrim(si$dog_symbols[i], 400)))
sink()

set_files <- paste0(out_prefix, c("_sets_primary.txt", "_sets_primary_reported.txt", "_sets_primary_p1e6.txt", "_sets_age_up.txt", "_sets_aux.txt"))
md5 <- tools::md5sum(set_files)
plan <- c(
  sprintf("06-5 분석 계획 (06-4가 %s에 씀; 06-5 결과를 보기 전에 고정, 사용자 결정 2026-11-12)", format(Sys.time(), "%Y-%m-%d %H:%M %Z")),
  "",
  "검정: MAGMA competitive gene-set (기본 공변량: 유전자 크기·밀도·inverse MAC 등), 유전자 결과 = 06-3 {run}.genes.raw (±10 kb)",
  sprintf("주 가설 세트 %d개 (_sets_primary.txt): %s", n_written["primary"], paste(names(variants$primary)[si$kept], collapse = ", ")),
  sprintf("  묶음: %s", paste(sprintf("%s %d", names(table(si$family[si$kept])), table(si$family[si$kept])), collapse = " · ")),
  "",
  "1단계 — 주 결론 (논문의 핵심 주장)",
  sprintf("  main 모형 (7,115마리) × %d개 세트 → %d개 안에서 BH FDR, q < 0.05 = '유의'", n_written["primary"], n_written["primary"]),
  "  유의한 세트는 다른 모형에서 같은 방향 (β > 0)·p < 0.05인지 확인 (재현 확인, 판정에 쓰지 않음)",
  "",
  "2단계 — 나이 가설 (사전 지정 부 가설): '이 경로의 유전 효과가 나이에 따라 달라지는가'",
  "  int_age (SNP × 나이 상호작용, 7,115마리) + het_age (young8 4,866 vs age8 2,249의 효과 차이)",
  sprintf("  2개 × %d개 세트 = %d개 검정을 한 family로 BH FDR, q < 0.05", n_written["primary"], 2L * n_written["primary"]),
  "  두 검정은 방향이 없는 p라 '나이 든 개에서 더 강하다'까지는 말하지 않음",
  "  참고 표 (판정에 쓰지 않음): age8·young8 세트 결과와 세트 β 차이 z = (β_age8 − β_young8) / √(se_age8² + se_young8²)",
  "",
  "3단계 — 탐색 (결론 아님, 표로만)",
  "  purelike · mixlike · sex_m · sex_f · wt_t1 ~ t3 · bin_age8 · age8 · young8: 모형마다 세트 안 BH",
  "  het_sex · het_ancestry · het_weight: 비교마다 세트 안 BH",
  "",
  "공통",
  "  null 보정 확인: 06-3 null1 ~ null20.genes.raw에 같은 세트 → 세트별 p < 0.05 비율, 세트 p의 λ",
  "  명목 p < 0.05: 탐색 표시만, null에서의 p < 0.05 비율과 함께",
  sprintf("  민감도 (1·2단계 결과가 같은 방향·p < 0.05로 남는지만): ±50 kb (06-3 {model}_w50, het_age_w50) / REPORTED GENE(S) (_sets_primary_reported.txt) / p ≤ 1e-6 (_sets_primary_p1e6.txt: %s)",
          paste(names(p1e6_out), collapse = ", ")),
  sprintf("             / 뇌 노화 세트를 증가 유전자만 (_sets_age_up.txt: AGE_brain_up_Palmer2021, %d개)", length(age_up_set[[1]])),
  "  보조 GO/KEGG (_sets_aux.txt): 모형마다 따로 BH, 결과 설명용",
  "  개 뇌에서의 방향 확인 (설명용): Palmer 세트 유전자가 GSE304125 개 피질 (어린 3 · 늙은 3 · CCD 3)에서 같은 방향인지 (_age_dog_check.csv)",
  "  묶음별 BH는 쓰지 않음 (family 2 ~ 4개라 보정이 가장 약함)",
  "",
  "알려진 한계: AD 세트에서 TREM2·CR1·MS4A4A·LILRB2 등 미세아교세포·면역 유전자가 CanFam3.1 1:1 대응이 없어 빠짐",
  "            뇌 노화 세트는 사람·쥐·랫드 메타분석 (뇌 자료 29개)을 옮긴 것이라 개 자체 자료가 아님 (개 자료는 마릿수가 적어 세트 정의에 쓰지 않음)",
  "",
  "세트 파일 md5",
  sprintf("  %s  %s", md5, basename(set_files)))
writeLines(plan, paste0(out_prefix, "_analysis_plan.txt"))
cat(sprintf("저장: %s_sets_primary.txt · _sets_primary_reported.txt · _sets_primary_p1e6.txt · _sets_age_up.txt · _sets_aux.txt · _age_genes.csv · _age_dog_check.csv · _set_info.csv · _aux_info.csv · _ortholog_map.csv · _gwas_genes.csv · _set_overlap.csv · _summary.txt · _analysis_plan.txt\n", out_prefix))

#        GO_0007623               5                 0
#        GO_0042752               7                 0
#        chronotype               3                 0
#        AD                       1                 1 (HLA-DQA1)
#        GO_0000398               10                0
#        GO_0048024               2                 0
#        SPLICE_neuronal          0                 0
