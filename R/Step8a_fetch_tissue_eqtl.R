###############################################################################
## Step 8a — 조직 cis-eQTL 자동 fetch (GTEx Portal API v2)  [eQTL Catalogue 500 → GTEx로 전환]
## Project: 01_GWAS_Asthma_MDD  (면역축: 조직특이 검정)
##
## 공유 면역유전자 7개 × GTEx v8 폐(Lung)·뇌(Brain_Cortex)·전혈(Whole_Blood) cis-eQTL fetch.
## GTEx singleTissueEqtl 필드: snpId(rsid) variantId(chr_pos_ref_alt_b38) pValue nes ...
##   → Step8b에서 nes/pValue로 beta/SE 유도, variantId로 EA(alt)/NEA(ref).
## 출력: raw/immune/tissue_eqtl/<TISSUE>__<GENE>.tsv  +  _fetch_log.csv
## 실행: source("Step8a_fetch_tissue_eqtl.R")  (인터넷 필요, GTEx API)
###############################################################################

suppressPackageStartupMessages({
  for (p in c("httr","jsonlite","data.table")) if (!requireNamespace(p,quietly=TRUE)) install.packages(p, repos="https://cloud.r-project.org")
  library(httr); library(jsonlite); library(data.table)
})

PROJ <- "PATH/TO/PROJECT"
OUT  <- file.path(PROJ,"raw","immune","tissue_eqtl")
dir.create(OUT, showWarnings=FALSE, recursive=TRUE)
GTEX <- "https://gtexportal.org/api/v2"

GENES   <- c("IL33","IL1RL1","IL4R","IL13","TSLP","TYK2","IL2RA")
TISSUES <- c("Lung","Brain_Cortex","Whole_Blood")

api_get <- function(path, query=list(), tries=4) {
  for (i in 1:tries) {
    r <- try(GET(paste0(GTEX,path), query=query, timeout(60)), silent=TRUE)
    if (!inherits(r,"try-error") && status_code(r)==200) {
      d <- try(fromJSON(content(r,"text",encoding="UTF-8"), flatten=TRUE), silent=TRUE)
      if (!inherits(d,"try-error")) return(d)
    }
    Sys.sleep(2*i)
  }
  NULL
}

## ---- 1) 심볼 → gencodeId (버전 포함) ------------------------------------ ##
message("=== gencodeId 해결 ===")
gencode <- setNames(rep(NA_character_, length(GENES)), GENES)
for (sym in GENES) {
  d <- api_get("/reference/gene", list(geneId=sym))
  gd <- d$data
  if (is.data.frame(gd) && nrow(gd)) {
    idx <- which(toupper(gd$geneSymbol)==toupper(sym))
    gencode[sym] <- gd$gencodeId[if (length(idx)) idx[1] else 1]
  }
  message("  ", sym, " → ", gencode[sym])
}

## ---- 2) gene × tissue eQTL fetch (페이지네이션) ------------------------- ##
fetch_eqtl <- function(gid, tissue) {
  out <- list(); page <- 0; PER <- 1000
  repeat {
    d <- api_get("/association/singleTissueEqtl",
                 list(gencodeId=gid, tissueSiteDetailId=tissue, datasetId="gtex_v8",
                      itemsPerPage=PER, page=page))
    dd <- d$data
    if (!is.data.frame(dd) || !nrow(dd)) break
    out[[length(out)+1]] <- as.data.table(dd)
    np <- d$paging_info$numberOfPages
    page <- page + 1
    if (is.null(np) || page >= np) break
    if (page > 50) break
  }
  if (!length(out)) return(NULL)
  rbindlist(out, fill=TRUE)
}

message("\n=== 유전자 x 조직 fetch (GTEx v8) ===")
log <- list()
for (t in TISSUES) for (g in GENES) {
  gid <- gencode[g]
  f <- file.path(OUT, paste0(t,"__",g,".tsv"))
  if (!is.na(gid) && file.exists(f)) { message("  [건너뜀] ",t," ",g);
    log[[length(log)+1]] <- data.table(tissue=t,gene=g,nSNP=as.integer(NA)); next }
  if (is.na(gid)) { message("  [",t," ",g,"] gencodeId 없음 → 건너뜀");
    log[[length(log)+1]] <- data.table(tissue=t,gene=g,nSNP=0L); next }
  dt <- try(fetch_eqtl(gid, t), silent=TRUE)
  n <- if (inherits(dt,"try-error")||is.null(dt)) 0L else nrow(dt)
  if (n>0) { dt[, `:=`(GENE=g, TISSUE=t)]; fwrite(dt, f, sep="\t") }
  message(sprintf("  %-13s %-7s → %6d SNP%s", t, g, n, if(n>0) "" else " (없음)"))
  log[[length(log)+1]] <- data.table(tissue=t, gene=g, nSNP=n)
  Sys.sleep(0.3)
}
logdt <- rbindlist(log)
fwrite(logdt, file.path(OUT,"_fetch_log.csv"))
message("\n=== fetch 요약 (유전자 x 조직 SNP수) ===")
print(dcast(logdt, gene~tissue, value.var="nSNP"))
message("\n저장: ", OUT)
message("관건: IL33/TSLP/IL13 이 폐(Lung)에서 나오나 (혈액엔 없었음).")
message("다음: 요약표 + 아무 파일 헤더 붙여주면 Step8b(cis-MR+coloc) 확정.")
