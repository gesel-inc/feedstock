# This defines equivalence classes for each gene - namely, all synonymous identifiers between Entrez and Ensembl. 
# This may group multiple identifiers of the same type together if they are all synonymous with a single identifier of another type.
# Indeed, we need to build a graph and isolate each component to resolve very complex synonym structures.
# (For obvious reasons, symbols are not used to identify synonyms as they are ambiguous and too often re-used.)

library(AnnotationHub)
ahub <- AnnotationHub(ask=FALSE)
source("../annotations.R")

library(BiocParallel)
output.dir <- "_built"
unlink(output.dir, recursive=TRUE)
dir.create(output.dir, showWarnings=FALSE)

dump <- function(x, out) {
    if (is.list(x)) {
        dump <- unlist(bplapply(x, paste, collapse="\t", BPPARAM=MulticoreParam()))
    } else {
        dump <- x
    }
    handle <- gzfile(file.path(output.dir, out), open="wb")
    writeLines(dump, con=handle)
    close(handle)
}

library(igraph)
for (species in names(annotations)) {
    ensdb <- ahub[[annotations[[species]]$ensdb]]
    orgdb <- ahub[[annotations[[species]]$orgdb]]

    from.o <- select(orgdb, keys=keys(orgdb), columns="ENSEMBL")
    from.e <- select(ensdb, keys=keys(ensdb), columns="ENTREZID")

    all.edges <- rbind(
        DataFrame(ensembl = from.o$ENSEMBL, entrez = from.o$ENTREZID),
        DataFrame(ensembl = from.e$GENEID, entrez = from.e$ENTREZID)
    )
    edges <- all.edges[!is.na(all.edges$ensembl) & !is.na(all.edges$entrez),]
    edges <- unique(edges)

    all.ensembl <- setdiff(unique(all.edges$ensembl), unique(edges$ensembl))
    all.entrez <- setdiff(unique(all.edges$entrez), unique(edges$entrez))
    g <- make_graph(rbind(edges$ensembl, edges$entrez), isolates=c(all.ensembl[!is.na(all.ensembl)], all.entrez[!is.na(all.entrez)]), directed=FALSE)
    classes <- components(g)
    pooled <- names(classes$membership)

    is.entrez <- grepl("^[0-9]", pooled)
    f <- factor(as.integer(unname(classes$membership)))
    stopifnot(identical(levels(f), as.character(seq_along(levels(f)))))

    names.entrez <- names(classes$membership)[is.entrez]
    f.entrez <- f[is.entrez]
    by.entrez <- split(names.entrez, f.entrez, drop=FALSE)
    dump(by.entrez, paste0(species, "_entrez.tsv.gz"))

    names.ensembl <- names(classes$membership)[!is.entrez]
    f.ensembl <- f[!is.entrez]
    by.ensembl <- split(names.ensembl, f.ensembl, drop=FALSE)
    dump(by.ensembl, paste0(species, "_ensembl.tsv.gz"))

    stopifnot(identical(names(by.entrez), names(by.ensembl)))
    by.entrez <- unname(by.entrez)
    by.ensembl <- unname(by.ensembl)

    entrez2sym <- select(orgdb, keys=unique(names.entrez), columns="SYMBOL")
    entrez2sym <- entrez2sym[!is.na(entrez2sym$SYMBOL),]
    sym.by.entrez <- split(entrez2sym$SYMBOL, entrez2sym$ENTREZID)

    ensembl2sym <- select(ensdb, keys=unique(names.ensembl), columns="SYMBOL")
    ensembl2sym <- ensembl2sym[!is.na(ensembl2sym$SYMBOL),]
    sym.by.ensembl <- split(ensembl2sym$SYMBOL, ensembl2sym$GENEID)

    m.entrez <- split(match(names.entrez, names(sym.by.entrez)), f.entrez, drop=FALSE)
    m.ensembl<- split(match(names.ensembl, names(sym.by.ensembl)), f.ensembl, drop=FALSE)
    m.entrez <- unname(m.entrez)
    m.ensembl <- unname(m.ensembl)

    # This bit isn't very fast, but whatever.
    by.sym <- bplapply(seq_along(by.entrez), FUN=function(i) {
        union(
            unlist(sym.by.entrez[m.entrez[[i]]], use.names=FALSE), 
            unlist(sym.by.ensembl[m.ensembl[[i]]], use.names=FALSE)
        )
    }, BPPARAM=MulticoreParam())
    dump(by.sym, paste0(species, "_symbols.tsv.gz"))
}

write("genes-v0.3.0", file="_tag")
payload <- capture.output(print(sessionInfo()))
write(c("<details>", "<summary>Session information</summary>", "", "```", payload, "```", "</details>"), file="_session")
