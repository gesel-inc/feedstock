# This defines equivalence classes for each gene - namely, all synonymous identifiers between Entrez and Ensembl. 
# This may group multiple identifiers of the same type together if they are all synonymous with a single identifier of another type.
# Indeed, we need to build a graph and isolate each component to resolve very complex synonym structures.
# (For obvious reasons, symbols are not used to identify synonyms as they are ambiguous and too often re-used.)

library(AnnotationHub)
ahub <- AnnotationHub(ask=FALSE)
source("../annotations.R")

output.dir <- "_built"
unlink(output.dir, recursive=TRUE)
dir.create(output.dir, showWarnings=FALSE)

dump <- function(x, out) {
    dump <- vapply(x, paste, collapse="\t", FUN.VALUE="")
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

    # Creating equivalent classes based on networks of Entrez/Ensembl synonyms. 
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

    # Dumping the Entrez and Ensembl IDs for each equivalence class.
    names.entrez <- names(classes$membership)[is.entrez]
    f.entrez <- f[is.entrez]
    by.entrez <- split(names.entrez, f.entrez, drop=FALSE)

    names.ensembl <- names(classes$membership)[!is.entrez]
    f.ensembl <- f[!is.entrez]
    by.ensembl <- split(names.ensembl, f.ensembl, drop=FALSE)

    stopifnot(identical(names(by.entrez), names(by.ensembl)))
    dump(by.entrez, paste0(species, "_entrez.tsv.gz"))
    dump(by.ensembl, paste0(species, "_ensembl.tsv.gz"))

    # Finding the symbols for each equivalence class.
    entrez2sym <- select(orgdb, keys=unique(names.entrez), columns="SYMBOL")
    entrez2sym <- entrez2sym[!is.na(entrez2sym$SYMBOL),]
    sym.by.entrez <- split(entrez2sym$SYMBOL, entrez2sym$ENTREZID)

    m.entrez <- match(names.entrez, names(sym.by.entrez))
    keep.entrez <- !is.na(m.entrez)
    remapped.sym.entrez <- unname(sym.by.entrez[m.entrez[keep.entrez]])
    remapped.f.entrez <- rep(f.entrez[keep.entrez], lengths(remapped.sym.entrez))
    remapped.sym.entrez <- unlist(remapped.sym.entrez) 

    ensembl2sym <- select(ensdb, keys=unique(names.ensembl), columns="SYMBOL")
    ensembl2sym <- ensembl2sym[!is.na(ensembl2sym$SYMBOL),]
    sym.by.ensembl <- split(ensembl2sym$SYMBOL, ensembl2sym$GENEID)

    m.ensembl <- match(names.ensembl, names(sym.by.ensembl))
    keep.ensembl <- !is.na(m.ensembl)
    remapped.sym.ensembl <- unname(sym.by.ensembl[m.ensembl[keep.ensembl]])
    remapped.f.ensembl <- rep(f.ensembl[keep.ensembl], lengths(remapped.sym.ensembl))
    remapped.sym.ensembl <- unlist(remapped.sym.ensembl) 

    by.sym <- split(c(remapped.sym.entrez, remapped.sym.ensembl), c(remapped.f.entrez, remapped.f.ensembl), drop=FALSE)
    by.sym <- lapply(by.sym, unique)
    stopifnot(identical(names(by.entrez), names(by.sym)))
    dump(by.sym, paste0(species, "_symbol.tsv.gz"))
}

payload <- capture.output(print(sessionInfo()))
write(c("<details>", "<summary>Session information</summary>", "", "```", payload, "```", "</details>"), file="_session")
