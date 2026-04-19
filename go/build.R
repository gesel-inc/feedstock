# Re-using the same annotations that were used to define the Entrez gene IDs.

library(AnnotationHub)
ahub <- AnnotationHub(ask=FALSE)
source("../annotations.R")

library(GO.db)
go.info <- GO.db::GO_dbInfo()
go.version <- go.info[go.info$name == "GOSOURCEDATE","value"]

output.tag <- "v1.1.0"
output.dir <- "_built"
unlink(output.dir, recursive=TRUE)
dir.create(output.dir, showWarnings=FALSE)

for (species in names(annotations)) {
    orgdb.name <- annotations[[species]]$orgdb
    orgdb <- ahub[[orgdb.name]]
    mappings <- select(orgdb, keytype="GO", keys=keys(orgdb, "GO"), columns="ENTREZID")

    library(BiocParallel)
    output <- split(mappings$ENTREZID, mappings$GO)
    output <- bplapply(output, function(x) {
        current <- unique(sort(x))
        paste(current, collapse="\t")
    }, BPPARAM=MulticoreParam())

    # Saving the names and descriptions.
    info <- select(GO.db, keys=names(output), column="TERM")
    payload <- sprintf("%s\t%s\t%s", names(output), info$TERM[match(names(output), info$GOID)], output)

    con <- gzfile(file.path(output.dir, paste0(species, ".gmt.gz")), open="wb")
    write(payload, file=con)
    close(con)

    write(
        file=file.path(output.dir, paste0(species, ".json")),
        jsonlite::toJSON(
            list(
                title="Gene ontology",
                description=paste0("Gene sets defined from the Gene Ontology (version ", go.version, "), sourced from ", orgdb.name, " in Bioconductor's AnnotationHub."),
                species=species,
                maintainer="Aaron Lun",
                id="entrez",
                source=paste0("https://github.com/gesel-inc/feedstock/blob/gene-ontology-", output.tag, "/go/build.R")
            ),
            auto_unbox=TRUE,
            pretty=4
        )
    )
}

write(paste0("gene-ontology-", output.tag), file="_tag")
payload <- capture.output(print(sessionInfo()))
write(c("<details>", "<summary>Session information</summary>", "", "```", payload, "```", "</details>"), file="_session")
