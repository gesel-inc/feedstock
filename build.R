library(jsonlite)
library(S4Vectors)

###########################################################
# Update these for new releases.

input.dirs <- c(
    file.path("go", "_built"),
    file.path("msigdb", "_built")
)

gene.dirs <- file.path("genes", "_built")

output.tag <- "v0.3.0"

###########################################################

stopifnot(file.exists(gene.dirs))
for (idir in input.dirs) {
    stopifnot(file.exists(idir))
}

species.set.info <- list()
species.collections <- list()
species.counter <- list()
species.mapping <- list()
species.set.membership <- list()
species.genes <- list()
species.ngenes <- list()

for (idir in input.dirs) {
    all.files <- list.files(idir, pattern="\\.json$")
    
    for (f in all.files) {
        gc()
        current <- jsonlite::fromJSON(file.path(idir, f), simplifyVector=FALSE)

        # Checking validity.
        stopifnot(isSingleString(current$title), !grepl("[\n\t]", current$title))
        stopifnot(isSingleString(current$description), !grepl("[\n\t]", current$description))
        stopifnot(isSingleString(current$maintainer), !grepl("[\n\t]", current$maintainer))
        stopifnot(isSingleString(current$species), !grepl("[\n\t]", current$species))
        stopifnot(isSingleString(current$source), !grepl("[\n\t]", current$source))
        stopifnot(current$id %in% c("entrez", "ensembl", "symbol"))

        # Setting up a species.
        cur.species <- current$species
        if (!(cur.species %in% names(species.genes))) {
            species.mapping[[cur.species]] <- list()
            species.counter[[cur.species]] <- 0L
            species.genes[[cur.species]] <- integer()
            species.set.membership[[cur.species]] <- list()
            species.set.info[[cur.species]] <- list()
            species.collections[[cur.species]] <- list()
        }

        # Loading the GMT file.
        gmt.path <- file.path(idir, sub("\\.json$", ".gmt.gz", f))
        fragments <- strsplit(readLines(gmt.path), "\t")

        set.names <- vapply(fragments, FUN=function(x) x[1], FUN.VALUE="")
        set.descriptions <- vapply(fragments, FUN=function(x) x[2], FUN.VALUE="")
        cursets <- lapply(fragments, FUN=tail, n=-2)
        gene.ids <- unlist(cursets, use.names=FALSE)
        set.ids <- rep(seq_along(cursets), lengths(cursets))

        # Mapping to the gene IDs.
        gene.mappings <- species.mapping[[cur.species]]
        if (!(current$id %in% names(gene.mappings))) {
            gene.path <- file.path(gene.dirs, paste0(cur.species, "_", current$id, ".tsv.gz"))
            all.lines <- readLines(gene.path)
            species.ngenes[[cur.species]] <- length(all.lines)

            collected <- vector("list", length(all.lines))
            keep <- all.lines != ""
            collected[keep] <- strsplit(all.lines[keep], "\t")

            gene.mappings[[current$id]] <- split(rep(seq_along(collected), lengths(collected)), unlist(collected))
            species.mapping[[cur.species]] <- gene.mappings
        }
        known.ids <- gene.mappings[[current$id]]

        m <- match(gene.ids, names(known.ids))
        keep <- !is.na(m)
        stopifnot(mean(keep) >= 0.95)

        gene.ids.mapped <- known.ids[m[keep]]
        set.ids.mapped <- rep(set.ids[keep], lengths(gene.ids.mapped))
        gene.ids.mapped <- unlist(gene.ids.mapped, use.names=FALSE)
        stopifnot(length(set.ids.mapped) == length(gene.ids.mapped))

        by.set <- split(gene.ids.mapped, factor(set.ids.mapped, seq_along(set.names)))
        by.set <- lapply(by.set, unique)

        # Filling the values.
        j <- length(species.set.info[[cur.species]]) + 1L
        species.set.info[[cur.species]][[j]] <- DataFrame(
            name=set.names,
            description=set.descriptions,
            size=unname(lengths(by.set)),
            collection=j,
            number=seq_along(set.names)
        )

        counter <- species.counter[[cur.species]]
        species.collections[[cur.species]][[j]] <- DataFrame(
            title=current$title,
            description=current$description,
            maintainer=current$maintainer,
            source=current$source,
            start=counter,
            size=length(set.names)
        )

        species.set.membership[[cur.species]][[j]] <- by.set
        species.counter[[cur.species]] <- counter + length(set.names)
    }
}

#################################
# Looping across species.

library(gesel)
output.dir <- "_built"
unlink(output.dir, recursive=TRUE)
dir.create(output.dir, showWarnings=FALSE)

for (species in names(species.set.info)) {
    set.info <- do.call(rbind, species.set.info[[species]])
    collections <- do.call(rbind, species.collections[[species]])
    set.membership <- do.call(c, species.set.membership[[species]])
    stopifnot(anyDuplicated(collections$title) == 0)

    prepareDatabaseFiles(
        species=species,
        collections=collections,
        set.info=set.info,
        set.membership=set.membership,
        num.genes=species.ngenes[[species]],
        path=output.dir
    )
}

payload <- capture.output(print(sessionInfo()))
write(c("<details>", "<summary>Session information</summary>", "", "```", payload, "```", "</details>"), file="_session")
