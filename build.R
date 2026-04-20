library(jsonlite)

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

isSingleString <- function(x) {
    (is.character(x) && length(x) == 1 && !is.na(x))
}

species.set.info <- list()
species.collections <- list()
species.mapping <- list()
species.set.membership <- list()
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
        if (!(cur.species %in% names(species.mapping))) {
            species.mapping[[cur.species]] <- list()
            species.set.membership[[cur.species]] <- list()
            species.set.info[[cur.species]] <- list()
            species.collections[[cur.species]] <- list()
        }

        # Loading the gene IDs.
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

        # Loading the GMT file.
        gmt.path <- file.path(idir, sub("\\.json$", ".gmt.gz", f))
        fragments <- strsplit(readLines(gmt.path), "\t")

        set.names <- vapply(fragments, FUN=function(x) x[1], FUN.VALUE="")
        set.descriptions <- vapply(fragments, FUN=function(x) x[2], FUN.VALUE="")
        set.memberships <- lapply(fragments, FUN=tail, n=-2)

        # Mapping to the gene identifiers, accounting for 1:many mapping to synonyms. 
        raw.genes <- unlist(set.memberships)
        raw.sets <- rep(seq_along(set.memberships), lengths(set.memberships))

        m <- match(raw.genes, names(known.ids))
        keep <- !is.na(m)
        stopifnot(mean(keep) >= 0.95) # check that most genes are found, to avoid obvious problems with the incorrect ID/species.
        remapped.genes <- unname(known.ids[m[keep]])
        remapped.sets <- rep(raw.sets[keep], lengths(remapped.genes))

        remapped.memberships <- split(
            unlist(remapped.genes),
            factor(remapped.sets, seq_along(set.memberships))
        )

        # Filling the values.
        j <- length(species.set.info[[cur.species]]) + 1L
        species.set.info[[cur.species]][[j]] <- data.frame(
            name=set.names,
            description=set.descriptions
        )

        species.collections[[cur.species]][[j]] <- data.frame(
            title=current$title,
            description=current$description,
            maintainer=current$maintainer,
            source=current$source
        )

        species.set.membership[[cur.species]][[j]] <- remapped.memberships
    }
}

#################################
# Looping across species.

library(gesel)
output.dir <- "_built"
unlink(output.dir, recursive=TRUE)
dir.create(output.dir, showWarnings=FALSE)

for (species in names(species.set.info)) {
    collections <- do.call(rbind, species.collections[[species]])
    stopifnot(anyDuplicated(collections$title) == 0)

    prepareDatabaseFiles(
        species=species,
        collections=collections,
        set.info=species.set.info[[species]],
        set.membership=species.set.membership[[species]],
        num.genes=species.ngenes[[species]],
        path=output.dir
    )
}

payload <- capture.output(print(sessionInfo()))
write(c("<details>", "<summary>Session information</summary>", "", "```", payload, "```", "</details>"), file="_session")
