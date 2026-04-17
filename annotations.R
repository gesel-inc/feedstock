# For Entrez, we use whatever the latest org.*.eg.db is available in AnnotationHub for each species.
# For Ensembl, we use whatever the latest EnsDb is available in AnnotationHub for each species.

annotations <- list(
    `10090` = list(
        ensdb = "AH119358", # query(ahub, c("EnsDb", "Mus musculus"))
        orgdb = "AH121954"  # query(ahub, "org.Mm.eg.db")
    ),
    `9606` = list(
        ensdb = "AH119325", # query(ahub, c("EnsDb", "Homo sapiens"))
        orgdb = "AH121953"  # query(ahub, "org.Hs.eg.db")
    ),
    `6239` = list(
        ensdb = "AH119242", # query(ahub, c("EnsDb", "Caenorhabditis elegans"))
        orgdb = "AH121958"  # query(ahub, "org.Ce.eg.db")
    ),
    `10116` = list(
        ensdb = "AH119437", # query(ahub, c("EnsDb", "Rattus norvegicus"))
        orgdb = "AH121956"  # query(ahub, "org.Rn.eg.db")
    ),
    `7227` = list(
        ensdb = "AH119285", # query(ahub, c("EnsDb", "Drosophila melanogaster"))
        orgdb = "AH121952"  # query(ahub, "org.Dm.eg.db")
    ),
    `7955` = list(
        ensdb = "AH119289", # query(ahub, c("EnsDb", "Danio rerio"))
        orgdb = "AH121961"  # query(ahub, "org.Dr.eg.db")
    ),
    `9598` = list(
        ensdb = "AH119429", # query(ahub, c("EnsDb", "Pan troglodytes"))
        orgdb = "AH121949"  # query(ahub, "org.Pt.eg.db")
    )
)
