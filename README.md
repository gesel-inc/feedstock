# Build gene sets to feed gesel

![Build genes](https://github.com/gesel-inc/feedstock/actions/workflows/build-genes.yaml/badge.svg)
![Build GO](https://github.com/gesel-inc/feedstock/actions/workflows/build-go.yaml/badge.svg)
![Build MSigDB](https://github.com/gesel-inc/feedstock/actions/workflows/build-msigdb.yaml/badge.svg)
![Build indices](https://github.com/gesel-inc/feedstock/actions/workflows/build-indices.yaml/badge.svg)

## Overview

This repository contains code to build the gene set database for **gesel**, the client-side gene set search interface.

- [`genes/build.R`](genes/build.R) will create gene identity files (Ensembl/Entrez identifiers, symbols) for various species,
  using `org.*.eg.db` and `EnsDb` resources from Bioconductor's AnnotationHub.
- [`go/build.R`](go/build.R) will extract gene ontology sets from Bioconductor's `org.*.eg.db` databases. 
- [`msigdb/build.R`](msigdb/build.R) will extract the MSigDB sets with permissive licenses.

All of these gene sets are assembled into Gesel database indices via the top-level [`build.R`](build.R) script.

Pre-computed database files are available on the [Releases page](https://github.com/gesel-inc/feedstock/releases).

## Instructions

Gesel maintainers can trigger a new build of any particular resource by simply updating the relevant `VERSION` file.
This will automatically run the appropriate workflow to create a new release for that resource.

Interested users can contribute their own gene sets by following [these instructions](custom/README.md).

## Further links

The specification for the database files is described [here](https://github.com/gesel-inc/gesel-spec).

Applications should use the [R](https://github.com/gesel-inc/gesel-R), [Python](https://github.com/gesel-inc/gesel-py)
or [Javascript](https://github.com/gesel-inc/gesel.js) clients to interface with the Gesel database indices.
