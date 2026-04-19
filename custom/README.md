# Contributing custom gene sets

Make a [pull request](https://github.com/gesel-inc/feedstock/pulls) to add a GMT file containing your gene sets to the `custom/` directory.

The GMT file should be uncompressed.
Its name should only have alphanumeric characters, underscores or dashes, and should end in a `*.gmt` suffix.
Each line should represent a gene set and contains at least 2 tab-separated values.
The first value should be the name of the gene set, the second value should be its description, and all subsequent values should be identifiers for genes in the set.
We suggest using stable identifiers like Entrez or Ensembl instead of symbols.
Use `gesel::mapGenesByName` to check that the provided gene identifiers are known to **gesel**.

The GMT file should be accompanied by a JSON file (same name but `.gmt` replaced by `.json`).
This should contain a JSON object with the following fields:

- `title`: the title of the collection.
  This should not contain tabs or newlines.
- `description`: the description of the collection.
  This should not contain tabs or newlines.
- `species`: the NCBI taxonomy ID for the species, e.g., `"9606"` for human.
  This should be a string, even though the IDs themselves are usually integer.
- `maintainer`: the name of the maintainer of the collection.
- `source`: the source of the collection.
  This may reference an article or the code used to generate the collection, and is intended for human readers.
- `id`: the type of identifier.
  This should be one of `"entrez"`, `"ensembl"` or `"symbol"`; the former two are more reliable and preferred.
