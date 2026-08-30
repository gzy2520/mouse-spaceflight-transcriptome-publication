# Checks before making the release public

The scientific release tree and its reproduction contract are complete. This directory has its own clean Git history and can be published independently of the exploratory research repository. Two author-level decisions remain outside the analysis code:

1. Add the license selected by the authors and the final citation or DOI metadata.
2. Decide whether the historical absolute catalog path in `results/tables/06_source_expression_read_audit.csv` may remain public. It is part of the approved table and is therefore preserved here. Sanitizing it would create a revised result table with a new checksum and should not be done silently.

The release-authored code and documentation were scanned for development-tool traces and credentials. No such markers, private keys, access tokens or project email addresses were found. Package-maintainer contact fields inside `renv.lock` are upstream CRAN metadata.
