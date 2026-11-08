# Public copy

This directory was copied separately from the experimental workspace. The original datasets, indexes, source trees, and recorded experiment files were preserved.

- The current query implementation is the final V4 source. All six `code/src/` file hashes are recorded in `source_hashes.json`.
- The original small models and query banks are copied without changes and checked by `data/manifest.json`.
- The original fitting/reference tools are isolated under `code/model/`; they are not substituted for the current query engine.
- Job scripts use relative public-directory paths and the current Python interpreter. Machine-specific run controllers, private server paths, temporary logs, and binaries are omitted.
- The case script accepts the already bundled MAP query file only when its content exactly matches the generated reachable targets.
- Main CSVs are copied from the V2, V3, and V4 deliveries. Original timeouts and statuses are retained. Large case pair/map outputs are reproducible with the included analysis script.
- The large original DIMACS inputs are a release attachment, with archive and per-file SHA-256 checks. Generated synthetic inputs and indexes are ignored by Git.

This organization changes paths, packaging, and documentation. It does not change the query or construction algorithms, the fitting procedure, or the experimental parameter values.
