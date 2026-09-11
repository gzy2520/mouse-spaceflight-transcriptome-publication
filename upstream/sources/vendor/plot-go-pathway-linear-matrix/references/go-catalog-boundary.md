# Seven-pathway GO catalog boundary

## Completeness definition

Treat the catalog as complete relative to this declared ontology contract,
not relative to proteins observed in any one experiment.

- Source official Gene Ontology `go-basic.obo`.
- Retain active canonical GO terms and their `alt_id` identifiers.
- Expand `exact` rules to the seed only.
- Expand `structural` rules through reverse `is_a` and `part_of` relations.
- Expand `pathway_related` rules through reverse `is_a`, `part_of`,
  `regulates`, `positively_regulates`, and `negatively_regulates` relations.
- Permit the same term to map to multiple pathways.
- Assign unmatched input terms to `Others`.

The NHEJ hierarchy requires an explicit boundary. Classical
`GO:0097680` and alternative `GO:0097681` are expanded independently.
Generic `GO:0006303`, its generic regulation branch, telomere-protection
terms, and generic gap filling `GO:0061674` map to both AEJ and NHEJ because
the GO terms do not resolve the two mechanisms.

## Frozen snapshot

The bundled `assets/seven_pathway_go_term_catalog.csv` was generated from:

- GO release: `releases/2026-07-26`
- OBO SHA-256:
  `b08d45b268b8c24ccb2513dbbbc7d4df9f6521c099b413f79eb31e06e0fa3bcc`
- Active canonical GO terms examined: 38,092
- Catalog: 229 canonical terms, 21 alternative IDs, 266 GO-ID–pathway rows

Canonical terms by pathway:

| Pathway | Total | BP | CC | MF |
|---|---:|---:|---:|---:|
| BER | 53 | 20 | 1 | 32 |
| NER | 33 | 16 | 15 | 2 |
| MMR | 39 | 8 | 8 | 23 |
| FA | 5 | 4 | 1 | 0 |
| HR | 88 | 77 | 3 | 8 |
| AEJ | 10 | 10 | 0 | 0 |
| NHEJ | 17 | 10 | 6 | 1 |

## Deliberate exclusions

Do not assign a broad term solely because a known pathway protein carries it.
Keep generic `DNA repair`, DNA-damage response, DNA binding, helicase,
polymerase, nuclease, ligase, chromatin, replication, and protein-complex
terms in `Others` unless the GO term or ontology branch is explicitly
pathway-specific.

Do not include telomeric D-loop subclasses through generic D-loop binding:
the bundled rule intentionally maps only `GO:0062037` exactly. Do not infer
pathway activation or inhibition from positive or negative GO regulation
names.

## Refresh procedure

Download a dated official `go-basic.obo`, run
`scripts/build_seven_pathway_go_catalog.R`, inspect the generated provenance
and pathway counts, and rerun the completeness regression tests. Do not
silently replace the bundled catalog without reporting the GO release and
hash.
