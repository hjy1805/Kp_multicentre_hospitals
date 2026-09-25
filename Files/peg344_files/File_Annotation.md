# peg-344 Annotation and Presence/Absence Files

This directory contains reference sequence and presence/absence data used for the identification and annotation of **peg-344 (`pagO`)**, a putative PhoPQ-activated integral membrane protein.

Three files are provided:

| File | Description |
|---|---|
| `peg344.fasta` | Nucleotide sequence used for peg-344 annotation/reference. |
| `BAH65947.1_protein.faa` | Reference Peg-344/PagO protein sequence. |
| `peg344_presence.tsv` | Binary presence/absence calls for peg-344 across bacterial isolates. |

## 1. `peg344.fasta`

This FASTA file contains the nucleotide sequence corresponding to **peg-344 / pagO**.

The sequence header is:

`>KP1_p023 pagO putative PhoPQ-activated integral membrane protein 12338:13240 forward`

The header contains the following information:

| Field | Description |
|---|---|
| `KP1_p023` | Feature/locus identifier. |
| `pagO` | Gene annotation. |
| `putative PhoPQ-activated integral membrane protein` | Annotated protein function. |
| `12338:13240` | Coordinates of the annotated sequence in the corresponding reference sequence. |
| `forward` | Coding-strand orientation. |

The nucleotide sequence can be used as a reference for sequence-based identification of peg-344 in bacterial genomes.

## 2. `BAH65947.1_protein.faa`

This FASTA file contains the reference amino-acid sequence for the Peg-344/PagO protein.

The sequence header is:

`>BAH65947.1 putative PhoPQ-activated integral membrane protein (plasmid) [Klebsiella pneumoniae subsp. pneumoniae NTUH-K2044]`

| Field | Description |
|---|---|
| `BAH65947.1` | Protein accession identifier. |
| Protein annotation | Putative PhoPQ-activated integral membrane protein. |
| Location | Annotated as plasmid encoded in the reference record. |
| Organism | *Klebsiella pneumoniae* subsp. *pneumoniae* NTUH-K2044. |

This protein sequence provides a reference for protein-level identification or confirmation of Peg-344 homologues.

## 3. `peg344_presence.tsv`

This tab-separated file reports the presence or absence of peg-344 across the analysed bacterial isolates.

### Columns

| Column | Description |
|---|---|
| `strain` | Unique identifier of the bacterial isolate/strain. |
| `peg344` | Binary indicator of peg-344 presence or absence. |

### Coding

- `1` = peg-344 detected/present
- `0` = peg-344 not detected/absent

Example:

```text
strain    peg344
DKPB001   0
DKPB002   0
DKPB003   0
DKPB009   1
```

In this example, peg-344 was detected in isolate `DKPB009` but not in `DKPB001`, `DKPB002`, or `DKPB003`.

## Intended use

These files provide the sequence references and isolate-level binary calls used for **peg-344 annotation and downstream genomic/phenotypic analyses**.

The nucleotide and protein FASTA files provide reference sequences for identifying peg-344, while `peg344_presence.tsv` provides the resulting isolate-level presence/absence information that can be integrated with phenotype or clinical metadata for downstream statistical analyses.

## File formats

- `.fasta` — nucleotide FASTA
- `.faa` — amino-acid FASTA
- `.tsv` — tab-separated presence/absence matrix
