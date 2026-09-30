#!/usr/bin/env python3
"""Targeted phylogenetic validation of MAG-encoded K10944 homologs.

The workflow extracts every local K10944 coding sequence from the original MAG
contigs, translates it, and places it with named archaeal AmoA, bacterial AmoA,
and bacterial PmoA references retrieved from UniProt. The resulting tree is an
annotation-validation layer, not evidence of transcription or process rate.
"""

from __future__ import annotations

import csv
import datetime as dt
import io
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import urllib.parse
import urllib.request


PROJECT = Path(os.environ.get("CHILE_PROJECT_ROOT", Path(__file__).resolve().parents[1]))
OUT_ROOT = Path(
    os.environ.get(
        "CHILE_V20_8_FIGURE_ROOT",
        PROJECT / "analysis/V20_8_candidate_Figure5_amoA_validation_2026-08-10",
    )
)
PHYLO_DIR = OUT_ROOT / "phylogeny"
SOURCE_DIR = OUT_ROOT / "source_data"
DOCS_DIR = OUT_ROOT / "docs"
for directory in (PHYLO_DIR, SOURCE_DIR, DOCS_DIR):
    directory.mkdir(parents=True, exist_ok=True)

ANNOTATIONS = (
    PROJECT
    / "data/processed/metagenome_output/MAGs/dram/annotations.tsv"
)
TAXONOMY = PROJECT / "data/processed/metagenome_output/MAGs/gtdb_taxonomy.tsv"
GENOMES = PROJECT / "data/processed/metagenome_output/MAGs/genomes"

BIN_DIR = Path(os.environ.get("CHILE_PHYLO_BIN", "/usr/local/bin"))
MAFFT = shutil.which("mafft") or str(BIN_DIR / "mafft")
TRIMAL = shutil.which("trimal") or str(BIN_DIR / "trimal")
IQTREE = shutil.which("iqtree3") or str(BIN_DIR / "iqtree3")


REFERENCE_SPECS = [
    # Named representatives spanning major ammonia-oxidising archaeal lineages.
    ("Archaeal AmoA", "Nitrososphaera viennensis", "ammonia monooxygenase subunit A"),
    ("Archaeal AmoA", "Nitrosopumilus maritimus", "ammonia monooxygenase subunit A"),
    ("Archaeal AmoA", "Nitrosotalea devaniterrae", "ammonia monooxygenase subunit A"),
    ("Archaeal AmoA", "Candidatus Nitrosocaldus cavascurensis", "ammonia monooxygenase subunit A"),
    ("Archaeal AmoA", "Cenarchaeum symbiosum", "ammonia monooxygenase subunit A"),
    ("Archaeal AmoA", "Candidatus Nitrosocosmicus franklandianus", "ammonia monooxygenase subunit A"),
    # Canonical bacterial and comammox AmoA representatives.
    ("Bacterial AmoA", "Nitrosomonas europaea", "ammonia monooxygenase subunit A"),
    ("Bacterial AmoA", "Nitrosospira multiformis", "ammonia monooxygenase subunit A"),
    ("Bacterial AmoA", "Nitrosococcus oceani", "ammonia monooxygenase subunit A"),
    ("Bacterial AmoA", "Nitrospira inopinata", "ammonia monooxygenase subunit A"),
    # PmoA references distinguish methane from ammonia monooxygenase homologs.
    ("PmoA", "Methylococcus capsulatus", "particulate methane monooxygenase subunit A"),
    ("PmoA", "Methylosinus trichosporium", "particulate methane monooxygenase subunit A"),
    ("PmoA", "Methylocystis parvus", "particulate methane monooxygenase subunit A"),
    ("PmoA", "Methylomicrobium album", "particulate methane monooxygenase subunit A"),
    ("PmoA", "Methylacidiphilum infernorum", "particulate methane monooxygenase subunit A"),
    ("PmoA", "Methylocaldum szegediense", "particulate methane monooxygenase subunit A"),
]


CODON_TABLE = {
    "TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
    "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
    "TAT": "Y", "TAC": "Y", "TAA": "*", "TAG": "*",
    "TGT": "C", "TGC": "C", "TGA": "*", "TGG": "W",
    "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
    "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
    "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
    "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
    "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
    "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
    "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
    "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
    "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
    "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
    "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
    "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
}


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, rows: list[dict[str, object]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def parse_fasta(path: Path) -> dict[str, str]:
    records: dict[str, list[str]] = {}
    current = None
    with path.open(encoding="utf-8") as handle:
        for raw in handle:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                current = line[1:].split()[0]
                if current in records:
                    raise ValueError(f"Duplicate FASTA identifier {current} in {path}")
                records[current] = []
            elif current is None:
                raise ValueError(f"Sequence before FASTA header in {path}")
            else:
                records[current].append(line.upper())
    return {name: "".join(parts) for name, parts in records.items()}


def reverse_complement(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def translate(sequence: str) -> str:
    usable = sequence[: len(sequence) - (len(sequence) % 3)]
    peptide = "".join(CODON_TABLE.get(usable[i : i + 3], "X") for i in range(0, len(usable), 3))
    if peptide.endswith("*"):
        peptide = peptide[:-1]
    return peptide


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    with path.open("w", encoding="utf-8") as handle:
        for name, sequence in records:
            handle.write(f">{name}\n")
            for start in range(0, len(sequence), 70):
                handle.write(sequence[start : start + 70] + "\n")


def ko_tokens(value: str) -> set[str]:
    return set(re.findall(r"K\d{5}", value or ""))


def safe_id(value: str) -> str:
    return re.sub(r"[^A-Za-z0-9_.-]+", "_", value).strip("_")


def extract_local_queries() -> tuple[list[dict[str, object]], list[tuple[str, str]], list[dict[str, object]]]:
    annotations = read_tsv(ANNOTATIONS)
    if annotations and "" in annotations[0]:
        for row in annotations:
            row["gene_id"] = row.pop("")
    elif annotations:
        first_key = next(iter(annotations[0]))
        if first_key.startswith("Unnamed"):
            for row in annotations:
                row["gene_id"] = row.pop(first_key)
    taxonomy = {row["user_genome"]: row for row in read_tsv(TAXONOMY)}

    targets = [row for row in annotations if "K10944" in ko_tokens(row.get("ko_id", ""))]
    if len(targets) != 5:
        raise ValueError(f"Expected 5 local K10944 genes, found {len(targets)}")

    fasta_cache: dict[str, dict[str, str]] = {}
    metadata: list[dict[str, object]] = []
    sequences: list[tuple[str, str]] = []
    context_rows: list[dict[str, object]] = []

    for row in sorted(targets, key=lambda item: item["fasta"]):
        mag = row["fasta"]
        if mag not in fasta_cache:
            fasta_cache[mag] = parse_fasta(GENOMES / f"{mag}.fasta")
        scaffold = row["scaffold"]
        contig = fasta_cache[mag].get(scaffold)
        if contig is None:
            raise ValueError(f"Scaffold {scaffold} not found in {mag}.fasta")
        start = int(float(row["start_position"]))
        end = int(float(row["end_position"]))
        strand = int(float(row["strandedness"]))
        if not (1 <= start <= end <= len(contig)):
            raise ValueError(f"Invalid coordinates for {row['gene_id']}: {start}-{end}")
        nucleotide = contig[start - 1 : end]
        if strand == -1:
            nucleotide = reverse_complement(nucleotide)
        elif strand != 1:
            raise ValueError(f"Unexpected strand {strand} for {row['gene_id']}")
        peptide = translate(nucleotide)
        internal_stops = peptide.count("*")
        if internal_stops or not 175 <= len(peptide) <= 320:
            raise ValueError(
                f"Translation QC failed for {row['gene_id']}: "
                f"aa={len(peptide)}, internal_stops={internal_stops}"
            )

        tax = taxonomy[mag]
        sequence_id = f"QUERY_{mag}"
        local_class = (
            "Local archaeal K10944"
            if tax["Domain"] == "Archaea"
            else "Local bacterial K10944"
        )
        metadata.append(
            {
                "sequence_id": sequence_id,
                "source": "local_MAG",
                "reference_class": local_class,
                "accession": "",
                "MAG": mag,
                "gene_id": row["gene_id"],
                "organism": mag,
                "Domain": tax["Domain"],
                "phylum": tax["phylum"],
                "class": tax["class"],
                "scaffold": scaffold,
                "start": start,
                "end": end,
                "strand": strand,
                "length_aa": len(peptide),
                "pfam_hits": row.get("pfam_hits", ""),
                "uniprot_url": "",
            }
        )
        sequences.append((sequence_id, peptide))

        mag_context = [
            item
            for item in annotations
            if item.get("fasta") == mag
            and ko_tokens(item.get("ko_id", "")) & {"K10944", "K10945", "K10946"}
        ]
        for context in mag_context:
            for ko in sorted(ko_tokens(context.get("ko_id", "")) & {"K10944", "K10945", "K10946"}):
                context_rows.append(
                    {
                        "MAG": mag,
                        "Domain": tax["Domain"],
                        "phylum": tax["phylum"],
                        "class": tax["class"],
                        "KO": ko,
                        "gene_id": context["gene_id"],
                        "scaffold": context["scaffold"],
                        "start": int(float(context["start_position"])),
                        "end": int(float(context["end_position"])),
                        "strand": int(float(context["strandedness"])),
                        "same_scaffold_as_K10944": context["scaffold"] == scaffold,
                        "kegg_hit": context.get("kegg_hit", ""),
                        "pfam_hits": context.get("pfam_hits", ""),
                    }
                )
    return metadata, sequences, context_rows


def uniprot_query(taxon: str, protein_name: str) -> list[dict[str, str]]:
    query = (
        f'(taxonomy_name:"{taxon}" AND protein_name:monooxygenase '
        "AND length:[175 TO 330])"
    )
    params = urllib.parse.urlencode(
        {
            "query": query,
            "format": "tsv",
            "fields": "accession,protein_name,organism_name,lineage,length,reviewed,sequence",
            "size": 50,
        }
    )
    request = urllib.request.Request(
        f"https://rest.uniprot.org/uniprotkb/search?{params}",
        headers={"User-Agent": "Chile-metagenomic-amoA-validation/1.0"},
    )
    with urllib.request.urlopen(request, timeout=90) as response:
        text = response.read().decode("utf-8")
    return list(csv.DictReader(io.StringIO(text), delimiter="\t"))


def choose_reference(
    reference_class: str, taxon: str, protein_name: str, candidates: list[dict[str, str]]
) -> dict[str, str]:
    usable = []
    taxon_core = taxon.lower().replace("candidatus ", "")
    for row in candidates:
        organism = row.get("Organism", "")
        sequence = re.sub(r"\s+", "", row.get("Sequence", "")).upper()
        protein = row.get("Protein names", "")
        organism_lower = organism.lower().replace("candidatus ", "")
        if not sequence or any(term in organism_lower for term in ("uncultured", "metagenome")):
            continue
        if taxon_core not in organism_lower:
            continue
        protein_lower = protein.lower()
        if reference_class == "PmoA":
            accepted_name = "methane monooxygenase" in protein_lower
            accepted_subunit = any(
                term in protein_lower for term in ("subunit a", "a subunit", "beta subunit")
            )
        else:
            accepted_name = "ammonia monooxygenase" in protein_lower
            accepted_subunit = "subunit a" in protein_lower or "alpha subunit" in protein_lower
        if not accepted_name or not accepted_subunit:
            continue
        if not 175 <= len(sequence) <= 330:
            continue
        row = dict(row)
        row["_sequence"] = sequence
        row["_score"] = (
            0 if row.get("Reviewed", "").lower() == "reviewed" else 1,
            0 if not re.search(r"\bputative\b|\bpartial\b", protein, re.I) else 1,
            abs(len(sequence) - 215),
            row.get("Entry", ""),
        )
        usable.append(row)
    if not usable:
        raise RuntimeError(
            f"No suitable UniProt reference for {reference_class}: {taxon} / {protein_name}"
        )
    return sorted(usable, key=lambda row: row["_score"])[0]


def fetch_references() -> tuple[list[dict[str, object]], list[tuple[str, str]], list[dict[str, object]]]:
    selected: list[dict[str, object]] = []
    sequences: list[tuple[str, str]] = []
    candidates_audit: list[dict[str, object]] = []
    seen_accessions: set[str] = set()

    for reference_class, taxon, protein_name in REFERENCE_SPECS:
        candidates = uniprot_query(taxon, protein_name)
        choice = choose_reference(reference_class, taxon, protein_name, candidates)
        accession = choice["Entry"]
        if accession in seen_accessions:
            raise RuntimeError(f"Reference accession selected twice: {accession}")
        seen_accessions.add(accession)
        sequence_id = f"REF_{safe_id(reference_class)}_{accession}"
        selected.append(
            {
                "sequence_id": sequence_id,
                "source": "UniProt",
                "reference_class": reference_class,
                "accession": accession,
                "MAG": "",
                "gene_id": "",
                "organism": choice["Organism"],
                "Domain": "Archaea" if reference_class == "Archaeal AmoA" else "Bacteria",
                "phylum": "",
                "class": "",
                "scaffold": "",
                "start": "",
                "end": "",
                "strand": "",
                "length_aa": len(choice["_sequence"]),
                "pfam_hits": "",
                "uniprot_url": f"https://www.uniprot.org/uniprotkb/{accession}/entry",
            }
        )
        sequences.append((sequence_id, choice["_sequence"]))
        for candidate in candidates:
            candidates_audit.append(
                {
                    "requested_class": reference_class,
                    "requested_taxon": taxon,
                    "requested_protein": protein_name,
                    "selected": candidate.get("Entry") == accession,
                    "accession": candidate.get("Entry", ""),
                    "protein_name": candidate.get("Protein names", ""),
                    "organism": candidate.get("Organism", ""),
                    "lineage": candidate.get("Taxonomic lineage", ""),
                    "length_aa": candidate.get("Length", ""),
                    "reviewed": candidate.get("Reviewed", ""),
                }
            )
    return selected, sequences, candidates_audit


def run_logged(command: list[str], log_path: Path, stdout_path: Path | None = None) -> None:
    with log_path.open("w", encoding="utf-8") as log_handle:
        log_handle.write("COMMAND\n" + " ".join(command) + "\n\n")
        if stdout_path is None:
            result = subprocess.run(
                command,
                cwd=PHYLO_DIR,
                stdout=log_handle,
                stderr=subprocess.STDOUT,
                text=True,
                check=False,
            )
        else:
            with stdout_path.open("w", encoding="utf-8") as stdout_handle:
                result = subprocess.run(
                    command,
                    cwd=PHYLO_DIR,
                    stdout=stdout_handle,
                    stderr=log_handle,
                    text=True,
                    check=False,
                )
        log_handle.write(f"\nRETURN_CODE\n{result.returncode}\n")
    if result.returncode != 0:
        raise RuntimeError(f"Command failed ({result.returncode}); see {log_path}")


def main() -> None:
    for executable in (MAFFT, TRIMAL, IQTREE):
        if not Path(executable).is_file() and shutil.which(executable) is None:
            raise FileNotFoundError(f"Required executable not found: {executable}")

    local_metadata, local_sequences, context_rows = extract_local_queries()
    reference_metadata, reference_sequences, candidates_audit = fetch_references()
    all_metadata = local_metadata + reference_metadata
    all_sequences = local_sequences + reference_sequences

    metadata_fields = [
        "sequence_id", "source", "reference_class", "accession", "MAG", "gene_id",
        "organism", "Domain", "phylum", "class", "scaffold", "start", "end",
        "strand", "length_aa", "pfam_hits", "uniprot_url",
    ]
    write_tsv(SOURCE_DIR / "amoA_pmoA_sequence_metadata.tsv", all_metadata, metadata_fields)
    write_tsv(
        SOURCE_DIR / "amoA_local_gene_context.tsv",
        context_rows,
        [
            "MAG", "Domain", "phylum", "class", "KO", "gene_id", "scaffold",
            "start", "end", "strand", "same_scaffold_as_K10944", "kegg_hit", "pfam_hits",
        ],
    )
    write_tsv(
        SOURCE_DIR / "amoA_uniprot_reference_candidates.tsv",
        candidates_audit,
        [
            "requested_class", "requested_taxon", "requested_protein", "selected",
            "accession", "protein_name", "organism", "lineage", "length_aa", "reviewed",
        ],
    )

    local_fasta = PHYLO_DIR / "amoA_local_K10944_queries.faa"
    reference_fasta = PHYLO_DIR / "amoA_pmoA_uniprot_references.faa"
    combined_fasta = PHYLO_DIR / "amoA_pmoA_combined_unaligned.faa"
    aligned_fasta = PHYLO_DIR / "amoA_pmoA_mafft.faa"
    trimmed_fasta = PHYLO_DIR / "amoA_pmoA_mafft_trimal.faa"
    write_fasta(local_fasta, local_sequences)
    write_fasta(reference_fasta, reference_sequences)
    write_fasta(combined_fasta, all_sequences)

    run_logged(
        [MAFFT, "--localpair", "--maxiterate", "1000", "--thread", "2", str(combined_fasta)],
        PHYLO_DIR / "mafft.log",
        aligned_fasta,
    )
    run_logged(
        [TRIMAL, "-in", str(aligned_fasta), "-out", str(trimmed_fasta), "-automated1"],
        PHYLO_DIR / "trimal.log",
    )
    prefix = PHYLO_DIR / "amoA_pmoA_ml"
    run_logged(
        [
            IQTREE, "-s", str(trimmed_fasta), "-m", "MFP", "-B", "1000",
            "--alrt", "1000", "-T", "AUTO", "-seed", "42", "-redo",
            "--prefix", str(prefix),
        ],
        PHYLO_DIR / "iqtree.log",
    )

    treefile = Path(str(prefix) + ".treefile")
    if not treefile.is_file() or treefile.stat().st_size == 0:
        raise RuntimeError(f"IQ-TREE output missing: {treefile}")

    today = dt.date.today().isoformat()
    selected_counts = {}
    for row in reference_metadata:
        selected_counts[row["reference_class"]] = selected_counts.get(row["reference_class"], 0) + 1
    readme = [
        "# Targeted amoA/pmoA validation",
        "",
        f"Generated: {today}",
        "",
        "## Scope",
        "",
        "- All five local DRAM K10944 homologs were extracted from the original MAG contigs.",
        "- Coordinates are interpreted as 1-based inclusive and strand-aware.",
        "- Named UniProt references represent archaeal AmoA, bacterial AmoA, and PmoA.",
        "- MAFFT L-INS-i, TrimAl automated1, and IQ-TREE ModelFinder + 1,000 UFBoot/SH-aLRT were used.",
        "- A phylogenetic placement can refine homolog identity; it does not demonstrate expression or activity.",
        "",
        "## Counts",
        "",
        f"- Local K10944 queries: {len(local_metadata)}",
        f"- Archaeal AmoA references: {selected_counts.get('Archaeal AmoA', 0)}",
        f"- Bacterial AmoA references: {selected_counts.get('Bacterial AmoA', 0)}",
        f"- PmoA references: {selected_counts.get('PmoA', 0)}",
        "",
        "## Reproducibility",
        "",
        "Selected accessions and all candidate records are frozen in source_data.",
        "Command lines and return codes are recorded in phylogeny/*.log.",
    ]
    (DOCS_DIR / "AMOA_TARGETED_VALIDATION_README.md").write_text(
        "\n".join(readme) + "\n", encoding="utf-8"
    )
    print(f"Local queries: {len(local_metadata)}")
    print(f"References: {len(reference_metadata)}")
    print(f"Tree: {treefile}")


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise
