#!/usr/bin/env python3
"""Collect the selected eleven-species Table 1 without modifying source data.

Standard library only. No new BUSCO/RepeatMasker/TE jobs are launched.
Outputs are new TSV/JSON files, not an Excel workbook. Missing evidence is NA,
never zero. Use --inputs with an edited selected_sources.json to resolve paths
or provide consistently classified transposon-protein ID lists.
"""
import argparse
import csv
import datetime as dt
import gzip
import hashlib
import json
import math
from pathlib import Path
import re
import sys
from collections import Counter, defaultdict
from urllib.parse import unquote


SPECIES = ["S_solidus", "S_cotti", "L_intestinalis", "D_latus",
           "Spirometra_Aus1", "S_proliferum", "E_multilocularis",
           "H_microstoma", "E_granulosus", "T_solium", "T_multiceps"]
VERSIONS = ["v1.0 pre-NCBI (20260831)", "v1.0 pre-NCBI (20261118)",
            "local published annotation", "GCA_900617775.1", "GCA_051168835.1",
            "GCA_053814205.1 / PRJDB8966", "GCA_000469725.3", "PRJEB124 / WBPS19",
            "GCA_021556725.1", "GCA_001870725.1", "PRJNA307624 / WBPS19"]
EXPECTED = dict(zip(SPECIES, [9796, 9234, 27256, 19930, 11064, 28993,
                             10656, 10139, 9985, 11854, 12502]))
EXON_METHOD = "transcript_exon_associations_v2"
METRICS = [
    ("version", "Assembly / annotation version", None),
    ("source", "Source", None),
    ("scope", "Selected assembly scope", None),
    ("genome_bp", "Genome size (Mbp)", 3),
    ("gc", "GC content (%)", 2),
    ("scaffolds", "Total scaffolds / sequences", 0),
    ("longest_bp", "Longest scaffold (Mbp)", 3),
    ("over_1mb", "Number of scaffolds >1 Mbp", 0),
    ("n_bp", "Total number of N bases (bp)", 0),
    ("n_percent", "N bases (%)", 3),
    ("n50_bp", "Scaffold N50 (Mbp)", 3),
    ("n75_bp", "Scaffold N75 (Mbp)", 3),
    ("l50", "Scaffold L50", 0), ("l75", "Scaffold L75", 0),
    ("genes", "Gene loci (unique GFF gene IDs)", 0),
    ("coding_genes", "Protein-coding gene loci (with CDS)", 0),
    ("transposon_proteins", "Number of transposon proteins (representative set)", 0),
    ("mean_gene_bp", "Average protein-coding gene span (bp)", 1),
    ("mean_all_gene_bp", "Average gene span (all GFF gene IDs; bp)", 1),
    ("exons", "Exon-transcript associations (all annotated transcripts)", 0),
    ("coding_exons", "Exon-transcript associations (protein-coding transcripts)", 0),
    ("exons_per_gene", "Mean exons per gene (representative isoform)", 2),
    ("mean_exon_bp", "Average protein-coding exon length (bp; transcript-level)", 1),
    ("isoforms", "Protein isoforms", 0),
    ("primary", "Representative proteins", 0),
    ("repeat_percent", "Bases masked by RepeatMasker (%)", 2),
]
for MODE in ["genome", "protein"]:
    for KEY, LABEL in [("C", "Complete"), ("S", "Single-copy"), ("D", "Duplicated"),
                       ("F", "Fragmented"), ("M", "Missing")]:
        METRICS.append((MODE+"_"+KEY, LABEL+" "+MODE+" BUSCO (%)", 1))
    METRICS.append((MODE+"_n", MODE.capitalize()+" BUSCO groups (metazoa_odb12)", 0))


def read_tsv(path):
    with path.open(newline="") as f:
        return list(csv.DictReader(f, delimiter="\t"))


def write_tsv(path, fields, rows):
    with path.open("x", newline="") as f:
        w = csv.DictWriter(f, fields, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)


def opengz(path):
    with path.open("rb") as f:
        compressed = f.read(2) == b"\x1f\x8b"
    return gzip.open(path, "rt") if compressed else path.open()


def fasta_stats(path, nucleotide=False):
    """Streaming sequence fingerprints match audit_final_busco_inputs.py."""
    lengths, digests = {}, {}
    bases = Counter()
    pid, digest, size = None, None, 0

    def finish():
        if pid is not None:
            if not size:
                raise ValueError("Empty FASTA record: "+pid)
            lengths[pid], digests[pid] = size, digest.hexdigest()

    with opengz(path) as f:
        for line in f:
            if line.startswith(">"):
                finish()
                parts = line[1:].split()
                if not parts:
                    raise ValueError("Empty FASTA header")
                pid = parts[0]
                if pid in lengths:
                    raise ValueError("Duplicate FASTA ID: "+pid)
                digest, size = hashlib.sha256(), 0
            else:
                # Process long unwrapped FASTA lines in bounded pieces.
                for start in range(0, len(line), 1024*1024):
                    seq = "".join(line[start:start+1024*1024].split()).upper()
                    if not seq:
                        continue
                    if pid is None:
                        raise ValueError("Sequence before FASTA header")
                    digest.update(seq.encode("ascii"))
                    size += len(seq)
                    if nucleotide:
                        bases.update({b: seq.count(b) for b in "ACGTN"})
        finish()
    if not lengths:
        raise ValueError("No FASTA sequences")
    canonical = hashlib.sha256()
    for name, value in sorted(digests.items()):
        canonical.update(name.encode()+b"\t"+value.encode()+b"\n")
    result = {"records": len(lengths), "residues": sum(lengths.values()),
              "id_sequence_sha256": canonical.hexdigest(), "lengths": lengths,
              "sequence_sha256": digests}
    if nucleotide:
        ordered = sorted(lengths.values(), reverse=True)
        total, accumulated = sum(ordered), 0
        result.update(genome_bp=total, scaffolds=len(ordered), longest_bp=ordered[0],
                      over_1mb=sum(x > 1000000 for x in ordered), n_bp=bases["N"],
                      n_percent=100*bases["N"]/total)
        valid = sum(bases[x] for x in "ACGT")
        result["gc"] = 100*(bases["G"]+bases["C"])/valid if valid else None
        for i, length in enumerate(ordered, 1):
            accumulated += length
            for pct in [50, 75]:
                if accumulated*100 >= total*pct and "l"+str(pct) not in result:
                    result["l"+str(pct)], result["n"+str(pct)+"_bp"] = i, length
    return result


def attributes(text):
    result = {}
    for item in text.split(";"):
        if "=" in item:
            key, value = item.split("=", 1)
            result[key] = value  # Split list values BEFORE URL decoding.
    return result


def aliases(identifier):
    vals = {identifier}
    for prefix in ["rna-", "cds-", "transcript:", "mRNA:", "protein:"]:
        if identifier.startswith(prefix):
            vals.add(identifier[len(prefix):])
    return vals


def is_transcript_feature(kind):
    return (kind in {"mRNA", "transcript"} or kind.endswith("RNA")
            or kind.endswith("_transcript"))


def is_pseudogene_feature(kind, attr):
    if kind == "pseudogene" or kind.startswith("pseudogenic_"):
        return True
    for key in ("gene_biotype", "transcript_biotype", "biotype"):
        if "pseudogene" in unquote(attr.get(key, "")).lower():
            return True
    return any(key in attr and unquote(attr[key]).lower() not in {"false", "0", "no"}
               for key in ("pseudo", "pseudogene"))


def annotation_stats(path, genome_lengths, primary_ids, diagnostics=None, exon_diagnostics=None):
    gene_spans, transcript_gene, exon_counts = {}, {}, Counter()
    loci, edges, cds_parents = {}, defaultdict(set), set()
    protein_to_transcripts = defaultdict(set)
    exon_records = []
    transcript_exon_coordinates = defaultdict(set)
    feature_types, pseudogene_nodes = {}, set()
    coordinates_bad = []
    with opengz(path) as f:
        for line_num, line in enumerate(f, 1):
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            cols = line.rstrip("\r\n").split("\t")
            if len(cols) != 9:
                raise ValueError(f"GFF line {line_num}: expected nine tab-separated fields")
            seq, _, kind, left, right, _, strand, _, attr = cols
            left, right = int(left), int(right)
            if (seq not in genome_lengths or not 1 <= left <= right <= genome_lengths[seq]):
                if len(coordinates_bad) < 5:
                    coordinates_bad.append(f"{seq}:{left}-{right}")
            a = attributes(attr)
            fid = unquote(a.get("ID", ""))
            parents = [unquote(x) for x in a.get("Parent", "").split(",") if x]
            if fid:
                if fid in feature_types and feature_types[fid] != kind:
                    raise ValueError("Conflicting feature types for GFF ID: "+fid)
                feature_types[fid] = kind
                edges[fid].update(parents)
                if is_pseudogene_feature(kind, a):
                    pseudogene_nodes.add(fid)
            # GFF3 permits a single discontinuous feature to occupy several rows.
            # Count its ID once and use its bounding genomic span. Reject changes
            # of reference sequence, strand, type, Parent or identity attributes.
            if kind in ["gene", "pseudogene", "mRNA", "transcript"] or (fid and is_transcript_feature(kind)):
                if not fid:
                    raise ValueError("Missing "+kind+" ID")
                signature = (seq, strand, kind, tuple(sorted(parents)))
                if fid in loci:
                    old = loci[fid]
                    if old["signature"] != signature:
                        raise ValueError("Incompatible multipart feature ID (sequence/strand/type/Parent): "+fid)
                    for tag in ["Name", "gene_biotype", "biotype", "locus_tag"]:
                        if tag in old["attributes"] and tag in a and old["attributes"][tag] != a[tag]:
                            raise ValueError("Incompatible multipart feature attribute "+tag+": "+fid)
                    old["left"], old["right"] = min(old["left"], left), max(old["right"], right)
                    old["parts"] += 1
                    old["attributes"].update(a)
                else:
                    loci[fid] = dict(signature=signature, left=left, right=right,
                                     parts=1, attributes=dict(a))
            if kind == "gene":
                gene_spans[fid] = loci[fid]["right"]-loci[fid]["left"]+1
            if fid and is_transcript_feature(kind):
                transcript_gene[fid] = parents
                for key in ["ID", "Name", "transcript_id", "protein_id"]:
                    for val in aliases(unquote(a.get(key, ""))):
                        if val:
                            protein_to_transcripts[val].add(fid)
            if kind == "CDS":
                cds_parents.update(parents)
                if is_pseudogene_feature(kind, a):
                    pseudogene_nodes.update(parents)
                for key in ["protein_id", "Name", "ID"]:
                    for val in aliases(unquote(a.get(key, ""))):
                        if val:
                            protein_to_transcripts[val].update(parents)
            if kind == "exon":
                exon_records.append((right-left+1, parents))
                for parent in parents:
                    transcript_exon_coordinates[parent].add((seq, left, right, strand))
    if coordinates_bad:
        raise ValueError("GFF/selected genome coordinate mismatch: "+", ".join(coordinates_bad))
    def ancestors(node):
        pending, seen = [node], set()
        while pending:
            current = pending.pop()
            if current in seen:
                continue
            seen.add(current)
            pending.extend(edges.get(current, ()))
        return seen

    coding_nodes = set()
    for parent in cds_parents:
        lineage = ancestors(parent)
        if not lineage & pseudogene_nodes:
            coding_nodes.update(lineage)
    coding_genes = set(gene_spans) & coding_nodes
    # Count associations with explicit transcript features only. A gene-level
    # exon mirrored at mRNA level is not another exon-transcript association.
    # Coordinates shared by different isoforms count once per transcript.
    all_transcripts = set(transcript_gene)
    coding_transcripts = {tid for tid in all_transcripts & coding_nodes
                          if ancestors(tid) & coding_genes
                          and not ancestors(tid) & pseudogene_nodes}
    pseudo_transcripts = {tid for tid in all_transcripts
                          if ancestors(tid) & pseudogene_nodes}
    all_exons = sum(len(transcript_exon_coordinates[tid]) for tid in all_transcripts)
    coding_exons = [end-start+1 for tid in coding_transcripts
                    for _, start, end, _ in transcript_exon_coordinates[tid]]
    missing_coding_exons = {tid for tid in coding_transcripts
                            if not transcript_exon_coordinates.get(tid)}
    coding_exons_complete = bool(coding_transcripts) and not missing_coding_exons
    exon_counts.update({tid: len(coords) for tid, coords in transcript_exon_coordinates.items()})
    result = {"genes": len(gene_spans) if gene_spans else None,
              "coding_genes": len(coding_genes) if gene_spans and cds_parents else None,
              "mean_gene_bp": sum(gene_spans[g] for g in coding_genes)/len(coding_genes) if coding_genes else None,
              "mean_all_gene_bp": sum(gene_spans.values())/len(gene_spans) if gene_spans else None,
              "exons": all_exons if all_exons else None,
              "coding_exons": len(coding_exons) if coding_exons and coding_exons_complete else None,
              "mean_exon_bp": sum(coding_exons)/len(coding_exons) if coding_exons and coding_exons_complete else None,
              "exons_per_gene": None}
    if exon_diagnostics is not None:
        parent_links = Counter(parent for _, parents in exon_records for parent in parents)
        raw_tx_links = sum(parent_links[t] for t in all_transcripts)
        exon_diagnostics.update(method=EXON_METHOD, status="PARSED", reason="",
            raw_exon_records=len(exon_records), raw_exon_parent_links=sum(parent_links.values()),
            raw_transcript_parent_links=raw_tx_links,
            raw_gene_parent_links=sum(n for p, n in parent_links.items() if feature_types.get(p) == "gene"),
            raw_pseudogene_parent_links=sum(n for p, n in parent_links.items() if feature_types.get(p) == "pseudogene"),
            raw_unresolved_parent_links=sum(n for p, n in parent_links.items() if p not in feature_types),
            raw_other_parent_links=sum(n for p, n in parent_links.items()
                if p in feature_types and p not in all_transcripts
                and feature_types[p] not in {"gene", "pseudogene"}),
            exon_records_without_Parent=sum(not parents for _, parents in exon_records),
            explicit_transcripts=len(all_transcripts), coding_transcripts=len(coding_transcripts),
            pseudogene_transcripts=len(pseudo_transcripts),
            all_transcript_exon_associations=all_exons,
            coding_transcript_exon_associations_observed=len(coding_exons),
            pseudogene_transcript_exon_associations=sum(len(transcript_exon_coordinates[t]) for t in pseudo_transcripts),
            duplicate_transcript_parent_links_removed=raw_tx_links-all_exons,
            coding_transcripts_without_exons=len(missing_coding_exons),
            coding_exon_metrics_complete=coding_exons_complete)
    if diagnostics is not None:
        for fid, locus in sorted(loci.items()):
            if locus["parts"] > 1 or (locus["signature"][2] == "gene" and fid not in coding_genes):
                diagnostics.append(dict(feature_id=fid, kind=locus["signature"][2],
                    sequence=locus["signature"][0], strand=locus["signature"][1],
                    rows=locus["parts"], start=locus["left"], end=locus["right"],
                    span_bp=locus["right"]-locus["left"]+1,
                    cds_supported_gene=fid in coding_genes,
                    note="multipart_feature" if locus["parts"] > 1 else "gene_without_CDS_support"))
    selected, used_genes = [], set()
    problems = []
    for pid in sorted(primary_ids):
        transcripts = set()
        for val in aliases(pid):
            transcripts.update(protein_to_transcripts[val])
        transcripts &= transcript_gene.keys()
        if len(transcripts) != 1:
            problems.append(pid+": missing/ambiguous transcript mapping")
            continue
        tid = next(iter(transcripts))
        parents = transcript_gene[tid]
        if len(parents) != 1 or parents[0] not in gene_spans or parents[0] in used_genes:
            problems.append(pid+": missing/ambiguous/duplicate gene mapping")
            continue
        if tid not in coding_transcripts:
            problems.append(pid+": transcript is not CDS-supported protein-coding or is pseudogene-associated")
            continue
        if exon_counts[tid] == 0:
            problems.append(pid+": transcript has no exon features")
            continue
        used_genes.add(parents[0])
        selected.append(exon_counts[tid])
    if primary_ids and not problems:
        result["exons_per_gene"] = sum(selected)/len(selected)
    return result, problems


def repeat_percent(path, lengths, exclude_ids=(), diagnostics=None):
    """RepeatMasker .out: union of inclusive intervals, not overlapping sums."""
    excluded = set(exclude_ids)
    if excluded - lengths.keys() or not lengths.keys() - excluded:
        raise ValueError("Invalid RepeatMasker sequence exclusions")
    intervals = defaultdict(list)
    records = Counter()
    with opengz(path) as f:
        for line in f:
            parts = line.split()
            if not parts or not parts[0].isdigit():
                continue
            if len(parts) < 14:
                raise ValueError("Malformed RepeatMasker .out record")
            name, start, end, remaining = parts[4], int(parts[5]), int(parts[6]), parts[7]
            if name not in lengths or not 1 <= start <= end <= lengths[name]:
                raise ValueError("RepeatMasker coordinates differ from selected genome: "+name)
            if remaining.startswith("(") and end+int(remaining.strip("()")) != lengths[name]:
                raise ValueError("RepeatMasker query length differs from selected genome: "+name)
            intervals[name].append((start, end))
            records["excluded" if name in excluded else "retained"] += 1
    covered = Counter()
    for name, ranges in intervals.items():
        category = "excluded" if name in excluded else "retained"
        left = right = None
        for start, end in sorted(ranges):
            if left is None:
                left, right = start, end
            elif start <= right+1:
                right = max(right, end)
            else:
                covered[category] += right-left+1
                left, right = start, end
        if left is not None:
            covered[category] += right-left+1
    if not intervals:
        raise ValueError("No repeat records; cannot establish a successful zero-repeat run")
    denominator = sum(size for name, size in lengths.items() if name not in excluded)
    if diagnostics is not None:
        diagnostics.update(masked_bp=covered["retained"], excluded_masked_bp=covered["excluded"],
                           retained_repeat_records=records["retained"],
                           excluded_repeat_records=records["excluded"], denominator_bp=denominator)
    return 100*covered["retained"]/denominator


def repeat_evidence(src, selected, original=None):
    """Project only explicitly excluded records; every retained sequence must match."""
    exclusions = src.get("repeatmasker_exclude_ids", [])
    if not isinstance(exclusions, list) or any(not isinstance(x, str) or not x for x in exclusions):
        raise ValueError("repeatmasker_exclude_ids must be a list of nonempty sequence IDs")
    if len(exclusions) != len(set(exclusions)):
        raise ValueError("Duplicate RepeatMasker exclusion IDs")
    if exclusions and original is None:
        raise ValueError("RepeatMasker exclusions require the original input FASTA")
    if original is not None:
        old_ids, new_ids = set(original["lengths"]), set(selected["lengths"])
        if old_ids - new_ids != set(exclusions) or new_ids - old_ids:
            raise ValueError("Original/selected FASTA ID difference does not equal explicit RepeatMasker exclusions")
        changed = [name for name in sorted(new_ids)
                   if original["sequence_sha256"][name] != selected["sequence_sha256"][name]]
        if changed:
            raise ValueError("Retained RepeatMasker input sequences changed: "+", ".join(changed[:5]))
    evidence = dict(source=src["repeatmasker_out"], selected_genome=src["genome"],
                    original_genome=src.get("repeatmasker_genome", ""),
                    excluded_ids=",".join(exclusions),
                    excluded_genome_bp=sum(original["lengths"][x] for x in exclusions) if original else 0,
                    selected_id_sequence_sha256=selected["id_sequence_sha256"],
                    original_id_sequence_sha256=original["id_sequence_sha256"] if original else "",
                    sequence_identity_status="all_retained_IDs_and_sequences_verified" if original else
                    "coordinates_and_lengths_only_original_FASTA_not_supplied")
    evidence["repeat_percent"] = repeat_percent(Path(src["repeatmasker_out"]),
        (original or selected)["lengths"], exclusions, evidence)
    return evidence


def default_sources(root):
    std = root/"02_standardized_data"
    sources = {}
    for i, sp in enumerate(SPECIES):
        sources[sp] = dict(version=VERSIONS[i],
                           source="Selected local release" if i < 2 else "Selected public dataset: "+VERSIONS[i],
                           scope="Nuclear only" if i < 2 else "mt exclusion not verified",
                           genome=str(std/"genomes"/(sp+".fasta")),
                           gff=str(std/"gff"/(sp+".gff3")),
                           isoforms=str(std/"proteins_raw"/(sp+".faa")),
                           primary=str(std/"proteins_primary"/(sp+".faa")),
                           repeatmasker_out="", repeatmasker_genome="", repeatmasker_exclude_ids=[],
                           transposon_ids="", transposon_method="")
    for sp, name, date in [("S_solidus", "Ssolidus", "20260831"),
                            ("S_cotti", "Scotti", "20261118")]:
        release = Path("/work/shared")/(name+"_release_v1.0_preNCBI_"+date)
        gff = ("Ssolidus_nuclear.annotation.v1.0.gff3" if sp == "S_solidus"
               else "Scotti_nuclear.protein_coding.v1.0.gff3")
        sources[sp].update(genome=str(release/"assembly"/(name+"_nuclear.fa")),
                           gff=str(release/"annotation"/gff),
                           isoforms=str(release/"proteins"/(name+"_proteins.all_isoforms.v1.0.fa")),
                           primary=str(release/"proteins"/(name+"_proteins.longest_isoform.v1.0.fa")))
    sources["S_cotti"]["repeatmasker_out"] = "/work/cyu/assembly/worm/scotti/repeats/02_repeatmasker/Scotti_nuclear.mt_filtered.ragtag.provisional.fa.out"
    sources["S_solidus"].update(
        repeatmasker_out="/work/cyu/annotation/ssol-annotation/results/repeats/worm.q10_r500k_no_contig_ec.nuclear.final.fa.out",
        repeatmasker_genome="/work/cyu/assembly/worm/hifiasm/purge_q10/worm.q10_r500k_no_contig_ec.nuclear.final.fa",
        repeatmasker_exclude_ids=["scaffold_61"])
    return sources


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=Path("/work/cyu/Schistocephalus_comparative_genomics"))
    p.add_argument("--out", type=Path, required=True)
    p.add_argument("--inputs", type=Path, help="Edited selected_sources.json from a previous collection")
    p.add_argument("--audit-dir", type=Path, help="Exact final_busco_audit_* directory; otherwise latest named audit")
    p.add_argument("--busco-dir", type=Path, help="run_missing_genome_busco.py output directory")
    args = p.parse_args()
    sources = json.loads(args.inputs.read_text()) if args.inputs else default_sources(args.root)
    if set(sources) != set(SPECIES):
        raise ValueError("Input JSON must specify exactly the eleven selected species")
    audit = args.audit_dir
    if audit is None:
        choices = sorted(x for x in (args.root/"03_qc").glob("final_busco_audit_*")
                         if (x/"selected_inputs.tsv").is_file())
        audit = choices[-1] if choices else None
    busco_dir = args.busco_dir or args.root/"03_qc/final_genome_busco_odb12_v1"
    # The recorded audit pins identities; do not silently select historical BUSCOs.
    identities, busco = {}, {}
    if audit and (audit/"selected_inputs.tsv").is_file():
        identities = {(r["species"], r["mode"]): r for r in read_tsv(audit/"selected_inputs.tsv")}
        review = audit/"eleven_species.genome_and_protein.review.tsv"
        if review.is_file():
            busco = {(r["species"], r["mode"]): r for r in read_tsv(review)
                     if r["review_status"] == "matched_current_input_check_run_time_provenance"}
    combined = busco_dir/"eleven_species.genome_and_protein.odb12.tsv"
    if combined.is_file():
        if not (busco_dir/"COMPLETE.json").is_file():
            raise ValueError("Combined BUSCO TSV has no COMPLETE.json: "+str(busco_dir))
        rows = read_tsv(combined)
        expected_keys = {(sp, mode) for sp in SPECIES for mode in ["genome", "proteins"]}
        if len(rows) != 22 or {(r["species"], r["mode"]) for r in rows} != expected_keys:
            raise ValueError("Combined BUSCO table must contain exactly 22 unique assessments")
        busco = {(r["species"], r["mode"]): r for r in rows}
        if not identities:
            spec_file = busco_dir/"run_spec.json"
            if spec_file.is_file():
                # Always request the original audit rather than guessing its schema.
                print("NOTE: supply --audit-dir to verify the selected FASTA identities", file=sys.stderr)
    args.out.mkdir(parents=True, exist_ok=False)
    (args.out/"selected_sources.json").write_text(json.dumps(sources, indent=2)+"\n")
    values, audit_rows, warnings, gff_audit, repeat_audit, exon_audit = {}, [], [], [], [], []
    cache = {}
    def get_fasta(path, nucleotide):
        resolved = Path(path).resolve()
        key = (str(resolved), nucleotide)
        if key not in cache:
            cache[key] = fasta_stats(resolved, nucleotide)
        return cache[key]

    for sp in SPECIES:
        print("COLLECT:", sp, flush=True)
        src = sources[sp]
        exon_evidence = dict(species=sp, source=src.get("gff", ""), method=EXON_METHOD,
                             status="UNAVAILABLE", reason="Selected genome unavailable")
        row = {k: None for k, _, _ in METRICS}
        origins, issues = {}, {}
        for key in ["version", "source", "scope"]:
            row[key], origins[key] = src.get(key, ""), "selected_sources.json"
        fasta = {}
        for key, nucleotide in [("genome", True), ("primary", False), ("isoforms", False)]:
            path = src.get(key, "")
            try:
                if not path:
                    raise ValueError("Source path not supplied")
                fasta[key] = get_fasta(path, nucleotide)
                if key == "genome":
                    for metric, _, _ in METRICS:
                        if metric in fasta[key]:
                            row[metric], origins[metric] = fasta[key][metric], path
                else:
                    metric = "primary" if key == "primary" else "isoforms"
                    row[metric], origins[metric] = fasta[key]["records"], path
            except (OSError, ValueError, UnicodeError) as exc:
                issues[key] = str(exc)
        if row["primary"] is not None and row["primary"] != EXPECTED[sp]:
            warnings.append(dict(species=sp, issue="Representative count differs from selected baseline", detail=str(row["primary"])))
        if "primary" in fasta and "isoforms" in fasta:
            absent = set(fasta["primary"]["lengths"])-set(fasta["isoforms"]["lengths"])
            if absent:
                warnings.append(dict(species=sp, issue="Primary IDs absent from all-isoform FASTA", detail=str(len(absent))))
        if "genome" in fasta:
            try:
                feature_diagnostics = []
                metrics, mapping_issues = annotation_stats(Path(src["gff"]), fasta["genome"]["lengths"],
                                                          set(fasta.get("primary", {}).get("lengths", {})),
                                                          feature_diagnostics, exon_evidence)
                gff_audit.extend(dict(species=sp, source=src["gff"], **r) for r in feature_diagnostics)
                for key, val in metrics.items():
                    row[key], origins[key] = val, src["gff"]
                if mapping_issues:
                    issues["exons_per_gene"] = "; ".join(mapping_issues[:5])+f"; total affected={len(mapping_issues)}"
                if exon_evidence.get("coding_transcripts_without_exons"):
                    for key in ("coding_exons", "mean_exon_bp"):
                        issues[key] = ("Exon features missing for "
                            +str(exon_evidence["coding_transcripts_without_exons"])
                            +" CDS-supported coding transcripts; partial totals not reported")
            except (OSError, ValueError, KeyError) as exc:
                issues["gff"] = str(exc)
                exon_evidence.update(status="UNAVAILABLE", reason=str(exc))
            repeat = src.get("repeatmasker_out", "")
            if repeat:
                try:
                    original = get_fasta(src["repeatmasker_genome"], True) if src.get("repeatmasker_genome") else None
                    evidence = repeat_evidence(src, fasta["genome"], original)
                    repeat_audit.append(dict(species=sp, **evidence))
                    row["repeat_percent"] = evidence["repeat_percent"]
                    origins["repeat_percent"] = repeat+"; excluded_ids="+(evidence["excluded_ids"] or "none")
                    warnings.append(dict(species=sp, issue="RepeatMasker input provenance",
                                         detail=evidence["sequence_identity_status"]))
                    print(f"REPEAT: {sp} {evidence['repeat_percent']:.2f}% "
                          f"masked_bp={evidence['masked_bp']} "
                          f"excluded_ids={evidence['excluded_ids'] or 'none'}", flush=True)
                except (OSError, ValueError) as exc:
                    issues["repeat_percent"] = str(exc)
        exon_audit.append(exon_evidence)
        te_path, method = src.get("transposon_ids", ""), src.get("transposon_method", "")
        if te_path and method and "primary" in fasta:
            try:
                ids = [x.strip() for x in Path(te_path).read_text().splitlines() if x.strip()]
                if len(ids) != len(set(ids)) or set(ids)-set(fasta["primary"]["lengths"]):
                    raise ValueError("TE ID list contains duplicates or IDs outside the representative set")
                row["transposon_proteins"], origins["transposon_proteins"] = len(ids), te_path+"; method="+method
            except (OSError, ValueError) as exc:
                issues["transposon_proteins"] = str(exc)
        for mode, file_key, prefix in [("genome", "genome", "genome"), ("proteins", "primary", "protein")]:
            record, identity = busco.get((sp, mode)), identities.get((sp, mode))
            problem = "No input-matched metazoa_odb12 report; complete/audit BUSCO first"
            if record and identity and file_key in fasta:
                if fasta[file_key]["id_sequence_sha256"] != identity["id_sequence_sha256"]:
                    problem = "Selected FASTA differs from the BUSCO-audited input"
                elif record["n"] != "672" or record["dataset_date"] != "2026-05-22" or record["busco_version"] != "6.0.0":
                    problem = "BUSCO dataset/version differs from pinned metazoa_odb12 / 672 / 2026-05-22 / BUSCO 6.0.0"
                else:
                    reported_input = record.get("selected_input")
                    if reported_input and get_fasta(reported_input, mode == "genome")["id_sequence_sha256"] != identity["id_sequence_sha256"]:
                        raise ValueError("Combined BUSCO table belongs to a different selected input: "+sp+" "+mode)
                    if record.get("input_provenance") == "frozen_input_run_hash_verified":
                        run_record = json.loads((busco_dir/"logs"/(sp+".run.json")).read_text())
                        if run_record["id_sequence_sha256"] != identity["id_sequence_sha256"]:
                            raise ValueError("Frozen BUSCO run identity differs from selected input: "+sp)
                    percentages = [float(record[key]) for key in ["C", "S", "D", "F", "M"]]
                    if any(not math.isfinite(x) or not 0 <= x <= 100 for x in percentages):
                        raise ValueError("Invalid BUSCO percentage: "+sp+" "+mode)
                    for key in ["C", "S", "D", "F", "M", "n"]:
                        metric = prefix+"_"+key
                        row[metric] = int(record[key]) if key == "n" else float(record[key])
                        origins[metric] = record.get("summary") or str(combined if combined.is_file() else audit/"eleven_species.genome_and_protein.review.tsv")
                    warnings.append(dict(species=sp, issue=mode+" BUSCO provenance", detail=record.get("input_provenance", "Current ID/sequence match; historical run-time checksum not recorded")))
            elif record and not identity:
                problem = "BUSCO table present but original --audit-dir identity table is unavailable"
            for key in ["C", "S", "D", "F", "M", "n"]:
                if row[prefix+"_"+key] is None:
                    issues[prefix+"_"+key] = problem
        for key, label, _ in METRICS:
            reason = issues.get(key, "")
            if row[key] is None and not reason:
                if key in ["genes", "coding_genes", "mean_gene_bp", "mean_all_gene_bp", "exons", "coding_exons", "mean_exon_bp", "exons_per_gene"]:
                    reason = issues.get("gff") or issues.get("genome") or "Required GFF feature/representative transcript mapping unavailable"
                elif key == "transposon_proteins":
                    reason = "Supply a representative-protein ID list AND a shared classification method; TE32 is sensitivity-only, not an eleven-species annotation"
                elif key == "repeat_percent":
                    reason = "No compatible RepeatMasker .out provided; soft-masked lowercase is not substituted"
                else:
                    reason = issues.get("primary" if key == "primary" else "isoforms" if key == "isoforms" else "genome", "Source unavailable")
            audit_rows.append(dict(species=sp, metric=label, value=row[key], status="MISSING" if row[key] is None else "AVAILABLE",
                                   source=origins.get(key, ""), reason=reason))
        values[sp] = row
        if src.get("scope") != "Nuclear only":
            warnings.append(dict(species=sp, issue="Assembly scope", detail="Mitochondrial exclusion is not established; verify and update the selected FASTA before calling these nuclear-only statistics."))
    fields = ["Assembly / annotation feature"]+SPECIES
    table = []
    for key, label, decimals in METRICS:
        out = {fields[0]: label}
        for sp in SPECIES:
            val = values[sp][key]
            if val is None:
                out[sp] = "n.a."
            elif decimals is None:
                out[sp] = str(val)
            else:
                if key in ["genome_bp", "longest_bp", "n50_bp", "n75_bp"]:
                    val /= 1000000
                out[sp] = f"{val:.{decimals}f}"
        table.append(out)
    write_tsv(args.out/"Table1_11species.tsv", fields, table)
    af = ["species", "metric", "value", "status", "source", "reason"]
    write_tsv(args.out/"metric_sources.tsv", af, audit_rows)
    missing = [r for r in audit_rows if r["status"] == "MISSING"]
    write_tsv(args.out/"missing_metrics.tsv", af, missing)
    write_tsv(args.out/"warnings.tsv", ["species", "issue", "detail"], warnings)
    write_tsv(args.out/"gff_feature_audit.tsv", ["species", "feature_id", "kind", "sequence", "strand",
        "rows", "start", "end", "span_bp", "cds_supported_gene", "note", "source"], gff_audit)
    write_tsv(args.out/"exon_transcript_audit.tsv", ["species", "method", "status", "source", "reason",
        "raw_exon_records", "raw_exon_parent_links", "raw_transcript_parent_links",
        "raw_gene_parent_links", "raw_pseudogene_parent_links", "raw_unresolved_parent_links",
        "raw_other_parent_links", "exon_records_without_Parent", "explicit_transcripts",
        "coding_transcripts", "pseudogene_transcripts", "all_transcript_exon_associations",
        "coding_transcript_exon_associations_observed", "pseudogene_transcript_exon_associations",
        "duplicate_transcript_parent_links_removed", "coding_transcripts_without_exons",
        "coding_exon_metrics_complete"], exon_audit)
    write_tsv(args.out/"repeatmasker_release_audit.tsv", ["species", "source", "original_genome",
        "selected_genome", "excluded_ids", "excluded_genome_bp", "sequence_identity_status",
        "original_id_sequence_sha256", "selected_id_sequence_sha256", "retained_repeat_records",
        "excluded_repeat_records", "masked_bp", "excluded_masked_bp", "denominator_bp", "repeat_percent"], repeat_audit)
    (args.out/"metric_definitions.txt").write_text(
        "All statistics come from selected files, not older paper tables.\n"
        "Species order: first six Diphyllobothriidea; last five Cyclophyllidea.\n"
        "Genome size and Nxx include Ns; Mbp=1,000,000 bp. N content uses all assembly bases.\n"
        "GC denominator is A+C+G+T only (ambiguous symbols excluded).\n"
        "Scaffold counts are FASTA record counts, not inferred chromosome counts.\n"
        "N50/N75 and L50/L75 use descending scaffold lengths and 50/75% cumulative assembly size.\n"
        "Mitochondrial exclusion is guaranteed only for explicitly verified nuclear inputs; see scope row.\n"
        "Gene loci count unique IDs of GFF feature type gene (not ncRNA_gene/pseudogene); composition may differ among annotations.\n"
        "Protein-coding gene loci require CDS ancestry linked to a gene ID without explicit pseudogene ancestry, flags or biotypes; this is not inferred from the representative FASTA count. Non-coding-classified gene features are reported separately in gff_feature_audit.tsv.\n"
        "Repeated gene/transcript IDs with compatible reference, strand, type, Parents and identity attributes are treated as multipart features. Conflicting IDs still fail validation.\n"
        "Gene length is the inclusive bounding genomic span (min start to max end for multipart genes), not CDS length or sum of parts. Coding and all-gene averages are separate.\n"
        "Exon method: "+EXON_METHOD+". Exon totals count distinct (transcript ID, sequence, start, end, strand) associations, not raw GFF rows or globally unique genomic intervals. Gene/pseudogene-level exon rows are excluded even when they mirror transcript-level rows.\n"
        "All-transcript associations include explicit mRNA/transcript and RNA transcript features, including pseudogene-associated transcripts. Coding associations and their length mean require CDS-to-transcript-to-gene ancestry without explicit pseudogene features, pseudo flags or pseudogene biotypes.\n"
        "Exon length is inclusive end-start+1, weighted once per distinct transcript association. Shared exons count separately for different transcripts. Repeated coordinates within the same transcript count once.\n"
        "Coding exon total/mean are NA if coding transcripts lack exon features; incomplete totals are not silently reported. Unrecognized or unresolved exon Parents are not promoted to transcripts. Raw counts and omitted parent levels are retained in exon_transcript_audit.tsv.\n"
        "Mean exons/gene uses only selected representative transcripts; every representative must map uniquely to a distinct gene, be CDS-supported and non-pseudogene-associated, and have exon features, otherwise NA.\n"
        "No CDS intervals are silently substituted for absent exon features.\n"
        "Transposon proteins count unique classified IDs in the representative set; supply the same classification procedure across species.\n"
        "RepeatMasker percentage is the union of .out intervals / selected assembly length; matching names/lengths do not by themselves prove sequence identity.\n"
        "Explicit RepeatMasker projection requires the original FASTA: removed IDs must exactly match the exclusion list, and every retained ID/sequence must match the selected assembly. All original .out coordinates are validated before excluded records are omitted.\n"
        "For solidus the recorded fish contaminant scaffold_61 is excluded only after this sequence-level verification. See repeatmasker_release_audit.tsv for exact masked bp, excluded bp and the release denominator; no RepeatMasker rerun is performed.\n"
        "The 32 cotti strong-TE candidates remain in the 9,234 baseline; do not call them a uniform transposon-protein census.\n"
        "BUSCO columns are two modes of metazoa_odb12 (genome and representative proteins), not odb10 versus odb12.\n"
        "BUSCO results must match the original audit's normalized FASTA IDs/sequences and pinned dataset/version. Historical results may lack run-time hashes.\n"
        "BUSCO rounded percentages need not sum exactly to 100. Protein BUSCO and genome BUSCO measure different things.\n"
        "NA means missing/unverified evidence, never zero. This collector launches no annotation or BUSCO jobs.\n")
    (args.out/"collection.json").write_text(json.dumps(dict(created_UTC=dt.datetime.now(dt.timezone.utc).isoformat(),
        root=str(args.root), audit_dir=str(audit) if audit else None, busco_dir=str(busco_dir),
        exon_method=EXON_METHOD,
        source_files_modified=False, missing_cells=len(missing), species=SPECIES,
        fingerprints={sp:{key:{k:v for k,v in cache[(str(Path(src[key]).resolve()), key == "genome")].items()
                              if k in ["records", "residues", "id_sequence_sha256"]}
                          for key in ["genome", "primary", "isoforms"]
                          if (str(Path(src[key]).resolve()), key == "genome") in cache}
                      for sp,src in sources.items()}), indent=2)+"\n")
    print("WROTE:", args.out/"Table1_11species.tsv")
    print("Missing/unverified cells:", len(missing), "(see missing_metrics.tsv)")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError) as exc:
        print("STOPPED:", exc, file=sys.stderr)
        sys.exit(1)
