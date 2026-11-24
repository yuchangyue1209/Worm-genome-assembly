#!/usr/bin/env python3
"""Read-only evidence review of focal CAFE candidate families; standard library only.

Uses the selected OrthoFinder run, its frozen FASTAs, matching InterPro TSVs
and the focal release GFFs. Does not run CAFE, alignments, RNA or repeat tests.
Existing output directories are never overwritten.
"""
import argparse
from collections import Counter, defaultdict
import csv
import datetime as dt
import hashlib
import json
from pathlib import Path
from urllib.parse import unquote


FOCAL = ["S_cotti", "S_solidus"]
DEFAULT_OGS = ["OG0000875", "OG0000880", "OG0000999", "OG0001937"]
EXPECTED_PFAM = {
    "OG0000875": {"PF03645"},
    "OG0000880": {"PF04968", "PF04969"},
    "OG0000999": {"PF00268"},
    "OG0001937": {"PF00650"},
}


def write_table(path, fields, rows):
    with path.open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({key: "NA" if value is None else value for key, value in row.items()})


def read_table(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fields = reader.fieldnames
        if not fields or len(fields) != len(set(fields)):
            raise ValueError(f"Missing/duplicate TSV header: {path}")
        return fields, list(reader)


def file_sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024*1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_fasta(path):
    records, name, parts = {}, None, []
    def save():
        if name is not None:
            seq = "".join(parts).upper()
            if name in records or not seq:
                raise ValueError(f"Duplicate/empty FASTA record {name}: {path}")
            records[name] = seq
    with path.open() as handle:
        for line in handle:
            if line.startswith(">"):
                save()
                header = line[1:].split()
                if not header:
                    raise ValueError(f"Empty FASTA header: {path}")
                name, parts = header[0], []
            elif line.strip():
                if name is None:
                    raise ValueError(f"Sequence before FASTA header: {path}")
                parts.append("".join(line.split()))
        save()
    if not records:
        raise ValueError(f"Empty FASTA: {path}")
    return records


def attributes(text):
    return dict(item.split("=", 1) for item in text.split(";") if "=" in item)


def aliases(value):
    result = {value} if value else set()
    for prefix in ("rna-", "cds-", "transcript:", "mRNA:", "protein:"):
        if value.startswith(prefix):
            result.add(value[len(prefix):])
    return result


def union_bp(intervals):
    total, left, right = 0, None, None
    for start, end in sorted(set(intervals)):
        if left is None:
            left, right = start, end
        elif start <= right+1:
            right = max(right, end)
        else:
            total += right-left+1
            left, right = start, end
    return total+(right-left+1 if left is not None else 0)


def locations(path, wanted):
    """Explicit gene->transcript->CDS mapping, never infer gene IDs by stripping .t."""
    genes, transcripts, protein_transcripts = {}, {}, defaultdict(set)
    cds = defaultdict(set)
    with path.open() as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith("##FASTA"):
                break
            if not line.strip() or line.startswith("#"):
                continue
            c = line.rstrip("\r\n").split("\t")
            if len(c) != 9:
                raise ValueError(f"Malformed GFF line {number}: {path}")
            if c[2] not in {"gene", "mRNA", "transcript", "CDS"}:
                continue
            a = attributes(c[8])
            fid = unquote(a.get("ID", ""))
            parents = tuple(sorted(unquote(x) for x in a.get("Parent", "").split(",") if x))
            start, end = int(c[3]), int(c[4])
            if not 1 <= start <= end:
                raise ValueError(f"Invalid GFF coordinates: {path}:{number}")
            record = (c[0], start, end, c[6], parents)
            if c[2] in {"gene", "mRNA", "transcript"}:
                if not fid:
                    raise ValueError(f"Missing GFF ID: {path}:{number}")
                store = genes if c[2] == "gene" else transcripts
                if fid in store:
                    old = store[fid]
                    if (old[0], old[3], old[4]) != (record[0], record[3], record[4]):
                        raise ValueError(f"Conflicting GFF ID {fid}: {path}")
                    record = (old[0], min(old[1], start), max(old[2], end), old[3], old[4])
                store[fid] = record
            if c[2] in {"mRNA", "transcript", "CDS"}:
                for tag in ("ID", "Name", "transcript_id", "protein_id"):
                    for pid in aliases(unquote(a.get(tag, ""))) & wanted:
                        protein_transcripts[pid].update([fid] if c[2] != "CDS" else parents)
            if c[2] == "CDS":
                for parent in parents:
                    cds[parent].add((c[0], start, end, c[6]))
    result, intervals = {}, {}
    for pid in sorted(wanted):
        hits = protein_transcripts[pid] & transcripts.keys()
        if len(hits) != 1:
            result[pid] = dict(Mapping_status="Missing_or_ambiguous_transcript")
            continue
        tid = next(iter(hits))
        transcript = transcripts[tid]
        if len(transcript[4]) != 1 or transcript[4][0] not in genes:
            result[pid] = dict(Mapping_status="Missing_or_ambiguous_gene")
            continue
        gid = transcript[4][0]
        gene, segments = genes[gid], sorted(cds[tid])
        if not segments:
            result[pid] = dict(Mapping_status="No_CDS_records")
            continue
        if (transcript[0] != gene[0] or transcript[3] != gene[3] or
            not gene[1] <= transcript[1] <= transcript[2] <= gene[2] or
            any(s[0] != transcript[0] or s[3] != transcript[3] or
                not transcript[1] <= s[1] <= s[2] <= transcript[2] for s in segments)):
            raise ValueError(f"Conflicting gene/transcript/CDS coordinates: {pid} in {path}")
        result[pid] = dict(Mapping_status="Explicit_transcript_gene_CDS_mapping",
            Gene_ID=gid, Transcript_ID=tid, Sequence=gene[0], Start=gene[1], End=gene[2],
            Strand=gene[3], Gene_span_bp=gene[2]-gene[1]+1, CDS_segments=len(segments),
            CDS_union_bp=union_bp((s[1], s[2]) for s in segments))
        intervals[pid] = segments
    return result, intervals


def amino_acid_sequence(sequence):
    """Remove one terminal stop marker only, never internal/multiple stops."""
    if len(sequence) > 1 and sequence.endswith("*") and "*" not in sequence[:-1]:
        return sequence[:-1]
    return sequence


def annotations(path, sequences):
    """Verify exact residues, permitting one explicitly recorded terminal stop."""
    hits, seen = defaultdict(list), set()
    normalized = {pid: amino_acid_sequence(seq) for pid, seq in sequences.items()}
    digests = {pid: hashlib.md5(seq.encode()).hexdigest() for pid, seq in normalized.items()}
    with path.open() as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            c = line.rstrip("\r\n").split("\t")
            if c[0] not in sequences:
                continue
            if len(c) < 11:
                raise ValueError(f"Malformed target InterPro row: {path}:{number}")
            pid, length = c[0], len(normalized[c[0]])
            if c[1].lower() != digests[pid] or int(c[2]) != length:
                raise ValueError(f"InterPro/selected FASTA sequence mismatch after allowing one terminal stop: {pid} in {path}:{number}")
            start, end = int(c[6]), int(c[7])
            if not 1 <= start <= end <= length:
                raise ValueError(f"InterPro hit outside protein: {pid} in {path}:{number}")
            row = dict(Protein_ID=pid, Analysis=c[3], Signature=c[4], Description=c[5],
                Start_aa=start, End_aa=end, Score_reported=c[8], Status_reported=c[9],
                InterPro=c[11] if len(c) > 11 else "-",
                InterPro_description=c[12] if len(c) > 12 else "-",
                Sequence_match_rule=("Single_terminal_stop_removed" if normalized[pid] != sequences[pid]
                                     else "Exact_input_sequence"), InterPro_sequence_length=length)
            key = tuple(row.items())
            if key not in seen:
                hits[pid].append(row)
                seen.add(key)
    return hits


def overrides(items):
    result = {}
    for item in items:
        sp, path = item.split("=", 1)
        if sp not in FOCAL or sp in result:
            raise ValueError("Override must specify each focal species at most once: "+item)
        result[sp] = Path(path)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("/work/cyu/Schistocephalus_comparative_genomics"))
    parser.add_argument("--results", type=Path)
    parser.add_argument("--input-dir", type=Path)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--orthogroups", nargs="+", default=DEFAULT_OGS)
    parser.add_argument("--expected-species", type=int, default=11)
    parser.add_argument("--gff", action="append", default=[], metavar="SPECIES=PATH")
    parser.add_argument("--interpro", action="append", default=[], metavar="SPECIES=PATH")
    args = parser.parse_args()
    if args.out.exists():
        raise FileExistsError("Output already exists; select a new version: "+str(args.out))
    results = args.results or args.root/(
        "04_orthofinder/runs/exploratory11_scotti_filtered7_v1/input/OrthoFinder/Results_Nov09")
    input_dir = args.input_dir or results.parent.parent
    source = results/"Orthogroups/Orthogroups.tsv"
    fields, table = read_table(source)
    if fields[0] != "Orthogroup" or len(fields)-1 != args.expected_species or not set(FOCAL) <= set(fields[1:]):
        raise ValueError("Unexpected species/header in Orthogroups.tsv")
    wanted = set(args.orthogroups)
    selected = {r["Orthogroup"]: r for r in table if r["Orthogroup"] in wanted}
    if len(selected) != len(wanted) or sum(r["Orthogroup"] in wanted for r in table) != len(wanted):
        raise ValueError("Missing/duplicate target orthogroups")
    paths = {
        "S_cotti": dict(gff=Path("/work/shared/Scotti_release_v1.0_preNCBI_20261118/annotation/Scotti_nuclear.protein_coding.v1.0.gff3"),
            interpro=Path("/work/cyu/annotation/scotti-annotation/results/functional/03_interproscan/scotti_braker.tsv")),
        "S_solidus": dict(gff=Path("/work/shared/Ssolidus_release_v1.0_preNCBI_20260831/annotation/Ssolidus_nuclear.annotation.v1.0.gff3"),
            interpro=args.root/"05_functional_annotation/interpro/S_solidus.tsv"),
    }
    for kind in ("gff", "interpro"):
        for sp, path in overrides(getattr(args, kind)).items():
            paths[sp][kind] = path
    members, proteins, species_sequences, identical, counts = {}, {}, {}, {}, []
    input_paths = {"orthogroups": source}
    for sp in fields[1:]:
        fasta = input_dir/(sp+".faa")
        input_paths["fasta:"+sp] = fasta
        seqs = read_fasta(fasta)
        all_ids = []
        for og in sorted(wanted):
            ids = [x.strip() for x in selected[og][sp].split(",") if x.strip()]
            if len(ids) != len(set(ids)) or set(ids)-seqs.keys():
                raise ValueError(f"Duplicate/missing FASTA member: {sp}/{og}")
            members[sp, og] = ids
            all_ids.extend(ids)
            for pid in ids:
                proteins[sp, pid] = og
        if len(all_ids) != len(set(all_ids)):
            raise ValueError("Protein assigned to multiple target families: "+sp)
        species_sequences[sp] = {pid: seqs[pid] for pid in all_ids}
        if sp in FOCAL:
            identical[sp] = Counter(seqs.values())
    for og in sorted(wanted):
        counts.append(dict(Orthogroup=og, **{sp: len(members[sp, og]) for sp in fields[1:]}))
    rows, hit_rows, beds, warnings = [], [], defaultdict(list), []
    for sp in FOCAL:
        seqs = species_sequences[sp]
        loc, cds = locations(paths[sp]["gff"], set(seqs))
        hits = annotations(paths[sp]["interpro"], seqs)
        for kind, path in paths[sp].items():
            input_paths[kind+":"+sp] = path
        for pid, seq in sorted(seqs.items(), key=lambda x: (proteins[sp, x[0]], x[0])):
            amino_acids = amino_acid_sequence(seq)
            length = len(amino_acids)
            og = proteins[sp, pid]
            matches = hits.get(pid, [])
            pfam = {h["Signature"] for h in matches if h["Analysis"] == "Pfam"}
            ipr = {h["InterPro"] for h in matches if h["InterPro"].startswith("IPR")}
            expected = EXPECTED_PFAM.get(og, set())
            coverage = union_bp((h["Start_aa"], h["End_aa"]) for h in matches if h["Analysis"] == "Pfam")
            status = "Sequence_verified_records_present" if matches else "No_record_in_InterPro_TSV"
            flags = []
            if loc[pid]["Mapping_status"] != "Explicit_transcript_gene_CDS_mapping":
                flags.append(loc[pid]["Mapping_status"])
            if length < 100:
                flags.append("LT100aa_review_not_automatic_rejection")
            if amino_acids != seq:
                flags.append("Single_terminal_stop_marker_not_an_amino_acid")
            if seq.endswith("**"):
                flags.append("Multiple_terminal_stop_markers")
            if "*" in seq.rstrip("*"):
                flags.append("Internal_stop")
            if "X" in seq:
                flags.append("Contains_X")
            if not matches:
                flags.append("Annotation_missing_not_evidence_of_domain_absence")
            elif expected and not expected <= pfam:
                flags.append("Expected_Pfam_not_all_detected")
            if identical[sp][seq] > 1:
                flags.append("Identical_protein_sequence_review_loci")
            if loc[pid].get("CDS_union_bp") not in {None, length*3, length*3+3}:
                flags.append("CDS_length_review_phase_partiality")
            rows.append(dict(Species=sp, Orthogroup=og, Protein_ID=pid, **loc[pid],
                Length_aa=length, Input_length_including_terminal_stop=len(seq),
                Terminal_stop_symbols=len(seq)-len(seq.rstrip("*")), X_aa=seq.count("X"),
                Input_sequence_MD5=hashlib.md5(seq.encode()).hexdigest(),
                Amino_acid_sequence_MD5=hashlib.md5(amino_acids.encode()).hexdigest(),
                InterPro_sequence_match=matches[0]["Sequence_match_rule"] if matches else None,
                Identical_sequence_count_in_species_input=identical[sp][seq], Annotation_status=status,
                Pfam=";".join(sorted(pfam)) or "NA", InterPro=";".join(sorted(ipr)) or "NA",
                Expected_Pfam=";".join(sorted(expected)) or "NA",
                All_expected_Pfam_detected=(int(expected <= pfam) if matches and expected else None),
                Pfam_covered_aa=coverage if matches else None,
                Pfam_coverage_percent=round(100*coverage/length, 2) if matches else None,
                RNR_active_site_entry_detected=(int("IPR030475" in ipr) if og == "OG0000999" and matches else None),
                Flags=";".join(flags) or "No_screening_flags_not_full_validation"))
            for match in matches:
                hit_rows.append(dict(Species=sp, Orthogroup=og, **match))
            for chrom, left, right, strand in cds.get(pid, []):
                beds[sp].append((chrom, left-1, right, pid, 0, strand))
            for flag in flags:
                warnings.append(dict(Species=sp, Orthogroup=og, Protein_ID=pid, Flag=flag))
    summary = []
    for sp in FOCAL:
        for og in sorted(wanted):
            rr = [r for r in rows if (r["Species"], r["Orthogroup"]) == (sp, og)]
            lengths = [r["Length_aa"] for r in rr]
            summary.append(dict(Species=sp, Orthogroup=og, Proteins=len(rr),
                Min_length_aa=min(lengths) if lengths else None,
                Max_length_aa=max(lengths) if lengths else None,
                With_InterProScan_records=sum(r["Annotation_status"] == "Sequence_verified_records_present" for r in rr),
                With_Pfam=sum(r["Pfam"] != "NA" for r in rr),
                With_InterPro_entry=sum(r["InterPro"] != "NA" for r in rr),
                With_single_terminal_stop_marker=sum(r["Length_aa"] != r["Input_length_including_terminal_stop"] for r in rr),
                With_all_expected_Pfam=sum(r["All_expected_Pfam_detected"] == 1 for r in rr),
                With_explicit_gene_CDS_mapping=sum(r["Mapping_status"] == "Explicit_transcript_gene_CDS_mapping" for r in rr),
                With_identical_sequence_elsewhere_in_input=sum(r["Identical_sequence_count_in_species_input"] > 1 for r in rr)))
    # No output directory is created until all selected inputs/hits are validated.
    provenance = {kind: dict(path=str(path.resolve()), sha256=file_sha256(path))
                  for kind, path in input_paths.items()}
    args.out.mkdir(parents=True)
    protein_fields = ["Species", "Orthogroup", "Protein_ID", "Gene_ID", "Transcript_ID", "Mapping_status",
        "Sequence", "Start", "End", "Strand", "Gene_span_bp", "CDS_segments", "CDS_union_bp",
        "Length_aa", "Input_length_including_terminal_stop", "Terminal_stop_symbols", "InterPro_sequence_match",
        "X_aa", "Identical_sequence_count_in_species_input", "Annotation_status", "Pfam",
        "InterPro", "Expected_Pfam", "All_expected_Pfam_detected", "Pfam_covered_aa", "Pfam_coverage_percent",
        "RNR_active_site_entry_detected", "Flags"]
    # Missing mapping fields must remain visibly unavailable, never blank/zero.
    write_table(args.out/"protein_sequence_match_audit.tsv", ["Species", "Orthogroup", "Protein_ID",
        "Input_length_including_terminal_stop", "Length_aa", "Terminal_stop_symbols", "Input_sequence_MD5",
        "Amino_acid_sequence_MD5", "InterPro_sequence_match"], rows)
    rows = [{k: row.get(k) for k in protein_fields} for row in rows]
    write_table(args.out/"priority_protein_review.tsv", protein_fields, rows)
    write_table(args.out/"family_copy_counts_allspecies.tsv", ["Orthogroup"]+fields[1:], counts)
    write_table(args.out/"family_review_summary.tsv", list(summary[0]), summary)
    write_table(args.out/"feature_hits.tsv", ["Species", "Orthogroup", "Protein_ID", "Analysis", "Signature",
        "Description", "Start_aa", "End_aa", "Score_reported", "Status_reported", "InterPro", "InterPro_description",
        "Sequence_match_rule", "InterPro_sequence_length"], hit_rows)
    write_table(args.out/"review_flags.tsv", ["Species", "Orthogroup", "Protein_ID", "Flag"], warnings)
    membership_rows = [dict(Species=sp, Orthogroup=og, Protein_count=len(members[sp, og]),
        Protein_IDs=",".join(members[sp, og])) for sp in fields[1:] for og in sorted(wanted)]
    write_table(args.out/"family_membership.tsv", ["Species", "Orthogroup", "Protein_count", "Protein_IDs"], membership_rows)
    write_table(args.out/"priority_family_selection.tsv",
        ["Species", "Orthogroup", "Use_for_RNA_repeat_validation"],
        [dict(Species=sp, Orthogroup=og, Use_for_RNA_repeat_validation=1)
         for sp in FOCAL for og in sorted(wanted) if members[sp, og]])
    fasta_dir = args.out/"all_species_fastas"
    fasta_dir.mkdir()
    for og in sorted(wanted):
        with (fasta_dir/(og+".faa")).open("x") as handle:
            for sp in fields[1:]:
                for pid in members[sp, og]:
                    handle.write(f">{sp}|{pid}\n{species_sequences[sp][pid]}\n")
    for sp in FOCAL:
        with (args.out/(sp+".priority_cds.bed")).open("x") as handle:
            for row in sorted(set(beds[sp])):
                handle.write("\t".join(map(str, row))+"\n")
    (args.out/"provenance.json").write_text(json.dumps(dict(
        created_UTC=dt.datetime.now(dt.timezone.utc).isoformat(), results=str(results.resolve()),
        orthogroups=sorted(wanted), expected_species=args.expected_species, input_files=provenance,
        source_files_modified=False, all_observed_target_InterPro_rows_MD5_and_length_verified=True,
        interpro_sequence_comparison_rule="exact_amino_acids_allow_one_terminal_stop_marker_only",
        focal_proteins_with_single_terminal_stop_marker=sum(r["Length_aa"] != r["Input_length_including_terminal_stop"] for r in rows),
        focal_proteins_without_InterPro_records=sum(r["Annotation_status"] == "No_record_in_InterPro_TSV" for r in rows),
        RNA_support_tested=False, repeat_overlap_tested=False, gene_trees_inferred=False), indent=2)+"\n")
    (args.out/"interpretation_notes.txt").write_text(
        "Exploratory family review, not a new significance test or a final functional annotation.\n"
        "OG identifiers and membership refer only to the supplied OrthoFinder run. No CAFE rerun occurs.\n"
        "InterPro target-row MD5, length and feature coordinates must match amino-acid residues in the frozen protein input; one terminal stop marker may be removed for comparison only. Internal/multiple stops and residue substitutions are not normalized away.\n"
        "Length_aa and Pfam coverage exclude that single terminal marker. Stored input lengths, MD5 values and comparison rules are recorded in protein_sequence_match_audit.tsv. Source files and exported family FASTAs retain their original sequences, including terminal markers.\n"
        "Missing InterPro rows mean unavailable evidence, not proven domain absence.\n"
        "Expected Pfam IDs are descriptive checks for the four selected baseline families, not universal completeness criteria.\n"
        "IPR030475 detection is a reported active-site feature, not proof of enzymatic activity; its absence is not proof of inactivity.\n"
        "Pfam coverage is the union of reported Pfam match intervals / protein length, not a full-length homology test.\n"
        "Feature hits include domain, family, superfamily and site entries; not every InterPro hit is a domain.\n"
        "Identical protein sequences can arise from genuine distinct loci; no automatic deduplication is performed.\n"
        "Gene coordinates are GFF 1-based inclusive; CDS BED intervals are 0-based half-open, deduplicated within each protein.\n"
        "GFF mapping requires an explicit gene/transcript/CDS relationship; missing/ambiguous mappings remain NA.\n"
        "RNA depth, uniquely mapped expression, TE overlap, assembly haplotigs, tandem duplication and gene-tree support remain to be tested.\n"
        "No RNA/repeat evidence from previous candidate sets is transferred automatically.\n")
    print("PASS:", args.out)
    print("Focal proteins reviewed:", len(rows))
    print("Family\tCotti\tSolidus")
    for og in sorted(wanted):
        print(og, len(members["S_cotti", og]), len(members["S_solidus", og]), sep="\t")


if __name__ == "__main__":
    try:
        main()
    except (OSError, ValueError, KeyError) as exc:
        raise SystemExit("STOPPED: "+str(exc))
