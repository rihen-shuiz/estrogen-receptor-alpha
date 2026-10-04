"""Where is each deposited residue in its canonical receptor? Prefer the
deposited UniProt segment, align the declared sequence (not a gappy CA-only
sequence), validate raw pairs, then join coordinates. Native ESR2 numbering
is related to ESR1 by a sequence equivalence, not copied numeric bounds.
"""
from Bio import Align
from Bio.Align import substitution_matrices
from Bio.SeqUtils import seq1

from .read_structure import category_rows


def validate_pairs(pairs):
    """Reject duplicates and backward correspondence before building a dictionary."""
    query = [q for q, c in pairs]
    canonical = [c for q, c in pairs]
    if len(query) != len(set(query)) or len(canonical) != len(set(canonical)):
        raise ValueError("Mapping is not one-to-one: repeated query/canonical position")
    for previous, current in zip(pairs, pairs[1:]):
        if current[0] <= previous[0] or current[1] <= previous[1]:
            raise ValueError("Mapping order reverses in original correspondence")


def alignment_pairs(reference, query, mode="local"):
    aligner = Align.PairwiseAligner()
    aligner.mode = mode
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    alignment = next(iter(aligner.align(reference, query)), None)
    if alignment is None:
        return []
    pairs = []
    for (r0, r1), (q0, q1) in zip(*alignment.aligned):
        if r1 - r0 != q1 - q0:
            raise ValueError("Alignment block has unequal correspondence lengths")
        for offset in range(int(r1 - r0)):
            pairs.append((int(q0 + offset + 1), int(r0 + offset + 1)))
    validate_pairs(pairs)
    return pairs


def receptor_equivalence(sequences):
    """Return native ESR2 -> ESR1 coordinates using the frozen full sequences."""
    pairs = alignment_pairs(sequences["esr1"], sequences["esr2"], "global")
    return dict(pairs)


def receptor_accessions(receptor, settings):
    """Return the canonical accession and explicitly evidenced secondary aliases."""
    return [settings[f"{receptor}_uniprot"],
            *settings.get("receptor_accession_aliases", {}).get(receptor, {})]


def identify_entity(entity_id, sequence, data, sequences, settings, polymer_type):
    refs = category_rows(data, "_struct_ref")
    accessions = [r.get("pdbx_db_accession", "") for r in refs if r.get("entity_id") == entity_id]
    for receptor in ("esr1", "esr2"):
        accepted = receptor_accessions(receptor, settings)
        deposited = [accession for accession in accessions if accession in accepted]
        if deposited:
            # A correct small segment cannot establish identity of the entire
            # declared construct. Count independent canonical matches against
            # its full length, retaining the configured mutant tolerance.
            full_pairs = alignment_pairs(sequences[receptor], sequence)
            full_agreement = (sum(sequence[q-1] == sequences[receptor][c-1] for q,c in full_pairs)
                              / len(sequence)) if sequence else 0.0
            if full_agreement < settings["annotated_min_identity"]:
                return "unresolved", f"deposited receptor accession contradicts sequence correspondence (full declared={full_agreement:.3f})", accessions
            pairs, source = map_entity(entity_id, sequence, data, sequences[receptor], accepted)
            agreement = sum(sequence[q-1] == sequences[receptor][c-1] for q,c in pairs) / len(pairs) if pairs else 0.0
            if agreement < settings["annotated_min_identity"]:
                return "unresolved", f"deposited receptor accession contradicts sequence correspondence (mapped={agreement:.3f}, full declared={full_agreement:.3f})", accessions
            provenance = "deposited UniProt accession"
            if accepted[0] not in deposited:
                alias = deposited[0]
                evidence = settings["receptor_accession_aliases"][receptor][alias]
                provenance = (f"deposited UniProt secondary accession {alias} -> {accepted[0]}; "
                              f"source={evidence['source_url']}; cache_sha256={evidence['cache_sha256']}")
            return receptor, f"{provenance}; declared correspondence agreement={agreement:.3f}; full declared agreement={full_agreement:.3f}", accessions
    if "ribonucleotide" in polymer_type:
        return ("dna" if "deoxy" in polymer_type else "rna"), "deposited polymer type", accessions
    if "polypeptide" not in polymer_type:
        return "unresolved", "non-protein polymer", accessions
    # A declared different accession takes precedence over a coincidental local match.
    if accessions:
        return ("peptide" if len(sequence) <= settings["peptide_max_length"] else "partner"), "other deposited accession", accessions
    matches = []
    for receptor in ("esr1", "esr2"):
        pairs = alignment_pairs(sequences[receptor], sequence)
        hits = sum(sequence[q-1] == sequences[receptor][c-1] for q, c in pairs)
        identity = hits / len(pairs) if pairs else 0
        if len(pairs) >= settings["fallback_min_positions"] and identity >= settings["fallback_min_identity"] and len(pairs)/len(sequence) >= .9:
            matches.append(receptor)
    if len(matches) == 1:
        return matches[0], "declared-sequence fallback; no accession; >=90% query aligned", accessions
    return ("peptide" if len(sequence) <= settings["peptide_max_length"] else "unresolved"), "no conclusive receptor evidence", accessions


def map_entity(entity_id, sequence, data, canonical, accession):
    """Map deposited segments for a primary accession or [primary, *aliases]."""
    accessions = [accession] if isinstance(accession, str) else list(accession)
    primary = accessions[0]
    reference_rows = [r for r in category_rows(data, "_struct_ref")
                      if r.get("entity_id") == entity_id and r.get("pdbx_db_accession") in accessions]
    aliases = sorted({r["pdbx_db_accession"] for r in reference_rows if r["pdbx_db_accession"] != primary})
    alias_source = f"; deposited secondary accession {','.join(aliases)} -> {primary}" if aliases else ""
    isoforms = sorted({r["pdbx_db_isoform"] for r in reference_rows
                       if r.get("pdbx_db_isoform") not in (None, "", ".", "?", *accessions)})
    if isoforms:
        # Deposited database bounds refer to the named isoform, not the frozen
        # canonical sequence. Align the full declared construct directly; do
        # not apply an assumed isoform offset or truncate it to those bounds.
        return alignment_pairs(canonical, sequence), (
            f"deposited isoform {','.join(isoforms)}; full declared-sequence alignment to canonical {primary}{alias_source}")
    refs = {r["id"] for r in reference_rows}
    segments = category_rows(data, "_struct_ref_seq")
    segment_pairs = []
    seen_segments = set()
    for row in segments:
        if row.get("ref_id") not in refs:
            continue
        bounds = tuple(int(row[name]) for name in ("seq_align_beg", "seq_align_end", "db_align_beg", "db_align_end"))
        if bounds in seen_segments:
            continue  # Same entity reference declared for several identical chain copies.
        seen_segments.add(bounds)
        q0, q1, c0, c1 = bounds
        segment = alignment_pairs(canonical[c0-1:c1], sequence[q0-1:q1], "global")
        query_length = q1-q0+1
        if query_length <= 0 or len(segment)/query_length < .9:
            raise ValueError(f"Deposited UniProt segment loses declared query positions: entity {entity_id}, "
                             f"query {q0}-{q1}, database {c0}-{c1}, mapped {len(segment)}/{query_length}")
        for q, c in segment:
            segment_pairs.append((q+q0-1, c+c0-1))
    if segment_pairs:
        validate_pairs(segment_pairs)
        return segment_pairs, "deposited UniProt segments plus global sequence alignment" + alias_source
    return alignment_pairs(canonical, sequence), "declared-sequence local alignment fallback" + alias_source


def residue_letter(name, data):
    parent = name
    for row in category_rows(data, "_chem_comp"):
        if row.get("id") == name and row.get("mon_nstd_parent_comp_id") not in (None, ".", "?"):
            parent = row["mon_nstd_parent_comp_id"]
            break
    return seq1(parent, custom_map={"MSE": "M", "CME": "C", "SEP": "S", "TPO": "T", "PTR": "Y"}, undef_code="X")


def missing_stretches(expected, observed):
    missing = sorted(set(expected) - set(observed))
    groups = []
    for position in missing:
        if groups and position == groups[-1][-1] + 1:
            groups[-1].append(position)
        else:
            groups.append([position])
    return ";".join(str(g[0]) if len(g)==1 else f"{g[0]}-{g[-1]}" for g in groups)
