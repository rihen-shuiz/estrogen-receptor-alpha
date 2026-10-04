"""Block 07: resolve receptor-specific action evidence without promoting provisional labels.

Read this after the chemical identity selection. These rules describe evidence,
not fingerprints, similarity or inferred binding. No files are written here.
"""
import hashlib
import json
import re
from rdkit import Chem

ACTION_COLUMNS = [
    "action_type", "action_status", "action_source", "action_evidence",
    "action_receptor_context", "action_doi", "action_reviewer",
    "action_record_hash", "action_notes", "supported_action_type",
    "evidence_supported", "supported_action_source", "supported_action_evidence",
    "supported_action_receptor_context", "supported_action_record_hash",
]


def record_hash(row: dict) -> str:
    return hashlib.sha256(json.dumps(row, sort_keys=True).encode()).hexdigest()


def _action(value: str) -> str:
    value = value.strip()
    if value.upper() in {"SERM", "SERD"}:
        return value.upper()
    return value.lower() or "unknown"


def receptor_contexts(context: str) -> set[str]:
    """Read explicit receptor names/identifiers; never infer from agonism."""
    receptors = set()
    for receptor, pattern in [
        ("ESR1", r"\bESR1\b|\bER\s*[-_ ]?alpha\b|ERα|\bCHEMBL206\b|\bP03372\b"),
        ("ESR2", r"\bESR2\b|\bER\s*[-_ ]?beta\b|ERβ|\bCHEMBL242\b|\bQ92731\b"),
    ]:
        if re.search(pattern, context, re.IGNORECASE):
            receptors.add(receptor)
    return receptors


def _human_review(row: dict, override: bool = False) -> bool:
    action = row.get("action_type" if override else "final_action_type", "")
    reviewer = row.get("reviewer", "").strip()
    doi = row.get("doi_or_url" if override else "doi", "")
    status = row.get("status", "").lower()
    confirmed = not override or status in {"confirmed", "human_confirmed", "human confirmed"}
    return bool(
        confirmed and _action(action) not in {"unknown", "conflicting"} and reviewer
        and status not in {"conflicting", "rejected"}
        and not reviewer.lower().startswith("auto")
        and re.search(r"10\.\d{4,9}/\S+", doi)
    )


def resolve_action(
    code: str, name: str, receptor: str, annotations: list[dict],
    overrides: list[dict], rules: list[dict], canonical_isomeric_smiles: str = "",
) -> dict:
    """Resolve labels without elevating provisional or other-receptor evidence.

    Human final labels with reviewer and DOI take precedence over versioned
    provisional overrides, then code/name seed rules. The historical review
    sheet was explicitly ESR1-scoped; a missing context there means ESR1 only.
    Supported action is a separate field: a human label or an exact cached
    ChEMBL mechanism with explicit receptor target, reference and high confidence.
    Generic MODULATOR/DEGRADER stay generic; they are not converted to SERM/SERD.
    """
    result = dict.fromkeys(ACTION_COLUMNS, "")
    result.update(action_type="unknown", action_status="unknown", evidence_supported=False)
    human, provisional, seeded, curated = [], [], [], []
    notes = []
    for row in annotations:
        if row.get("ligand_code", "").upper() != code:
            continue
        if canonical_isomeric_smiles and row.get("smiles"):
            annotated_mol = Chem.MolFromSmiles(row["smiles"])
            if annotated_mol is None or Chem.MolToSmiles(annotated_mol, isomericSmiles=True) != canonical_isomeric_smiles:
                notes.append("annotation_chemical_identity_mismatch")
                continue
        context = row.get("receptor_context", "") or "ESR1 (annotation workflow CHEMBL206)"
        if _human_review(row):
            if receptor in receptor_contexts(context):
                human.append((_action(row["final_action_type"]), row, context, "human_annotation"))
            else:
                notes.append("human_annotation_receptor_mismatch")
        elif row.get("final_action_type"):
            notes.append("incomplete_or_conflicting_human_review")
        mechanism = row.get("chembl_mechanism_action_types", "")
        evidence = row.get("chembl_mechanism_evidence", "")
        confidence = row.get("confidence", "") or row.get("auto_confidence", "")
        # Read the target attached to each actual mechanism, not a target name
        # mentioned elsewhere in an annotation or a chemical-similarity hint.
        targets = re.findall(r"->\s*(CHEMBL\d+)\s*:", evidence, re.IGNORECASE)
        target = "CHEMBL206" if receptor == "ESR1" else "CHEMBL242"
        actions = {_action(value) for value in mechanism.split(";") if value.strip()}
        evidence_actions = {
            _action(value) for value in re.findall(r"->\s*CHEMBL\d+\s*:\s*([A-Za-z ]+?)(?:\s*\(|\s*\||$)", evidence)
        }
        if (len(actions) == 1 and actions != {"unknown"}
                and evidence_actions == actions and targets and set(targets) == {target}
                and confidence.lower() == "high" and row.get("chembl_mechanism_refs")
                and (row.get("candidate_source") == "chembl_mechanism" or row.get("auto_source") == "chembl_mechanism")
                and (not canonical_isomeric_smiles or row.get("smiles"))):
            curated.append((next(iter(actions)), row, f"{receptor} ({target})", "cached_chembl_mechanism"))

    for row in overrides:
        if row.get("ligand_code", "").upper() != code or _action(row.get("action_type", "")) == "unknown":
            continue
        context = row.get("receptor_context", "")
        if not context:
            notes.append("override_missing_receptor_context")
            continue
        receptors = receptor_contexts(context)
        if receptors and receptor not in receptors:
            notes.append("override_receptor_mismatch: " + context)
            continue
        candidate = (_action(row["action_type"]), row, context or "unspecified", "versioned_override")
        if _human_review(row, override=True) and receptor in receptors:
            human.append(candidate)
        elif row.get("status", "").lower() == "provisional":
            provisional.append(candidate)
    for row in rules:
        pattern = row.get("pattern", "").upper()
        match = row.get("match_type")
        if pattern and ((match == "code_exact" and pattern == code)
                        or (match == "name_contains" and pattern in name.upper())):
            context = row.get("receptor_context", "")
            if context and receptor_contexts(context) and receptor not in receptor_contexts(context):
                continue
            seeded.append((_action(row.get("action_type", "")), row, context or "unspecified (seed rule)", "seeded_rule"))

    selected = []
    for tier, status in [(human, "human_confirmed"), (provisional, "provisional"), (seeded, "provisional")]:
        if tier:
            selected = tier
            labels = {candidate[0] for candidate in tier}
            if len(labels) != 1:
                result.update(action_type="conflicting", action_status="conflicting")
                notes.append("conflicting_labels_at_highest_precedence")
            else:
                label, row, context, source = sorted(tier, key=lambda item: record_hash(item[1]))[0]
                result.update(
                    action_type=label, action_status=status, action_source=source,
                    action_evidence=row.get("evidence_source", "") or row.get("notes", "") or row.get("note", ""),
                    action_receptor_context=context, action_doi=row.get("doi", "") or row.get("doi_or_url", ""),
                    action_reviewer=row.get("reviewer", ""), action_record_hash=record_hash(row),
                )
            break

    supported = human if human else curated
    if supported and len({candidate[0] for candidate in supported}) == 1:
        label, row, context, source = sorted(supported, key=lambda item: record_hash(item[1]))[0]
        evidence = row.get("chembl_mechanism_evidence", "") if source == "cached_chembl_mechanism" else result["action_evidence"]
        references = row.get("chembl_mechanism_refs", "") if source == "cached_chembl_mechanism" else result["action_doi"]
        result.update(
            supported_action_type=label, evidence_supported=True,
            supported_action_source=source, supported_action_evidence=f"{evidence}; {references}",
            supported_action_receptor_context=context, supported_action_record_hash=record_hash(row),
        )
        if result["action_type"] not in {"unknown", label}:
            notes.append("supported_mechanism_differs_from_provisional_label")
    elif supported:
        notes.append("conflicting_supported_evidence")
    if selected and result["action_status"] == "conflicting":
        result["action_evidence"] = json.dumps(
            [item[1] for item in sorted(selected, key=lambda item: record_hash(item[1]))], sort_keys=True,
        )
    result["action_notes"] = "; ".join(sorted(set(notes)))
    return result



def main():
    import pandas as pd
    from ._task_io import parser, prepare, records
    args = parser('Block 07: exact receptor/identity action evidence; unknowns remain unknown').parse_args()
    root, source, out, settings = prepare(args)
    annotations = records(root / 'data/annotations/ligand_annotation_candidates.csv')
    overrides = records(root / 'config/ligand_action_overrides.csv')
    rules = records(root / 'config/ligand_class_rules.csv') + records(root / 'config/ligand_class_rules_ext.csv')
    rows = records(source / 'ligand_identities.csv')
    for row in rows:
        row.update(resolve_action(row['ligand_code'], row['ligand_name'], row['receptor_context'], annotations, overrides, rules, row['canonical_isomeric_smiles']))
    result = pd.DataFrame(rows)
    result.to_csv(out / 'ligand_actions.csv', index=False)
    print('Block 07: ' + str(result.supported_action_type.replace('', 'unknown').value_counts().to_dict()))


if __name__ == "__main__":
    main()
