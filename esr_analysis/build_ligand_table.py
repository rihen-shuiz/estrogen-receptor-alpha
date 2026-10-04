"""Block 06: exact local CCD chemical identities for contacting ligands.
Input: local CCD path. Method: canonical isomeric SMILES; no standardisation.
Output: identity/source/error fields and RDKit molecule; no network or I/O writes."""
from __future__ import annotations
import csv
import json
import hashlib
from pathlib import Path
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from rdkit import Chem
from .assign_action_classes import ACTION_COLUMNS, resolve_action, receptor_contexts


def _file_hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_ccd(path: Path) -> tuple[dict, object | None]:
    """Read local CCD descriptors, preferring canonical CACTVS SMILES.

    Canonical isomeric SMILES defines identity, retaining salts, charges and
    specified stereochemistry. Fingerprint chirality is a separate option.
    """
    row = dict.fromkeys([
        "ligand_name", "smiles", "canonical_isomeric_smiles", "chemical_identity_key",
        "smiles_source", "ccd_sha256", "chemistry_error",
    ], "")
    row.update(ccd_path=str(path), chemistry_status="missing_ccd")
    if not path.exists():
        row["chemistry_error"] = "local CCD file missing; no network retrieval attempted"
        return row, None
    row["ccd_sha256"] = _file_hash(path)
    try:
        ccd = MMCIF2Dict(str(path))
    except (ValueError, StopIteration) as error:
        row.update(chemistry_status="invalid_ccd", chemistry_error=f"{type(error).__name__}: {error}")
        return row, None
    row["ligand_name"] = ccd.get("_chem_comp.name", [""])[0]
    declared_code = ccd.get("_chem_comp.id", [""])[0]
    if declared_code and declared_code.upper() != path.stem.upper():
        row.update(chemistry_status="invalid_ccd", chemistry_error=f"CCD id {declared_code} does not match file code {path.stem}")
        return row, None
    descriptors = []
    for kind, program, version, descriptor in zip(
        ccd.get("_pdbx_chem_comp_descriptor.type", []),
        ccd.get("_pdbx_chem_comp_descriptor.program", []),
        ccd.get("_pdbx_chem_comp_descriptor.program_version", []),
        ccd.get("_pdbx_chem_comp_descriptor.descriptor", []),
    ):
        if kind.upper() in {"SMILES", "SMILES_CANONICAL"} and descriptor not in {"", ".", "?"}:
            priority = (kind.upper() != "SMILES_CANONICAL", program.upper() != "CACTVS", program, version, descriptor)
            descriptors.append((priority, kind, program, version, descriptor))
    if not descriptors:
        row.update(chemistry_status="missing_smiles", chemistry_error="CCD has no SMILES descriptor")
        return row, None
    _, kind, program, version, smiles = sorted(descriptors)[0]
    row.update(smiles=smiles, smiles_source=f"CCD {kind} {program} {version}")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or mol.GetNumAtoms() == 0:
        row.update(chemistry_status="invalid_smiles", chemistry_error="RDKit could not parse selected CCD SMILES")
        return row, None
    canonical = Chem.MolToSmiles(mol, canonical=True, isomericSmiles=True)
    row.update(
        chemistry_status="valid", canonical_isomeric_smiles=canonical,
        chemical_identity_key="smiles_sha256:" + hashlib.sha256(canonical.encode()).hexdigest(),
    )
    return row, mol


CHEMISTRY_COLUMNS = [
    "ligand_code", "receptor_context", "ligand_name", "cohort",
    "ligand_instance_keys", "chain_keys", "pdb_ids", "n_instances",
    "n_receptor_chains", "site_classification", "site_classifications",
    "n_pocket_compatible_contacts", "n_other_or_uncertain_contacts",
    "historical_receptors", "smiles", "canonical_isomeric_smiles",
    "chemical_identity_key", "smiles_source", "ccd_path", "ccd_sha256",
    "chemistry_status", "chemistry_error", "fingerprint_index",
    "ecfp_on_bits", "maccs_on_bits", "chem_cluster", "pca1", "pca2",
] + ACTION_COLUMNS
REVIEW_COLUMNS = [
    "ligand_code", "receptor_context", "ligand_name", "chemical_identity_key",
    "chemistry_status", "review_required",
] + ACTION_COLUMNS


def _read_csv(path: Path, required: bool = False) -> list[dict]:
    if not path.exists():
        if required:
            raise FileNotFoundError(f"Chemistry prerequisite missing: {path}; run inventory first")
        return []
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def build_ligand_table(root: Path, out: Path, settings: dict):
    """Return contacting identity/action rows, molecules and source provenance.

    Read saved ligand_instances/contacts and local CCD/annotation tables.
    Retain failed identities and uncertain contacts. No fingerprints or files written.
    """
    # Load observed ligand instances and receptor contacts.
    root, out = Path(root), Path(out)
    instances_path, contacts_path = out / "ligand_instances.csv", out / "contacts.csv"
    instances = _read_csv(instances_path, required=True)
    contacts = _read_csv(contacts_path, required=True)
    instance_index = {}
    for row in instances:
        key = row["ligand_instance_key"]
        if key in instance_index:
            raise ValueError(f"Duplicate ligand_instance_key in {instances_path}: {key}")
        instance_index[key] = row
    # Select contacting candidate molecules; retain uncertain sites.
    grouped = {}
    default_receptor = settings.get("chemistry_receptor", "ESR1")
    if default_receptor not in {"ESR1", "ESR2"}:
        raise ValueError("chemistry_receptor must be ESR1 or ESR2")
    for contact in contacts:
        if contact.get("contact_kind") != "receptor_ligand":
            continue
        key = contact["ligand_instance_key"]
        if key not in instance_index:
            raise ValueError(f"Contact references absent ligand_instance_key: {key}")
        instance = instance_index[key]
        if instance.get("category") != "candidate_small_molecule":
            continue
        context = contact.get("receptor_context", "") or contact.get("receptor", "") or instance.get("receptor", "") or default_receptor
        receptors = receptor_contexts(context)
        if len(receptors) != 1:
            raise ValueError(f"Contact needs one explicit receptor context: {key}: {context}")
        receptor = next(iter(receptors))
        code = instance["ligand_code"].upper()
        grouped.setdefault((code, receptor), []).append(contact)

    paths = [instances_path, contacts_path]
    historical = {}
    historical_names = {}
    for receptor in ("ESR1", "ESR2"):
        path = root / "data" / f"ligands_{receptor.lower()}.csv"
        paths.append(path)
        for row in _read_csv(path):
            code = row["ligand_code"].upper()
            historical.setdefault(code, set()).add(receptor)
            historical_names.setdefault(code, set()).add(row.get("ligand_name", ""))
    annotation_path = root / "data" / "annotations" / "ligand_annotation_candidates.csv"
    override_path = root / "config" / "ligand_action_overrides.csv"
    annotations, overrides = _read_csv(annotation_path), _read_csv(override_path)
    paths.extend([annotation_path, override_path])
    rules = []
    for name in ("ligand_class_rules.csv", "ligand_class_rules_ext.csv"):
        path = root / "config" / name
        paths.append(path)
        rules.extend(_read_csv(path))

    # Prepare exact chemical identities and attach separate evidence labels.
    molecules, chemistry_rows = {}, []
    for (code, receptor), assignments in sorted(grouped.items()):
        path = root / "data" / "ccd" / f"{code}.cif"
        paths.append(path)
        chemistry, mol = read_ccd(path)
        if not chemistry["ligand_name"]:
            chemistry["ligand_name"] = next(iter(sorted(historical_names.get(code, {""}))))
        if mol is not None:
            molecules[chemistry["chemical_identity_key"]] = mol
        sites = []
        for assignment in assignments:
            site = assignment.get("site_classification", "") or assignment.get("site", "")
            if site not in {"pocket_compatible", "other_or_uncertain"}:
                raise ValueError(f"Unknown ligand site classification for {code}: {site!r}")
            sites.append(site)
        instance_keys = sorted({row["ligand_instance_key"] for row in assignments})
        chain_keys = sorted({row["chain_key"] for row in assignments})
        unique_sites = sorted(set(sites))
        row = dict.fromkeys(CHEMISTRY_COLUMNS, "")
        row.update(
            ligand_code=code, receptor_context=receptor, cohort="receptor_contacting_candidate_small_molecule",
            ligand_instance_keys=json.dumps(instance_keys), chain_keys=json.dumps(chain_keys),
            pdb_ids=json.dumps(sorted({row["pdb_id"] for row in assignments})),
            n_instances=len(instance_keys), n_receptor_chains=len(chain_keys),
            site_classification=unique_sites[0] if len(unique_sites) == 1 else "mixed",
            site_classifications=json.dumps(unique_sites),
            n_pocket_compatible_contacts=sites.count("pocket_compatible"),
            n_other_or_uncertain_contacts=sites.count("other_or_uncertain"),
            historical_receptors=json.dumps(sorted(historical.get(code, set()))), **chemistry,
        )
        row.update(resolve_action(code, row["ligand_name"], receptor, annotations, overrides, rules,
                                  row["canonical_isomeric_smiles"]))
        chemistry_rows.append(row)

    return chemistry_rows, molecules, paths, historical


def main():
    import pandas as pd
    from ._task_io import parser, prepare, json_write
    args = parser('Block 06: local exact CCD identities for saved contacts; no fingerprints').parse_args()
    root, source, out, settings = prepare(args)
    rows, molecules, paths, historical = build_ligand_table(root, source, settings)
    result = pd.DataFrame(rows).drop(columns=['fingerprint_index', 'ecfp_on_bits', 'maccs_on_bits', 'chem_cluster', 'pca1', 'pca2'], errors='ignore')
    result.to_csv(out / 'ligand_identities.csv', index=False)
    json_write(out / 'identity_sources.json', {str(path.relative_to(root)) if path.is_relative_to(root) else str(path): _file_hash(path) for path in paths if path.is_file()})
    print(f'Block 06: {len(result)} contacting CCD/context rows; {len(molecules)} exact identities')


if __name__ == "__main__":
    main()
