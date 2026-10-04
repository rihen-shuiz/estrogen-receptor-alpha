"""Block 02: deposited mmCIF downloads, reusing the completed cache policy.
Input: explicit PDB IDs and output folder. Output: mmCIF files and manifest.
Existing valid files remain cached. Network runs only when this function is called.
"""
from pathlib import Path
import time
import requests
import pandas as pd
from Bio.PDB import MMCIFParser
from Bio.PDB.MMCIF2Dict import MMCIF2Dict

RCSB_CIF_URL = "https://files.rcsb.org/download/{pdb_id}.cif"


def validate_mmcif(path: Path, expected_pdb_id=None) -> tuple[bool, str]:
    """Check atoms and the deposited entry identity where declared."""
    try:
        data = MMCIF2Dict(str(path))
        declared = [str(value).upper() for field in ("_entry.id", "_struct.entry_id")
                    for value in data.get(field, []) if value not in {".", "?", ""}]
        if expected_pdb_id and declared and any(value != expected_pdb_id.upper() for value in declared):
            raise ValueError(f"deposited entry ID {declared} does not match requested {expected_pdb_id.upper()}")
        structure = MMCIFParser(QUIET=True).get_structure(path.stem, str(path))
        if not any(structure.get_atoms()):
            raise ValueError("no atoms")
    except Exception as exc:
        return False, f"invalid mmCIF: {exc}"
    return True, "valid"


def download_one(pdb_id: str, out_path: Path, retries: int = 3) -> tuple[bool, str]:
    if out_path.exists():
        valid, _ = validate_mmcif(out_path, pdb_id)
        if valid:
            return True, "cached"

    url = RCSB_CIF_URL.format(pdb_id=pdb_id)
    last_error = ""
    for attempt in range(1, retries + 1):
        temporary_path = out_path.with_suffix(f"{out_path.suffix}.tmp")
        try:
            r = requests.get(url, timeout=90)
            r.raise_for_status()
            out_path.parent.mkdir(parents=True, exist_ok=True)
            temporary_path.write_bytes(r.content)
            valid, status = validate_mmcif(temporary_path, pdb_id)
            if not valid:
                raise RuntimeError(status)
            temporary_path.replace(out_path)
            return True, "downloaded"
        except Exception as exc:
            last_error = str(exc)
            temporary_path.unlink(missing_ok=True)
            if attempt < retries:
                time.sleep(attempt * 1.5)
    return False, last_error


def download_structures(pdb_ids, output_folder):
    """Return and save the explicit-ID download manifest; never query or change IDs."""
    folder = Path(output_folder)
    folder.mkdir(parents=True, exist_ok=True)
    rows = []
    for pdb_id in pdb_ids:
        pdb_id = str(pdb_id).upper()
        path = folder / f"{pdb_id}.cif"
        ok, status = download_one(pdb_id, path)
        rows.append({"pdb_id": pdb_id, "download_ok": ok, "status": status,
                     "path": str(path), "source_url": RCSB_CIF_URL.format(pdb_id=pdb_id)})
    manifest = pd.DataFrame(rows, columns=["pdb_id", "download_ok", "status", "path", "source_url"])
    manifest.to_csv(folder / "download_manifest.csv", index=False)
    return manifest


def main():
    import shutil
    from ._task_io import parser, prepare
    from .run_io import sha256
    p = parser('Block 02: validate/copy bundled raw pilot; --live explicitly downloads listed IDs')
    p.add_argument('--ids', nargs='+', default=['3Q95', '7NFB', '7BAA', '6CHZ'])
    p.add_argument('--raw-folder', type=Path)
    p.add_argument('--live', action='store_true')
    args = p.parse_args()
    root, source, out, settings = prepare(args)
    if args.live:
        result = download_structures(args.ids, out)
        if not result.download_ok.all():
            raise RuntimeError('Some explicit downloads failed; see download_manifest.csv')
    else:
        folder = args.raw_folder or root / 'data/raw/esr1'
        rows = []
        for pdb in args.ids:
            path = folder / (pdb.upper() + '.cif')
            ok, reason = validate_mmcif(path, pdb)
            if not ok:
                raise ValueError(f'Offline prerequisite {path}: {reason}; restore it or choose --live explicitly')
            target = out / path.name
            shutil.copy2(path, target)
            rows.append({'pdb_id': pdb.upper(), 'download_ok': True, 'status': 'bundled_cached', 'path': str(target), 'sha256': sha256(target)})
        pd.DataFrame(rows).to_csv(out / 'download_manifest.csv', index=False)
    print(f'Block 02: {len(args.ids)} explicit structures, live={args.live}')


if __name__ == "__main__":
    main()
