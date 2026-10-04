"""Block 01: explicitly refresh RCSB accession/X-ray ID queries.
Input: frozen project settings, separate output directory. Output: new ID lists.
Never replace the 487-entry ESR1 master list. Network only when explicitly run."""
import json
from pathlib import Path


def query_structures(accession, method="X-RAY DIFFRACTION"):
    """Return sorted RCSB entry IDs for the existing accession/method query."""
    import requests

    query = {"query": {"type": "group", "logical_operator": "and", "nodes": [
        {"type": "terminal", "service": "text", "parameters": {
            "attribute": "rcsb_polymer_entity_container_identifiers.reference_sequence_identifiers.database_accession",
            "operator": "exact_match", "value": accession}},
        {"type": "terminal", "service": "text", "parameters": {
            "attribute": "exptl.method", "operator": "exact_match", "value": method}}]},
             "return_type": "entry", "request_options": {"return_all_hits": True}}
    response = requests.post("https://search.rcsb.org/rcsbsearch/v2/query", json=query, timeout=60)
    response.raise_for_status()
    ids = sorted(hit["identifier"] for hit in response.json()["result_set"])
    return ids


def refresh_cohort(root, output):
    """Save a newly queried cohort separately. Does not change the frozen lists."""
    settings = json.loads((Path(root) / "config" / "analysis_settings.json").read_text())
    output = Path(output)
    output.mkdir(parents=True, exist_ok=True)
    for receptor in ("esr1", "esr2"):
        ids = query_structures(settings[f"{receptor}_uniprot"])
        (output / f"pdb_ids_{receptor}.txt").write_text("\n".join(ids) + "\n")
        print(f"Refresh: {receptor} {len(ids)} IDs saved to {output}; raw download is separate.")



def main():
    from ._task_io import parser, prepare, json_write
    args_parser = parser('Block 01: frozen offline cohort by default; live query is explicit and separate')
    args_parser.add_argument('--live', action='store_true')
    args = args_parser.parse_args()
    root, source, out, settings = prepare(args)
    if args.live:
        ids = query_structures(settings['esr1_uniprot'])
        name, provenance = 'pdb_ids_esr1_live.txt', {'scope': 'new live cohort; not the frozen report baseline'}
    else:
        cached = json.loads((source / 'query_snapshot.json').read_text())
        ids = cached['pdb_ids']
        name, provenance = 'pdb_ids_esr1.txt', cached
        if ids != (root / 'data/pdb_ids_esr1.txt').read_text().splitlines():
            raise ValueError('Frozen query IDs disagree with the master cohort')
    (out / name).write_text('\n'.join(ids) + '\n')
    json_write(out / 'query_provenance.json', provenance)
    print(f'Block 01: {len(ids)} IDs, live={args.live}')


if __name__ == "__main__":
    main()
