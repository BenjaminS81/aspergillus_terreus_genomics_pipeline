#!/usr/bin/env python3
"""Read-only deliverable catalogue; writes only to a separate report directory."""
from __future__ import annotations

import argparse
import collections
import datetime as dt
import hashlib
import html
import json
import os
import re
from pathlib import Path

RESULT_TREES = ("results", "results1", "results2")
SKIP_NAMES = {
    "train_misc", "predict_misc", "annotate_misc", "genemark-es", "augustus",
    "species", "hmm_db", "busco_downloads", "hmmer_output", "miniprot_output",
    "metaeuk_output", "busco_sequences", "single_copy_busco_sequences",
    "multi_copy_busco_sequences", "fragmented_busco_sequences", "corrected",
    "tmp", "temp", "pipeline_state", "spades_corrected", "orthologues", "translated_proteins",
    "orthogroup_sequences", "workingdirectory", "gene_trees", "resolved_gene_trees",
    "multiplesequencealignments", "phylogenetic_hierarchical_orthogroups",
    "orthologues_stats", "html_content", "js", "css", "images", "img", "fonts",
}
EXTENSIONS = {".fa", ".fasta", ".fna", ".faa", ".gff", ".gff3", ".gtf", ".gbk",
              ".gbff", ".tsv", ".csv", ".txt", ".json", ".yaml", ".yml", ".md",
              ".pdf", ".png", ".svg", ".html", ".treefile", ".nwk", ".newick",
              ".iqtree", ".log", ".phy", ".aln", ".nex", ".xlsx", ".ods", ".ipynb",
              ".sh", ".bash", ".py", ".r", ".vcf", ".bed", ".out", ".tab", ".stats", ".signalp5", ""}
PRIORITY = {".treefile": 0, ".iqtree": 1, ".nwk": 2, ".pdf": 3, ".svg": 4,
            ".png": 5, ".html": 6, ".tsv": 7, ".csv": 8, ".fa": 9}


def stamp(t):
    return dt.datetime.fromtimestamp(t, dt.timezone.utc).astimezone().isoformat(timespec="seconds")


def category(stage):
    if "Houbraken" in stage:
        return "Taxonomic marker phylogenies"
    if stage.startswith("1_") or stage.startswith("3_0"):
        return "Read mapping, variants and structural variation"
    if stage.startswith("3_1"):
        return "Lovastatin cluster analyses"
    if "orthofinder" in stage or "iqtree" in stage or "mst" in stage or "clustering" in stage:
        return "Comparative genomics and phylogenies"
    if "braker" in stage or "funannotate" in stage or "busco" in stage:
        return "Annotation and completeness"
    if "fastp" in stage or "spades" in stage:
        return "Reads and assemblies"
    return "Functional analyses"


def candidates(root, branch, max_depth, exclusions, errors, limit=6000):
    records = []
    stack = [(branch, 0)]
    truncated = False
    dirs_seen = 0
    while stack:
        folder, depth = stack.pop()
        dirs_seen += 1
        if dirs_seen > 12000:
            exclusions.append({"path": str(branch), "reason": "Directory visit cap reached"})
            truncated = True
            break
        try:
            with os.scandir(folder) as entries:
                for n, entry in enumerate(entries):
                    if n >= 12000:
                        exclusions.append({"path": str(folder), "reason": "Per-directory entry cap reached"})
                        truncated = True
                        break
                    try:
                        if entry.is_symlink() or (getattr(entry.stat(follow_symlinks=False), "st_file_attributes", 0) & 0x400):
                            exclusions.append({"path": entry.path, "reason": "Reparse point not followed"})
                            continue
                        if entry.is_dir(follow_symlinks=False):
                            name = entry.name.lower()
                            if name in SKIP_NAMES or re.fullmatch(r"k\d+", name) or name.startswith("run_eurotiales") or ("busco" in branch.name.lower() and name.startswith("run_")):
                                exclusions.append({"path": entry.path, "reason": "Intermediate/database/large per-gene directory excluded"})
                            elif depth < max_depth:
                                stack.append((Path(entry.path), depth + 1))
                            else:
                                exclusions.append({"path": entry.path, "reason": f"Depth limit {max_depth} reached"})
                            continue
                        path = Path(entry.path)
                        suffix = path.suffix.lower()
                        if re.match(r"OG\d+.*\.(fa|fasta|faa|aln)$", path.name, re.I) or re.search(r"\.region\d+\.gbk$", path.name, re.I):
                            exclusions.append({"path": str(path), "reason": "Per-gene alignment or per-region cluster file omitted from release-candidate ranking"})
                            continue
                        if suffix not in EXTENSIONS and not (suffix == ".gz" and path.name.lower().endswith((".fna.gz", ".fasta.gz", ".vcf.gz"))):
                            continue
                        stat = entry.stat(follow_symlinks=False)
                        rel = path.relative_to(root)
                        parts = rel.parts
                        tree = parts[0]
                        stage = parts[1] if len(parts) > 2 else tree
                        records.append({"path": str(path), "relative_path": str(rel), "tree": tree,
                                        "stage": stage, "category": category(stage) if tree in RESULT_TREES else tree,
                                        "within_stage": str(Path(*parts[2:])) if tree in RESULT_TREES else str(Path(*parts[1:])),
                                        "bytes": stat.st_size, "modified": stamp(stat.st_mtime),
                                        "created_or_ctime": stamp(stat.st_ctime), "mtime_epoch": stat.st_mtime,
                                        "suffix": suffix, "archive_named": any(p.lower() in ("old", "archive", "archived") for p in parts),
                                        "zero_bytes": stat.st_size == 0})
                        if len(records) >= limit:
                            exclusions.append({"path": str(branch), "reason": f"Deliverable cap {limit} reached"})
                            truncated = True
                            stack.clear()
                            break
                    except OSError as exc:
                        errors.append({"path": entry.path, "error": str(exc)})
        except OSError as exc:
            errors.append({"path": str(folder), "error": str(exc)})
    return records, truncated


def sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("root", type=Path)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    output = args.output.resolve()
    if output.is_relative_to(root):
        raise SystemExit("Report directory must be outside the inventoried project")
    output.mkdir(parents=True, exist_ok=True)
    errors, exclusions, files, stages = [], [], [], []
    for tree in RESULT_TREES:
        base = root / tree
        for stage in sorted(base.iterdir()):
            if not stage.is_dir():
                continue
            recs, limited = candidates(root, stage, 4, exclusions, errors)
            files.extend(recs)
            children = sorted(p.name for p in stage.iterdir() if p.is_dir())
            stages.append({"tree": tree, "stage": stage.name, "category": category(stage.name),
                           "path": str(stage), "folder_modified": stamp(stage.stat().st_mtime),
                           "folder_created_or_ctime": stamp(stage.stat().st_ctime), "children": children,
                           "catalogued_files": len(recs), "limited": limited})
            print(f"{tree}/{stage.name}: {len(recs)} selected files", flush=True)
    for name, depth in (("data", 3), ("config", 3), ("docs", 2), ("notebooks", 0),
                        ("workflow", 3), ("terreus_qc_check", 4)):
        recs, _ = candidates(root, root / name, depth, exclusions, errors, limit=12000)
        files.extend(recs)
    # Input/archive material is not treated as an automatically superseding analysis release.
    exclusions.append({"path": str(root / "inbox"), "reason": "Import/archive tree excluded from deliverable version ranking; inventoried separately in prior structural review"})

    file_by_path = {f["path"]: f for f in files}
    groups = collections.defaultdict(list)
    for f in files:
        if f["tree"] in RESULT_TREES:
            groups[(f["stage"], f["within_stage"])].append(f)
    comparisons = []
    hash_budget = 256 * 1024 * 1024
    for (stage, within), group in sorted(groups.items()):
        if len({f["tree"] for f in group}) < 2:
            continue
        ordered = sorted(group, key=lambda f: (f["mtime_epoch"], {"results": 2, "results1": 1, "results2": 0}[f["tree"]]), reverse=True)
        hashes = []
        total = sum(f["bytes"] for f in ordered)
        if total <= hash_budget and max(f["bytes"] for f in ordered) <= 4 * 1024 * 1024:
            for f in ordered:
                try:
                    f["sha256"] = sha256(f["path"])
                    hashes.append(f["sha256"])
                except OSError as exc:
                    errors.append({"path": f["path"], "error": str(exc)})
            hash_budget -= total
        identical = len(hashes) == len(ordered) and len(set(hashes)) == 1
        tie = len({f["mtime_epoch"] for f in ordered}) == 1
        state = "identical content (SHA-256)" if identical else (("same timestamp; different content (SHA-256)" if len(hashes) == len(ordered) else "same timestamp; unresolved content") if tie else "latest timestamp candidate; workflow lineage unverified")
        comparisons.append({"stage": stage, "relative_identity": within, "status": state,
                            "candidate_path": ordered[0]["path"], "candidate_modified": ordered[0]["modified"],
                            "members": [{k: f[k] for k in ("path", "tree", "modified", "bytes")} | ({"sha256": f["sha256"]} if "sha256" in f else {}) for f in ordered]})

    for stage in stages:
        selected = [f for f in files if f["tree"] == stage["tree"] and f["stage"] == stage["stage"]]
        stage["newest_catalogued_modified"] = max((f["modified"] for f in selected), default=None)
        ordered = sorted(selected, key=lambda f: (-f["mtime_epoch"], PRIORITY.get(f["suffix"], 30)))
        stage["newest_files"] = [{k: f[k] for k in ("path", "modified", "bytes", "archive_named")} for f in ordered[:6]]

    for f in files:
        f["observed_version"] = "observed-" + f["modified"][:19].replace(":", "").replace("-", "")
        if "sha256" in f:
            f["content_id"] = "sha256:" + f["sha256"]

    busco_runs = []
    for f in files:
        if not (Path(f["path"]).name.startswith("short_summary.specific.") and f["suffix"] == ".json"):
            continue
        try:
            obj = json.loads(Path(f["path"]).read_text(encoding="utf-8"))
            params = obj.get("parameters", {})
            lineage = obj.get("lineage_dataset", {}).get("name", params.get("lineage_dataset", "unknown"))
            busco_runs.append({"path": f["path"], "tree": f["tree"], "stage": f["stage"],
                               "modified": f["modified"], "mode": params.get("mode", "unknown"),
                               "lineage": lineage, "auto_lineage": params.get("auto-lineage", "unknown"),
                               "recorded_input": params.get("in", "unknown")})
        except (OSError, ValueError, AttributeError) as exc:
            errors.append({"path": f["path"], "error": "BUSCO metadata could not be parsed: " + str(exc)})

    # Associate IQ-TREE outputs only by an identical filename prefix in the same directory.
    bundles = []
    for f in files:
        if f["suffix"] != ".treefile":
            continue
        stem = f["path"][:-len(".treefile")]
        members = [stem + ext for ext in ("", ".treefile", ".iqtree", ".log", ".mldist", ".model.gz", ".ckp.gz") if Path(stem + ext).is_file()]
        log_status = "not checked"
        logfile = Path(stem + ".log")
        if logfile.is_file():
            try:
                with logfile.open("rb") as handle:
                    handle.seek(max(0, logfile.stat().st_size - 16384))
                    tail = handle.read().decode("utf-8", errors="replace")
                log_status = "completion wording present" if "Total wall-clock time used:" in tail else "completion wording not found in log tail"
            except OSError:
                log_status = "log unreadable"
        bundles.append({"treefile": f["path"], "tree": f["tree"], "stage": f["stage"],
                        "modified": f["modified"], "associated_paths": members,
                        "association_basis": "Same directory and exact IQ-TREE output prefix", "log_status": log_status})

    payload = {"root": str(root), "generated": dt.datetime.now().astimezone().isoformat(),
               "scope": "Selected deliverables, max depth 4 below result stage; not an exhaustive file inventory",
               "rules": ["Dates are filesystem metadata, not a validated analysis history.",
                         "Creation dates can reflect copying; prefer file modification dates and run evidence.",
                         "Same relative path across results trees defines a candidate comparison, not proof of identical parameters or inputs.",
                         "Different samples, mod/unmod, good/low_quality and named run folders remain separate variants.",
                         "Files in old folders are historical-labelled; dates are retained.",
                         "No semantic version numbers or authoritative CURRENT labels assigned."],
               "stage_inventory": stages, "files": files, "cross_tree_comparisons": comparisons,
               "iqtree_bundles": bundles, "busco_run_metadata": busco_runs, "exclusions": exclusions, "errors": errors}
    (output / "chronological_catalogue.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    e = lambda v: html.escape(str(v))
    link = lambda p: '<a href="' + e(Path(p).as_uri()) + '">' + e(p) + '</a>'
    body = ["<!doctype html><html><head><meta charset='utf-8'><title>Terreus chronological catalogue</title>",
            "<style>body{font:15px system-ui;margin:32px;color:#172638;max-width:1600px}table{border-collapse:collapse;width:100%;margin:20px 0}th,td{border:1px solid #ccd5df;padding:8px;vertical-align:top;text-align:left}th{background:#eaf0f7}small{color:#506078}a{overflow-wrap:anywhere}details{margin:14px 0}h2{margin-top:36px}.note{background:#fff3cd;padding:16px}</style></head><body>",
            "<h1>Terreus chronological deliverable catalogue</h1>", f"<p>Source: {e(root)}<br>Generated: {e(payload['generated'])}</p>",
            "<div class='note'><strong>Review catalogue, not an authoritative release manifest.</strong> Source files were only read. Selected final outputs were catalogued; intermediate files, imported archives and deeper paths are excluded. Latest timestamp does not prove scientific validity or completion.</div>",
            "<h2>Interpretation rules</h2><ul>" + "".join(f"<li>{e(rule)}</li>" for rule in payload["rules"]) + "</ul>",
            f"<p>{len(files):,} selected files; {len(comparisons):,} matching-path comparisons; {len(bundles):,} IQ-TREE bundles; {len(errors)} read errors. Full paths and exclusions are in chronological_catalogue.json.</p>",
            "<h2>Analysis stages and candidate recent outputs</h2><table><tr><th>Category / stage</th><th>Tree</th><th>Newest catalogued file date</th><th>Recent paths</th></tr>"]
    for s in sorted(stages, key=lambda s: (s["category"], s["stage"], s["tree"])):
        paths = "<br>".join(link(f["path"]) + "<br><small>" + e(f["modified"]) + (" · old-folder label" if f["archive_named"] else "") + "</small>" for f in s["newest_files"][:3])
        body.append(f"<tr><td>{e(s['category'])}<br><strong>{e(s['stage'])}</strong></td><td>{e(s['tree'])}</td><td>{e(s['newest_catalogued_modified'])}<br><small>{s['catalogued_files']} catalogued files</small></td><td>{paths}</td></tr>")
    body.append("</table><h2>Matching relative paths across result trees</h2><p>Candidate latest is chosen by file modification time. Equal timestamps are not resolved by folder numbering. Checksums confirm content equality only for selected small files.</p>")
    for c in comparisons:
        body.append(f"<details><summary>{e(c['stage'])} / {e(c['relative_identity'])} — {e(c['status'])}</summary><ul>" + "".join(f"<li>{link(m['path'])}<br>{e(m['modified'])} · {m['bytes']:,} bytes</li>" for m in c['members']) + "</ul></details>")
    body.append("<h2>Associated phylogeny files</h2>")
    for b in sorted(bundles, key=lambda b: b["modified"], reverse=True):
        body.append(f"<details><summary>{e(b['modified'])} · {e(Path(b['treefile']).relative_to(root))}</summary><p>{e(b['log_status'])}; association by exact output prefix.</p><ul>" + "".join(f"<li>{link(p)}</li>" for p in b['associated_paths']) + "</ul></details>")
    body.append("<h2>BUSCO method categories</h2><p>Genome and transcriptome modes and different lineage datasets are distinct analyses, not automatically superseding versions.</p><table><tr><th>Tree / stage</th><th>Mode</th><th>Lineage</th><th>Auto-lineage</th><th>Specific summary files catalogued</th></tr>")
    busco_counts = collections.Counter((r["tree"], r["stage"], r["mode"], r["lineage"], str(r["auto_lineage"])) for r in busco_runs)
    for key, count in sorted(busco_counts.items()):
        body.append("<tr>" + "".join("<td>" + e(v) + "</td>" for v in (key[0] + "/" + key[1], *key[2:], count)) + "</tr>")
    body.append("</table></body></html>")
    (output / "chronological_catalogue.html").write_text("\n".join(body), encoding="utf-8")

    # A compact entry point with specific final-output links and explicitly bounded claims.
    chosen = [
        ("Genome phylogeny", "results/0_4_mafft_iqtree/Sep28/concatenated.fa.treefile", "Latest observed main concatenated tree; IQ-TREE completion wording checked. Sep28 is a run-folder label, not the final tree date."),
        ("Genome phylogeny: modified Newick variant", "results/0_4_mafft_iqtree/Sep28/concatenated.fa.mod.nwk", "Later modified Newick file; exact edits not verified. Retained as a separate derivative, not evidence of a later inference."),
        ("Taxonomic marker phylogeny", "results/Houbraken_genes/merged/concatenated.fa.treefile", "Latest observed concatenated marker tree in this family; merged/2 remains a separate run."),
        ("Taxonomic marker phylogeny: modified Newick variant", "results/Houbraken_genes/merged/concatenated.fa.mod.nwk", "Modified Newick derivative; exact edits not verified. Retain association with the original tree and alignment."),
        ("RPB2 individual marker", "results/Houbraken_genes/mafft_iqtree/RPB2_aligned.fa.treefile", "RPB2 tree; sample-specific SRR33623979 run remains separate."),
        ("Expanded genome QC", "terreus_qc_check/complete_data/complete_genome_qc_summary.tsv", "Expanded 66-genome summary snapshot. Does not establish coverage of assemblies added later."),
        ("Expanded transcript QC", "terreus_qc_check/complete_data/complete_transcript_qc_summary.tsv", "59-transcript summary snapshot."),
        ("QC data availability", "terreus_qc_check/complete_data/data_availability.tsv", "Associates the expanded QC samples with applicable data types."),
        ("QC graphics", "terreus_qc_check/figures/figure4_integrated_qc.png", "Existing 28-sample graphics. No expanded complete_data/figures folder observed."),
        ("Minimum spanning tree", "results/0_5_mst/mst.png", "Latest observed non-old-folder MST graphic."),
        ("BiG-SCAPE graphic", "results/0_6_1_bigscape/bigscape.png", "Current-tree graphic; supporting tables remain linked by folder, not inferred to be interchangeable."),
        ("BiG-SCAPE presence/absence", "results/0_6_1_bigscape/absence_presence(3).tsv", "Numbered filename is a candidate revision; relationship to earlier tables unverified."),
        ("Recent annotation example: MCA10", "results/0_2_repeat_mask_funannotate/MCA10/funannotate/predict_results/Aspergillus_hortae.summary.json", "A recent sample-specific annotation output, not a replacement for other samples' annotations."),
    ]
    picks = []
    for title, rel, note in chosen:
        path = root / Path(rel)
        if not path.is_file():
            continue
        stat = path.stat()
        picks.append({"analysis": title, "path": str(path), "modified": stamp(stat.st_mtime),
                      "observed_version": "observed-" + stamp(stat.st_mtime)[:19].replace(":", "").replace("-", ""),
                      "interpretation": note, "associated_paths": next((b["associated_paths"] for b in bundles if b["treefile"] == str(path)), [])})
    recent_busco = sorted([r for r in busco_runs if r["tree"] == "results" and r["stage"] == "0_2_2_busco"], key=lambda r: r["modified"], reverse=True)
    if recent_busco:
        r = recent_busco[0]
        picks.append({"analysis": "Additional BUSCO series: recent sample example", "path": r["path"], "modified": r["modified"],
                      "observed_version": "observed-" + r["modified"][:19].replace(":", "").replace("-", ""),
                      "interpretation": f"Recorded mode {r['mode']}; lineage {r['lineage']}; auto-lineage {r['auto_lineage']}. Different from the September 10 genome-QC series; not automatically interchangeable.", "associated_paths": []})
    (output / "latest_output_index.json").write_text(json.dumps(picks, indent=2), encoding="utf-8")
    index = ["<!doctype html><html><head><meta charset='utf-8'><title>Terreus latest output index</title><style>body{font:15px system-ui;margin:32px;max-width:1450px;color:#172638}table{border-collapse:collapse;width:100%}td,th{border:1px solid #cbd5df;padding:10px;text-align:left;vertical-align:top}th{background:#eaf0f7}a{overflow-wrap:anywhere}.note{background:#fff3cd;padding:16px}</style></head><body>",
             "<h1>Terreus: chronological categories and recent output index</h1>",
             "<p class='note'>Source project unchanged. These are observed-date version labels for review, not formal releases or validated workflow provenance. Creation dates around 30 October suggest copying; do not use them as analysis dates.</p>",
             "<h2>Chronological organisation</h2><ul><li>April–August 2026: baseline and earlier read, annotation, functional and variant analyses.</li><li>9–10 September: original and expanded QC snapshots; the figure set still covers 28 samples.</li><li>Mid/late September: additional assembly/QC runs, OrthoFinder run folders, functional outputs and clustering.</li><li>October: genome phylogeny, marker phylogenies, selected variant outputs and additional annotations.</li></ul>",
             "<p>results, results1 and results2 are overlapping analysis collections, not an ordered v1/v2/v3 sequence. Different samples, run folders, quality groups and mod/unmod outputs remain separate variants.</p>",
             "<h2>Associated candidate recent outputs</h2><table><tr><th>Analysis</th><th>File modification / observed version</th><th>File and associated outputs</th><th>Interpretation</th></tr>"]
    for p in picks:
        paths = link(p["path"]) + "".join("<br>" + link(q) for q in p["associated_paths"] if q != p["path"])
        index.append(f"<tr><td>{e(p['analysis'])}</td><td>{e(p['modified'])}<br>{e(p['observed_version'])}</td><td>{paths}</td><td>{e(p['interpretation'])}</td></tr>")
    index.append("</table><h2>Confidence and scope</h2><p>Confirmed identical files have matching SHA-256 hashes. Latest-by-date candidates have unverified input/parameter lineage. Timestamp ties and renamed paths remain unresolved. The full catalogue uses a depth limit and excludes intermediate directories; it is not an exhaustive current-release manifest.</p><p><a href='chronological_catalogue.html'>Full chronological catalogue and cross-tree comparisons</a> · <a href='chronological_catalogue.json'>Machine-readable catalogue, exclusions and file associations</a></p></body></html>")
    (output / "latest_output_index.html").write_text("\n".join(index), encoding="utf-8")
    summary = {"selected_files": len(files), "comparison_groups": len(comparisons),
               "comparison_status_counts": dict(collections.Counter(c["status"] for c in comparisons)),
               "iqtree_bundles": len(bundles), "read_errors": len(errors), "excluded_directories": len(exclusions),
               "report_html": str(output / "chronological_catalogue.html")}
    (output / "catalogue_summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")
    print(json.dumps(summary), flush=True)


if __name__ == "__main__":
    main()
