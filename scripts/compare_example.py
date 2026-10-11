#!/usr/bin/env python3
"""Optional, fresh MiXCR 3.0.13 versus Python example comparison.

Run from an installed environment or a source checkout. The user supplies the
original JAR; its SHA-256 must match the verified artifact. Outputs are written
to a new directory. Exit 0 means all declared comparisons passed, 1 means a
result mismatch, and 2 means an input, execution, or comparison error.
"""
from __future__ import annotations

import argparse
from collections import Counter
import csv
from decimal import Decimal, InvalidOperation
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import subprocess
import sys
import time

ORIGINAL_SHA256 = "1526a6d36221e86b775e39831b5ebaab6d5877639212c672bef2049a3fd37149"
ROOT = Path(__file__).resolve().parents[1]
FIELDS = [f"all{kind}HitsWithScore" for kind in "VDJC"] + ["nSeqCDR3", "aaSeqCDR3", "qualCDR3"]


def sha256(path):
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def count_fastq(path):
    opener = gzip.open if path.suffix.lower() == ".gz" else open
    count = 0
    with opener(path, "rt", encoding="ascii", newline="") as stream:
        while True:
            name = stream.readline()
            if not name:
                break
            sequence, plus, quality = (stream.readline() for _ in range(3))
            if any("\r" in line for line in (name, sequence, plus, quality)):
                raise ValueError("MiXCR 3.0.13 requires FASTQ files with LF line endings")
            if not name.startswith("@") or not plus.startswith("+") or not quality:
                raise ValueError(f"Invalid four-line FASTQ record {count}")
            if len(sequence.rstrip("\n")) != len(quality.rstrip("\n")):
                raise ValueError(f"Sequence/quality length mismatch at read {count}")
            count += 1
    if not count:
        raise ValueError("FASTQ is empty")
    return count


def rows(path):
    with path.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream, delimiter="\t"))


def keyed(records, field, label, integer=False):
    result = {}
    for record in records:
        key = int(record[field]) if integer else record[field]
        if key in result:
            raise ValueError(f"Duplicate {field} in {label}: {key!r}; this verifier requires unique CDR3 clone keys")
        result[key] = record
    return result


def number(text):
    value = Decimal(text)
    if not value.is_finite():
        raise ValueError(f"Non-finite result value: {text}")
    return value


def normalize(field, value):
    if field.endswith("HitsWithScore"):
        result = []
        for item in value.split(",") if value else []:
            match = re.fullmatch(r"([^,()]+)\(([^()]+)\)", item)
            if match is None:
                raise ValueError(f"Malformed gene hit: {item!r}")
            result.append((match[1], number(match[2])))
        return sorted(result)
    if field in ("cloneCount", "cloneFraction"):
        return number(value)
    return value


def compare(original, prefix, input_reads):
    counts, differences = Counter(), []

    def check(scope, key, field, expected, observed):
        counts[scope + "_fields_checked"] += 1
        if normalize(field, expected) != normalize(field, observed):
            differences.append(dict(scope=scope, key=key, field=field, original=expected, python=observed))

    oa = keyed(rows(original / "alignments.tsv"), "readId", "original alignments", True)
    pa = keyed(rows(Path(str(prefix) + ".alignments.tsv")), "readId", "Python alignments", True)
    if set(pa) != set(range(input_reads)):
        raise ValueError("Python alignment output must include every input read exactly once")
    if not set(oa) <= set(pa):
        raise ValueError("Original output contains an out-of-range read ID")
    for rid, py in sorted(pa.items()):
        expected = oa.get(rid)
        check("alignment", rid, "aligned", expected is not None, py["status"] == "aligned")
        if expected is not None and py["status"] == "aligned":
            for field in FIELDS:
                check("alignment", rid, field, expected[field], py[field])

    oc = keyed(rows(original / "clones.tsv"), "nSeqCDR3", "original clones")
    pc = keyed(rows(Path(str(prefix) + ".clones.tsv")), "nSeqCDR3", "Python clones")
    for key in sorted(set(oc) | set(pc)):
        expected, observed = oc.get(key), pc.get(key)
        check("clone", key, "present", expected is not None, observed is not None)
        if expected is not None and observed is not None:
            for field in FIELDS + ["cloneCount", "cloneFraction"]:
                check("clone", key, field, expected[field], observed[field])
    original_ids = {r["cloneId"]: r["nSeqCDR3"] for r in oc.values()}
    python_ids = {r["cloneId"]: r["nSeqCDR3"] for r in pc.values()}
    om = keyed(rows(original / "read-assignments.tsv"), "readId", "original assignments", True)
    pm = keyed(rows(Path(str(prefix) + ".read-assignments.tsv")), "readId", "Python assignments", True)
    if set(pm) != set(range(input_reads)) or not set(om) <= set(pm):
        raise ValueError("Assignment identifiers do not cover the expected input reads")
    for rid in range(input_reads):
        original_key, original_type = None, "dropped"
        if rid in om:
            clone_id, original_type = om[rid]["cloneMapping"].split(":", 1)
            original_key = original_ids.get(clone_id)
        python_key = python_ids.get(pm[rid]["cloneId"])
        check("assignment", rid, "clone_cdr3", original_key, python_key)
        if original_key is not None:
            python_type = pm[rid]["mappingType"]
            if python_type == "core,clustered":
                python_type = "clustered"
            check("assignment", rid, "mappingType", original_type, python_type)
    counts.update(input_reads=input_reads, original_aligned=len(oa),
                  python_aligned=sum(row["status"] == "aligned" for row in pa.values()),
                  original_clones=len(oc), python_clones=len(pc), mismatch_count=len(differences))
    return dict(passed=not differences, counts=dict(counts), differences=differences)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--jar", required=True, type=Path, help="User-supplied, verified MiXCR 3.0.13 JAR")
    parser.add_argument("--java", default="java", help="Java executable or path (default: java on PATH)")
    parser.add_argument("--fastq", type=Path, default=ROOT / "examples" / "trb.fastq")
    parser.add_argument("--output", required=True, type=Path, help="New directory for logs, tables and report.json")
    args = parser.parse_args(argv)
    output = args.output.resolve()
    try:
        output.mkdir(parents=True, exist_ok=False)
    except OSError as error:
        print(f"Cannot create a new output directory: {error}", file=sys.stderr)
        return 2
    report = dict(passed=False, scope="human single-end default alignment and CDR3 assembly",
                  baseline_version="MiXCR 3.0.13", expected_jar_sha256=ORIGINAL_SHA256,
                  normalization=["gene hits compared as sorted gene/score pairs", "numeric text compared as decimal values",
                                 "clone identity is unique CDR3 nucleotide sequence", "clone IDs and row order are not comparison keys",
                                 "assigned-read mapping types are compared; unassigned reasons are grouped as dropped"],
                  commands=[])

    def run(label, command, env):
        start = time.perf_counter()
        log = output / (label + ".log")
        entry = dict(label=label, command=list(map(str, command)), log=log.name)
        report["commands"].append(entry)
        try:
            with log.open("wb") as stream:
                process = subprocess.run(entry["command"], stdout=stream, stderr=subprocess.STDOUT,
                                         cwd=ROOT, env=env, timeout=300, check=False)
            entry["returncode"] = process.returncode
        finally:
            entry["seconds"] = round(time.perf_counter() - start, 3)
            save(output / "report.json", report)
        if process.returncode:
            raise RuntimeError(f"{label} exited {process.returncode}; see {log.name}")
        print(label, "passed", flush=True)

    try:
        # Absolute lexical paths also work in restricted Windows environments
        # where querying a file's final resolved name can be denied.
        jar, fastq = args.jar.absolute(), args.fastq.absolute()
        if not jar.is_file() or not fastq.is_file():
            raise ValueError("Both --jar and --fastq must identify readable files")
        report["jar_sha256"] = sha256(jar)
        if report["jar_sha256"] != ORIGINAL_SHA256:
            raise ValueError("JAR SHA-256 differs from the verified MiXCR 3.0.13 artifact")
        report["fastq_sha256"] = sha256(fastq)
        input_reads = count_fastq(fastq)
        original, python = output / "original", output / "python"
        original.mkdir(); python.mkdir()
        vdjca, clna = original / "alignments.vdjca", original / "clones.clna"
        base = [args.java, "-Xmx1g", "-jar", jar]
        env = os.environ.copy()
        local_package = ROOT / "python"
        if (local_package / "mixcrpy").is_dir():
            env["PYTHONPATH"] = str(local_package) + (os.pathsep + env["PYTHONPATH"] if env.get("PYTHONPATH") else "")
        operations = [
            ("original-align", ["align", "-s", "hs", "-t", "1", "-OsaveOriginalReads=true", "-j", original / "align-report.json", fastq, vdjca]),
            ("original-assemble", ["assemble", "-a", "-t", "1", "-j", original / "assemble-report.json", vdjca, clna]),
            ("original-export-alignments", ["exportAlignments", "-readIds", "-vHitsWithScore", "-dHitsWithScore", "-jHitsWithScore", "-cHitsWithScore", "-nFeature", "CDR3", "-aaFeature", "CDR3", "-qFeature", "CDR3", vdjca, original / "alignments.tsv"]),
            ("original-export-clones", ["exportClones", "-cloneId", "-count", "-fraction", "-vHitsWithScore", "-dHitsWithScore", "-jHitsWithScore", "-cHitsWithScore", "-nFeature", "CDR3", "-aaFeature", "CDR3", "-qFeature", "CDR3", clna, original / "clones.tsv"]),
            ("original-export-assignments", ["exportAlignments", "-readIds", "-cloneIdWithMappingType", "-nFeature", "CDR3", clna, original / "read-assignments.tsv"]),
        ]
        for label, operation in operations:
            run(label, base + operation, env)
        prefix = python / "run"
        run("python-analyze", [sys.executable, "-m", "mixcrpy", "analyze", fastq, prefix], env)
        report.update(compare(original, prefix, input_reads))
        save(output / "report.json", report)
        print(json.dumps({key: report[key] for key in ("passed", "counts")}, indent=2))
        return 0 if report["passed"] else 1
    except (OSError, ValueError, KeyError, InvalidOperation, RuntimeError, subprocess.TimeoutExpired) as error:
        report["error"] = f"{type(error).__name__}: {error}"
        save(output / "report.json", report)
        print(report["error"], file=sys.stderr)
        return 2


if __name__ == "__main__":
    raise SystemExit(main())
