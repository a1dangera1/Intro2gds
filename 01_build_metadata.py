#!/usr/bin/env python3
"""
01_build_metadata.py  <metadata_dir>

Builds metadata/sample_metadata.tsv (QIIME2 format) from the NCBI run table
and BioSample XML, and logs every correction to metadata_corrections.tsv.

INPUT : SraRunInfo.csv, biosample.xml, ena_filereport.tsv
OUTPUT: sample_metadata.tsv, metadata_corrections.tsv

Corrections made (the raw files are never edited, only this derived table):
  * "Unknown cultiovar"  -> "Unknown cultivar"   (typo in the submitter's
    SampleName and BioSample host field). Without the fix the six-cultivar
    design would contain a seventh, misspelled level wherever the raw string is
    used, and labels would not match the paper (Figs. 2-6 use "Unknown cultivar").
  * host "Apple_<cv>" vs "Apple-<cv>" -> the separator is inconsistent between
    canker and healthy BioSamples; the cultivar is parsed out so both map to
    the same level.
  * "Cankered wood"/"Healthy wood" -> Health_status "Canker"/"Healthy".
"""
import csv, sys, re, os
import xml.etree.ElementTree as ET

d = sys.argv[1]
runs = list(csv.DictReader(open(os.path.join(d, "SraRunInfo.csv"))))

# BioSample attributes
attrs = {}
for bs in ET.parse(os.path.join(d, "biosample.xml")).getroot().findall("BioSample"):
    a = {x.get("attribute_name"): (x.text or "").strip() for x in bs.iter("Attribute")}
    attrs[bs.get("accession")] = a

ena = {r["run_accession"]: r for r in csv.DictReader(open(os.path.join(d, "ena_filereport.tsv")), delimiter="\t")}

TYPO = {"Unknown cultiovar": "Unknown cultivar"}
ABBR = {"Aztec Fuji": "Fuji", "Cripps Pink": "CP", "Gala": "Gala", "Granniwinkle": "Nic",
        "Hewe's Crab": "HW", "Unknown cultivar": "UN"}       # labels used in paper Fig. 5C
# Canker-forming wood susceptibility, as listed in the paper's Methods
SUSC = {"Gala": "highly susceptible", "Aztec Fuji": "tolerant/medium",
        "Cripps Pink": "intermediate/moderately susceptible", "Granniwinkle": "susceptible",
        "Hewe's Crab": "highly susceptible", "Unknown cultivar": "unknown"}
ROOT = {"Gala": "B.10", "Aztec Fuji": "B.9", "Cripps Pink": "G.11", "Granniwinkle": "Nic.29 (M.9)",
        "Hewe's Crab": "G.11", "Unknown cultivar": "unknown"}

rows, fixes = [], []
for r in runs:
    raw_name = r["SampleName"]                       # e.g. "Aztec Fuji_1_Canker_HiFi.fastq"
    m = re.match(r"^(.+?)_(\d+)_(Canker|Healthy)_HiFi\.fastq$", raw_name)
    if not m:
        sys.exit(f"Unexpected SampleName format: {raw_name}")
    cv_raw, rep, status = m.groups()
    cv = TYPO.get(cv_raw, cv_raw)
    if cv != cv_raw:
        fixes.append((r["Run"], "SampleName/cultivar", cv_raw, cv, "spelling typo in SRA submission"))

    a = attrs.get(r["BioSample"], {})
    host_raw = a.get("host", "")
    host_cv = TYPO.get(re.sub(r"^Apple[_-]", "", host_raw), re.sub(r"^Apple[_-]", "", host_raw))
    if host_raw and not host_raw.startswith("Apple_"):
        fixes.append((r["Run"], "BioSample host", host_raw, f"Apple_{host_cv}",
                      "inconsistent separator ('-' vs '_') in BioSample host"))
    if host_cv != cv:
        sys.exit(f"{r['Run']}: SampleName cultivar '{cv}' != BioSample host '{host_cv}'")
    iso = a.get("isolation_source", "")
    if {"Cankered wood": "Canker", "Healthy wood": "Healthy"}.get(iso) != status:
        sys.exit(f"{r['Run']}: isolation_source '{iso}' disagrees with name status '{status}'")

    sid = f"{ABBR[cv]}{rep}_{status}"              # short, unique, QIIME2-safe sample id
    rows.append({
        "sample-id": sid, "Run": r["Run"], "BioSample": r["BioSample"],
        "Cultivar": cv, "Cultivar_abbr": ABBR[cv], "Replicate": rep,
        "Health_status": status, "Group": f"{cv}_{status}",
        "Rootstock": ROOT[cv], "Canker_susceptibility": SUSC[cv],
        "Collection_date": a.get("collection_date", ""), "Spots_NCBI": r["spots"],
        "Bases_NCBI": r["bases"], "ENA_fastq_md5": ena.get(r["Run"], {}).get("fastq_md5", ""),
        "SRA_SampleName_raw": raw_name,
    })

rows.sort(key=lambda x: (x["Health_status"], x["Cultivar"], int(x["Replicate"])))
cols = list(rows[0].keys())
with open(os.path.join(d, "sample_metadata.tsv"), "w") as fh:
    fh.write("\t".join(cols) + "\n")
    fh.write("#q2:types\t" + "\t".join("numeric" if c in ("Spots_NCBI", "Bases_NCBI") else "categorical"
                                        for c in cols[1:]) + "\n")
    for x in rows:
        fh.write("\t".join(str(x[c]) for c in cols) + "\n")

with open(os.path.join(d, "metadata_corrections.tsv"), "w") as fh:
    fh.write("Run\tField\tOriginal_value\tCorrected_value\tReason\n")
    for f in fixes:
        fh.write("\t".join(f) + "\n")

n_c = sum(x["Health_status"] == "Canker" for x in rows)
print(f"[metadata] {len(rows)} samples: {n_c} Canker, {len(rows)-n_c} Healthy; "
      f"{len({x['Cultivar'] for x in rows})} cultivars; "
      f"{sum(int(x['Spots_NCBI']) for x in rows):,} reads in NCBI; {len(fixes)} corrections logged")
