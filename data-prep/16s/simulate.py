#!/usr/bin/env python3
"""Build the simulated 16S dataset for the 16S lesson (site/lessons/16s/data).

    curl -LO https://raw.githubusercontent.com/nccr-microbiomes-org/Advanced_methods_in_microbial_community_analysis/main/6_reproducibility/16s_example/ASV_table.tsv
    python3 data-prep/16s/simulate.py ASV_table.tsv

Input: ASV_table.tsv from the NCCR course repository
(nccr-microbiomes-org/Advanced_methods_in_microbial_community_analysis,
6_reproducibility/16s_example/), the public zAMP output of the cystic fibrosis airway study
(Nassirnia et al. 2025, doi:10.1099/jmm.0.002079): 6,355 ASVs x 371 samples, rarefied to 20,000
reads, with taxonomy (EzBioCloud). Its sample IDs can't be matched to the study's metadata, so
everything about the samples below is simulated; only the ASVs, their taxonomy and the shape of
real airway communities are real.

What is simulated (all with a fixed seed):
- 11 children x 2 visits x (sputum, throat swab) = 44 samples. Each child's sputum profile is a
  real sample from the table; visit 2 is a perturbed copy of visit 1. The throat swab of a visit
  mixes the sputum profile with an "oral" profile and scales CF pathogens down, so sputum and
  throat differ moderately while each pair stays similar.
- Metadata: clinical status (some visits are pulmonary exacerbations, with more pathogens), age,
  FEV1, and DNA concentration (lower in throat swabs).
- Reagent contamination: 15 typical kit-contaminant ASVs, added to every sample in inverse
  proportion to its DNA concentration (the model behind decontam's frequency method).
- 4 extraction blanks and 2 PCR blanks: mostly contaminants, extraction blanks with some
  cross-contamination from real samples.
- 1 mock community: ZymoBIOMICS D6300, reads following its theoretical 16S composition
  (equal genomic DNA, 16S copy numbers 4-10 per genome).

The script checks its own output and exits with an error if a check fails.
"""
import gzip
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SEED = 2026
OUT = Path("site/lessons/16s/data")

N_PATIENTS = 11
EXACERBATION_VISITS = 6           # of 22 visits
PATHOGEN_GENERA = ["Pseudomonas", "Staphylococcus", "Haemophilus", "Achromobacter",
                   "Stenotrophomonas", "Moraxella"]
THROAT_ORAL_SHARE = 0.45          # throat swab = 55% sputum profile + 45% oral profile
THROAT_PATHOGEN_FACTOR = 0.35     # pathogens are scaled down in throat swabs
EXACERBATION_PATHOGEN_FACTOR = 2.0
CONC_VISIT = 2000                 # Dirichlet concentration: higher = closer to the base profile
CONC_VISIT2 = 600
CONC_THROAT = 1500
DEPTH = {"Sputum": (60_000, 0.35), "Throat swab": (35_000, 0.35)}   # median, log-sd
DNA = {"Sputum": (25.0, 0.6), "Throat swab": (2.5, 0.6)}            # ng/µl, median, log-sd
CONTAM_K = 0.25                   # contaminant fraction = k / (k + DNA concentration) ...
CONTAM_NOISE = 0.5                # ... times log-normal noise (log-sd), so the relation isn't exact
CONTAMINANTS = {                  # ASV: relative weight in the reagent "kitome"
    "ASV_2887": 6, "ASV_7200": 5, "ASV_633": 5, "ASV_5247": 2, "ASV_5580": 1, "ASV_2832": 1,
    "ASV_5517": 3, "ASV_1423": 2, "ASV_5242": 2, "ASV_2610": 1, "ASV_6156": 1,
    "ASV_3016": 2, "ASV_4020": 1, "ASV_5629": 1, "ASV_2769": 1,
}
BLANKS = {"EB1": "Extraction blank", "EB2": "Extraction blank", "EB3": "Extraction blank",
          "EB4": "Extraction blank", "PCR1": "PCR blank", "PCR2": "PCR blank"}
BLANK_CROSS = {"Extraction blank": 0.08, "PCR blank": 0.02}   # share of reads from real samples
# ZymoBIOMICS Microbial Community Standard D6300 (protocol, Tables 1 and 2):
# species, 16S copies per genome, theoretical composition (%) by genomic DNA, 16S genes, cells
ZYMO = [
    ("Pseudomonas aeruginosa", 4, 12, 4.2, 6.1),
    ("Escherichia coli", 7, 12, 10.1, 8.5),
    ("Salmonella enterica", 7, 12, 10.4, 8.8),
    ("Lactobacillus fermentum", 5, 12, 18.4, 21.9),
    ("Enterococcus faecalis", 4, 12, 9.9, 14.6),
    ("Staphylococcus aureus", 6, 12, 15.5, 15.3),
    ("Listeria monocytogenes", 6, 12, 14.1, 13.9),
    ("Bacillus subtilis", 10, 12, 17.4, 10.3),
]
ZYMO_TAX = {  # Phylum, Class, Order, Family, Genus (EzBioCloud names, as in the table)
    "Pseudomonas aeruginosa": ("Proteobacteria", "Gammaproteobacteria", "Pseudomonadales", "Pseudomonadaceae", "Pseudomonas"),
    "Escherichia coli": ("Proteobacteria", "Gammaproteobacteria", "Enterobacterales", "Enterobacteriaceae", "Escherichia"),
    "Salmonella enterica": ("Proteobacteria", "Gammaproteobacteria", "Enterobacterales", "Enterobacteriaceae", "Salmonella"),
    "Lactobacillus fermentum": ("Firmicutes", "Bacilli", "Lactobacillales", "Lactobacillaceae", "Lactobacillus"),
    "Enterococcus faecalis": ("Firmicutes", "Bacilli", "Lactobacillales", "Enterococcaceae", "Enterococcus"),
    "Staphylococcus aureus": ("Firmicutes", "Bacilli", "Bacillales", "Staphylococcaceae", "Staphylococcus"),
    "Listeria monocytogenes": ("Firmicutes", "Bacilli", "Bacillales", "Listeriaceae", "Listeria"),
    "Bacillus subtilis": ("Firmicutes", "Bacilli", "Bacillales", "Bacillaceae", "Bacillus"),
}
TAX_COLS = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]

rng = np.random.default_rng(SEED)


def dirichlet_around(p, conc):
    """A random profile around p (zeros stay zero)."""
    out = np.zeros_like(p)
    nz = p > 0
    out[nz] = rng.dirichlet(p[nz] * conc + 1e-3)
    return out


def scale(p, mask, factor):
    q = p.copy()
    q[mask] *= factor
    return q / q.sum()


def draw(p, n):
    return rng.multinomial(int(n), p / p.sum())


def lognormal(median, sd):
    return float(np.exp(rng.normal(np.log(median), sd)))


def main(table_path):
    table = pd.read_csv(table_path, sep="\t").set_index("ASV")
    meta_cols = ["ASV_length", "ASV_sequence"] + TAX_COLS + ["Confidence"]
    counts = table.drop(columns=meta_cols)
    tax = table[meta_cols].rename(columns={"ASV_sequence": "sequence", "ASV_length": "length"})
    asvs = counts.index.to_numpy()
    rel = (counts / counts.sum()).to_numpy()          # ASV x sample profiles

    pathogen = tax["Genus"].fillna("").str.contains("|".join(PATHOGEN_GENERA)).to_numpy()
    contam_idx = np.array([np.where(asvs == a)[0][0] for a in CONTAMINANTS])
    kitome = np.zeros(len(asvs))
    kitome[contam_idx] = list(CONTAMINANTS.values())
    kitome /= kitome.sum()

    # Base profiles: real samples with a clear CF-pathogen signal and a rich community
    path_share = rel[pathogen].sum(0)
    richness = (rel > 0).sum(0)
    candidates = np.where((path_share > 0.15) & (path_share < 0.8) & (richness > np.median(richness)))[0]
    bases = rng.choice(candidates, N_PATIENTS, replace=False)
    oral = rel[:, np.argsort(path_share)[:100]].mean(1)
    oral[pathogen] = 0
    oral /= oral.sum()

    # Samples and metadata
    visits = [(p, v) for p in range(N_PATIENTS) for v in (1, 2)]
    exacerbation = set(rng.choice(len(visits), EXACERBATION_VISITS, replace=False).tolist())
    ages = rng.uniform(2, 17, N_PATIENTS)
    rows, cols, truth_profiles = [], {}, {}
    for i, (p, v) in enumerate(visits):
        pid = f"P{p + 1:02d}"
        base = rel[:, bases[p]]
        visit1 = dirichlet_around(base, CONC_VISIT) if v == 1 else None
        if v == 1:
            truth_profiles[pid] = visit1
            sputum = visit1
        else:
            sputum = dirichlet_around(truth_profiles[pid], CONC_VISIT2)
        status = "Pulmonary exacerbation" if i in exacerbation else "Regular follow-up"
        if i in exacerbation:
            sputum = scale(sputum, pathogen, EXACERBATION_PATHOGEN_FACTOR)
        throat = (1 - THROAT_ORAL_SHARE) * scale(sputum, pathogen, THROAT_PATHOGEN_FACTOR) \
            + THROAT_ORAL_SHARE * oral
        throat = dirichlet_around(throat / throat.sum(), CONC_THROAT)
        age = round(ages[p] + (0.5 if v == 2 else 0), 1)
        fev1 = np.nan if age < 6 else round(float(np.clip(
            rng.normal(70 if i in exacerbation else 90, 12), 30, 130)))
        for stype, profile, code in (("Sputum", sputum, "SPU"), ("Throat swab", throat, "THR")):
            sid = f"{pid}_V{v}_{code}"
            depth = lognormal(*DEPTH[stype])
            dna = lognormal(*DNA[stype])
            f = min(0.9, CONTAM_K / (CONTAM_K + dna) * lognormal(1, CONTAM_NOISE))
            cols[sid] = draw(profile, depth * (1 - f)) + draw(dirichlet_around(kitome, 200), depth * f)
            rows.append(dict(sample_id=sid, patient=pid, visit=f"V{v}", sample_type=stype,
                             sample_group="sample", clinical_status=status, age=age, fev1=fev1,
                             dna_ng_ul=round(dna, 3)))

    real_ids = list(cols)
    for sid, stype in BLANKS.items():
        depth = lognormal(3000, 0.4)
        cross = BLANK_CROSS[stype]
        source = cols[rng.choice(real_ids)]
        cols[sid] = draw(dirichlet_around(kitome, 200), depth * (1 - cross)) + draw(source / source.sum(), depth * cross)
        rows.append(dict(sample_id=sid, patient=None, visit=None, sample_type=stype,
                         sample_group="negative_control", clinical_status=None, age=None, fev1=None,
                         dna_ng_ul=round(lognormal(0.05 if stype == "Extraction blank" else 0.02, 0.5), 3)))

    # Mock community: real table ASVs where the species is in the table, new IDs otherwise
    mock_ids, mock_new = [], []
    for species, *_ in ZYMO:
        hit = tax.index[tax["Species"].fillna("").str.startswith(species)]
        hit = [a for a in hit if a not in CONTAMINANTS]
        if hit:
            mock_ids.append(max(hit, key=lambda a: counts.loc[a].sum()))
        else:
            new = "MOCK_" + species.replace(" ", "_")
            mock_ids.append(new)
            mock_new.append((new, species))
    for new, species in mock_new:
        tax.loc[new] = [np.nan, np.nan, "Bacteria", *ZYMO_TAX[species], species, np.nan]
    asvs = tax.index.to_numpy()
    for sid in cols:                      # pad count vectors for the new mock ASVs
        cols[sid] = np.concatenate([cols[sid], np.zeros(len(mock_new), dtype=int)])
    kitome = np.concatenate([kitome, np.zeros(len(mock_new))])
    expected_16s = np.array([z[3] for z in ZYMO]) / sum(z[3] for z in ZYMO)
    mock_profile = np.zeros(len(asvs))
    mock_profile[[np.where(asvs == a)[0][0] for a in mock_ids]] = expected_16s
    mock_profile = dirichlet_around(mock_profile, 3000)
    depth, f = 40_000, 0.005
    cols["MOCK"] = draw(mock_profile, depth * (1 - f)) + draw(kitome, depth * f)
    rows.append(dict(sample_id="MOCK", patient=None, visit=None, sample_type="Mock community",
                     sample_group="positive_control", clinical_status=None, age=None, fev1=None,
                     dna_ng_ul=10.0))

    # Tables: keep ASVs seen at least once
    counts_out = pd.DataFrame(cols, index=asvs)
    counts_out = counts_out[counts_out.sum(axis=1) > 0]
    samples = pd.DataFrame(rows)
    taxonomy = tax.loc[counts_out.index, TAX_COLS + ["Confidence", "sequence"]]
    role = pd.Series("community", index=counts_out.index)
    role[role.index.isin(CONTAMINANTS)] = "reagent contaminant"
    role[role.index.isin(mock_ids)] = "mock community"
    taxonomy = taxonomy.assign(simulated_role=role)
    mock = pd.DataFrame(ZYMO, columns=["species", "copies_16s", "pct_genomic_dna", "pct_16s", "pct_cells"])
    mock["asv"] = mock_ids

    checks(counts_out, samples, taxonomy, mock)

    OUT.mkdir(parents=True, exist_ok=True)
    counts_out.index.name = "asv"
    counts_out.to_csv(OUT / "asv_counts.csv.gz")
    taxonomy.rename_axis("asv").to_csv(OUT / "taxonomy.csv.gz")
    samples.to_csv(OUT / "samples.csv", index=False)
    mock.to_csv(OUT / "mock_community.csv", index=False)
    (OUT / "provenance.json").write_text(json.dumps({
        "simulated": True,
        "source_table": "nccr-microbiomes-org/Advanced_methods_in_microbial_community_analysis "
                        "6_reproducibility/16s_example/ASV_table.tsv (zAMP output of Nassirnia et al. 2025, "
                        "doi:10.1099/jmm.0.002079; ENA PRJEB41059)",
        "real": "ASV sequences, taxonomy (EzBioCloud) and the composition of the base sputum profiles",
        "simulated_parts": "samples, pairing, metadata, depths, DNA concentrations, contamination, "
                           "negative controls, mock community reads",
        "base_samples": [counts.columns[b] for b in bases],
        "seed": SEED,
        "parameters": {k: v for k, v in globals().items() if k.isupper() and k not in ("OUT", "TAX_COLS")},
        "mock_reference": "ZymoBIOMICS Microbial Community Standard D6300 protocol, Tables 1-2",
    }, indent=2, default=str) + "\n")
    print(f"Wrote {OUT}: {counts_out.shape[0]} ASVs x {counts_out.shape[1]} samples")


def clr(x):
    lx = np.log(x + 0.5)
    return lx - lx.mean(0)


def checks(counts, samples, taxonomy, mock):
    s = samples.set_index("sample_id")
    real = s.index[s.sample_group == "sample"]
    contam = taxonomy.simulated_role == "reagent contaminant"
    frac = counts.loc[contam, real].sum() / counts[real].sum()
    rho = pd.Series(frac).rank().corr(s.loc[real, "dna_ng_ul"].rank())

    z = clr(counts[real].to_numpy())
    d = np.sqrt(((z[:, :, None] - z[:, None, :]) ** 2).sum(0))
    ids = list(real)
    pair = [d[ids.index(f"{p}_{v}_SPU"), ids.index(f"{p}_{v}_THR")]
            for p, v in s.loc[real, ["patient", "visit"]].drop_duplicates().itertuples(index=False)]
    patient = s.loc[real, "patient"].to_numpy()
    between = d[patient[:, None] != patient[None, :]]

    obs = counts.loc[mock.asv, "MOCK"]
    obs = 100 * obs / obs.sum()
    blank_contam = counts.loc[contam, list(BLANKS)].sum() / counts[list(BLANKS)].sum()

    print("Depth of real samples:", counts[real].sum().describe().round().to_dict())
    print("Contaminant fraction: sputum median %.3f, throat median %.3f; Spearman with DNA %.2f" % (
        frac[s.loc[real, "sample_type"] == "Sputum"].median(),
        frac[s.loc[real, "sample_type"] == "Throat swab"].median(), rho))
    print("Aitchison distance: within visit pairs median %.1f, between patients median %.1f" % (
        np.median(pair), np.median(between)))
    print("Contaminant share in blanks:", blank_contam.round(2).to_dict())
    print("Mock community, observed vs expected 16S %:")
    print(pd.DataFrame({"observed": obs.round(1).to_numpy(), "expected": mock.pct_16s.to_numpy(),
                        "cells": mock.pct_cells.to_numpy()}, index=mock.species).to_string())

    problems = []
    if counts[real].sum().min() < 10_000:
        problems.append("a sample has fewer than 10,000 reads")
    if rho > -0.6:
        problems.append(f"contaminant fraction doesn't fall with DNA concentration (rho = {rho:.2f})")
    if np.median(pair) >= np.median(between):
        problems.append("visit pairs aren't closer than samples from different patients")
    if (obs.to_numpy() - mock.pct_16s.to_numpy()).__abs__().max() > 2:
        problems.append("mock community is more than 2 points off its expected composition")
    if problems:
        sys.exit("Checks failed: " + "; ".join(problems))
    print("All checks passed.")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])
