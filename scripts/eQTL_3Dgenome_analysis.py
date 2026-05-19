"""
Trans-eQTL Hi-C Contact Analysis
=================================
Analyses Hi-C contact enrichment for trans-eQTL SNP–gene pairs (gill tissue).

Steps
-----
1. Map SNP and gene TSS positions to Hi-C bins.
2. Extract Hi-C contact values (raw + GW_SCALE balanced).
3. Annotate each pair with:
   - same_tad  (intra-chromosomal pairs only)
   - same_compartment
   - same_loop (intra-chromosomal pairs only)
4. Compare observed vs. background contacts (Mann-Whitney test) and plot.
"""

import hictkpy as htk
import pandas as pd
import bioframe as bf
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt

HIC_PATH         = "HiC_gill_Ss.mcool"
TAD_PATH         = "HiC_liver_Ss_tad_50_domains.bed"
LOOP_PATH        = "gill_Ss_loop.0_05.tsv"
COMPARTMENT_PATH = "ab_compartment_gill_Ss.tsv"
GENE_POS_PATH    = "Ssal_genePos.tsv"
EQTL_PATH        = "hic_random_pairs_unique_for_fillin_with_GeneID_strand.tsv"
# EQTL_PATH        = "trans_eqtl_pairs.txt"
EQTL_RANDOM_PATH = "trans_eqtl_pairs_shuffled.txt"

RESOLUTION    = 50_000

def load_annotations():
    tad_df = bf.read_table(TAD_PATH, schema="bed12")[["chrom", "start", "end", "score"]]
    tad_df["chrom"] = "chr" + tad_df["chrom"].astype(str)

    loops = pd.read_table(LOOP_PATH)[
        ["BIN1_CHR", "BIN1_START", "BIN1_END", "BIN2_CHROMOSOME", "BIN2_START", "BIN2_END"]
    ].rename(columns={"BIN2_CHROMOSOME": "BIN2_CHR"})

    compartment = pd.read_table(COMPARTMENT_PATH)
    compartment["chrom"] = "chr" + compartment["chrom"].astype(str)

    return tad_df, loops, compartment

def load_eqtls():

    def _tss_start(row):
        return row["start"] if row["strand"] == "+" else row["end"] - 1

    def _tss_end(row):
        return row["start"] + 1 if row["strand"] == "+" else row["end"]

    gene_pos = pd.read_table(GENE_POS_PATH)
    gene_pos["TSS_start"] = gene_pos.apply(_tss_start, axis=1)
    gene_pos["TSS_end"]   = gene_pos.apply(_tss_end,   axis=1)

    # tss_start_lookup = gene_pos.set_index("geneID")["TSS_start"]
    # tss_end_lookup   = gene_pos.set_index("geneID")["TSS_end"]
    # trans_eqtls        = pd.read_csv(EQTL_PATH, sep=" ")
    # trans_eqtls["Gene_start"] = trans_eqtls["GeneID"].map(tss_start_lookup)
    # trans_eqtls["Gene_end"] = trans_eqtls["GeneID"].map(tss_end_lookup)


    trans_eqtls        = pd.read_table(EQTL_PATH)
    mask_plus = trans_eqtls['Strand'] == 1
    trans_eqtls.loc[mask_plus, 'Gene_end'] = trans_eqtls.loc[mask_plus, 'Gene_start'] + 1
    mask_minus = trans_eqtls['Strand'] == -1
    trans_eqtls.loc[mask_minus, 'Gene_start'] = trans_eqtls.loc[mask_minus, 'Gene_end'] - 1

    trans_eqtls_random = pd.read_csv(EQTL_RANDOM_PATH, sep=" ")

    return trans_eqtls, trans_eqtls_random

def _fetch_contact(f, gene_chrom, gene_start, gene_end, snp_chrom, snp_pos):
    region1 = f"{gene_chrom}:{gene_start}-{gene_end}"
    region2 = f"{snp_chrom}:{snp_pos}-{snp_pos + 1}"
    try:
        value = f.fetch(region1, region2, normalization="GW_SCALE").sum()
        return value if value is not None else np.nan
    except Exception:
        return np.nan

def build_contact_vectors(hic_path, trans_eqtls, trans_eqtls_random):
    inter = trans_eqtls[trans_eqtls["Snp_Chr"] != trans_eqtls["Gene_Chr"]]
    random_inter = (
        trans_eqtls_random[
            trans_eqtls_random["Snp_Chr"] != trans_eqtls_random["Gene_Chr"]
        ]
        .sample(frac=1, random_state=42)
        .reset_index(drop=True)
    )

    obs_contacts = []
    bg_contacts  = []

    with htk.File(hic_path, RESOLUTION) as f:
        for i, row in inter.iterrows():
            obs = _fetch_contact(
                f,
                str(row["Gene_Chr"]), int(row["Gene_start"]), int(row["Gene_end"]),
                str(row["Snp_Chr"]),  int(row["Snp_Pos"]),
            )
            rnd_row = random_inter.iloc[i]
            rdm = _fetch_contact(
                f,
                str(rnd_row["Gene_Chr"]), int(rnd_row["Gene_start"]), int(rnd_row["Gene_end"]),
                str(row["Snp_Chr"]),      int(row["Snp_Pos"]),
            )

            obs_contacts.append({
                "GeneID":     row.get("GeneID"),
                "Gene_Chr":   row["Gene_Chr"],
                "Gene_start": row["Gene_start"],
                "Snp_Chr":    row["Snp_Chr"],
                "Snp_Pos":    row["Snp_Pos"],
                "contact":    obs,
            })
            bg_contacts.append({
                "GeneID":     rnd_row.get("GeneID"),
                "Gene_Chr":   rnd_row["Gene_Chr"],
                "Gene_start": rnd_row["Gene_start"],
                "Snp_Chr":    rnd_row["Snp_Chr"],
                "Snp_Pos":    rnd_row["Snp_Pos"],
                "contact":    rdm,
            })

    return pd.DataFrame(obs_contacts), pd.DataFrame(bg_contacts)

def run_stats_and_plot(obs_df, bg_df):
    obs_vals = obs_df["contact"].dropna()
    bg_vals  = bg_df["contact"].dropna()

    mw_stat, mw_pval = stats.mannwhitneyu(obs_vals, bg_vals, alternative="greater")
    print(f"Mann-Whitney U:  stat={mw_stat:.2f},  p={mw_pval:.4e}")

    obs_nz = obs_vals[obs_vals > 0]
    bg_nz  = bg_vals[bg_vals  > 0]

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    axes[0].hist(np.log10(bg_nz),  bins=50, alpha=0.6, density=True,
                 label=f"Background",  color="steelblue")
    axes[0].hist(np.log10(obs_nz), bins=50, alpha=0.7, density=True,
                 label=f"Observed", color="tomato")
    axes[0].set_xlabel("log10(Hi-C contact)")
    axes[0].set_ylabel("Density")
    axes[0].set_title("Non-zero contacts (log scale)")
    axes[0].legend()

    for vals, label, color in [
        (obs_nz, "Observed",   "tomato"),
        (bg_nz,  "Background", "steelblue"),
    ]:
        sorted_vals = np.sort(vals)
        ecdf = np.arange(1, len(sorted_vals) + 1) / len(sorted_vals)
        axes[1].plot(np.log10(sorted_vals), ecdf, label=label, color=color, linewidth=2)
    axes[1].set_xlabel("log10(Hi-C contact)")
    axes[1].set_ylabel("Cumulative proportion")
    axes[1].set_title(f"ECDF — MW p={mw_pval:.2e}")
    axes[1].legend()

    plt.suptitle("Trans-eQTL Hi-C Contact Analysis", fontsize=13, fontweight="bold", y=1.02)
    plt.tight_layout()
    plt.show()

def get_contacts(f, row):
    chrom1 = str(row.Gene_Chr)
    chrom2 = str(row.Snp_Chr)
    gene_region = f"{chrom1}:{row.Gene_start}-{row.Gene_end}"
    snp_region  = f"{chrom2}:{row.Snp_Pos}-{row.Snp_Pos + 1}"

    if chrom1 == chrom2:
        if row.Gene_start > row.Snp_Pos:
            gene_region, snp_region = snp_region, gene_region
    else:
        chrom_order = list(f.chromosomes().keys())
        if chrom_order.index(chrom1) > chrom_order.index(chrom2):
            gene_region, snp_region = snp_region, gene_region

    norm     = "GW_SCALE" if "GW_SCALE" in f.avail_normalizations() else f.avail_normalizations()[0]
    raw      = f.fetch(gene_region, snp_region).sum()
    balanced = f.fetch(gene_region, snp_region, normalization=norm).sum()
    return raw, balanced

def in_same_tad(row, tad_df):
    if str(row["Gene_Chr"]) != str(row["Snp_Chr"]):
        return "Not same chromosome", np.nan, np.nan

    tad_df = tad_df.copy()
    tad_df["tad_id"] = np.arange(len(tad_df))

    gene_df = pd.DataFrame({
        "chrom":  [f"chr{row['Gene_Chr']}"],
        "start":  [row["Gene_start"]],
        "end":    [row["Gene_start"]],
        "GeneID": [row["GeneID"]],
    })
    snp_df = pd.DataFrame({
        "chrom":  [f"chr{row['Snp_Chr']}"],
        "start":  [row["Snp_Pos"]],
        "end":    [row["Snp_Pos"]],
        "GeneID": [row["GeneID"]],
    })

    gene_overlap = bf.overlap(gene_df, tad_df, suffixes=("_gene", "_tad"))
    snp_overlap  = bf.overlap(snp_df,  tad_df, suffixes=("_snp",  "_tad"))

    gene_tad = (
        gene_overlap[["GeneID_gene", "tad_id_tad", "start_tad", "end_tad"]]
        .rename(columns={
            "GeneID_gene": "GeneID",
            "tad_id_tad":  "tad_id_gene",
            "start_tad":   "tad_start_gene",
            "end_tad":     "tad_end_gene",
        })
        .dropna(subset=["tad_start_gene", "tad_end_gene"])
    )
    snp_tad = (
        snp_overlap[["GeneID_snp", "tad_id_tad", "start_tad", "end_tad"]]
        .rename(columns={
            "GeneID_snp": "GeneID",
            "tad_id_tad": "tad_id_snp",
            "start_tad":  "tad_start_snp",
            "end_tad":    "tad_end_snp",
        })
        .dropna(subset=["tad_start_snp", "tad_end_snp"])
    )

    if gene_tad.empty and snp_tad.empty:
        return False, np.nan, np.nan
    if gene_tad.empty:
        snp_coord = (f"chr{row['Snp_Chr']}:"
                     f"{int(snp_tad['tad_start_snp'].values[0])}-"
                     f"{int(snp_tad['tad_end_snp'].values[0])}")
        return False, np.nan, snp_coord
    if snp_tad.empty:
        gene_coord = (f"chr{row['Gene_Chr']}:"
                      f"{int(gene_tad['tad_start_gene'].values[0])}-"
                      f"{int(gene_tad['tad_end_gene'].values[0])}")
        return False, gene_coord, np.nan

    same = gene_tad["tad_id_gene"].values[0] == snp_tad["tad_id_snp"].values[0]
    gene_coord = (f"chr{row['Gene_Chr']}:"
                  f"{int(gene_tad['tad_start_gene'].values[0])}-"
                  f"{int(gene_tad['tad_end_gene'].values[0])}")
    snp_coord  = (f"chr{row['Snp_Chr']}:"
                  f"{int(snp_tad['tad_start_snp'].values[0])}-"
                  f"{int(snp_tad['tad_end_snp'].values[0])}")
    return same, gene_coord, snp_coord

def in_same_comp(row, comp_df):
    if str(row["Gene_Chr"]) != str(row["Snp_Chr"]):
        return "Not same chromosome", np.nan, np.nan

    gene_df = pd.DataFrame({
        "chrom":  [f"chr{row['Gene_Chr']}"],
        "start":  [row["Gene_start"]],
        "end":    [row["Gene_end"]],
        "GeneID": [row["GeneID"]],
    })
    snp_df = pd.DataFrame({
        "chrom":  [f"chr{row['Snp_Chr']}"],
        "start":  [row["Snp_Pos"]],
        "end":    [row["Snp_Pos"]],
        "GeneID": [row["GeneID"]],
    })

    gene_overlap = bf.overlap(gene_df, comp_df, suffixes=("_gene", "_comp"))
    snp_overlap  = bf.overlap(snp_df,  comp_df, suffixes=("_snp",  "_comp"))

    if gene_overlap["compartment_comp"].values[0] == snp_overlap["compartment_comp"].values[0]:
        gene_coord = (f"chr{row['Gene_Chr']}:"
                      f"{gene_overlap['start_comp'].values[0]}-"
                      f"{gene_overlap['end_comp'].values[0]}")
        snp_coord  = (f"chr{row['Snp_Chr']}:"
                      f"{snp_overlap['start_comp'].values[0]}-"
                      f"{snp_overlap['end_comp'].values[0]}")
        return True, gene_coord, snp_coord
    return False, np.nan, np.nan

def in_same_loop(row, loops):
    gene_chr = f"chr{row['Gene_Chr']}"
    snp_chr  = f"chr{row['Snp_Chr']}"
    gene_tss = (row["Gene_start"] + row["Gene_end"]) // 2
    snp_pos  = int(row["Snp_Pos"])

    if gene_chr != snp_chr:
        return "Not same chromosome", np.nan

    def _fmt(r):
        return (f"{r['BIN1_CHR']}:{r['BIN1_START']}-{r['BIN1_END']}"
                f"<>{r['BIN2_CHR']}:{r['BIN2_START']}-{r['BIN2_END']}")

    fwd = loops[
        (loops["BIN1_CHR"] == gene_chr) & (loops["BIN2_CHR"] == snp_chr)
        & (loops["BIN1_START"] <= gene_tss) & (loops["BIN1_END"] >= gene_tss)
        & (loops["BIN2_START"] <= snp_pos)  & (loops["BIN2_END"] >= snp_pos)
    ]
    if not fwd.empty:
        return True, _fmt(fwd.iloc[0])

    rev = loops[
        (loops["BIN1_CHR"] == snp_chr)  & (loops["BIN2_CHR"] == gene_chr)
        & (loops["BIN1_START"] <= snp_pos)  & (loops["BIN1_END"] >= snp_pos)
        & (loops["BIN2_START"] <= gene_tss) & (loops["BIN2_END"] >= gene_tss)
    ]
    if not rev.empty:
        return True, _fmt(rev.iloc[0])

    return False, np.nan

def expected_contacts_intra(f, row, normalization="GW_SCALE"):
    norm = normalization if normalization in f.avail_normalizations() else "NONE"
    chrom  = f"{row.Gene_Chr}"
    pos1   = row.Gene_start
    pos2   = row.Snp_Pos
    start  = min(pos1, pos2)
    end    = max(pos1, pos2)
    region = f"{chrom}:{start}-{end}"
    pixels = f.fetch(region, normalization=norm).describe()
    return pixels["mean"]

def expected_contacts_inter(f, chr1, chr2, normalization="GW_SCALE"):
    norm        = normalization if normalization in f.avail_normalizations() else "NONE"
    chroms      = f.chromosomes()
    chrom_order = list(chroms.keys())

    if chrom_order.index(str(chr1)) > chrom_order.index(str(chr2)):
        chr1, chr2 = chr2, chr1

    region1 = f"{chr1}:0-{chroms[str(chr1)]}"
    region2 = f"{chr2}:0-{chroms[str(chr2)]}"
    pixels  = f.fetch(region1, region2, normalization=norm).describe()
    return pixels["mean"]

def annotate_eqtls(hic_path, trans_eqtls, tad_df, compartment, loops):
    raw_list         = []
    balanced_list    = []
    expected         = []
    same_tad_flags   = []
    gene_tad_coords  = []
    snp_tad_coords   = []
    same_comp_flags  = []
    gene_comp_coords = []
    snp_comp_coords  = []
    same_loop_flags  = []
    loop_coords      = []

    with htk.File(hic_path, RESOLUTION) as f:
        for _, row in trans_eqtls.iterrows():
            raw_contact, balanced_contact = get_contacts(f, row)
            raw_list.append(raw_contact)
            balanced_list.append(balanced_contact)

            if row.Gene_Chr == row.Snp_Chr:
                exp_contacts = expected_contacts_intra(f, row)
                expected.append(exp_contacts)
                same_tad,  gene_tad_coord,  snp_tad_coord  = in_same_tad(row, tad_df)
                same_comp, gene_comp_coord, snp_comp_coord  = in_same_comp(row, compartment)
                same_loop, loop_coord                       = in_same_loop(row, loops)
            else:
                exp_contacts = expected_contacts_inter(f, row.Gene_Chr, row.Snp_Chr)
                expected.append(exp_contacts)
                same_tad,  gene_tad_coord,  snp_tad_coord  = "Not same chr", np.nan, np.nan
                same_comp, gene_comp_coord, snp_comp_coord  = "Not same chr", np.nan, np.nan
                same_loop, loop_coord                       = "Not same chr", np.nan

            same_tad_flags.append(same_tad)
            gene_tad_coords.append(gene_tad_coord)
            snp_tad_coords.append(snp_tad_coord)
            same_comp_flags.append(same_comp)
            gene_comp_coords.append(gene_comp_coord)
            snp_comp_coords.append(snp_comp_coord)
            same_loop_flags.append(same_loop)
            loop_coords.append(loop_coord)

    trans_eqtls = trans_eqtls.copy()
    trans_eqtls["raw_contact"]      = raw_list
    trans_eqtls["balanced_contact"] = balanced_list
    trans_eqtls["expected_contact"] = expected
    trans_eqtls["same_tad"]         = same_tad_flags
    # trans_eqtls["gene_tad_coord"]   = gene_tad_coords
    # trans_eqtls["snp_tad_coord"]    = snp_tad_coords
    trans_eqtls["same_comp"]        = same_comp_flags
    # trans_eqtls["gene_comp_coord"]  = gene_comp_coords
    # trans_eqtls["snp_comp_coord"]   = snp_comp_coords
    trans_eqtls["same_loop"]        = same_loop_flags
    # trans_eqtls["loop_coord"]       = loop_coords

    return trans_eqtls


def print_summaries(trans_eqtls):

    print("\n--- Pairs in same compartment ---")
    print(len(trans_eqtls[trans_eqtls["same_comp"] == True]))

    print("\n--- Pairs in same TAD ---")
    print(len(trans_eqtls[trans_eqtls["same_tad"] == True]))

    intra = trans_eqtls[trans_eqtls["Snp_Chr"] == trans_eqtls["Gene_Chr"]]
    print("\n--- Intra-chromosomal pairs sharing a loop ---")
    print(len(intra[intra["same_loop"] == True]))

def main():
    print("Loading annotation data ...")
    tad_df, loops, compartment = load_annotations()

    print("Loading eQTL tables ...")
    trans_eqtls, trans_eqtls_random = load_eqtls()

    obs_df, bg_df = build_contact_vectors(HIC_PATH, trans_eqtls, trans_eqtls_random)
    run_stats_and_plot(obs_df, bg_df)

    trans_eqtls = annotate_eqtls(HIC_PATH, trans_eqtls, tad_df, compartment, loops)
    trans_eqtls.to_csv('trans_eqtls.3D.rand.exp.csv', index=False)

    print_summaries(trans_eqtls)

if __name__ == "__main__":
    main()