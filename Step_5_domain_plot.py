# ── IMPORTS ───────────────────────────────────────────────────────────────────

import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as pdf_backend
import numpy as np
import os
import glob

# ── CONFIGURATION ─────────────────────────────────────────────────────────────

RESULTS_FOLDER = "domain_analysis_results"
OUTPUT_PDF = "domain_analysis_results/all_24_plots_lost_vs_preserved.pdf"

COLOUR_A = "#4C72B0"   # Blue   → prop_A
COLOUR_B = "#DD8452"   # Orange → prop_B
COLOUR_C = "#55A868"   # Green  → prop_C

# The exact strings used in your Interaction_Outcome column.
# We'll detect these automatically — but define fallbacks just in case.
# NEW:
LOST_KEYWORDS      = ["lost", "loss", "Lost", "Loss", "LOST", "negative", "Negative", "NEGATIVE"]
PRESERVED_KEYWORDS = ["preserved", "maintained", "Preserved", "Maintained", "PRESERVED", "positive", "Positive", "POSITIVE"]

# ── STEP A: Helper — plot one grouped bar chart into one axes object ──────────

def plot_domain_bars(ax, domain_stats, title, n_total_isoforms):
    # This function draws one grouped bar chart into an axes object (ax).
    # We define it as a function so we can call it twice per page
    # (once for lost, once for preserved) without repeating code.
    #
    # ax               = the axes object (the "canvas") to draw on
    # domain_stats     = a DataFrame with columns: Domain_name, prop_A, prop_B, prop_C
    # title            = the title for this subplot
    # n_total_isoforms = number of isoforms in this group (shown in title)

    if domain_stats.empty:
        # If there are no domains to plot (e.g. no lost isoforms have domain data),
        # show a message instead of a blank plot.
        ax.text(0.5, 0.5, "No data available",
                ha="center", va="center", fontsize=12, color="gray",
                transform=ax.transAxes)
        # transform=ax.transAxes means coordinates are 0-1 relative to the axes,
        # so 0.5, 0.5 is always the centre regardless of axis scale.
        ax.set_title(title, fontsize=11, fontweight="bold")
        return
        # Exit the function early — nothing more to draw.

    n_domains = len(domain_stats)
    x = np.arange(n_domains)
    # np.arange creates [0, 1, 2, ...] — one position per domain on the x-axis.

    bar_width = 0.25
    # Each of the 3 bars is 0.25 wide, total group width = 0.75.

    ax.bar(x - bar_width, domain_stats["prop_A"], bar_width,
           label="prop_A (interface in domain)",
           color=COLOUR_A, alpha=0.85)
    # Left bar — prop_A. Shifted left by one bar_width from centre.

    ax.bar(x, domain_stats["prop_B"], bar_width,
           label="prop_B (spliced in domain)",
           color=COLOUR_B, alpha=0.85)
    # Middle bar — prop_B. Centred at x.

    ax.bar(x + bar_width, domain_stats["prop_C"], bar_width,
           label="prop_C (spliced + interface) ★",
           color=COLOUR_C, alpha=0.85)
    # Right bar — prop_C. Shifted right.

    ax.set_xticks(x)
    ax.set_xticklabels(domain_stats["Domain_name"],
                       rotation=35, ha="right", fontsize=8)
    # Rotate domain name labels 35 degrees so they don't overlap.
    # ha="right" aligns the right end of each label to its tick.

    ax.set_ylabel("Mean Proportion of Domain", fontsize=9)
    ax.set_ylim(0, 1.1)
    # Y-axis from 0 to 1.1 — slightly above 1 for breathing room.

    ax.set_title(f"{title}\n(n = {n_total_isoforms} isoforms)", fontsize=10, fontweight="bold")
    # Title includes the number of isoforms in this outcome group.
    # \n puts the isoform count on a second line.
    # n_total_isoforms tells your supervisor how many isoforms each bar is based on.

    ax.legend(fontsize=7, loc="upper right")
    # Show the legend in the top-right corner with small font.

    ax.yaxis.grid(True, linestyle="--", alpha=0.5)
    ax.set_axisbelow(True)
    # Dashed grid lines behind the bars.

# ── STEP B: Find all 24 proportion CSV files ──────────────────────────────────

all_files = glob.glob(f"{RESULTS_FOLDER}/**/*_proportions.csv", recursive=True)
all_files = sorted(all_files)
print(f"Found {len(all_files)} proportion files.")

# ── STEP C: Open PDF and loop over all 24 files ───────────────────────────────

with pdf_backend.PdfPages(OUTPUT_PDF) as pdf:

    for filepath in all_files:

        # ── Identify group and model ──────────────────────────────────────────
        parts = filepath.replace("\\", "/").split("/")
        group = parts[-2]
        model = parts[-1].replace("_proportions.csv", "")
        print(f"\nProcessing: {group} / {model}")

        # ── Load the proportion data ──────────────────────────────────────────
        df = pd.read_csv(filepath)

        if df.empty:
            print(f"  SKIPPED — empty file")
            continue

        # ── Check what Interaction_Outcome values exist in this file ──────────
        outcomes = df["Interaction_Outcome"].unique()
        print(f"  Interaction outcomes found: {outcomes.tolist()}")
        # This prints what values are in the column — e.g. ["lost", "preserved"]
        # so we can verify the exact strings used in your data.

        # ── Detect lost and preserved labels automatically ────────────────────
        lost_label = None
        preserved_label = None

        for outcome in outcomes:
            outcome_lower = str(outcome).lower()
            if any(k.lower() in outcome_lower for k in LOST_KEYWORDS):
                lost_label = outcome
                # "any(...)" checks if ANY of our keywords appear in this outcome string.
                # We use .lower() on both sides so matching is case-insensitive.
            elif any(k.lower() in outcome_lower for k in PRESERVED_KEYWORDS):
                preserved_label = outcome

        print(f"  Lost label: '{lost_label}' | Preserved label: '{preserved_label}'")

        # ── Filter data into lost and preserved subsets ───────────────────────
        df_lost = df[df["Interaction_Outcome"] == lost_label] if lost_label else pd.DataFrame()
        df_pres = df[df["Interaction_Outcome"] == preserved_label] if preserved_label else pd.DataFrame()
        # df[condition] filters a DataFrame to only rows where condition is True.
        # If no label was detected, we use an empty DataFrame so the plot shows "No data".

        # ── Aggregate prop_A, prop_B, prop_C per domain for each subset ───────
        def aggregate(subset_df):
            # Takes a subset (lost OR preserved) and returns mean props per domain.
            if subset_df.empty:
                return pd.DataFrame()
                # Return empty DataFrame if no data — handled gracefully in plot function.
            stats = subset_df.groupby("Domain_name")[["prop_A", "prop_B", "prop_C"]].mean().reset_index()
            stats = stats.sort_values("prop_C", ascending=False)
            # Sort by prop_C descending so most interesting domains appear first (leftmost).
            return stats

        stats_lost = aggregate(df_lost)
        stats_pres = aggregate(df_pres)

        n_lost = df_lost["Isoform_ID"].nunique() if not df_lost.empty else 0
        n_pres = df_pres["Isoform_ID"].nunique() if not df_pres.empty else 0
        # .nunique() counts the number of UNIQUE values.
        # This tells us how many unique isoforms are in each group.

        # ── Create one figure with TWO side-by-side subplots ─────────────────
        fig, (ax_lost, ax_pres) = plt.subplots(1, 2, figsize=(20, 6))
        # plt.subplots(1, 2) creates 1 row, 2 columns of subplots.
        # figsize=(20, 6) — wide enough to show both plots clearly.
        # (ax_lost, ax_pres) unpacks the two axes objects so we can use them separately.

        # ── Draw the lost plot (left side) ────────────────────────────────────
        plot_domain_bars(
            ax_lost,
            stats_lost,
            title=f"LOST interactions\n{group} | {model}",
            n_total_isoforms=n_lost
        )

        # ── Draw the preserved plot (right side) ──────────────────────────────
        plot_domain_bars(
            ax_pres,
            stats_pres,
            title=f"PRESERVED interactions\n{group} | {model}",
            n_total_isoforms=n_pres
        )

        # ── Add a shared super-title for the whole page ───────────────────────
        fig.suptitle(f"Domain Overlap Proportions — {group}  |  Model: {model}",
                     fontsize=13, fontweight="bold", y=1.02)
        # suptitle() adds a title above both subplots.
        # y=1.02 pushes it slightly above the plots so it doesn't overlap.

        plt.tight_layout()
        # Automatically adjusts spacing between the two subplots and their labels.

        pdf.savefig(fig, bbox_inches="tight")
        plt.close(fig)

        print(f"  → Page added (lost: {n_lost} isoforms | preserved: {n_pres} isoforms)")

print(f"\nDone! PDF saved to: {OUTPUT_PDF}")
print("Each page shows LOST (left) vs PRESERVED (right) for one group × model combination.")