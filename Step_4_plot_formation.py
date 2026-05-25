
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf as pdf_backend
import os
import glob
RESULTS_FOLDER = "domain_analysis_results"
OUTPUT_PDF = "domain_analysis_results/all_24_plots.pdf"
# Colours for the 3 bars per domain — chosen to be distinct and thesis-friendly.
COLOUR_A = "#4C72B0"   # Blue   → prop_A (interface residues in domain)
COLOUR_B = "#DD8452"   # Orange → prop_B (spliced residues in domain)
COLOUR_C = "#55A868"   # Green  → prop_C (both spliced AND interface — key metric)

#  Finds all 24 proportion CSV files 

all_files = glob.glob(f"{RESULTS_FOLDER}/**/*_proportions.csv", recursive=True)
all_files = sorted(all_files)
print(f"Found {len(all_files)} proportion files.")


#  Open the PDF and loop over all 24 files
with pdf_backend.PdfPages(OUTPUT_PDF) as pdf:
    for filepath in all_files:

        # Identify group and model from the file path 

        parts = filepath.replace("\\", "/").split("/")
        group = parts[-2]
        model = parts[-1].replace("_proportions.csv", "")
        
        print(f"Plotting: {group} / {model}")

        # ── Load the proportion data ──────────────────────────────────────────

        df = pd.read_csv(filepath)
        # Load this model's proportion results into a DataFrame.

        if df.empty:
            print(f"  SKIPPED — no data in {filepath}")
            continue
            # If the file has no rows, skip it — nothing to plot.
            # This shouldn't happen but is a safe check.

        # ── Aggregate prop_A, prop_B, prop_C per domain ───────────────────────

        domain_stats = df.groupby("Domain_name")[["prop_A", "prop_B", "prop_C"]].mean().reset_index()
        # .groupby("Domain_name") groups all rows by domain name.
        # [["prop_A", "prop_B", "prop_C"]] selects just these 3 columns.
        # .mean() computes the average of each proportion across all isoforms
        #   that have that domain — so each domain gets one mean value per proportion.
        # .reset_index() converts the group labels back into a regular column
        #   instead of being the row index.

        domain_stats = domain_stats.sort_values("prop_C", ascending=False)
        # Sort domains by mean prop_C (highest first).
        # This puts the most biologically interesting domains at the left of the plot.
        # ascending=False means largest value first.

        # ── Set up the figure ─────────────────────────────────────────────────

        n_domains = len(domain_stats)
        # How many domains are in this plot.

        fig_width = max(10, n_domains * 1.2)
        # Calculate figure width based on number of domains.
        # max(10, ...) ensures a minimum width of 10 inches even if there are few domains.
        # n_domains * 1.2 gives each domain about 1.2 inches of space.

        fig, ax = plt.subplots(figsize=(fig_width, 6))
        # plt.subplots() creates a figure (the whole image) and axes (the plot area).
        # fig = the whole figure object
        # ax  = the axes object where we draw the bars
        # figsize=(width, height) sets the size in inches.

        # ── Draw the grouped bars ─────────────────────────────────────────────

        import numpy as np
        # numpy is a library for numerical operations.
        # We use it here to calculate bar positions.
        # Install if needed: pip install numpy

        x = np.arange(n_domains)
        # np.arange() creates an array of evenly spaced numbers: [0, 1, 2, 3, ...]
        # One number per domain — these are the x-axis positions.

        bar_width = 0.25
        # Width of each individual bar.
        # With 3 bars per group and width 0.25, the total group width = 0.75,
        # leaving 0.25 space between groups.

        bars_A = ax.bar(x - bar_width, domain_stats["prop_A"], bar_width,
                        label="prop_A (interface in domain)",
                        color=COLOUR_A, alpha=0.85)
        # ax.bar() draws a bar chart.
        # x - bar_width shifts prop_A bars to the LEFT of centre.
        # domain_stats["prop_A"] = the heights of the bars.
        # bar_width = width of each bar.
        # label = text shown in the legend.
        # color = bar colour.
        # alpha = transparency (0=invisible, 1=fully opaque). 0.85 looks clean.

        bars_B = ax.bar(x, domain_stats["prop_B"], bar_width,
                        label="prop_B (spliced in domain)",
                        color=COLOUR_B, alpha=0.85)
        # prop_B bars go in the CENTRE (no x offset).

        bars_C = ax.bar(x + bar_width, domain_stats["prop_C"], bar_width,
                        label="prop_C (spliced + interface in domain) ★",
                        color=COLOUR_C, alpha=0.85)
        # prop_C bars go to the RIGHT of centre.
        # ★ in the label flags this as the key metric for your supervisor.

        # ── Labels and formatting ─────────────────────────────────────────────

        ax.set_xticks(x)
        # Place tick marks at positions 0, 1, 2, ... (one per domain).

        ax.set_xticklabels(domain_stats["Domain_name"], rotation=35, ha="right", fontsize=9)
        # Label each tick with the domain name.
        # rotation=35 tilts the labels 35 degrees so they don't overlap.
        # ha="right" aligns the right end of each label to its tick mark.
        # fontsize=9 makes them slightly smaller to fit more domains.

        ax.set_ylabel("Mean Proportion of Domain", fontsize=11)
        # Y-axis label.

        ax.set_xlabel("Pfam Domain", fontsize=11)
        # X-axis label.

        ax.set_ylim(0, 1.05)
        # Set Y-axis range from 0 to 1.05.
        # 1.05 (slightly above 1) gives a little breathing room at the top.

        ax.set_title(f"Domain Overlap Proportions\nGroup: {group}  |  Model: {model}",
                     fontsize=13, fontweight="bold")
        # Title with group and model name.
        # \n inserts a line break — puts group/model info on a second line.
        # fontweight="bold" makes the title stand out.

        ax.legend(loc="upper right", fontsize=9)
        # Show the legend (prop_A / prop_B / prop_C labels with colours).
        # loc="upper right" places it in the top-right corner.

        ax.yaxis.grid(True, linestyle="--", alpha=0.5)
        # Add horizontal dashed grid lines for easier reading.
        # linestyle="--" = dashed lines.
        # alpha=0.5 = semi-transparent so they don't overpower the bars.

        ax.set_axisbelow(True)
        # Draw grid lines BEHIND the bars, not on top of them.

        plt.tight_layout()
        # Automatically adjusts spacing so nothing gets cut off
        # (e.g. long x-axis labels, title, legend).

        # ── Save this plot as one page in the PDF ─────────────────────────────

        pdf.savefig(fig, bbox_inches="tight")
        # pdf.savefig() adds this figure as a new page in the PDF.
        # bbox_inches="tight" ensures nothing is clipped at the edges.

        plt.close(fig)
        # Close the figure to free up memory.
        # Without this, all 24 figures stay open in memory — on a Mac this
        # can slow things down or cause warnings.

        print(f"  → Added to PDF ({n_domains} domains)")

print(f"\nDone! PDF saved to: {OUTPUT_PDF}")
print(f"Open it to review all 24 plots.")