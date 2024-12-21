import matplotlib.pyplot as plt
import matplotlib.cm as cm
import seaborn as sns


import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

from pathlib import Path
from typing import List, Literal, Tuple


from matplotlib.patches import Ellipse
import matplotlib.pyplot as plt
from scipy import linalg
from skbio import OrdinationResults


def confidence_ellipse(x, y, ax, n_std=2.0, facecolor="none", **kwargs):
    """
    Create a plot of the covariance confidence ellipse of *x* and *y*.

    Parameters
    ----------
    x, y : array-like, shape (n, )
        Input data.
    ax : matplotlib.axes.Axes
        The axes object to draw the ellipse into.
    n_std : float
        The number of standard deviations to determine the ellipse's radius.
    **kwargs
        Forwarded to `~matplotlib.patches.Ellipse`
    """
    if x.size != y.size:
        raise ValueError("x and y must be the same size")

    # Calculate the eigenvectors and eigenvalues
    cov = np.cov(x, y)

    # Using eigenvalue decomposition
    eigvals, eigvecs = np.linalg.eigh(cov)

    # Angle in radians
    angle = np.degrees(np.arctan2(eigvecs[1, 0], eigvecs[0, 0]))

    # Width and height are "full" widths, not radius
    width, height = 2 * n_std * np.sqrt(eigvals)

    # Create the ellipse
    ellip = Ellipse(
        (np.mean(x), np.mean(y)),
        width,
        height,
        angle=angle,
        facecolor=facecolor,
        **kwargs,
    )

    # add a "cross" at the center of the ellipse
    ax.plot(np.mean(x), np.mean(y), "+", markersize=15, color=facecolor)

    ax.add_patch(ellip)


def plot_pcoa_with_confidence_ellipse(
    pcs: OrdinationResults,
    metadata: pd.DataFrame,
    category: str = "Cohort_Number",
    metric: str = "Jaccard",
    palette: str = "Set1",
    output_file: Path = None,
):
    """
    Visualizes PCoA ordination with confidence ellipses for each cohort.

    Parameters:
    - pcs: PCoA results from skbio.
    - category: Column in the metadata file to group samples by.
    - palette: Color palette to use for distinct cohorts (default is 'Set1').

    Raises:
    - ValueError: If the specified category is not found in the
        metadata DataFrame.
    """
    # get the first two axes
    pcs_data = pcs.samples.iloc[:, :2]
    pcs_data.columns = ["Axis 1", "Axis 2"]

    # merge the metadata with the PCoA data
    pcs_data_with_md = pcs_data.merge(
        metadata, left_index=True, right_index=True, how="left"
    )

    # Get unique categories
    categories = pcs_data_with_md[category].unique()

    palette = sns.color_palette(palette, n_colors=len(categories))

    g = sns.scatterplot(
        data=pcs_data_with_md, x="Axis 1", y="Axis 2", hue=category, palette=palette
    )

    # Add confidence ellipses
    for i, cohort in enumerate(categories):
        selected = pcs_data_with_md[category] == cohort
        confidence_ellipse(
            pcs_data_with_md[selected]["Axis 1"],
            pcs_data_with_md[selected]["Axis 2"],
            ax=g,
            n_std=1.0,  # 95% confidence
            facecolor=palette[i],
            alpha=0.2,
            edgecolor=palette[i],
        )

    g.set_title(f"PCoA Ordination of {category} ({metric})")
    g.set_xlabel(f"PC1 ({pcs.proportion_explained[0]*100:.2f}% variance)")
    g.set_ylabel(f"PC2 ({pcs.proportion_explained[1]*100:.2f}% variance)")

    if output_file:
        g.get_figure().savefig(output_file)


def plot_rarefaction_curves(
    data_dir: Path,
    category: str,
    selected_depth: int,
    output_file: Path = None,
    figsize: Tuple[int, int] = None,
):
    """
    Visualizes rarefaction curves for observed features and Shannon diversity at different sequencing depths.

    Parameters:
    - data_dir: Path to the directory containing the data files.
    - category: Column in the metadata file to group samples by.
    - selected_depth: Sequencing depth to highlight in the plot.
    - output_file: Path to save the plot as an image file (e.g., 'rarefaction_curves.png').
    """
    obs_feat = pd.read_csv(data_dir / "observed_features.csv")
    shannon = pd.read_csv(data_dir / "shannon.csv")
    data_cols = [c for c in obs_feat.columns if "depth" in c]
    depths = list(set([int(c.split("_")[0].split("-")[1]) for c in data_cols]))
    depths.sort()

    palette = sns.color_palette("Set2")

    boxplot_styles = {
        "showfliers": False,
        "width": 0.5,
        "linewidth": 0.2,
        "zorder": 114514,
    }

    if figsize is None:
        figsize = (8, 15)
    fig, (ax1, ax2, ax3) = plt.subplots(3, 1, figsize=(8, 15), sharex=True)
    conditions = obs_feat[category].unique()
    obs_feat_per_condition = {}
    shannon_per_condition = {}
    num_samples_per_condition = {}
    for i, condition in enumerate(conditions):
        obs_feat_sub = obs_feat[obs_feat[category] == condition]
        shannon_sub = shannon[shannon[category] == condition]
        obs_feat_box_data = []
        shannon_box_data = []
        num_samples_all = []
        for d in depths:
            # get columns for this depth (e.g., depth-10000) and drop NA rows
            selected = [c for c in data_cols if f"depth-{d}_" in c]
            depth_data = obs_feat_sub[selected].dropna(axis=0)
            data = depth_data.mean(axis=0).values
            # get number of samples at this depth
            num_samples = len(depth_data)
            num_samples_all.append(num_samples)

            # similarly, get shannon diversity data
            data_shannon = shannon_sub[selected].mean(axis=0).values

            obs_feat_box_data.append(data)
            shannon_box_data.append(data_shannon)

        obs_feat_per_condition[condition] = [np.mean(x) for x in obs_feat_box_data]
        g1 = sns.boxplot(
            data=obs_feat_box_data,
            ax=ax1,
            color=palette[i],
            **boxplot_styles,
        )
        g1.plot(
            range(len(depths)),
            obs_feat_per_condition[condition],
            color=palette[i],
            label=condition,
            zorder=-1,
        )

        shannon_per_condition[condition] = [np.mean(x) for x in shannon_box_data]
        g2 = sns.boxplot(
            data=shannon_box_data,
            ax=ax2,
            color=palette[i],
            **boxplot_styles,
        )
        g2.plot(
            range(len(depths)),
            shannon_per_condition[condition],
            color=palette[i],
            label=condition,
            zorder=-1,
        )

        num_samples_per_condition[condition] = num_samples_all
        g3 = sns.lineplot(
            x=range(len(depths)),
            y=num_samples_all,
            marker="o",
            markersize=6,
            markeredgecolor=palette[i],
            markerfacecolor=palette[i],
            ax=ax3,
            color=palette[i],
            label=condition,
        )

    g3.set_ylim(g3.get_ylim()[0] - g3.get_ylim()[1] * 0.1, g3.get_ylim()[1] * 1.1)

    # add a vertical line at the selected depth
    vline_styles = {
        "color": "red",
        "linestyle": "--",
    }
    g1.axvline(depths.index(selected_depth), **vline_styles)
    g2.axvline(depths.index(selected_depth), **vline_styles)
    g3.axvline(depths.index(selected_depth), **vline_styles)

    # annotate the selected depth
    annotation_styles = {
        "ha": "left",
        "va": "bottom",
        "fontsize": 10,
        "color": "red",
    }
    g1.text(
        depths.index(selected_depth) + 0.2,
        g1.get_ylim()[1] / 2,
        f"Selected Depth: {selected_depth}",
        **annotation_styles,
    )
    g2.text(
        depths.index(selected_depth) + 0.2,
        g2.get_ylim()[1] / 2,
        f"Selected Depth: {selected_depth}",
        **annotation_styles,
    )
    g3.text(
        depths.index(selected_depth) + 0.2,
        (g3.get_ylim()[0] + g3.get_ylim()[1]) // 2,
        f"Selected Depth: {selected_depth}",
        **annotation_styles,
    )

    # add annotations of the observed features at the selected depth
    annotation_styles2 = {
        "ha": "right",
        "va": "bottom",
        "fontsize": 10,
    }
    bbox_styles = {
        "alpha": 0.5,
        "edgecolor": "none",
    }

    for i, d in enumerate(depths):
        if d == selected_depth:
            for j, condition in enumerate(conditions):

                g1.text(
                    i - 0.2,
                    0 + j * g1.get_ylim()[1] / 10,
                    f"{condition}: {obs_feat_per_condition[condition][i] / max(obs_feat_per_condition[condition]) * 100:.1f}%",
                    bbox=dict(facecolor=palette[j], **bbox_styles),
                    **annotation_styles2,
                )
                g2.text(
                    i - 0.2,
                    0 + j * g2.get_ylim()[1] / 10,
                    f"{condition}: {shannon_per_condition[condition][i] / max(shannon_per_condition[condition]) * 100:.1f}%",
                    bbox=dict(facecolor=palette[j], **bbox_styles),
                    **annotation_styles2,
                )
                g3.text(
                    i - 0.2,
                    g3.get_ylim()[0] + j * (g3.get_ylim()[1] - g3.get_ylim()[0]) / 10,
                    f"{condition}: lose {max(num_samples_per_condition[condition]) - num_samples_per_condition[condition][i]} samples",
                    bbox=dict(facecolor=palette[j], **bbox_styles),
                    **annotation_styles2,
                )

    g1.set_xticklabels(depths, rotation=90)
    g1.set_xlabel("Sequencing Depth")
    g1.legend(title=category)
    g1.set_ylabel("Observed Features")
    g1.set_title("Observed Features vs. Sequencing Depth")

    g2.set_xticklabels(depths, rotation=90)
    g2.set_xlabel("Sequencing Depth")
    g2.legend(title=category)
    g2.set_ylabel("Shannon Diversity")
    g2.set_title("Shannon Diversity vs. Sequencing Depth")

    g3.set_xticklabels(depths, rotation=90)
    g3.set_xlabel("Sequencing Depth")
    g3.set_ylabel("Number of Samples")
    g3.set_title("Number of Samples vs. Sequencing Depth")

    fig.tight_layout()
    if output_file:
        fig.savefig(output_file)


def plot_taxonomic_abundance(
    df: pd.DataFrame,
    ranks: List[Literal["Phylum", "Class", "Order", "Family", "Genus", "Species"]],
    color_map: str = "Accent",
    sort_target: str = "Firmicutes",
    output_file: Path = None,
    figsize: Tuple[int, int] = None,
):
    """
    Plots relative abundance stacked bar plots for each taxonomic rank in ranks.

    Parameters:
    - df: DataFrame containing taxonomic data with a 'Taxon' column and abundance columns for each sample.
    - ranks: List of taxonomic ranks to extract and plot (e.g., ["Phylum", "Class", "Order", "Family", "Genus", "Species"]).
    - color_map: Colormap to use for distinct colors (default is 'Accent').
    - sort_target: The taxonomic entry to use for sample sorting (default is 'Firmicutes').
    - output_file: Path to save the plot as an image file (e.g., 'taxa-bar-plots.png').
    - figsize: Tuple specifying the width and height of the plot (default is (15, 8)).

    Raises:
    - ValueError: If the DataFrame does not contain a 'Taxon' column.
    """
    # If df does not have a 'Taxon' column, raise an error
    if "Taxon" not in df.columns:
        raise ValueError("DataFrame must contain a 'Taxon' column.")

    # Extract each taxonomic rank from the 'Taxon' column
    for rank in ranks:
        df[rank] = df["Taxon"].str.extract(rf"{rank[:1].lower()}__([^;]+)")

    # Define a color palette using the specified color map
    colors = cm.get_cmap(color_map).colors

    # Loop over each taxonomic rank and generate the plot
    for rank in ranks:
        # Prepare data for the specific rank by dropping other rank columns
        data_for_plot = df.reset_index(drop=True).drop(
            columns=["Taxon", "Confidence"] + [r for r in ranks if r != rank]
        )

        # Sum values by the current rank level for each sample
        rank_sums = data_for_plot.groupby(rank).sum()

        # Convert to relative abundance for each sample
        rank_relative = rank_sums.div(rank_sums.sum(axis=0), axis=1) * 100

        # Sort samples by the abundance of the specified sort_target (e.g., 'Firmicutes') if available
        if sort_target in rank_relative.index:
            sorted_samples = (
                rank_relative.loc[sort_target].sort_values(ascending=True).index
            )
        else:
            sort_target = (
                rank_relative.sum(axis=1).sort_values(ascending=False).index[0]
            )
            sorted_samples = (
                rank_relative.loc[sort_target].sort_values(ascending=True).index
            )
        rank_relative_sorted = rank_relative[sorted_samples]

        # Order rank indices based on total abundance across all samples
        rank_total_abundance = rank_relative_sorted.sum(axis=1)
        rank_order = rank_total_abundance.sort_values(ascending=True).index
        rank_relative_sorted = rank_relative_sorted.reindex(rank_order)

        # Map colors to the ordered indices (reverse order for a visual gradient)
        rank_colors = {
            entry: colors[i % len(colors)] for i, entry in enumerate(rank_order[::-1])
        }

        # Plot with specified colormap colors for the current rank
        fig, ax = plt.subplots(figsize=figsize if figsize else (15, 8))
        rank_relative_sorted.T.plot(
            kind="bar",
            stacked=True,
            figsize=(20, 10),
            width=0.8,
            color=[rank_colors[entry] for entry in rank_relative_sorted.index],
            ax=ax,
        )
        ax.set_xlabel("Sample", fontsize=14)
        ax.set_ylabel("Relative Abundance (%)", fontsize=14)
        ax.set_title(f"Relative Abundance of {rank} Across Samples", fontsize=18)
        ax.legend(title=rank, bbox_to_anchor=(1.05, 1), loc="upper left")
        fig.tight_layout()
        if output_file:
            fig.savefig(output_file)
