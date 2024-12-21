import pandas as pd
from scipy.stats import spearmanr
from statsmodels.stats.multitest import multipletests

def correct_for_multiple_testing(
    df: pd.DataFrame, p_value_column: str = "p_value", method: str = "fdr_bh"
) -> pd.DataFrame:
    """
    Correct p-values for multiple testing using user specified method (default: Ben
    jamini-Hochberg).

    Parameters:
    df (pd.DataFrame): DataFrame containing p-values.
    p_value_column (str): Name of the column containing p-values.
    method (str): Method to use for multiple testing correction.

    Returns:
    df (pd.DataFrame): DataFrame with corrected p-values and rejection decision.
    """

    # Apply Benjamini-Hochberg (FDR) correction
    p_values = df[p_value_column]
    reject, p_value_bh, _, _ = multipletests(p_values, method=method)

    # Add a column for the adjusted p-values and rejection decision
    df['q_value'] = p_value_bh
    df['reject_null'] = reject  # True if null hypothesis is rejected (significant)

    return df


# Define function to calculate Spearman correlation
def calculate_spearman_correlation(
    taxonomic_composition: pd.DataFrame,
    exclude_columns: list = ["Sample_ID", "Cohort_Number", "Recovery_Day"],
) -> pd.DataFrame:
    """
    Calculate Spearman correlation between Recovery_Day and each column in taxonomic_composition.

    Parameters:
    taxonomic_composition (pd.DataFrame): DataFrame containing taxonomic composition data.
    exclude_columns (list): List of columns to exclude from correlation calculation.

    Returns:
    res (pd.DataFrame): DataFrame containing correlation results.
    """
    # Filter the rows where Cohort_Number equals 1 and drop NaN in 'Recovery_Day'
    tax_comp_pretransplant = taxonomic_composition[
        taxonomic_composition["Cohort_Number"] == 1
    ]
    tax_comp_pretransplant_clean = tax_comp_pretransplant.dropna(
        subset=["Recovery_Day"]
    )

    # Get columns to correlate
    columns_to_correlate = [
        col
        for col in tax_comp_pretransplant_clean.columns
        if col not in exclude_columns
    ]

    # Calculate Spearman correlation for each column
    correlation_results = []
    for column in columns_to_correlate:
        if tax_comp_pretransplant_clean[column].nunique() == 1:
            corr, p_value = float("nan"), float("nan")
            sample_size = 0
        else:
            corr, p_value = spearmanr(
                tax_comp_pretransplant_clean["Recovery_Day"],
                tax_comp_pretransplant_clean[column],
            )
            sample_size = (
                tax_comp_pretransplant_clean[["Recovery_Day", column]].dropna().shape[0]
            )

        correlation_results.append(
            {
                "correlated_column": column,
                "correlation": corr,
                "p_value": p_value,
                "sample_size": sample_size,
            }
        )

    res = pd.DataFrame(correlation_results)

    return res
