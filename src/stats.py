"""
Statistical testing functions for responder vs non-responder comparisons.
"""

import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests


def compare_response_groups(frequency_df: pd.DataFrame) -> pd.DataFrame:
    """
    Compare cell population percentages between responders and non-responders.

    Uses a two-sided Mann-Whitney U test for each immune cell population.
    Also applies Benjamini-Hochberg FDR correction across the five tests.
    """

    results = []

    for population, group_df in frequency_df.groupby("population"):
        responders = group_df.loc[
            group_df["response"].str.lower() == "yes", "percentage"
        ].dropna()

        non_responders = group_df.loc[
            group_df["response"].str.lower() == "no", "percentage"
        ].dropna()

        if len(responders) == 0 or len(non_responders) == 0:
            statistic = None
            p_value = None
        else:
            test = mannwhitneyu(
                responders,
                non_responders,
                alternative="two-sided",
            )
            statistic = test.statistic
            p_value = test.pvalue

        results.append(
            {
                "population": population,
                "n_responders": len(responders),
                "n_non_responders": len(non_responders),
                "responder_mean_percentage": responders.mean(),
                "non_responder_mean_percentage": non_responders.mean(),
                "mean_difference_responder_minus_non_responder": (
                    responders.mean() - non_responders.mean()
                ),
                "test": "Mann-Whitney U",
                "statistic": statistic,
                "p_value": p_value,
            }
        )

    results_df = pd.DataFrame(results)

    valid_p_values = results_df["p_value"].notna()

    results_df["p_value_fdr_bh"] = None
    results_df["significant_fdr_0_05"] = False

    if valid_p_values.any():
        reject, corrected_p_values, _, _ = multipletests(
            results_df.loc[valid_p_values, "p_value"],
            alpha=0.05,
            method="fdr_bh",
        )

        results_df.loc[valid_p_values, "p_value_fdr_bh"] = corrected_p_values
        results_df.loc[valid_p_values, "significant_fdr_0_05"] = reject

    numeric_columns = [
        "responder_mean_percentage",
        "non_responder_mean_percentage",
        "mean_difference_responder_minus_non_responder",
        "statistic",
        "p_value",
        "p_value_fdr_bh",
    ]

    for column in numeric_columns:
        results_df[column] = pd.to_numeric(results_df[column], errors="coerce").round(6)

    return results_df.sort_values("p_value")