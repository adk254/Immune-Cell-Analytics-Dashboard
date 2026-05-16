"""
Plotting functions for the Teiko clinical trial analysis.
"""

from pathlib import Path

import pandas as pd
import plotly.express as px


def create_response_boxplot(
    frequency_df: pd.DataFrame,
    output_path: Path,
) -> None:
    """
    Create and save an interactive boxplot comparing responders and non-responders.
    """

    fig = px.box(
        frequency_df,
        x="population",
        y="percentage",
        color="response",
        points="all",
        title=(
            "Immune Cell Population Frequencies in Melanoma PBMC Samples "
            "Treated with Miraclib"
        ),
        labels={
            "population": "Immune Cell Population",
            "percentage": "Relative Frequency (%)",
            "response": "Treatment Response",
        },
    )

    fig.update_layout(
        xaxis_title="Immune Cell Population",
        yaxis_title="Relative Frequency (%)",
        legend_title="Response",
        boxmode="group",
    )

    fig.write_html(output_path)