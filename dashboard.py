"""
dashboard.py

Interactive Streamlit dashboard for the Teiko technical assessment.

Run from the repository root with:

    streamlit run dashboard.py
"""

import pandas as pd
import plotly.express as px
import streamlit as st

from src.analysis import (
    get_baseline_melanoma_miraclib_pbmc_samples,
    get_baseline_project_counts,
    get_baseline_response_subject_counts,
    get_baseline_sex_subject_counts,
    get_cell_frequency_table,
    get_male_responder_baseline_bcell_average,
    get_melanoma_miraclib_pbmc_frequencies,
)
from src.db import DB_PATH
from src.stats import compare_response_groups


st.set_page_config(
    page_title="Teiko Immune Cell Analysis",
    page_icon="🧬",
    layout="wide",
)


@st.cache_data
def load_frequency_table() -> pd.DataFrame:
    """Load the full sample-level immune cell frequency table."""
    return get_cell_frequency_table()


@st.cache_data
def load_response_frequency_table() -> pd.DataFrame:
    """Load melanoma PBMC miraclib frequency data for response analysis."""
    return get_melanoma_miraclib_pbmc_frequencies()


@st.cache_data
def load_stats_table() -> pd.DataFrame:
    """Load responder versus non-responder statistical comparison results."""
    frequency_df = get_melanoma_miraclib_pbmc_frequencies()
    return compare_response_groups(frequency_df)


@st.cache_data
def load_baseline_samples() -> pd.DataFrame:
    """Load baseline melanoma PBMC miraclib samples."""
    return get_baseline_melanoma_miraclib_pbmc_samples()


@st.cache_data
def load_baseline_project_counts() -> pd.DataFrame:
    """Load baseline sample counts by project."""
    return get_baseline_project_counts()


@st.cache_data
def load_baseline_response_counts() -> pd.DataFrame:
    """Load baseline subject counts by response."""
    return get_baseline_response_subject_counts()


@st.cache_data
def load_baseline_sex_counts() -> pd.DataFrame:
    """Load baseline subject counts by sex."""
    return get_baseline_sex_subject_counts()


@st.cache_data
def load_bcell_average() -> float:
    """Load average B-cell count for baseline male melanoma responders."""
    return get_male_responder_baseline_bcell_average()


def format_population_name(population: str) -> str:
    """Convert population names into readable dashboard labels."""
    labels = {
        "b_cell": "B cells",
        "cd8_t_cell": "CD8 T cells",
        "cd4_t_cell": "CD4 T cells",
        "nk_cell": "NK cells",
        "monocyte": "Monocytes",
    }
    return labels.get(population, population)


def add_readable_population_labels(df: pd.DataFrame) -> pd.DataFrame:
    """Add readable immune-cell population labels for plotting."""
    output = df.copy()
    output["population_label"] = output["population"].map(format_population_name)
    return output


def main() -> None:
    st.title("Teiko Clinical Trial Immune Cell Analysis")

    st.write(
        "This dashboard analyzes immune-cell count data from a clinical trial. "
        "The main question is whether melanoma patients who responded to the drug "
        "candidate miraclib had different immune-cell profiles than patients who "
        "did not respond."
    )

    if not DB_PATH.exists():
        st.error(
            "Database not found. Run `python load_data.py` or `make pipeline` "
            "from the repository root before starting the dashboard."
        )
        st.stop()

    frequency_df = load_frequency_table()
    response_frequency_df = add_readable_population_labels(
        load_response_frequency_table()
    )
    stats_df = load_stats_table()
    baseline_df = load_baseline_samples()

    st.sidebar.header("Dashboard Sections")

    page = st.sidebar.radio(
        "Choose a section",
        [
            "Project Summary",
            "Cell Frequencies",
            "Responder Analysis",
            "Baseline Subset",
        ],
    )

    if page == "Project Summary":
        st.header("Project Summary")

        st.info(
            "For each sample, the five immune-cell counts are converted into "
            "percentages of the total measured immune cells in that sample. This "
            "makes it possible to compare immune-cell composition across patients "
            "even when total cell counts differ."
        )

        st.subheader("Main Takeaway")

        st.write(
            "Responders had a slightly higher average CD4 T-cell frequency than "
            "non-responders, and CD4 T cells showed the strongest unadjusted "
            "difference between response groups. However, after correcting for "
            "multiple testing across the five immune-cell populations, no cell "
            "population remained statistically significant at FDR-adjusted "
            "alpha = 0.05."
        )

        st.write(
            "In plain English, this dataset suggests a possible CD4 T-cell pattern "
            "worth follow-up, but it does not provide strong corrected statistical "
            "evidence that any single measured immune-cell population clearly "
            "predicts response to miraclib."
        )

        col1, col2, col3, col4 = st.columns(4)

        col1.metric("Total samples", f"{frequency_df['sample'].nunique():,}")
        col2.metric("Cell populations", f"{frequency_df['population'].nunique():,}")
        col3.metric(
            "Miraclib melanoma PBMC samples",
            f"{response_frequency_df['sample'].nunique():,}",
        )
        col4.metric(
            "Miraclib melanoma PBMC subjects",
            f"{response_frequency_df['subject'].nunique():,}",
        )

        st.subheader("Analysis Subset Used for Response Comparison")

        analysis_subset = response_frequency_df[
            [
                "sample",
                "subject",
                "condition",
                "treatment",
                "response",
                "sample_type",
                "time_from_treatment_start",
                "sex",
            ]
        ].drop_duplicates()

        response_counts = (
            analysis_subset.groupby("response", as_index=False)
            .agg(n_samples=("sample", "nunique"), n_subjects=("subject", "nunique"))
            .sort_values("response")
        )

        fig_response_counts = px.bar(
            response_counts,
            x="response",
            y="n_samples",
            text="n_samples",
            title="Miraclib Melanoma PBMC Samples by Treatment Response",
            labels={
                "response": "Treatment Response",
                "n_samples": "Number of Samples",
            },
        )
        fig_response_counts.update_traces(textposition="outside")
        st.plotly_chart(fig_response_counts, use_container_width=True)

        st.dataframe(response_counts, use_container_width=True)

        st.caption(
            "The response analysis is restricted to melanoma PBMC samples from "
            "patients treated with miraclib, matching the task requirements."
        )

    elif page == "Cell Frequencies":
        st.header("Cell Population Frequencies by Sample")

        st.write(
            "Bob's first question asks for the frequency of each cell type in each "
            "sample. The table below reports each immune-cell population as a "
            "percentage of the total cells counted in that sample."
        )

        sample_options = sorted(frequency_df["sample"].unique())
        selected_sample = st.sidebar.selectbox("Select sample", sample_options)

        selected_df = add_readable_population_labels(
            frequency_df[frequency_df["sample"] == selected_sample]
        )

        st.subheader(f"Selected Sample: {selected_sample}")

        total_count = selected_df["total_count"].iloc[0]
        st.metric("Total cell count", f"{total_count:,}")

        fig_sample = px.bar(
            selected_df,
            x="population_label",
            y="percentage",
            text="percentage",
            hover_data=["count", "total_count"],
            title=f"Immune-Cell Relative Frequencies for {selected_sample}",
            labels={
                "population_label": "Cell Population",
                "percentage": "Relative Frequency (%)",
                "count": "Cell Count",
                "total_count": "Total Cell Count",
            },
        )
        fig_sample.update_traces(texttemplate="%{text:.2f}%", textposition="outside")
        st.plotly_chart(fig_sample, use_container_width=True)

        display_df = selected_df[
            ["sample", "total_count", "population", "count", "percentage"]
        ].copy()
        display_df["percentage"] = display_df["percentage"].round(2)

        st.dataframe(display_df, use_container_width=True)

        st.subheader("Required Summary Table: Cell-Type Frequencies for All Samples")

        st.write(
            "This table answers Bob's first question directly. Each row represents one "
            "immune-cell population from one sample and includes the sample ID, total "
            "cell count, population name, raw count, and relative frequency percentage."
        )

        full_display_df = frequency_df.copy()
        full_display_df["percentage"] = full_display_df["percentage"].round(2)

        st.dataframe(full_display_df, use_container_width=True)

    elif page == "Responder Analysis":
        st.header("Responder vs Non-Responder Analysis")

        st.write(
            "This is the main treatment-response analysis. It compares relative "
            "immune-cell frequencies between melanoma PBMC samples from miraclib "
            "responders and non-responders."
        )

        col1, col2, col3 = st.columns(3)

        col1.metric("Filtered samples", f"{response_frequency_df['sample'].nunique():,}")
        col2.metric("Subjects", f"{response_frequency_df['subject'].nunique():,}")
        col3.metric("Cell-frequency rows", f"{len(response_frequency_df):,}")

        st.subheader("Boxplot: Immune-Cell Frequencies by Response")

        fig_box = px.box(
            response_frequency_df,
            x="population_label",
            y="percentage",
            color="response",
            points="outliers",
            title=(
                "Relative Immune-Cell Frequencies in Miraclib-Treated "
                "Melanoma PBMC Samples"
            ),
            labels={
                "population_label": "Immune-Cell Population",
                "percentage": "Relative Frequency (%)",
                "response": "Treatment Response",
            },
        )
        fig_box.update_layout(boxmode="group")
        st.plotly_chart(fig_box, use_container_width=True)

        st.subheader("Mean Difference Between Responders and Non-Responders")

        plot_stats_df = stats_df.copy()
        plot_stats_df["population_label"] = plot_stats_df["population"].map(
            format_population_name
        )
        plot_stats_df = plot_stats_df.sort_values(
            "mean_difference_responder_minus_non_responder",
            ascending=False,
        )

        fig_diff = px.bar(
            plot_stats_df,
            x="population_label",
            y="mean_difference_responder_minus_non_responder",
            text="mean_difference_responder_minus_non_responder",
            title=(
                "Average Percentage-Point Difference "
                "(Responders Minus Non-Responders)"
            ),
            labels={
                "population_label": "Immune-Cell Population",
                "mean_difference_responder_minus_non_responder": (
                    "Mean Difference in Relative Frequency "
                    "(percentage points)"
                ),
            },
        )
        fig_diff.update_traces(texttemplate="%{text:.2f}", textposition="outside")
        st.plotly_chart(fig_diff, use_container_width=True)

        st.subheader("Statistical Results")

        st.write(
            "Each cell population was compared between responders and "
            "non-responders using a two-sided Mann–Whitney U test. Because five "
            "cell populations were tested, p-values were also adjusted using "
            "Benjamini-Hochberg FDR correction."
        )

        stats_display_df = stats_df.copy()
        stats_display_df["population"] = stats_display_df["population"].map(
            format_population_name
        )

        ordered_columns = [
            "population",
            "n_responders",
            "n_non_responders",
            "responder_mean_percentage",
            "non_responder_mean_percentage",
            "mean_difference_responder_minus_non_responder",
            "p_value",
            "p_value_fdr_bh",
            "significant_fdr_0_05",
        ]

        stats_display_df = stats_display_df[ordered_columns]

        st.dataframe(stats_display_df, use_container_width=True)

        st.subheader("Plain-English Interpretation")

        cd4_row = stats_df[stats_df["population"] == "cd4_t_cell"].iloc[0]

        st.write(
            f"CD4 T cells showed the strongest unadjusted difference between "
            f"responders and non-responders. On average, responders had "
            f"{cd4_row['mean_difference_responder_minus_non_responder']:.2f} "
            f"percentage points higher CD4 T-cell frequency than non-responders "
            f"(unadjusted p = {cd4_row['p_value']:.4f})."
        )

        significant_df = stats_df[stats_df["significant_fdr_0_05"] == True]

        if significant_df.empty:
            st.info(
                "After FDR correction, none of the five immune-cell populations "
                "were statistically significant at alpha = 0.05. This means the "
                "analysis does not provide strong corrected evidence that any "
                "single measured immune-cell population clearly distinguishes "
                "miraclib responders from non-responders."
            )
        else:
            significant_populations = ", ".join(
                significant_df["population"].map(format_population_name).tolist()
            )

            st.success(
                "After FDR correction, the following immune-cell population(s) "
                f"remained statistically significant: {significant_populations}."
            )

        st.caption(
            "These results should be interpreted as exploratory associations, not "
            "proof that a cell population causes treatment response. The CD4 T-cell "
            "signal may be useful for follow-up analysis, especially with additional "
            "data or a predictive model that considers multiple features together."
        )

    elif page == "Baseline Subset":
        st.header("Baseline Melanoma PBMC Subset")

        st.write(
            "This section focuses on melanoma PBMC samples from miraclib-treated "
            "patients at baseline, where time_from_treatment_start = 0."
        )

        st.info(
            "Looking at baseline samples helps Bob understand immune-cell profiles "
            "at the start of treatment, before mixing baseline measurements with "
            "later post-treatment samples."
        )

        project_counts_df = load_baseline_project_counts()
        response_counts_df = load_baseline_response_counts()
        sex_counts_df = load_baseline_sex_counts()
        avg_bcell = load_bcell_average()

        col1, col2, col3 = st.columns(3)

        col1.metric("Baseline samples", f"{baseline_df['sample'].nunique():,}")
        col2.metric("Baseline subjects", f"{baseline_df['subject'].nunique():,}")
        col3.metric(
            "Avg B cells, male responders",
            f"{avg_bcell:.2f}",
        )

        st.subheader("Baseline Samples by Project")

        fig_project = px.bar(
            project_counts_df,
            x="project",
            y="n_samples",
            text="n_samples",
            title="Baseline Melanoma PBMC Miraclib Samples by Project",
            labels={
                "project": "Project",
                "n_samples": "Number of Samples",
            },
        )
        fig_project.update_traces(textposition="outside")
        st.plotly_chart(fig_project, use_container_width=True)

        st.dataframe(project_counts_df, use_container_width=True)

        st.subheader("Baseline Subjects by Response and Sex")

        col_left, col_right = st.columns(2)

        with col_left:
            fig_response = px.bar(
                response_counts_df,
                x="response",
                y="n_subjects",
                text="n_subjects",
                title="Baseline Subjects by Response",
                labels={
                    "response": "Response",
                    "n_subjects": "Number of Subjects",
                },
            )
            fig_response.update_traces(textposition="outside")
            st.plotly_chart(fig_response, use_container_width=True)
            st.dataframe(response_counts_df, use_container_width=True)

        with col_right:
            fig_sex = px.bar(
                sex_counts_df,
                x="sex",
                y="n_subjects",
                text="n_subjects",
                title="Baseline Subjects by Sex",
                labels={
                    "sex": "Sex",
                    "n_subjects": "Number of Subjects",
                },
            )
            fig_sex.update_traces(textposition="outside")
            st.plotly_chart(fig_sex, use_container_width=True)
            st.dataframe(sex_counts_df, use_container_width=True)

        st.subheader("Baseline Subset Interpretation")

        st.write(
            f"The baseline subset contains {baseline_df['sample'].nunique():,} "
            f"samples from {baseline_df['subject'].nunique():,} subjects. These "
            "samples come from projects prj1 and prj3. The responder and "
            "non-responder groups are similarly sized at baseline, with "
            "331 responder subjects and 325 non-responder subjects."
        )

        st.write(
            f"For melanoma male responders at time 0, the average B-cell count is "
            f"{avg_bcell:.2f}. This value is reported as a raw count because the "
            "task specifically asks for the average number of B cells, not the "
            "relative B-cell frequency."
        )

        with st.expander("View all baseline melanoma PBMC miraclib samples"):
            st.dataframe(baseline_df, use_container_width=True)


if __name__ == "__main__":
    main()