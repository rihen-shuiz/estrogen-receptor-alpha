"""Saved-input report figure: regional_B."""
from ._common import np, read_table, frame, save

def plot(source, output, resources):
    annotations = read_table(source, 'helix_annotations.csv')
    summary = read_table(source, 'helix_distribution_and_coverage.csv')
    summary = summary.loc[summary.coverage_threshold.eq(0.8)]
    order = annotations.interval.tolist()
    labels = [f"{('H12*' if r.interval == 'H12_reference' else r.interval)} {r.start}–{r.end}" for r in annotations.itertuples()]
    fig, axes = frame('Coverage-aware ESR1 regional B z — within-chain normalization', size=(14, 6))
    axis = axes[0]
    x = np.arange(len(order))
    for shift, cohort, colour in [(-0.18, 'full', '#4477aa'), (0.18, 'combined_QC', '#228833')]:
        table = summary.loc[summary.cohort.eq(cohort)].set_index('interval').loc[order]
        axis.bar(x + shift, table.median_chain_median_z, width=0.36, color=colour, label=cohort)
    axis.axhline(0, color='#999999', lw=0.8)
    axis.set_xticks(x, labels, rotation=35, ha='right')
    axis.set(ylabel='Median chain regional median B z (dimensionless)')
    axis.legend()
    save(fig, output, '10_regional_B', 'Valid observed LBD CA B>=0; population SD ddof=0 within the same selected chain. At least 80% valid regional coverage. Full and QC-pass summaries have explicit chain counts in the saved table. Provisional transferred boundaries, including H12* 536–545, are not universal assignments or evidence of molecular motion.')
