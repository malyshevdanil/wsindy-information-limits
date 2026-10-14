### S23.2. Boundary-valid intervals for paired frequency differences

Let $p_{10}=P(W=1,M=0)$ and $p_{01}=P(W=0,M=1)$ for the targeted Wald
outcome $W$ and MSTLS inclusion $M$ on the same path. Their frequency
difference is $\Delta=p_{10}-p_{01}$. For $n$ independent paths, the two
discordant counts have binomial marginals; they are not independent of each
other. We construct two-sided Clopper–Pearson intervals [@clopper1934]
$[L_{10},U_{10}]$ and $[L_{01},U_{01}]$, each with confidence
$1-\alpha/2$. Bonferroni then gives simultaneous coverage at least
$1-\alpha$ for these two probabilities, without assuming independence.
Projection yields the conservative pointwise interval

$$
\Delta\in[L_{10}-U_{01},\ U_{10}-L_{01}]\cap[-1,1].
$$

For $\alpha=0.05$, each individual interval has confidence 97.5% and each
tail probability is 0.0125. At $n_{10}=100,n_{01}=0$ the result is
$[0.914252,1]$; at $n_{10}=n_{01}=0$ it is
$[-0.042874,0.042874]$. Thus extreme samples retain uncertainty. These
intervals can be wider than paired score intervals [@newcombe1998]; we use
the conservative construction for an explicit finite-sample guarantee.
Coverage is pointwise for each grid cell, not simultaneous across the
whole alternative grid.

This is a post-study correction of uncertainty reporting prompted by
degenerate percentile bootstrap intervals. All 1600 paired outcomes,
the rejection/inclusion counts, nominal thresholds and Wilson compatibility
gate are unchanged. The original bootstrap results remain archived;
`reports/support_wald_paired_intervals.json` stores both versions and all
four joint-cell counts. The local one-unit difference remains 91 percentage
points, now with interval 78.2–96.3 points. It compares different operating
outcomes, not tests with common size or general structure-selection accuracy.
