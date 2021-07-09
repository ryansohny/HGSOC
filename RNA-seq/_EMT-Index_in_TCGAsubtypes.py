import pandas as pd
import seaborn as sns
from statsmodels.stats.anova import anova_lm
from statsmodels.formula.api import ols
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import matplotlib.pyplot as plt

dat = pd.read_csv("TCGA_OV_EMThighlow_TCGAsubtypes.csv", index_col="ID")

# One-way ANOVA test
model = ols('EMT_score ~ C(gene_expression_subtype)', dat).fit()
anova_lm(model)
#                               df       sum_sq     mean_sq          F        PR(>F)
#C(gene_expression_subtype)    3.0  1584.032355  528.010785  40.226055  1.986196e-21
#Residual                    263.0  3452.161437   13.126089        NaN           NaN

# Post-hoc (Tukey test)
dat2 = dat.reset_index()
posthoc = pairwise_tukeyhsd(dat2['EMT_score'], dat2['gene_expression_subtype'], alpha=0.05)
print(posthoc)
#        Multiple Comparison of Means - Tukey HSD, FWER=0.05         
#====================================================================
#    group1         group2     meandiff p-adj   lower   upper  reject
#--------------------------------------------------------------------
#differentiated immunoreactive  -0.5184 0.8189 -2.1313  1.0944  False
#differentiated    mesenchymal   5.9032  0.001  4.1994   7.607   True
#differentiated  proliferative   1.7049 0.0353  0.0819  3.3279   True
#immunoreactive    mesenchymal   6.4216  0.001  4.7866  8.0567   True
#immunoreactive  proliferative   2.2233 0.0015  0.6726   3.774   True
#   mesenchymal  proliferative  -4.1983  0.001 -5.8434 -2.5533   True
#--------------------------------------------------------------------

# Violinplot

hue_palette ={"EMT_high": "r", "EMT_low": "b"}
sns.violinplot(x="gene_expression_subtype", y="EMT_score", data=dat, order=["proliferative", "differentiated", "immunoreactive", "mesenchymal"], inner=None, palette="muted")
#sns.boxplot(x="gene_expression_subtype", y="EMT_score", data=dat, order=["proliferative", "differentiated", "immunoreactive", "mesenchymal"], fliersize=2)
sns.swarmplot(x="gene_expression_subtype", y="EMT_score", data=dat, order=["proliferative", "differentiated", "immunoreactive", "mesenchymal"], color=".25", size=3, hue="EMT_class", palette=hue_palette)
sns.despine()
#plt.savefig("Boxplot_EMT-Index_TCGAsubtypes.pdf")
plt.savefig("Violinplot_EMT-Index_TCGAsubtypes.pdf")
plt.close()

# Chi-squared test
from scipy.stats import chi2_contingency
from itertools import combinations
from statsmodels.sandbox.stats.multicomp import multipletests

contingency_dat = dat.groupby(['gene_expression_subtype', 'EMT_class']).size().unstack()

# https://neuhofmo.github.io/chi-square-and-post-hoc-in-python/
def get_asterisks_for_pval(p_val):
    """Receives the p-value and returns asterisks string."""
    if p_val > 0.05:
        p_text = "ns"  # above threshold => not significant
    elif p_val < 1e-4:  
        p_text = '****'
    elif p_val < 1e-3:
        p_text = '***'
    elif p_val < 1e-2:
        p_text = '**'
    else:
        p_text = '*'
    
    return p_text

def chisq_and_posthoc_corrected(df):
    """Receives a dataframe and performs chi2 test and then post hoc.
    Prints the p-values and corrected p-values (after FDR correction)"""
    # start by running chi2 test on the matrix
    chi2, p, dof, ex = chi2_contingency(df, correction=True)
    print(f"Chi2 result of the contingency table: {chi2}, p-value: {p}")
    
    # post-hoc
    all_combinations = list(combinations(df.index, 2))  # gathering all combinations for post-hoc chi2
    p_vals = []
    print("Significance results:")
    for comb in all_combinations:
        new_df = df[(df.index == comb[0]) | (df.index == comb[1])]
        chi2, p, dof, ex = chi2_contingency(new_df, correction=True)
        p_vals.append(p)

    # checking significance
    # correction for multiple testing
    reject_list, corrected_p_vals = multipletests(p_vals, method='fdr_bh')[:2]
    for p_val, corr_p_val, reject, comb in zip(p_vals, corrected_p_vals, reject_list, all_combinations):
        print(f"{comb}: p_value: {p_val}; corrected: {corr_p_val} ({get_asterisks_for_pval(p_val)}) reject: {reject}")



