import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt


dat = pd.read_table("/clinix1/Analysis/mongol/phenomata/09.HGSOC/Genes_Revision/TCGA_EMT-index.txt", index_col="ID")

sns.boxplot(x="EMT_score", data=dat)
sns.swarmplot(x="EMT_score", data=dat, color=".25", size=3)
sns.despine()
plt.savefig("Boxplot_EMT-Index_TCGA-HGSOC.pdf")
plt.close()
