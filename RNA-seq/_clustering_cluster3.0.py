import os
import sys

# sys.argv[1] ==> txt file for clustering
# sys.argv[2] ==> K means (k) or Hierarchical clustering (h)

# output ==> cluster -u <-- here
#options:
#  -f filename   File loading
#  -cg a|m       Specifies whether to center each row (gene)
#                in the data
#                a: Subtract the mean of each row
#                m: Subtract the median of each row
#  -ng           Specifies to normalize each row (gene) in the data
#  -ca a|m       Specifies whether to center each column (microarray)
#                in the data
#                a: Subtract the mean of each column
#                m: Subtract the median of each column
#  -na           Specifies to normalize each column (microarray) in the data
#  -u jobname    Allows you to specify a different name for the output files
#  -g [0..8]     Specifies the distance measure for gene clustering
#                0: No gene clustering
#                1: Uncentered correlation
#                2: Pearson correlation
#                3: Uncentered correlation, absolute value
#                4: Pearson correlation, absolute value
#                5: Spearman's rank correlation
#                6: Kendall's tau
#                7: Euclidean distance
#                8: City-block distance
#  -e [0..8]     Specifies the distance measure for microarray clustering
#                0: No clustering
#                1: Uncentered correlation
#                2: Pearson correlation
#                3: Uncentered correlation, absolute value
#                4: Pearson correlation, absolute value
#                5: Spearman's rank correlation
#                6: Kendall's tau
#                7: Euclidean distance
#                8: City-block distance
#  -m [msca]     Specifies which hierarchical clustering method to use
#                m: Pairwise complete-linkage
#                s: Pairwise single-linkage
#                c: Pairwise centroid-linkage
#                a: Pairwise average-linkage
#  -k number     Specifies whether to run k-means clustering
#                instead of hierarchical clustering, and the number
#                of clusters k to use
#  -r number     For k-means clustering, the number of times the
#                k-means clustering algorithm is run
try:
  if sys.argv[2] == 'k':
	  os.system("cluster-1.59/src/cluster -k " + sys.argv[3] + " -cg a -e 7 -r 100 -u " + sys.argv[1][:-4] + '.e7.' + " -f " + sys.argv[1])
  elif sys.argv[2] == 'h':
	  os.system("cluster-1.59/src/cluster -cg a -g 7 -e 7 -m m -u " + sys.argv[1][:-4] + '_euclidean'+ " -f " + sys.argv[1])
except IndexError:
    print("You must specify an Input file to be clustered and which clustering method you would use\n")
