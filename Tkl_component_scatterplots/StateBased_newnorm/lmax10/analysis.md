
The effect of module interaction

Questions: 
1. Which module has the highest effect on correlations?
2. what are the properties of individual modules - the team strength, density within and density across? Which of these properties measures the effect on correlation best?
3. Which module has the most inconsistent clustering of networks based on their phenotypic scores?
4. Which property of the module best correlates with this inconsistency?

We first calculate the correlations of phenotypic score and network TS metrics for entire corpus of networks. Then, we pick one module at a time and consider the networks that do not have that module in them, and calculate the correlations for this corpus. For each correlation measure, we find the difference before and after removal, and plot this difference against module labels and module properties separately. 

We do the same thing for clustering. Cluster the networks based on their PS proerties (all combined together). Then, reduce the corpus for each module and run the clustering again. Compare the relative positions in the dendrogram for original and reduced corpus, and create a hamming distance. plot that against module metrics