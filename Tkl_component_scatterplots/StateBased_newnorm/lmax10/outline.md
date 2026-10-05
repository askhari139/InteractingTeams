Outline: 
Two-team networks give rise to robust, canalizing systems. Irrespective of the size of number of nodes, the phenotypic landscape emergent from two-team networks is often bimodal, with the two steady states exhibiting strong resilience to perturbations of various kind. Previously, we have treated two team networks as isolated entities while interrogating their response to stochasticity and perturbations. However, these networks, each representing an individual cell-fate decision module, often interact with other modules. In this project, we aim to evaluate the effect of such interactions on the emergent phenotypic landscapes of two-team networks. Using a set of large-scale GRNs constructed in the context of endothelial cells (Sullivan et al, 2023, Sizek et al., 2024, Greene et al., 2025), modelling the interaction of various cell-fate decision modules such as EMT, cell-cycle, and apoptosis, we first characterize the nature of inter-module interactions based on how well they align with the team configurations within each module. Then, using artificial networks and the biological examples, we study how these interactions affect the phenotypic landscapes of the individual modules and the characteristics of the collective network’s landscape. We try to answer the question of how important the inter-module interactions are to predict the emergent landscape.



Network analysis PPT: MechanosensingNet
Figures PPT: Figure_plan

Figure outline:
Figure 1 reiterates the findings from earlier papers on teams using the modules from Sullivan et al. Figure 2 demonstrates the diversity in the inter-module interactions, and the corresponding diversity in how the phenotypic landscapes of individual modules are effected. Figure 3 aims to characterize the inter-module interaction strengths using various interaction metrics. Figure 4 connects these interaction metrics with the properties of the emergent phenotypic landscape, to score the metrics. Figure 5 verifies these scores using artificial networks created to preserve the metrics like the biological network. I expect that the metrics I proposed in the layout below will not be good predictors and correspondingly will lead to a significant diversity in the artificial network analysis, but once we see these results, we will be in a better position to figure out appropriate metrics. 

Simulation layout:
We will primarily use Ising formalism to carry out our analysis, with RACIPE to support our results wherever possible. We will use the logical rules to validate our findings.

Figures:

Properties of individual modules and demonstration of the effect of two teams on the phenotypic landscape
Module properties
Influence matrices of strong and weak team modules (Apoptotic switch vs phase switch) in the network - demonstrating the diversity of network structures
Heatmap showing the important properties of individual network modules - team strength, feedback loops, density and impurity for AS, PS, RS, EMT, CIP and Migration. Same heatmap for other modules in supplementary, along with network science-based analysis of the larger network
Properties of the phenotypic landscape
Heatmap of steady states (ising) demonstrating the nature of teams, use of phenotypic scores as a means to simplify steady states. Smaller adjacent heatmap on the right showing the phenotypic scores of the steady states, demonstrating the distribution of the phenotypic score
Violin plot with the phenotypic score on y axis, network name with the team strength on the x axis and colored by the simulation formalism (we’ll simulate the modules using ising, RACIPE and logical formalisms). The networks will be arranged by the team strength.
Similar violin plots as ii, but frustration on the y axis (probably in supplementary)
Similar plots as ii, but total frequency of hybrid states (Fhy) on y axis (hybrid states defined strictly on the phenotypic score to start with. The definition can be revisited based on the results)
Similar plots as ii, but for the average extent of perturbation (Pmean) required to transition between terminal states. 
Conclusive trends against team strength
Team strength on the x axis and the phenotypic score on y-axis, with two colored lines demonstrating the phenotypic score of the most frequent state (PStop) and the weighted sum of the phenotypic scores of states (PSw) that make up the top 60% of the state space
PC1 variance vs team strength
Fhy vs team strength
Demonstrating the effect of module connection on the phenotypic landscape using two module combinations
Demonstration of the module-interaction networks
Three influence matrices of the interactions - one case of no interaction between the modules, another case of interactions that align closely with the team configuration of the original module and the third with the weaker, possibly conflicting nature of team interactions.
Effect of different interactions on the properties of phenotypic landscape - ising and RACIPE
Heatmap with modules on x and y axes, and tile color PStop. The modules are arranged by team strength and labelled the same as 2bii. 
Same as i but for PSw.
Same as i but for PC1Var
Same as i but for FHy
Characterizing the nature of inter-module interactions in an effort to identify metrics to describe them meaningfully
Properties of inter-module interactions
Voilin plot with the number of modules on the x-axis and edge density on y axis - normalized against the densities within modules (n*density across / sum(density within)) if makes better sense. 
Same as i, but with the y axis as the probability of edges originating from the same team in one module to different teams in another module will have opposite signs (This can serve as a interaction-matrix level measure of team strength, and I predict would not be as informative as the influence matrix based team strength. Side note - we should compare this metric with team strength for single modules too)
Same as i, but with the net-intermodule interaction strength captured using the influence matrix
Any other metrics?


Relationship between the changes in phenotypic landscape and the metrics measuring inter-module interactions 
Same as 4, but for artificial networks

Metrics calculated
Consider two module (A and B) combinations - each with two teams (A1, A2, B1 and B2). Let E(T) be the fraction of nodes (out of the set of nodes T) active in a given steady state
Measures of the phenotypic landscape
Phenotypic score of individual modules - top steady state, weighted average over top 30% and top 60%
Module phenotypic score per steady state = |E(A1) - E(A2)|
n% phenotypic score = Sum(Module phenotypic score of steady state * SSF) 
Phenotypic score to all modules (sum of the phenotypic score of individual modules)
Total phenotypic score = |E(A1) - E(A2)| + |E(B1) - E(B2)|
% explained variance explained by the first “n” PC axes - calculated for n = 1,2,3
Number of PC axes needed to explain 90% variance
Net expression (inspired from the ideas of conditional probability) - |E(A1) - E(A2)| * |E(B1) - E(B2)|

Measures of the interactions


