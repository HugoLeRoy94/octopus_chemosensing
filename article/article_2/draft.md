# information theoric constraints on the octopus chemosensitive organ

## clear in
- at first, heteromerization may seem interesting because it expands the diversity of receptors.
- however, heteromers are chemically related to homomers, which induces correlation between receptors activity.
- additionally, it is commonly accepted that multiple genes expression in a single cell is detrimental from a MI point of view.
- finally, there are thousands of possible pattern of cell, the downstream wiring is probably not one-to-one with each cell genetic expression pattern.
- then, is heteromerization still usefull ?

- we model the channel of information flux from ligands identity, and concentration, through individual ion channel receptor opening, cellular activation, to downstream wired cells.
- we find that heteromerization increase the MI at the receptor level, at the cell level, and require a specific wiring to be effective.
- additionally, we show that there exist an optimal average number of genes per cell expressed to maximize MI, and the octopus sits nearby.

## unsure
### analytic computation
1) one family case
could we use a simplified version of the activation channel, where the $c^*$ simply grows quadratically, or linearly (idk), single familly: a gaussian distribution curve, and compute the optimal position of the receptors to maximize MI ?
This represent one extreme of the chemical geometry. A case, where receptors actually work together to generate a complex fingerprint of each ligand and concentration. A situation where maximizing MI make sense because ligands are more blurred.[I don't know exactly what is the argument, but what I'm trying to say, is that to finely distinguish several ligands that share chemical characteristics, one needs receptors to act collaboratively].

2) infinitely separated families
this part only make sense if the 1 family case as been well studied, the idea is to study the extreme inverse problem. we consider several families "infinitely distant". assuming that the task is to be able to distinguish the different families, I think we can show that if the families are too chemically distincts, receptors cannot activate collaboratively, and that the only way to solve these tasks is to place one receptor in the center of each family. Which means, if one open, the others are closed and reciprocally. In that case MI scales as log(R) and not R as in the first case.

[Remaining question is do I use several ligand per sniff or not ? Screening the value of mu (# ligands per sniff) seems like a lot of work and blurring the main message.]

Knowing the two previous situation, one could ask: how do heteromers perform in each case.

### concentration channel

I know already that when optimizing for single ligand per sniff, heteromers naturally align to extract more information about concentration than what homomers do naturally. When several ligands per sniff, this natural optimization strategy disappear.

can I recover some of that effect in the infinitely separated families case ? Because heteromerizing a receptor design for a specific family with a receptor made for an absolute orthogonal family probably doesn't mean much else than shifting c* the same way for everyone: making it a little higher. Because it would perform the same way for every ligand, it cannot really discriminate different chemistry.

second topic, in each case : single family, distinct families, I can tune the concentration channel to reduce more or less the std, or even switch between log-normal to normal distribution and see how it impact the heteromer MI in both case. If the heteromers rely a lot on the concentration channel to generate information, I should see MI decay as I reduce the width of the concentraiton channel.

### genes expression

the first thing we have to mention is that the measurement of the genes expression might biased.
1) unbiased measurement. writing $h_i$ the measured histogram of number of cells expressing gene $i$, and $p_i$ the real probability that a cell express $i$. Assuming the bias works as $h_i = p_i e_i$ where $e_i$ is the error. then we can write $C_{ij} = p_{ij}/(p_i p_j) = h_{ij}/(h_i h_j)$ and that's in theory unbiased. However, I don't know how to relate that measurement to anything. it looks like a form of coexpression measurement, 0 correlation gives $C_{ij} = 0$, but complete correlation gives $C_{ij} = 1/p$, not knowing p, I cannot interpret the value of $C_{ij}$.

2) Accept the bias, and just say that we know it's just experiment and it's worth what it's worth. In that case, I have access to a couple of interesting measurements: $p_i$ the proba of expressing a gene $i$, $p(k)$ the distribution of number of genes expressed, and the coexpression or correlation.

I've tried to recapitulate $p_i$ and $p(k)$ using a max-ent model on the pattern diversity, while fixing the average $<k>$, and assuming the rank of individual genes expression. Problem, it doesn't recapitulate the correlation distribution.

The correlation seems to me like the most important feature. intuitively, one would imagine that some specific genes work together, whereas some other tend to not associate. This is not clear, aside of few blocks, the correlation seems to be unstructured, but non-negligeable. That's where the idea of joiner and loner comes from. Some genes tend to be expressed only one another genes is already expressed : joiner, whereas other can be expressed alone : loner. The problem is: how do I relate this to any MI computation ?

### pattern distribution

when looking at the ranked frequency of pattern distribution for individual sucker, or for global expression, we observe a power law with an exponent -1. This is called a Zipf law, and it's not sure whether it's informative or not. I've found an article that show that zipf law emerges naturally over max entropy maximization of pattern diversity + fixing the logarithm of the mean .... I am not sure to understand this demonstration, and not sure to understand what is the difference between fixing the mean and fixing the log of the mean. Either way, if we can relate the zipf law to fixing a mean, that can make sense in term of MI, and that seems valuable as a small brick in a broader paper.