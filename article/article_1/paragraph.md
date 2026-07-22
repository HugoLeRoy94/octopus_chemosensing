# suggestion color-coding
-> blue portions can be removed, they are not technically needed, but helps the understanding
-> red portions, should probably be introduced earlier in the paper

# paragraph explaining our result : to be intergrated in the main text

\color{red}{
-> The biological role of an array of sensory receptors is to  be able to feel a wide variety of odorant, and concentration.
-> To achieve this, chemosensitive systems, such as olfaction rely on genes expansion to generate a wide diversity of receptors, each sensitive to different chemical.
-> The combined activation profile of the overall array provide a coding between the input (odorant's identity and concentration, respectively denoted $\ell$ and $c$) and the encoding signal (the array's response denoted $A$)}
-> We measure the encoding capacity of an array using the mutual information (MI):
$$
MI(A;(c,\ell)) = H(A) - H(A | c,\ell),
$$
-> where H(A) is the entropy of the encoding signal: $H(A) = \sum_A p(A) \log(A)$, and $MI$ is measured in bits
-> We consider the limit of noiseless receptor, which means $A$ is obtained deterministically from $(c,\ell)$.
-> As a result, the contional entropy becomes $0$, as the mutual information reduces to the empirical entropy of $A$, submitted to a given environment.
\color{blue}{
-> Intuitively, the MI (and in that case, the entropy $H(A)$) measures how much an organism can learn about its environment, by only looking at the encoding signal.
-> It increases with the number of encoding patterns, while decreasing with the number of redundant coding (two environmental inputs leading to the same output signal), which are indinstiguishable from the organism point of view.
-> In the best case, there exist a single code-words for every input, which maximize the activation entropy H(A)=R : the number of receptors.
}
-> It is well known that for such sensory array to work optimally, receptors must be uncorrelated [zwicker paper].
-> If heteromers massively expand the receptors' combinatorial family of receptors, they are also chemically coupled.
-> It is thus not clear whether heteromerization brings about additional information, or if it provides redundant information.
-> To answer this question, we simulate the response of an array to random sensing events.
-> To do so, we first derive a physical model of the ion-channel expand a widely used model of ion-gated channel : the Moyen-Wyaman-Changeux model to heteromers, and derive a relation between the opening probability of hetero-pentamer and homo-pentamers as a function of the affinity of the units constituing the channel with the ligand. See details in the supplementary material.
-> We consider the limit of binary response, in which case, we show that a single dissociation constant per-genes-ligand couple controls the threshold concentration at which the channel opens:
$$
EC_{50} = ...
$$
-> We next design an environment to simulate random encounter of ligands.
-> Unlike previous study [cite the world of studies] that relies on uncorrelated chemical environment, we need to include correlation otherwise we would loose the exact correlation that makes heteromers special / different from homomers.
-> To model a realistic environment, we propose that the combination of morphological and chemical properties of ligands, and binding pocket can be described by a vector embedded in a euclidian space respectively denotes $\mathbf{v}_\ell$ and $\mathbf{v}_u$.
-> We propose that the activation threshold vary between a fixed lowest value when ligand-unit couples share similar chemical properties, and a fixed largest value for ligand-unit couples that are morphochemically very different.
-> *can we find a graphic representing this idea ?*
-> We implement this idea as a gaussian radial basis function, by setting the activation threshold as : $\log(EC_{50}) \propto \exp(\|\mathbf{v}_u - \mathbf{v}_\ell\|^2)$, see details in supplementary materials.
-> To simulate how evolution optimized the chemical characteristics of the receptors to the environment, we optimize the vector $\mathbf{v}_u$ to maximize the mutual information of the array.

-> We this framework, we can generate random environments, arbitrarely complex, in fact, we define a regime of parameters where our results are independant of the exact environmental parameters see supplementary material for the details about $\mathbf{v}_\ell$ generation. We then and compare genes expansion strategy with heteromerization.
-> To make sure that the difference in MI between homomers and heteromers are characteristics of the array, independant of the exact environmental conditions, we place ourself in a regime where homomers behave like a perfect array.
-> That way, we make sure that the environment itself doesn't limit the entropy of the array, and whatever difference comes out is a true array characteristic.
-> We compute the averaged mutual information for a growing number of encoding genes, and degrees of heteromerization over many sample of environment in Fig.1.
-> Error bars are small for small arrays, and grows over large arrays, showing that most of the uncertainties come from entropy estimation rather than environmental variability.

-> Fig.1 Shows how MI increase with the heteromerization index, showcasing that heteromerization is a viable strategy to expand the coding capacity of the array.
-> There are, however two limitations : 1) the gain decay with the heteromerization index 2) at fixed number of receptors, heteromerization is less effective than genes expansion, see Supplementary material.
\color{blue}{
-> The first limitation comes from the increasing correlation between receptors as the heteromerization index increases.
-> Hence, the response of a heteromers made of 3-genes, can be almost obtained through a combination of all the 2-genes combination receptors.
-> The second point, is un-surprising as introduced, the heteromers are correlated, which by definition cannot behave "perfectly"
}


# legend

-> Mutual information of an array of receptors as a function of the number of encoding genes for growing heteromerization index.
-> R/G = 1 corresponds to the situation where each receptor is encoded by an independant gene.
-> R is the number of receptors.
-> We observe that increasing heteromerization level (raising $R/n_\text{genes}$) increase the mutual information.
-> However, the gain at each heteromerization level is decay. 
-> Notice that the growth of error bars comes from the exponentially increasing difficulty to numerically estimate entropy, see supplementary material for details.
-> The phase space of $A$ grows exponentially with $R$, which limits our ability to compute the Shannon entropy for large R. instead we compute a lower and an upper bound to the entropy, respectively the blocked shannon entropy and the Rényi entropy, see Supplementary material for details.

# vocabulary :
List of vocabulary that needs to be reduce, instead of using words that are similar, or related, adapt the phrasing to minimize the number of technical words.

- input / output signal | activation pattern (for the output)
- activation threshold / EC50 / affinity
- homomers = genes expansion strategy = us \neq heteromerization
- ligand / odors | units / protein units / channel /sub-units
- environmental conditions : what are they, 