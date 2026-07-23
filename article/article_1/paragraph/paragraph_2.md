# suggestion color-coding
-> blue portions can be removed, they are not technically needed, but help the understanding
-> red portions should probably be introduced earlier in the paper

# paragraph explaining our result : to be integrated in the main text

\color{red}{
-> The biological role of an array of sensory receptors is to sense a wide variety of odorants over a wide range of concentrations.
-> To achieve this, chemosensitive systems such as olfaction rely on gene expansion to generate a diverse repertoire of receptors, each sensitive to different chemicals \cite{valencia-montoya_evolution_2024}.
-> The combined activation profile of the array defines a code between the input (an odorant's identity and concentration, denoted $\ell$ and $c$) and the encoding signal (the array's response, denoted $A$).
-> To study the encoding capacity of an array of CR receptors, we place ourselves upstream of any post processing. 
-> Experimental probing of a complete olfactome would require to test thousands of additional ligands \cite{mainland_human_2015, munch_door_2016}.
-> For this reason, we will simulate an array of CR receptors, including heteromers, and comput their coding capacity.
}

-> Following the mainstream approach to evaluate the encoding capacity of an array \cite{fonollosa_quality_2012, del_castillo_convergent_2026}, we compute the mutual information (MI) between the response and the input:
$$
MI(A;(c,\ell)) = H(A) - H(A \mid c,\ell),
$$
-> where $H(A) = -\sum_A p(A)\log p(A)$ is the entropy of the encoding signal, and MI is expressed in bits.
-> We consider the limit of noiseless receptors, so that $A$ is obtained deterministically from $(c,\ell)$.
-> The conditional entropy then vanishes, $H(A \mid c,\ell)=0$, and the mutual information reduces to the empirical entropy $H(A)$ of the array's response under a given environment.

\color{blue}{
-> Maximizing $H(A)$ is an instance of the infomax principle \cite{bell_information-maximization_1995}: information transfer through a noiseless nonlinear channel is maximized by maximizing the output entropy.
-> Intuitively, $H(A)$ measures how much an organism can learn about its environment by looking only at the encoding signal.
-> It grows with the number of distinct activation patterns and shrinks with the number of redundant codes (two inputs mapping to the same output), which are indistinguishable from the organism's point of view.
-> In the best case there is a single codeword per input, which maximizes the activation entropy up to $H(A)=R$, the number of receptors.
}

-> It is well established that for such an array to operate optimally its receptors must be broadly tuned (each active for about half of the odorants) and mutually uncorrelated \cite{zwicker_receptor_2016, fonollosa_quality_2012}.
-> Heteromers massively expand the combinatorial family of receptors, but the subunits they share also make them chemically coupled, hence correlated.
-> It is therefore not clear whether heteromerization brings additional information or merely redundant information.
-> To answer this question, we simulate the response of an array to random sensing events.

-> We first derive a physical model of the receptor by extending the Monod-Wyman-Changeux (MWC) model of ligand-gated ion channels to heteromers, and obtain a relation between the opening probability of a hetero-pentamer and that of homo-pentamers as a function of the affinity of its constituent units for the ligand \cite{einav_monod-wyman-changeux_2017}. See the supplementary material for details.
-> Because an allosteric channel acts as a logarithmic sensor \cite{olsman_allosteric_2016}, this threshold naturally lives in log-concentration space, which is the variable our activation curve depends on.


-> We next design an environment to simulate random encounters with ligands.
-> Previous information-theoretic models of olfaction do include environmental correlations, but only in the concentration and co-occurrence of odorants \cite{zwicker_receptor_2016, tesileanu_adaptation_2019} leaving the underlying chemistry is left unstructured, with receptor affinities drawn independently.
-> This is insufficient here, for a reason specific to heteromers: receptors that share a subunit respond in a coupled way, and this coupling is the very signature of heteromerization, would only produces systematic redundancy when the chemistry of the environment is itself correlated.
-> We therefore need a realistic morpho-chemical environment, in which chemically similar ligands occupy nearby regions of chemical space.
-> We represent the combined morphological and chemical properties of ligands and binding pockets as vectors embedded in a Euclidean space represented in Fig.~\ref{fig:env}~\textbf{A}, denoted $\mathbf{v}_\ell$ and $\mathbf{v}_u$ respectively.
-> The idea of such a Euclidean morpho-chemical space was first introduced in the context of immunology \cite{perelson_theoretical_1979}, and is supported empirically \-> We represent how $c^*$ rises smoothly with $d$ in Fig.~\ref{fig:env}~\textbf{B}, saturating between a best-match floor $c^*_\text{min}$​ and a full-mismatch ceiling $c^*_\text{max}$​. This saturating shape is the natural form for a distribution of binding constants set by subsite complementarity \cite{lancet_probability_1993}, implemented as a Gaussian radial-basis kernel (supplementary material).
-> For each simulated sniff we draw a mixture (a Bernoulli presence mask and log-normal concentrations $\mathbf{c}$) and pass it through the activation curve to read out a binary activation pattern $A$.
-> Iterating this Monte-Carlo loop yields the empirical distribution $p(A)$ (Fig.~\ref{fig:env}~\textbf{C}), from which we estimate the mutual information.
-> To model how evolution tuned the chemistry of the receptors to the environment, we optimize the unit vectors vu\mathbf{v}_u
vu​ to maximize the array's mutual information \cite{tesileanu_adaptation_2019, zwicker_receptor_2016}.


-> Within this framework we can generate arbitrarily complex random environments; in fact, we identify a regime of parameters in which our results are independent of the exact environmental parameters (see supplementary material for details about the generation of $\mathbf{v}_\ell$).
-> To ensure that the MI difference between homomers and heteromers is a property of the array rather than of the environment, we place ourselves in a regime where homomers behave like a perfect array \cite{zwicker_receptor_2016}.
-> That way the environment itself does not limit the entropy of the array, and whatever residual difference comes out is a genuine array characteristic.
-> We compute the averaged mutual information for a growing number of encoding genes and degrees of heteromerization, over many environment samples, in Fig.~\ref{fig:MI}.

-> Fig.~\ref{fig:MI} shows that MI increases with the heteromerization index, demonstrating that heteromerization is a viable strategy to expand the coding capacity of the array.
-> There are, however, two limitations: (1) the gain decays with the heteromerization index. (2) at a fixed number of receptors, heteromerization is less effective than gene expansion (see supplementary material).

\color{blue}{
-> The first limitation reflects the growing redundancy between receptors as the heteromerization index increases. In partial-information-decomposition terms \cite{williams_nonnegative_2010}, the information contributed by a new high-order heteromer is increasingly redundant with the lower-order combinations. Concretely, the response of a 3-gene heteromer is almost reconstructible from the set of all 2-gene combinations, so its unique/synergistic contribution (the part that actually raises $H(A)$) shrinks.
-> The second limitation is unsurprising: as introduced above, heteromers are correlated, and a correlated array cannot by definition behave "perfectly".
}

# legend

-> Mutual information of an array of receptors as a function of the number of encoding genes, for a growing heteromerization index.
-> $R/n_\text{genes}=1$ corresponds to the case where each receptor is encoded by an independent gene; $R$ is the number of receptors.
-> Increasing the heteromerization level (raising $R/n_\text{genes}$) increases the mutual information, but the gain at each successive level decays.
-> The state space of $A$ grows exponentially with $R$ ($2^R$ in the binary case), which limits our ability to compute the exact Shannon entropy for large $R$.
-> Since $P(A)$ is a mixture over sensing events, we instead bracket $H(A)$ between an upper and a lower bound obtained from pairwise distances between the per-event activation distributions \cite{kolchinsky_estimating_2017}; see supplementary material for details.
-> Error bars are small for small arrays and grow for large arrays; this growth reflects the widening gap between these entropy bounds (numerical estimation difficulty) rather than environmental variability.
