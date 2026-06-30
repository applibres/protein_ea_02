# Protein engineering via Bayesian optimization-guided evolutionary algorithm and robotic experiments

**Ruyun Hu†, Lihao Fu†, Yongcan Chen, Junyu Chen, Yu Qiao, Tong Si**

†Equal contribution.

Corresponding authors: Tong Si (tong.si@siat.ac.cn), Yu Qiao (yu.qiao@siat.ac.cn), Ruyun Hu (ry.hu@siat.ac.cn)

*Briefings in Bioinformatics*, 2023, 24(1), 1–9
DOI: https://doi.org/10.1093/bib/bbac570
Problem Solving Protocol

Received: September 30, 2022. Revised: November 14, 2022. Accepted: November 22, 2022
© The Author(s) 2022. Published by Oxford University Press.

---

## Abstract

Directed protein evolution applies repeated rounds of genetic mutagenesis and phenotypic screening and is often limited by experimental throughput. Through *in silico* prioritization of mutant sequences, machine learning has been applied to reduce wet lab burden to a level practical for human researchers. On the other hand, robotics permits large batches and rapid iterations for protein engineering cycles, but such capacities have not been well exploited in existing machine learning-assisted directed evolution approaches. Here, we report a scalable and batched method, Bayesian Optimization-guided EVOlutionary (BO-EVO) algorithm, to guide multiple rounds of robotic experiments to explore protein fitness landscapes of combinatorial mutagenesis libraries. We first examined various design specifications based on an empirical landscape of protein G domain B1. Then, BO-EVO was successfully generalized to another empirical landscape of an *Escherichia coli* kinase PhoQ, as well as simulated NK landscapes with up to moderate epistasis. This approach was then applied to guide robotic library creation and screening to engineer enzyme specificity of RhlA, a key biosynthetic enzyme for rhamnolipid biosurfactants. A 4.8-fold improvement in producing a target rhamnolipid congener was achieved after examining less than 1% of all possible mutants after four iterations. Overall, BO-EVO proves to be an efficient and general approach to guide combinatorial protein engineering without prior knowledge.

**Keywords:** adaptive experimental design; Bayesian optimization; evolutionary algorithm; machine learning; directed evolution

---

## Introduction

A protein fitness landscape describes the metaphoric, high-dimensional surface relating a select property ("fitness") to amino acid sequences [1]. Exploring this landscape via protein engineering is challenging for the following reasons: (i) the search space grows exponentially with the number of amino acid positions considered; (ii) functional proteins are extremely rare, and high-performance sequences decrease exponentially with increasing fitness [2, 3]; (iii) the fitness landscape is rugged attributing to epistasis [4]; (iv) experimental profiling is laborious, expensive and slow.

Pioneered by David Baker and others [5–9], rational protein redesign evaluates and recommends engineering strategies based on structural conformations, energy functions and electrostatic interactions. However, the success rate is highly dependent on the prior biophysical and biochemical understanding of a target protein, which is not always available. In contrast, directed evolution does not require prior knowledge and performs multiple rounds of random library creation and screening [1, 10], which often fixes one top mutation in each round. Although efficient, such greedy exploitation may be trapped in local optima, especially in a rugged fitness landscape due to epistasis among mutations.

Machine learning (ML) algorithms are increasingly applied to both model fitness landscapes and guide protein engineering [11]. To learn sequence-activity relationships from labeled data, supervised models have been created to predict various properties including thermostability [12, 13], fluorescence [14], ligand-binding affinity [15] and catalytic performance [16]. It has been noted that ML algorithms may present overfitting when designing stability-enhancing mutations [12]. On the other hand, gene and protein sequences accumulate at an unprecedented speed in public databases (e.g. UniProt [17]). Unsupervised models, such as UniRep [18], TAPE [19], ESM-1v [20] and ProtT5-XL-U50 [21], have been developed to learn representation from vast unlabeled data and discover latent patterns in protein sequences. Although achieving proof-of-concept successes, the performance of current ML fitness models is impaired by data scarcity and inadequate features.

To guide protein engineering, both supervised and unsupervised ML models have been utilized to improve sample efficiency. For example, ML-assisted directed evolution (MLDE) employed both zero-shot, physics-based models and experimental measurements to train fitness models [22]. The evolutionary context-integrated neural network (ECNet) exploited homologous sequence information in a supervised way to predict high-order mutational effects [23]. Moreover, utilizing learned representation from >20 million natural sequences by unsupervised models, the low-N protein engineering method required as few as 24 sequence-function data to train a supervised fitness model and exploited the model using a Monte-Carlo method [24]. Notably, current ML-guided approaches seek to keep experimental batch sizes and iteration numbers as low as possible, which are confined by the physical limit of human researchers.

On the other hand, biofoundries accelerate design-build-test-learn cycles in biological engineering through physical and informatic automation [25–27]. When applied to protein engineering, automation permits large-scale library creation and screening in a short time, so that new sequence-function data of sufficient amount and quality can be iteratively collected [28, 29], leading to continuous improvement of both model prediction and sequence design. Iterative data acquisition, fitness modeling and sequence proposal are necessary to achieve efficient feedback between ML algorithms and robotic experiments. In this regard, Bayesian optimization (BO) [30] is well suited for stepwise optimization on complex fitness landscapes, whereas model uncertainty can be explicitly considered to negotiate exploration and exploitation using a defined acquisition function. Indeed, BO has been successfully applied to engineer proteins [31], pathways [32, 33] and fermentation strains [34].

Here, we developed a batched, BO-guided EVOlutionary algorithm (BO-EVO) for protein engineering via iterations of ML models and robotic experiments. Two main customizations were applied to BO for combinatorial protein optimization, including batched experiment design and search space evolution. Using an empirical GB1 dataset of a four-residue combinatorial fitness landscape (20⁴ = 160,000) [35], we explored and optimized many design specifications, including an experimental budget, sequence representation and iteration initialization, and identified surrogate model quality and moderate batch size as key considerations. The generalization capability of BO-EVO was evaluated using another empirical landscape of PhoQ [36] and a reformulated version of the mathematical NK landscape [37], which reveals landscape ruggedness as a limiting factor. Excitingly, BO-EVO successfully guided iterative robotic engineering of RhlA mutants, achieving a 4.8-fold improvement in selective production of a rhamnolipid (RL) congener after only four iterations.

---

## Methods

### Bayesian optimization

BO [30] is a derivative-free approach to optimize an objective function that is expensive to evaluate and lacks known structures like concavity or linearity. BO builds a surrogate model of the objective, estimates model uncertainty and utilizes an acquisition function to negotiate exploration and exploitation and decide where to sample next.

We used Gaussian process regression (GPR) [38] to model protein fitness as the objective, which is defined in the Supplementary Information. A Gaussian process (GP) is a generalization of the Gaussian probability distribution. Leaving mathematical sophistication aside, one can loosely think of a function as a very long vector, each entry in the vector specifying the function value f(x) ~ Normal(μ(x), σ²(x)), at a particular input x. A GPR uses GP for regression and is fully described by a mean function μ(x) and a covariance function or kernel Σ(x, x′). The kernel is chosen so that points x, x′ that are closer in the input space have a large positive correlation, encoding the belief that they should have more similar function values than points that are far apart. In this article, zero mean function prior was used, and a scaled radial basis function (RBF) kernel was applied:

Σ(x, x′) = k² exp(−γ‖x − x′‖²₂)

where γ can be interpreted as a similarity measure of the input space and k is related to prior uncertainty.

The problem of learning in GP is exactly the problem of finding suitable properties for the kernel, specifically, determining the hyperparameters k and γ. We used a marginal log-likelihood as a loss function to estimate the kernel parameters using training sets and inferenced for a new query x_q according to the learned posterior to obtain the mean μ(x_q) and variance σ²(x_q). The readers are referred to [38] for more details. To accelerate the training and inference procedures, we used GPyTorch [39] with CUDA acceleration as an efficient GP implementation. Regarding training, maximum likelihood estimation was performed with the Adam optimizer.

We used the upper confidence bound (UCB) [40] as an acquisition function, which is to be optimistic in the face of uncertainty. And the Gaussian process upper confidence bound [41] was proposed as a Bayesian optimistic algorithm with provable cumulative regret bounds in the bandit setting. In the GP case, since the posterior at any arbitrary point x is a Gaussian, any quantile of the distribution of f(x) is computed with its corresponding value β as follows:

α_UCB(x) = μ(x) + βσ(x)

The hyperparameter β was tuned to achieve optimal performance.

### BO batching strategy and search space evolution

Classical BO takes only one new query each time. However, a batch of mutant sequences can be readily created and evaluated in parallel by robotics. In this aspect, batched BO has been practiced to design biological sequences [42]. Here, we built batches iteratively for batched BO. Within a batch, the surrogate model was kept constant for each query. Between batches, while it was updated by the newly proposed and measured batch of sequence-fitness data, specifically, we re-trained the model with all the measured data available until the current round.

To make BO scalable for each single-sequence proposal, we restricted the search space to a subspace of all the possible sequences. The subspace was evolved by the integration of the evolutionary algorithm and BO, balancing the exploration and exploitation of the whole search space. Specifically, in terms of the evolutionary algorithm, we used random mutation on only one parent sequence to generate a child population. The per-sequence mutation rate is set to introduce a single mutation in one sequence on average. One parent sequence is mutated into 6–10 child sequences suggested in Supplementary Figure S1. Then, a new sequence was selected by applying BO to the child population. Finally, a new state sequence was set according to the uncertainty of the surrogate model, that is, the newly selected sequence by BO is set as a new parent sequence if its predicted uncertainty is lower than twice that of the first parent sequence in the current round, and otherwise, a new parent sequence is sampled from the last measured batch of sequences following a fitness-based sampling strategy. With the sampling strategy, a higher fitness sequence is sampled with a higher probability proportional to the exponential function value of its fitness. The first parent sequence in each round is sampled with the same strategy.

> Method details on protein fitness landscapes, simulation setup and algorithm evaluation, protein sequence encoding strategies and wet-lab measurements can be found in the Supplementary Information.

---

## Results

### Design principles of BO-EVO and software implementation

BO-EVO integrates BO and evolutionary algorithms with complementary advantages. The computational time of BO increases exponentially with the number of targeted residues (N) when evaluating the whole protein combinatorial space (20^N). To improve algorithm scalability, BO-EVO restricts the search to subspaces generated by random mutation of a parent sequence (evolution). The parent sequence is either sampled from measured sequences by fitness-based sampling or set as a newly proposed sequence depending on the confidence level of the surrogate model (Figure 1A). On the other hand, BO allows enhanced exploration than greedy evolutionary algorithms. We used GPR as a surrogate model to quantify uncertainty. GPR is trained on measured sequences encoded by numerical protein sequence representations, such as ESM-1v [20]. To prioritize mutant sequences for experiments, we utilized UCB as the acquisition function to negotiate exploration and exploitation. Once reaching a set budget for each round (i.e. 24, 96, 384 sequences), a mutant batch is experimentally created and analyzed. The resulting sequence-function data are used to refine the surrogate model for guiding further engineering rounds (Figure 1B).

To facilitate fine-tuning of BO-EVO design specifications, we implemented FAST-HIT as an in-house, model-based optimization software framework (Figure 1B). FAST-HIT consists of four modules:

- **Landscape module** — provides sequence queries of different types, such as empirical landscapes like GB1 [35] and PhoQ [36], a mathematical landscape generated using the NK model [37], and an interactive landscape coupled with wet-lab measurement.
- **Encoder module** — transforms amino acid sequences into numerical encoding. In addition to the categorical Onehot encoding, several learned representations and one physiochemical encoding are also included.
- **Top model module** — executes training and evaluating diverse ML models. Here, the GPR model was primarily implemented with GPyTorch [39] for accelerated model training and evaluation.
- **Generator module** — proposes sequences through interaction with the other three modules.

*[Figure 1: Scheme of BO-EVO. (A) Conceptual diagram — candidate sequences (children) are generated by random mutation of a parent sequence, and one new sequence is proposed by BO on the candidates. (B) Software framework (FAST-HIT) consisting of the four modules described above.]*

### Moderate batch size is beneficial for low-round measurement

BO-EVO was first developed on an empirical fitness landscape of GB1 [35], which consists of the sequence-affinity data of a four-site SSM library (20⁴ = 160,000 possible sequences). To run a simulation, 20 starting sequences were selected to be diverse in Hamming distance to the global optima (fitness value set as 1.0) and with similar fitness levels as the wild-type (WT) sequence (fitness value ~0.1). For each of the 20 starting sequences, five random seeds were used to initiate BO-EVO, and hence in total 100 simulations were conducted when evaluating each algorithm design.

Five metrics were used to evaluate algorithm performance: for round-wise evaluation, the cumulative maximum fitness and the success ratio to reach the global optima after each iteration were used; for the outcome, the maximum and mean fitness values of all the proposed sequences were used; and for model performance evaluation, the Pearson correlation coefficient was used.

To optimize the selection strategy for batch size and iteration number, key determinants for experimental budgets, two types of BO-EVO simulations were performed using the simple Onehot encoding:

1. **Fixed iteration number (4)** — completing the protein engineering campaign in about 1 month. The maximum fitness increased with larger batch sizes, and the global optima could be consistently reached when batch sizes were ≥384, or total sample budgets were ≥1536 (Figure 2A). Interestingly, mean fitness of all proposed variants first increased then decreased with greater batch sizes, possibly due to the limited number of high-fitness variants in the sequence space (only 2.31% of variants exhibit fitness higher than WT [35]), as well as enhanced exploration in less-fitted regions with larger batches by BO-EVO.

2. **Fixed total sample budget (1536, i.e. 0.096% of the whole design space 160,000)** — the minimum sample number required to consistently reach the global optima within four BO-EVO iterations. The maximum and mean fitness of 1536 proposed mutants deteriorated with increasing batch sizes and decreasing iteration numbers (Figure 2B), indicating that algorithm-experiment feedback is necessary to improve BO-EVO performance.

Although favored in the second simulation, large iteration numbers are more time-consuming, as each experimental round generally takes a fixed turnabout of 5–7 days to execute. Overall, moderate iteration numbers (4–5) and batch sizes (384) were chosen for BO-EVO experiments.

*[Figure 2: BO-EVO simulation on different experimental batch sizes. (A) Fixed iteration number (4). (B) Fixed total sample budget (1536/160,000 = 0.096%). Maximum (red) and mean (green) fitness of proposed variants across 100 simulations, shown as violin/box plots. Gray dashed line marks WT fitness (0.1).]*

### Mutant data improve surrogate model performance

As elucidated in [22], informative encodings can improve MLDE outcomes by enhancing surrogate model quality. Because landscape information also enhanced surrogate model quality (Supplementary Figure S2A) and hence BO-EVO performance (Supplementary Figure S2C), the authors explored whether single-residue SSM data — easy to measure and small in scale (20 × N, e.g. N = 4 for GB1) — could be beneficial for model initialization ("warm start"). Relative to a cool-start strategy that fed only one sequence-fitness data point to surrogate models, the warm-start strategy resulted in a substantial gain in success ratio (~10%) throughout different rounds of simulated iterations (Figure 3). The warm-start strategy was therefore recommended for BO-EVO implementation.

*[Figure 3: A warm start enhanced BO-EVO success ratios. After the fifth round, BO-EVO reached success ratios of 78% and 68% for warm-start and cool-start strategies, respectively.]*

### Comparison of fitness exploration algorithms

BO-EVO was benchmarked against pure evolutionary (AdaLead [43]) and pure BO algorithms to examine the necessity of combining these two exploration strategies; random mutagenesis (Random) and a Metropolis-Hastings Markov chain Monte Carlo algorithm (MCMC) were also evaluated as baselines (see Supplementary Table S2 for setup details).

Random performed worst of the five algorithms, as expected (Figure 4). In terms of round-wise success ratios (Figure 4A) and maximum/mean fitness after five rounds (Figure 4B), BO-EVO outperformed MCMC and AdaLead by 11% and 21% increases in success ratio at Round 5, respectively. Utilizing the UCB acquisition function, BO-EVO was substantially better than AdaLead, which takes greedy exploitation on top variants with high predicted fitness. Ignoring model uncertainty, MCMC avoids local optima by accepting worse proposals probabilistically. These results demonstrate the importance of considering both variant fitness and model uncertainty when exploring a rugged fitness landscape.

On the other hand, by exhaustively exploring the whole design space (160,000 sequences) in each iteration, pure BO achieved better performance than BO-EVO (Figure 4), which only evaluated 3072 sequences (1.92% of the whole design space) *in silico* per round. Although BO-EVO's performance was not as good as pure BO, the computational time of BO-EVO was almost constant for exploring combinatorial mutagenesis landscapes, whereas the computational time of pure BO scaled exponentially with targeted residue numbers and rapidly became intractable (Supplementary Figure S3).

*[Figure 4: Fitness landscape exploration algorithms. (A) Success ratio reached until each round. (B) Maximum (top) and mean (bottom) fitness obtained in all proposed variants.]*

### The ruggedness of fitness landscapes challenges the performance of BO-EVO

The general applicability of BO-EVO was further assessed using two additional fitness landscapes: one empirical, one simulated. The empirical landscape reflects the protein–protein interaction (PPI) between *E. coli* protein kinase PhoQ mutants and its substrate PhoP [36]. GB1 and PhoQ are, to the authors' knowledge, the only two comprehensive datasets on four-site combinatorial libraries in the literature. The simulated fitness landscapes were generated based on the NK model [37], where N denotes protein sequence length and K denotes a tuneable degree of epistasis.

To better mimic empirical fitness landscapes, the original NK model [37] was reformulated in three steps:

1. The fitness distribution of empirical landscapes often obeys exponential distribution [3], so the "fitness table" of a select residue was sampled from an exponential rather than uniform distribution.
2. Only a small number of key residues exert substantial impacts on a target function, so a weighted average of "fitness contribution" was proposed, sampled from the Zipf distribution.
3. Because functional proteins constitute only a tiny fraction of the whole sequence space [2], low-fitness variants were set to "non-functional" with a cut-off threshold, creating a black hole-filled, exponentially distributed NK landscape (Supplementary Figure S4).

The final reformulated NK landscape was applied for BO-EVO simulation with N = 4 for four-site combinatorial landscapes and various degrees of epistasis (K = 0, 1, 2, 3).

Using the ruggedness-to-slope ratio [45] (r/s) to measure landscape ruggedness, NK landscape ruggedness scaled exponentially with K (Figure 5A). The ruggedness of GB1 and PhoQ landscapes was comparable to the K=0 and K=1 landscapes, respectively (Figure 5B), indicating PhoQ exhibits stronger epistasis than GB1. For the simulated NK landscape, the success ratio of reaching the global maximum decreased logarithmically with ruggedness (Figure 5B, green dots). Interestingly, although the empirical GB1 landscape and simulated K=0 landscape had similar epistasis levels, GB1 proved a more difficult task for BO-EVO than the simulated K=0 landscape. Furthermore, BO-EVO achieved similar success ratios for GB1 and PhoQ, despite PhoQ being more rugged than GB1. These results indicate BO-EVO generalizes reasonably well for fitness landscapes with moderate epistasis, with landscape ruggedness as a major challenge to locating global optima.

*[Figure 5: Generalization of BO-EVO to empirical and simulated landscapes with diverse ruggedness. (A) Ruggedness of NK landscapes. (B) Success ratio after the fifth BO-EVO round — NK (green dots), GB1 (purple star), PhoQ (blue triangle).]*

### Enzyme engineering via algorithm-experiment iterations by BO-EVO

BO-EVO was applied to a real-world protein engineering task via iterative feedback between ML models and robotic experiments. RhlA is a key enzyme responsible for synthesizing the lipid moieties of rhamnolipids (RLs), an important biosurfactant. The enzyme specificity of RhlA determines the chemical structures of the lipid moieties (Figure 6A), affecting the physiochemical and biological activities of the corresponding RL molecules. However, altering RhlA's enzyme specificity through (semi-)rational design or directed evolution has proven difficult [46, 47].

A robotic protocol was developed to build and test any 384 out of the 160,000 members in the four-residue combinatorial SSM libraries of RhlA (Supplementary Figure S6), with each experimental round taking less than 1 week to complete (Supplementary Figure S7). Enzyme specificity was measured using a previously reported MALDI-ToF mass spectrometry assay [28, 48] to quantify two RL products, Rha-(C8-C10) and Rha-C10-C10 (Figure 6B), in liquid cultures of recombinant *E. coli*. WT RhlA produced Rha-C10-C10 as the main product; the goal was to shift product specificity toward the smaller Rha-C18 product in RhlA mutants. Normalized production of Rha-(C8-C10) relative to WT was used as fitness (WT fitness = 1).

Arg74, Ala101, Leu148 and Ser173 (RALS) were selected as the four target residues for combinatorial mutagenesis, since many mutations at these ligand-binding residues substantially enhanced Rha-(C8-C10) production. Single-residue SSM data for these residues were used for a "warm start" (Round 1). Across BO-EVO iterations, the cumulative maximum fitness increased round-wise, reaching 7.35 (AACA sequence) at Round 4. Secondary confirmation via plasmid re-transformation and repeated MALDI-ToF MS analysis of the top 5 hits from Round 4 showed consistent improvement in Rha-(C8-C10) percentiles, though overall RL production levels varied substantially (Supplementary Table S4) — possibly due to limited quantitation capability of MALDI-ToF MS for mono-RLs [48]. Alternative high-throughput MS modalities (e.g. RapidFire ESI-MS) may help reduce this discrepancy in the future [49].

Ultimately, one RhlA mutant was identified conferring a **4.8-fold improvement** in RL-(C8-C10) production relative to WT (Supplementary Table S4 and Figure 6C).

For comparison, simulated engineering iterations also explored the GB1 and PhoQ landscapes: starting from WT sequences, 10 BO-EVO rounds were executed per simulation, repeated five times for each protein (Supplementary Figure S8). The round-wise pattern was similar between RhlA (Figure 6) and GB1 (Supplementary Figure S8A) — an initial increase followed by a decrease in median and top fitness, with GB1's turning point near the round that identified the global optima. PhoQ, by contrast, proved an easier target than RhlA and GB1, identifying the best-performing variant within two rounds in four of five simulations (Supplementary Figure S8B), with no obvious round-wise changes in median/overall fitness — likely due to differences in WT fitness values, fitness distributions, and landscape ruggedness between these landscapes (Supplementary Figures S4 and S5).

*[Figure 6: BO-EVO guided four-residue combinatorial SSM engineering of RhlA. (A) Molecular structure of mono-RL Rha-C10-C10. (B) MALDI mass spectra of mono-RLs produced from WT and representative RhlA mutants (residues 74, 101, 148, 173 labeled). (C) Normalized production of Rha-(C8-C10) relative to WT across rounds; cumulative maximum fitness labeled and plotted as a solid line (0→1.00, Round1→2.99, Round2→4.35, Round3→6.26, Round4→7.35, Round5→7.35).]*

---

## Discussion

In this work, a scalable and batched BO algorithm, BO-EVO, was developed for guiding multiple rounds of robotic experiments to explore combinatorial protein fitness landscapes. For four-site combinatorial mutagenesis, the total experimental budget of BO-EVO was reduced to 1536 mutants or less than 1% of the theoretical library size of 160,000, with four iterations of wet-lab measurement, model refinement and mutant design executed at a moderate batch size of 384. Within these sample budgets, BO-EVO achieved decent success rates (75%) to reach global optima on a range of simulated and empirical protein fitness landscapes. A real protein engineering task was performed to modify the product specificity of RhlA, with substantially improved mutants rapidly designed by BO-EVO within a month. To the authors' knowledge, this is the first report of algorithm-guided, automated protein engineering via multiple feedback rounds between ML models and robotic experiments.

Previously, a BO algorithm-guided pathway engineering platform, BioAutomata [32], robotically constructed and evaluated 136 pathway variants out of 13,824 possible library members, achieving a 1.77-fold higher lycopene titer. While both BO-EVO and BioAutomata share similar biofoundry setups, BO-EVO's problem dimension in this study is one order of magnitude higher, urging advances in both algorithms and robotics. For BO-EVO, an evolutionary algorithm was applied to restrict BO's search space for better computational efficiency; for BioAutomata, new robotic protocols for synthetic biology and MS analysis were developed to create and profile a batch of 384 enzyme mutants within a week for rapid experiment-model feedback.

Several adaptive experimental design approaches — MLDE [22], cluster learning-assisted directed evolution (CLADE) [50] and ODBO [51] — have also utilized ML models for sampling-efficient protein engineering. Two key differences distinguish BO-EVO:

1. **No prior knowledge requirement**: BO-EVO does not require prior knowledge of structure or homolog sequences of a target protein, which may be unavailable or scarce. By contrast, MLDE [22] used protein structure information (e.g. ΔΔG) and CLADE [50] used homolog sequences as zero-shot predictors to exclude low-fitness mutants.

2. **No full-space evaluation**: BO-EVO does not need to evaluate the whole design space *in silico*, restricting sequence design via evolutionary algorithms instead — only 3072 candidates (1.92% of all 160,000 possible sequences) were evaluated per iteration. In contrast, MLDE, CLADE and ODBO employed fitness models for *in silico* evaluation of all possible sequences; CLADE additionally applied deep hierarchical clustering of the whole design space, and ODBO used an outlier mining algorithm to identify low-fitness candidates. However, full-space prediction time increases exponentially and becomes intractable with high-dimensional problems (Supplementary Figure S3). With ever-decreasing experimental time and cost for robotic sequence-function profiling, the sampling efficiencies of both experiment and computation need to be carefully balanced.

### Future directions

- **Uncertainty quantification**: more accurate and computationally tractable models — ensemble methods, Bayesian neural networks [52], evidential deep learning [53] — could replace GPR.
- **Mutant residue selection**: the four target residues of RhlA were manually picked based on previous experiments; algorithm-guided residue selection [54] is desirable. MutCompute [55], for example, has trained a self-supervised CNN on protein structure data to identify hot spots in PETases for stabilizing mutations [56].
- **Experimental capabilities**: this study explored a four-site combinatorial library (theoretical size 160,000) at 384 mutants/week; further robotics improvements may enable exploration of more complex variant libraries. Comprehensive profiling of every RhlA mutant in this library, given enhanced throughput, could allow retrospective evaluation of BO-EVO performance.

Overall, a new approach for protein engineering via algorithm-guided robotic experiments for combinatorial mutagenesis was successfully developed.

### Key Points

- BO-EVO is developed as a scalable and batched exploration algorithm for protein engineering.
- BO-EVO guides efficient iterations between machine learning models and robotic experiments.
- A 4.8-fold improvement in enzyme specificity is achieved by wet-lab measurement of <1% of the whole design space.

---

## Supplementary Data

Supplementary data are available online at https://academic.oup.com/bib

## Acknowledgement

We thank Professor Lei Dai for the helpful discussion on NK landscapes. We thank the Shenzhen Synthetic Biology Infrastructure for robotic workflow development.

## Funding

National Key Research and Development Program of China (2021YFA0910800 and 2020YFA0908500); National Natural Science Foundation of China (32071428).

## Data availability

The software FAST-HIT is publicly available in a GitHub repository: https://github.com/hury07/fasthit
The code for running experiments and generating figures in the manuscript is publicly available: https://github.com/hury07/bo-evo_paper_data

---

## References

1. Romero PA, Arnold FH. Exploring protein fitness landscapes by directed evolution. *Nat Rev Mol Cell Biol* 2009;10:866–76.
2. Keefe AD, Szostak JW. Functional proteins from a random-sequence library. *Nature* 2001;410:715–8.
3. Orr HA. The distribution of fitness effects among beneficial mutations in Fisher's geometric model of adaptation. *J Theor Biol* 2006;238:279–85.
4. Nishikawa KK, Hoppe N, Smith R, et al. Epistasis shapes the fitness landscape of an allosteric specificity switch. *Nat Commun* 2021;12:5562.
5. Röthlisberger D, Khersonsky O, Wollacott AM, et al. Kemp elimination catalysts by computational enzyme design. *Nature* 2008;453:190–5.
6. Anishchenko I, Pellock SJ, Chidyausiku TM, et al. De novo protein design by deep network hallucination. *Nature* 2021;600:547–52.
7. Gribenko AV, Patel MM, Liu J, et al. Rational stabilization of enzymes by computational redesign of surface charge-charge interactions. *Proc Natl Acad Sci U S A* 2009;106:2601–6.
8. Contessoto VG, de Oliveira VM, Fernandes BR, et al. TKSA-MC: a web server for rational mutation through the optimization of protein charge interactions. *Proteins Struct Funct Bioinforma* 2018;86:1184–8.
9. Gopi S, Devanshu D, Krishna P, et al. PStab: prediction of stable mutants, unfolding curves, stability maps and protein electrostatic frustration. *Bioinformatics* 2018;34:875–7.
10. Arnold FH. Design by directed evolution. *Acc Chem Res* 1998;31:125–31.
11. Sinai S, Kelsic ED. A primer on model-guided exploration of fitness landscapes for biological sequence design. arXiv preprint arXiv:2010.10614. 2020.
12. Fang J. A critical review of five machine learning-based algorithms for predicting protein stability changes upon mutation. *Brief Bioinform* 2020;21:1285–92.
13. Marabotti A, Scafuri B, Facchiano A. Predicting the stability of mutant proteins by computational approaches: an overview. *Brief Bioinform* 2021;22:bbaa074.
14. Gelman S, Fahlberg SA, Heinzelman P, et al. Neural networks to learn protein sequence-function relationships from deep mutational scanning data. *Proc Natl Acad Sci USA* 2021;118:2104878118.
15. Hie B, Bryson BD, Berger B. Leveraging uncertainty in machine learning accelerates biological discovery and design. *Cell Syst* 2020;11:461–477.e9.
16. Wu Z, Kan SBJ, Lewis RD, et al. Machine learning-assisted directed protein evolution with combinatorial libraries. *Proc Natl Acad Sci USA* 2019;116:8852–8.
17. Bateman A, Martin MJ, Orchard S, et al. UniProt: The universal protein knowledgebase in 2021. *Nucleic Acids Res* 2021;49:D480–9.
18. Alley EC, Khimulya G, Biswas S, et al. Unified rational protein engineering with sequence-based deep representation learning. *Nat Methods* 2019;16:1315–22.
19. Rao R, Bhattacharya N, Thomas N, et al. Evaluating protein transfer learning with TAPE. *Advances in neural information processing systems* 2019;32:9689–9701.
20. Meier J, Rao R, Verkuil R, et al. Language models enable zero-shot prediction of the effects of mutations on protein function. bioRxiv 2021;2021.07.09.450648.
21. Elnaggar A, Heinzinger M, Dallago C, et al. ProtTrans: towards cracking the language of life's code through self-supervised deep learning and high performance computing. *IEEE Trans Pattern Anal Mach Intell* 2021;14:1–16.
22. Wittmann BJ, Yue Y, Arnold FH. Informed training set design enables efficient machine learning-assisted directed protein evolution. *Cell Syst* 2021;12:1026–1045.e7.
23. Luo Y, Jiang G, Yu T, et al. ECNet is an evolutionary context-integrated deep learning framework for protein engineering. *Nat Commun* 2021;12:1–14.
24. Biswas S, Khimulya G, Alley EC, et al. Low-N protein engineering with data-efficient deep learning. *Nat Methods* 2021;18:389–96.
25. Hillson N, Caddick M, Cai Y, et al. Building a global alliance of biofoundries. *Nat Commun* 2019;10:1038–41.
26. Chao R, Mishra S, Si T, et al. Engineering biological systems using automated biofoundries. *Metab Eng* 2017;42:98–108.
27. Zhang J, Chen Y, Fu L, et al. Accelerating strain engineering in biofuel research via build and test automation of synthetic biology. *Curr Opin Biotechnol* 2021;67:88–98.
28. Zhang S, Zhu J, Fan S, et al. Directed evolution of a cyclodipeptide synthase with new activities via label-free mass spectrometric screening. *Chem Sci* 2022;13:7581–6.
29. Dörr M, Fibinger MPC, Last D, et al. Fully automatized high-throughput enzyme library screening using a robotic platform. *Biotechnol Bioeng* 2016;113:1421–32.
30. Shahriari B, Swersky K, Wang Z, et al. Taking the human out of the loop: a review of Bayesian optimization. *Proc IEEE* 2016;104:148–75.
31. Greenhalgh JC, Fahlberg SA, Pfleger BF, et al. Machine learning-guided acyl-ACP reductase engineering for improved in vivo fatty alcohol production. *Nat Commun* 2021;12:1–10.
32. HamediRad M, Chao R, Weisberg S, et al. Towards a fully automated algorithm driven platform for biosystems design. *Nat Commun* 2019;10:1–10.
33. Radivojević T, Costello Z, Workman K, et al. A machine learning automated recommendation tool for synthetic biology. *Nat Commun* 2020;11:4879.
34. Zhang J, Petersen SD, Radivojevic T, et al. Combining mechanistic and machine learning models for predictive engineering and optimization of tryptophan metabolism. *Nat Commun* 2020;11:4880.
35. Wu NC, Dai L, Olson CA, et al. Adaptation in protein fitness landscapes is facilitated by indirect paths. *Elife* 2016;5:e16965.
36. Podgornaia AI, Laub MT. Pervasive degeneracy and epistasis in a protein-protein interface. *Science* 2015;347:673–7.
37. Kauffman SA, Weinberger ED. The NK model of rugged fitness landscapes and its application to maturation of the immune response. *J Theor Biol* 1989;141:211–45.
38. Rasmussen CE, Williams CKI. Gaussian processes for machine learning. *Adapt Comput Mach Learn* 2006;7:32–46.
39. Gardner JR, Pleiss G, Bindel D, et al. GPyTorch: blackbox matrix-matrix Gaussian process inference with GPU acceleration. *Adv Neural Inf Process Syst* 2018;31:7576–86.
40. Lai TL, Robbins H. Asymptotically efficient adaptive allocation rules. *Adv Appl Math* 1985;6:4–22.
41. Srinivas N, Krause A, Kakade SM, et al. Information-theoretic regret bounds for Gaussian process optimization in the bandit setting. *IEEE Trans Inf Theory* 2012;58:3250–65.
42. Belanger D, Vora S, Mariet Z, et al. Biological Sequence Design using Batched Bayesian Optimization. 2019;1–8.
43. Sinai S, Wang R, Whatley A, et al. AdaLead: a simple and robust adaptive greedy search algorithm for sequence design. arXiv preprint arXiv:2010.02141. 2020.
44. Szendro IG, Schenk MF, Franke J, et al. Quantitative analyses of empirical fitness landscapes. *J Stat Mech Theory Exp* 2013;2013:P01005.
45. Aita T, Iwakura M, Husimi Y. A cross-section of the fitness landscape of dihydrofolate reductase. *Protein Eng* 2001;14:633–8.
46. Han L, Liu P, Peng Y, et al. Engineering the biosynthesis of novel rhamnolipids in Escherichia coli for enhanced oil recovery. *J Appl Microbiol* 2014;117:139–50.
47. Dulcey CE, López de los Santos Y, Létourneau M, et al. Semi-rational evolution of the 3-(3-hydroxyalkanoyloxy)alkanoate (HAA) synthase RhlA to improve rhamnolipid production in Pseudomonas aeruginosa and Burkholderia glumae. *FEBS J* 2019;286:4036–59.
48. Si T, Li B, Comi TJ, et al. Profiling of microbial colonies for high-throughput engineering of multistep enzymatic reactions via optically guided matrix-assisted laser desorption/ionization mass spectrometry. *J Am Chem Soc* 2017;139:12466–73.
49. Fu L, Guo E, Zhang J, et al. Towards one sample per second for mass spectrometric screening of engineered microbial strains. *Curr Opin Biotechnol* 2022;76:102725.
50. Qiu Y, Hu J, Wei G-W. Cluster learning-assisted directed evolution. *Nat Comput Sci* 2021;1:809–18.
51. Cheng L, Yang Z, Liao B, et al. ODBO: Bayesian optimization with search space prescreening for directed protein evolution. arXiv preprint arXiv:2205.09548 2022;1–25.
52. Wang H, Yeung D-YY. A survey on Bayesian deep learning. *ACM Comput Surv* 2020;53:1–37.
53. Soleimany AP, Amini A, Goldman S, et al. Evidential deep learning for guided molecular property prediction and discovery. *ACS Cent Sci* 2021;7:1356–67.
54. Yu H, Ma S, Li Y, et al. Hot spots-making directed evolution easier. *Biotechnol Adv* 2022;56:107926.
55. Shroff R, Cole AW, Diaz DJ, et al. Discovery of novel gain-of-function mutations guided by structure-based deep learning. *ACS Synth Biol* 2020;9:2927–35.
56. Lu H, Diaz DJ, Czarnecki NJ, et al. Machine learning-aided engineering of hydrolases for PET depolymerization. *Nature* 2022;604:662–7.

---

## Author biographies

**Ruyun Hu** is an assistant professor at the Shenzhen Institute of Advanced Technology, Chinese Academy of Sciences. His research interest focuses on AI for Life Sciences, especially for protein mutation effects prediction and protein engineering.

**Lihao Fu** is a PhD student at the Shenzhen Institute of Advanced Technology, Chinese Academy of Sciences. His research interest focuses on biofoundry automation.

**Yongcan Chen** is an assistant professor at the Shenzhen Institute of Advanced Technology, Chinese Academy of Sciences. His research interests include protein mutation effects profiling and protein engineering.

**Junyu Chen** is a master student at the Shenzhen Institute of Advanced Technology, Chinese Academy of Sciences. His research interest focuses on protein mutation effects prediction.

**Yu Qiao** is a professor at the Shenzhen Institute of Advanced Technology, Chinese Academy of Sciences. His research interests include computer vision, deep learning and bioinformatics.

**Tong Si** is a professor at the Shenzhen Institute of Advanced Technology, Chinese Academy of Sciences. His research interests include synthetic biology and biofoundry automation.