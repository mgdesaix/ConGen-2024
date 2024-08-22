Population assignment with WGSassign
================
Matt DeSaix
2024-08-21

### Introduction to WGSassign

Let’s start using `WGSassign` with a subset of American Redstart data (a
Neotropical migratory songbird), stored with the prefix
`amre.western.*`. This dataset has 47 individuals from two populations
on the breeding ground. These data are published in
<https://doi.org/10.1111/mec.17137>, and Figure 1 of the paper shows the
population structure. We have a `.beagle.gz` file and a tab-delimited
text file (`.wgsassign.ids.txt`) with two columns: individual and
population IDs (BR = Basin Rockies; WB = Western Boreal). The pairwise
Fst between these populations is around `0.006`…let’s see how well we
can assign these individuals back to their population!

The two files described above are the only files needed to initially
perform leave-one-out cross-validation assignment and calculate the
effective sample sizes (of the populations and individuals).
**Importantly, the ID text file for WGSassign must have individuals in
the same order as in the Beagle file!**

Let’s run our first command in `WGSassign`. We’ll load a Conda
environment (provided with the installation of WGSassign containing the
dependencies) and create an out directory. We’ll specify options below
to calculate reference allele frequencies (`--get_reference_af`),
calculate effective sample sizes (`--ne_obs`), and perform leave-one-out
cross-validation (`--loo`):

``` sh
conda activate WGSassign
# change file paths as needed
beagle=./data/amre.western.beagle.gz
popIDs=./data/amre.western.wgsassign.ids.txt
outdir=./out/western
mkdir -p ${outdir}
out=${outdir}/amre.western.reference

WGSassign --beagle ${beagle} --pop_af_IDs ${popIDs} --get_reference_af --ne_obs --loo --out ${out}
```

``` sh
## WGSassign
## Matt DeSaix.
## Using 1 thread(s).
## 
## Parsing Beagle file.
## Loaded 99999 sites and 47 individuals.
## Parsing reference population ID file.
## EM (MAF) converged at iteration: 11
## EM (MAF) converged at iteration: 14
## Saved reference population allele frequencies as ./out/western/amre.western.reference.pop_af.npy (Binary - np.float32)
## 
## Column order of populations is: ['BR' 'WB']
## Saved reference population names as ./out/western/amre.western.reference.pop_names.txt (String: Order of pops for .pop_af.npy, .ne_obs.npy, and fisher_obs.npy files)
...
...
```

The output is verbose and provides some handy information - such as
checking there’s the expected number of sites and individuals in the
file you provide.

Let’s switch into R and look at the file with the log likelihoods of
assignment (`*.pop_like_LOO.txt`):

``` r
library(tidyverse)
# read in initial IDs file for individuals
western.ids <- read_table2("./data/amre.western.wgsassign.ids.txt",
                          col_names = c("Sample", "Breeding_pop"))
# read in wgsassign output
# population names
western.pops <- readLines("./out/western/amre.western.reference.pop_names.txt")
# LOO assignment log likelihoods
western.loo <- cbind(western.ids, read_table2("./out/western/amre.western.reference.pop_like_LOO.txt",
                          col_names = western.pops))
```

``` r
head(western.loo)
```

| Sample | Breeding_pop |        BR |        WB |
|:-------|:-------------|----------:|----------:|
| Ind0   | WB           | -84764.77 | -84599.44 |
| Ind1   | WB           | -82439.19 | -81845.05 |
| Ind2   | WB           | -80902.37 | -80748.42 |
| Ind3   | WB           | -79201.48 | -78675.23 |
| Ind40  | BR           | -77307.96 | -79393.15 |
| Ind41  | BR           | -78018.68 | -79657.28 |

Note that the `*.pop_like_LOO.txt` does not have the sample names, but
is in the same order as the initial file we provided to `WGSassign` so
we `cbind()` em.

We then want to convert the log likelihoods to posterior probabilities
by exponentiating them and dividing by the sum. To deal with underflow
issues of values getting too small for the computer when we exponentiate
them, we’ll also subtract the max value.

Many ways to do this, but below I pivot the table so I can easily group
by the Sample IDs:

``` r
western.loo.probs <- western.loo %>%
  pivot_longer(cols = BR:WB, # pivot by the assigned populations
               names_to = "AssignedPop",
               values_to = "AssignedLike") %>%
  group_by(Sample) %>%
  mutate(AssignedProb = round(exp(AssignedLike - max(AssignedLike)) / sum(exp(AssignedLike - max(AssignedLike))),2 )) # summarize posterior probabilities, grouped by individual
```

``` r
head(western.loo.probs)
```

| Sample | Breeding_pop | AssignedPop | AssignedLike | AssignedProb |
|:-------|:-------------|:------------|-------------:|-------------:|
| Ind0   | WB           | BR          |    -84764.77 |            0 |
| Ind0   | WB           | WB          |    -84599.44 |            1 |
| Ind1   | WB           | BR          |    -82439.19 |            0 |
| Ind1   | WB           | WB          |    -81845.05 |            1 |
| Ind2   | WB           | BR          |    -80902.37 |            0 |
| Ind2   | WB           | WB          |    -80748.42 |            1 |

Ideally, we’d capture some uncertainty in the assignments with the
probabilities - but it is interesting to note that they are all 1 and 0,
even when wrong! We’ll get into this later. For now, we just want to
filter out the population the individual is assigned to (`AssignedPop`)
and determine the accuracy of assignment back to their known population
(`Breeding_pop`).

``` r
western.loo.summary <- western.loo.probs %>%
  filter(AssignedLike == max(AssignedLike)) %>% # keep only the maximum likelihood
  ungroup() %>%
  mutate(Correct = ifelse(Breeding_pop == AssignedPop, 1, 0)) %>% # create binary column of correct assignment
  group_by(Breeding_pop) %>%
  summarize(Accuracy = sum(Correct)/n()) # summarize accuracy by population
```

``` r
western.loo.summary
```

| Breeding_pop | Accuracy |
|:-------------|---------:|
| BR           |     1.00 |
| WB           |     0.65 |
