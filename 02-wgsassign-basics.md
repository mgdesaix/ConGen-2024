Population assignment with WGSassign
================
Matt DeSaix
2024-08-21

## Introduction to WGSassign

### Leave-one-out Assignment

Let’s start using `WGSassign` with a subset of American Redstart data (a
Neotropical migratory songbird), stored with the prefix
`amre.western.*`. This dataset has 47 individuals from two populations
on the Western portion of the breeding ground. These data are published
in <https://doi.org/10.1111/mec.17137>, and Figure 1 of the paper shows
the population structure. We have a `.beagle.gz` file and a
tab-delimited text file (`.wgsassign.ids.txt`) with two columns:
individual and population IDs (BR = Basin Rockies; WB = Western Boreal).
The pairwise Fst between these populations is around `0.006`…let’s see
how well we can assign these individuals back to their population!

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

    WGSassign
    Matt DeSaix.
    Using 1 thread(s).

    Parsing Beagle file.
    Loaded 99999 sites and 47 individuals.
    Parsing reference population ID file.
    EM (MAF) converged at iteration: 11
    EM (MAF) converged at iteration: 14
    Saved reference population allele frequencies as ./out/western/amre.western.reference.pop_af.npy (Binary - np.float32)

    Column order of populations is: ['BR' 'WB']
    Saved reference population names as ./out/western/amre.western.reference.pop_names.txt (String: Order of pops for .pop_af.npy, .ne_obs.npy, and fisher_obs.npy files)
    ...

The output is verbose and provides some handy information - such as
checking there’s the expected number of sites and individuals in the
file you provide. The `EM (MAF) converged at iteration` note comes up
each time the population allele frequency (minor) is calculated. It
occurs twice at the beginning (one for each pop) and also occurs during
LOO as the allele frequency for a given population has to be
recalculated without the individual that is being assigned…in our
example, this occured 47 times. This is the slowest step, but
fortunately this is fast in `WGSassign` and can be sped up by increasing
the number of threads used for parallelization (using the `--threads`
argument).

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

Well, mixed success. All individuals from Basin Rockies are accurately
assigned back to the population, but only 65% of Western Boreal
individuals are correctly assigned.

Let’s take a look at the effective sample sizes of these populations

``` r
western.ess <- read_table2("./out/western/amre.western.reference.ne_obs.txt")
```

``` r
western.ess
```

|       BR |       WB |
|---------:|---------:|
| 16.23298 | 10.35069 |

Basin Rockies has a notably higher effective sample size than Western
Boreal and this may be the root of the problem.

### Standard Assignment

The “standard” assignment test in `WGSassign` takes a Beagle file of
individuals and population allele frequencies of the populations in
question. Notably, it is assumed that the individuals in the Beagle file
were not used in the allele frequency calculations, in order to avoid
assignment bias…i.e. increased assignment accuracy will be perceived
when an individual’s genotype (likelihood) is included in the population
it is being assigned to.

Let’s look at this. Using the reference population allele frequencies
for these American Redstarts of the West, let’s determine the assignment
accuracy of these reference individuals with the standard test. Below
we’ll specify the Beagle file (same as before) and the allele frequency
file (created previously), along with the `--get_pop_like` argument to
get the likelihood of assignment from the standard assignment (in
contrast to `--loo` which we used before).

------------------------------------------------------------------------

Note: the allele frequencies are produced in a Numpy file. This is from
Python (which WGSassign is coded in) as a fast and easy to load file for
the calculations but if you want to look at it in your terminal you will
need to use Python. We won’t be getting into that here.
\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*\*

``` sh
# change file paths as needed
beagle=./data/amre.western.beagle.gz
af=./out/western/amre.western.reference.pop_af.npy
outdir=./out/western
out=${outdir}/amre.western.assign

WGSassign --beagle ${beagle} --pop_af_file ${af} --get_pop_like --out ${out}
```

    WGSassign
    Matt DeSaix.
    Using 1 thread(s).

    Parsing Beagle file.
    Loaded 99999 sites and 47 individuals.
    Parsing population allele frequency file.
    Calculating likelihood of population assignment
    47 individuals to assign to 2 populations
    Saved population assignment log likelihoods as ./out/western/amre.western.assign.pop_like.txt (text)

Super fast not having to recalculate allele frequencies for assigning
each individual! Why don’t we do this all the time for assessing
reference population assignment?? Let’s see…

Let’s load in the likelihoods of assignment and see what the assignment
accuracy looks like. I use the same R code as above, just pipe a little
more of it together:

``` r
western.assign <- cbind(western.ids, read_table2("./out/western/amre.western.assign.pop_like.txt",
                          col_names = western.pops))
```

``` r
head(western.assign)
```

| Sample | Breeding_pop |        BR |        WB |
|:-------|:-------------|----------:|----------:|
| Ind0   | WB           | -84764.77 | -79879.77 |
| Ind1   | WB           | -82439.19 | -76798.41 |
| Ind2   | WB           | -80902.37 | -75353.43 |
| Ind3   | WB           | -79201.48 | -72923.70 |
| Ind40  | BR           | -73404.52 | -79261.27 |
| Ind41  | BR           | -74071.25 | -79557.99 |

``` r
western.assign.summary <- western.assign %>%
  pivot_longer(cols = BR:WB,
               names_to = "AssignedPop",
               values_to = "AssignedLike") %>%
  group_by(Sample) %>%
  mutate(AssignedProb = round(exp(AssignedLike - max(AssignedLike)) / sum(exp(AssignedLike - max(AssignedLike))),2 )) %>%
  filter(AssignedLike == max(AssignedLike)) %>%
  ungroup() %>%
  mutate(Correct = ifelse(Breeding_pop == AssignedPop, 1, 0)) %>%
  group_by(Breeding_pop) %>%
  summarize(Accuracy = sum(Correct)/n(),
            .groups = "drop")
```

``` r
western.assign.summary
```

| Breeding_pop | Accuracy |
|:-------------|---------:|
| BR           |        1 |
| WB           |        1 |

Wow, 100% accuracy for both pops! In contrast, we had much lower
assignment accuracy before for the Western Boreal group. Clearly, we
upwardly bias the results when including the individual’s genotype
likelihood in the allele frequency calculation…and this is why we use
leave-one-out cross-validation for determining our assignment accuracy!
