Full WGSassign workflow
================
Matt DeSaix
2024-08-21

# WGSassign start to finish: Reference to unknown individuals

## Evaluating assignment uncertainty

In the previous section, we looked at the posterior probabilities of
assignment and they were all `1` for an assigned population and `0` for
the other population, even when the assignment was wrong. When using the
composite likelihoods from our assignment test in `WGSassign` across
thousands/millions of markers we’ve noticed they tend to underestimate
the uncertainty of assignment. One way we’ve developed to capture this
uncertainty, is to partition the Beagle file into `k` partitions of
equal numbers of markers (all with the same individuals) across the
genome and then performing assignment, to evaluate the proportion of
partitions that are assigned to each of the populations.

### Create Beagle partitions

We’ll try this out with a new data set of American Redstarts from the
eastern portion of the breeding range (`./data/eastern`). In this data
set we have three populations: MP = Maritime Provinces, NT = Northern
Temperate, ST = Southern Temperate. Again, details on the population
structure can be found in <https://doi.org/10.1111/mec.17137>. I provide
a script (`./data/scripts/partition-k-beagle.sh`) which creates $k$ new
beagles based on the following arguments:

`-k` Number of equal-sized partitions to make `-b` Input Beagle file
`-o` Prefix of the output $k$ number of Beagle files

We’ll create 5 subsets with the following code:

``` sh
k=5
beagle=./data/eastern/amre.eastern.beagle.gz
prefix=./data/eastern/amre.eastern.partition
./data/scripts/partition-k-beagle.sh -k ${k} -b ${beagle} -o ${prefix}
```

The original file was 99,999 SNPs - so not completely divisible by 5 for
each of the subsets, but we can see they are all roughly equal (the new
files are proper Beagle files with headers!) - 4 of the files have
20,000 markers and 1 has 19,999 which is what we want to see:

``` sh
for i in `ls ./data/eastern/*.partition.*`; do zcat $i | wc -l; done
## 20001
## 20000
## 20001
## 20001
## 20001
```

Now that we created these 5 different partitions, we’ll loop through
them and run leave-one-out assignment in WGSassign:

``` sh
# change file paths as needed
popIDs=./data/eastern/amre.eastern.wgsassign.ids.txt
outdir=./out/eastern
mkdir -p ${outdir}

for i in {1..5}
do
  beagle=./data/eastern/amre.eastern.partition.${i}.beagle.gz
  out=${outdir}/amre.eastern.partition.${i}
  WGSassign --beagle ${beagle} --pop_af_IDs ${popIDs} --get_reference_af --loo --out ${out}
done
```

This takes a minute or two…

### Evaluate assignment consistency of leave-one-out partitions

Now to evaluate likelihoods in R

``` r
library(tidyverse)
# number of partitions
k_partitions <- 5

# same IDs for each 
eastern.ids <- read_table2("./data/eastern/amre.eastern.wgsassign.ids.txt",
                          col_names = c("Sample", "Breeding_pop"))
loo.list <- list()
for(i in 1:k_partitions){
  prefix <- paste0("./out/eastern/amre.eastern.partition.", i)
  pop_names_file <- paste0(prefix, ".pop_names.txt")
  loo_file <- paste0(prefix, ".pop_like_LOO.txt")
  
  pop_names <- readLines(pop_names_file) # pop names should all be in the same order but just to be careful we'll read them in each time
  loo.list[[i]] <- eastern.ids %>%
    cbind(read_table2(file = loo_file,
                               col_names = pop_names)) %>%
    add_column("Partition" = i) # add partition signifier
}

loo.df <- do.call("rbind", loo.list)
```

``` r
head(loo.df)
```

| Sample | Breeding_pop |        MP |        NT |        ST | Partition |
|:-------|:-------------|----------:|----------:|----------:|----------:|
| Ind4   | ST           | -18063.29 | -17842.81 | -17810.49 |         1 |
| Ind5   | ST           | -16109.58 | -15852.50 | -15669.88 |         1 |
| Ind6   | ST           | -15807.76 | -15587.83 | -15349.10 |         1 |
| Ind7   | ST           | -16341.76 | -16160.37 | -15893.84 |         1 |
| Ind8   | ST           | -16416.22 | -16078.91 | -15943.50 |         1 |
| Ind9   | ST           | -17293.37 | -17069.62 | -16910.92 |         1 |

We’ll calculate posterior probabilities same as before, just need to
group the data by Sample *and* partition this time.

``` r
eastern.k.loo.probs <- loo.df %>%
  pivot_longer(cols = MP:ST, # pivot by the assigned populations
               names_to = "AssignedPop",
               values_to = "AssignedLike") %>%
  group_by(Sample, Partition) %>% # group by sample and partition here!
  mutate(AssignedProb = round(exp(AssignedLike - max(AssignedLike)) / sum(exp(AssignedLike - max(AssignedLike))),2 ))
```

``` r
head(eastern.k.loo.probs)
```

| Sample | Breeding_pop | Partition | AssignedPop | AssignedLike | AssignedProb |
|:-------|:-------------|----------:|:------------|-------------:|-------------:|
| Ind4   | ST           |         1 | MP          |    -18063.29 |            0 |
| Ind4   | ST           |         1 | NT          |    -17842.81 |            0 |
| Ind4   | ST           |         1 | ST          |    -17810.49 |            1 |
| Ind5   | ST           |         1 | MP          |    -16109.58 |            0 |
| Ind5   | ST           |         1 | NT          |    -15852.50 |            0 |
| Ind5   | ST           |         1 | ST          |    -15669.88 |            1 |

Now we can summarize assignment consistency - i.e. the proportion of the
genomic partitions assigned to each population:

``` r
k_partitions <- 5
eastern.k.consistency <- eastern.k.loo.probs %>%
  group_by(Sample, Partition) %>%
  filter(AssignedLike == max(AssignedLike)) %>%
  group_by(Sample, Breeding_pop, AssignedPop) %>%
  summarize(N = n(),
            .groups = "drop") %>%
  pivot_wider(names_from = AssignedPop,
              values_from = N) %>%
  mutate(across(where(is.numeric), ~replace_na(.x, 0)/k_partitions))
```

``` r
head(eastern.k.consistency, n = 8)
```

| Sample | Breeding_pop |  ST |  NT |  MP |
|:-------|:-------------|----:|----:|----:|
| Ind10  | ST           | 1.0 | 0.0 |   0 |
| Ind11  | ST           | 1.0 | 0.0 |   0 |
| Ind12  | ST           | 1.0 | 0.0 |   0 |
| Ind13  | ST           | 1.0 | 0.0 |   0 |
| Ind14  | ST           | 1.0 | 0.0 |   0 |
| Ind15  | ST           | 1.0 | 0.0 |   0 |
| Ind16  | ST           | 0.8 | 0.2 |   0 |
| Ind17  | ST           | 0.6 | 0.4 |   0 |

We can see that the partitions in some individuals are not getting
consistently assigned to the same population. For example, Ind17 has 3
partitions assigned to ST (0.6) and 2 partitions assigned to NT (0.4).

As in the previous section with the Western individuals, we can check
accuracy of assignment by breeding population. For now, we’ll treat the
most consistently assigned population as *the* assigned population

``` r
eastern.k.consistency.max <- eastern.k.consistency %>%
  pivot_longer(cols = ST:MP, names_to = "AssignedPop", values_to = "Consistency") %>%
  group_by(Sample) %>%
  filter(Consistency == max(Consistency)) # filter to most consistent population

eastern.k.acc <- eastern.k.consistency.max %>%
  mutate(Correct = ifelse(Breeding_pop == AssignedPop, 1, 0)) %>% # create binary column of correct assignment
  group_by(Breeding_pop) %>%
  summarize(Accuracy = sum(Correct)/n())
```

``` r
eastern.k.acc
```

| Breeding_pop |  Accuracy |
|:-------------|----------:|
| MP           | 0.6785714 |
| NT           | 0.9230769 |
| ST           | 0.9787234 |

Other useful ways of summarizing these results is with a confusion
matrix, with the first column listing the breeding population the
individuals were sampled from, and the other columns providing the
counts of individuals assigned to each population

``` r
eastern.k.conf <- eastern.k.consistency.max %>%
  group_by(Breeding_pop, AssignedPop) %>%
  summarize(N = n(), .groups = "drop") %>%
  pivot_wider(names_from = AssignedPop, values_from = N) %>%
  mutate(across(where(is.numeric), ~replace_na(.x, 0))) %>%
  select(Breeding_pop, MP, NT, ST)
```

``` r
eastern.k.conf
```

| Breeding_pop |  MP |  NT |  ST |
|:-------------|----:|----:|----:|
| MP           |  19 |   9 |   0 |
| NT           |   0 |  36 |   3 |
| ST           |   0 |   1 |  46 |

Still, we haven’t really examined the consistency results. Say we only
want to consider “successfully assigned individuals” as having all 5 out
of 5 genomic partitions assigned to the correct population - we can also
add a category in the confusion matrix for individuals with an
assignment consistency of \< 1.0

``` r
threshold <- 1
eastern.k.conf2 <- eastern.k.consistency.max %>%
  mutate(AssignedPop = ifelse(Consistency < threshold, "Uncertain", AssignedPop)) %>% # filter by consistency threshold
  group_by(Breeding_pop, AssignedPop) %>%
  summarize(N = n(), .groups = "drop") %>%
  pivot_wider(names_from = AssignedPop, values_from = N) %>%
  mutate(across(where(is.numeric), ~replace_na(.x, 0))) %>%
  select(Breeding_pop, MP, NT, ST, Uncertain)
eastern.k.conf2
```

    ## # A tibble: 3 × 5
    ##   Breeding_pop    MP    NT    ST Uncertain
    ##   <chr>        <int> <int> <int>     <int>
    ## 1 MP              10     1     0        17
    ## 2 NT               0    29     1         9
    ## 3 ST               0     1    44         2

# Break-out session (30 min.)

For this second break-out session, you’ll work together in groups to
evaluate assignment accuracy of the eastern data set of American
Redstarts. You will then use the standard assignment test on a Beagle
file of individuals sampled from the nonbreeding range to determine
their populations of origin (i.e. their migratory connection)

1.  Evaluate assignment accuracy of the American Redstart eastern
    breeding range individuals (i.e. follow steps 1-6 of the last
    exercise for this data set)…*but now make sure to consider
    assignment consistency*. Select individuals to create balanced
    effective samples for the three populations (“LOO reference
    individuals”). Create a new `*.wgsassign.ids.txt` file and use
    `subset-beagle-ind.sh` to create the new Beagle file. *Break this
    file into $k$ partitions of your choice* as we did at the start of
    this session. Then perform LOO assignment on the partitions and
    evaluate assignment accuracy of this “reference set” as well as for
    the remaining breeding individuals.

2.  Determine the eastern breeding populations of origin of nonbreeding
    sampled individuals with the standard assignment test in
    `WGSassign`. Characterize the assignment of the nonbreeding
    individuals with assignment consistency as well (i.e. partition the
    nonbreeding individuals genome into $k$ partitions and do $k$
    assignments). Summarize the migratory assignments by nonbreeding
    site, which are coded by country in the 2nd column of the
    `./data/amre.nonbreeding.ids.sites.txt` file:

CO = Colombia CR1 = Costa Rica site 1 CR2 = Costa Rica site 2 CU1 = Cuba
JM1 = Jamaica TR = Trinidad and Tobago

The files you’ll need for this are:

- `./data/eastern/amre.eastern.beagle.gz`
- `./data/eastern/amre.eastern.wgsassign.ids.txt`
- `./data/unknown/amre.nonbreeding.beagle.gz`
- `./data/unknown/amre.nonbreeding.ids.sites.txt`
- `./data/scripts/subset-beagle-ind.sh`
- `./data/scripts/partition-k-beagle.sh`
- `./data/full.amre.beagle.meta.csv`

\`\`\`

We’ll reconvene shortly and discuss these results as a group!

------------------------------------------------------------------------

These wraps up the session!
