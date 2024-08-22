Genotype likelihood data
================
Matt DeSaix
2024-08-21

One aspect that is different when working with low-coverage data is the
file types of the variant information. Typically you don’t work with the
standard variant call format (VCF) files that are commonplace with
called genotypes. Unfortunately, what this means, is you don’t have
access to all the handy VCF handling software (ex. bcftools) that you
may be used to. Fortunately,
[ANGSD](http://www.popgen.dk/angsd/index.php/ANGSD) is a software that
was specifically developed to analyze low-coverage sequencing data, and
notably has many methods designed for accounting for genotype
uncertainty. Here, we’ll be making use of Beagle files produced in
ANGSD.

## What is a Beagle file?

<figure id="id">
<img src="./images/beagle_puppy_small.jpg" class="class" width="300"
alt="Sadly, not the beagles we’re talking about" />
<figcaption aria-hidden="true">Sadly, not the beagles we’re talking
about</figcaption>
</figure>

Beagle files provide genotype likelihoods for each marker (row) and
individual (columns). The layout is as follows:

- Column 1 (marker) = chromosome and position
- Column 2 (allele1) = the major allele coded as an integer (0=A, 1=C,
  2=G, 3=T)
- Column 3 (allele2) = the minor allele coded as an integer

After that, there are 3 columns for each individual (with the header
being the individual ID) with the genotype likelihoods in the order of:

- major/major genotype
- major/minor genotype
- minor/minor genotype

Let’s use the Terminal to look at an example **tiny** Beagle file in the
data directory that has genotype likelihoods for two individuals across
five markers:

``` sh
cat ./data/tiny.beagle.txt
```

    ## marker   allele1 allele2 Ind0    Ind0    Ind0    Ind1    Ind1    Ind1
    ## scaffold1|size5275185_1104   0   2   0.000000    0.999788    0.000212    0.941168    0.058832    0.000000
    ## scaffold1|size5275185_4870   1   3   0.333333    0.333333    0.333333    0.799889    0.200111    0.000000
    ## scaffold1|size5275185_10574  1   2   0.666649    0.333333    0.000018    0.333333    0.333333    0.333333
    ## scaffold1|size5275185_11007  0   2   0.992246    0.007754    0.000000    0.996108    0.003892    0.000000
    ## scaffold1|size5275185_11167  0   2   0.000058    0.999942    0.000000    0.633983    0.366017    0.000000

Note that the genotype likelihoods are scaled to sum to 1 for each
individual, let’s look at just these values in the first individual:

``` sh
cut -f4-6 ./data/tiny.beagle.txt
```

    ## Ind0 Ind0    Ind0
    ## 0.000000 0.999788    0.000212
    ## 0.333333 0.333333    0.333333
    ## 0.666649 0.333333    0.000018
    ## 0.992246 0.007754    0.000000
    ## 0.000058 0.999942    0.000000

Standard naming of individuals is *Ind0* to Ind “n-1” for $n$
individuals.

------------------------------------------------------------------------

*Question 1.1:* For Ind0, what does it mean that the 2nd marker has a
genotype likelihood of $0.333$ for all 3 genotypes?

------------------------------------------------------------------------

## Working with Beagle files

Beagle files are typically very large and cumbersome to work with, and
as far as I know there aren’t common tools for manipulating them - such
as subsetting the files to specific markers or individuals. I typically
do this with command line tools such as `awk`, and when I do I always
double-check that the produced file is a proper Beagle file. Certain
things I check are:

- has a header!
- line count! A simple check with `wc -l` should provide $k+1$ for a
  file with $k$ markers
- all 3 genotype likelihoods for each individual. At the very least, I
  check the start of the header has the 3 initial columns and the 3
  replicate names of the first individual, as well as checking the end
  of the header has the final individual’s header 3 times.

Also, because Beagle files are large they are typically gzipped. Let’s
check out another Beagle file in the `./data/` directory

``` sh
# zcat ./data/amre.western.beagle.gz | wc -l
# Macs don't 'zcat' the same as on Linux, on a Mac I use the command below
zcat < ./data/amre.western.beagle.gz | wc -l
```

    ##   100000

------------------------------------------------------------------------

*Question 1.2:* How many markers are in the `amre.western.beagle.gz`
file?

*Bonus Question 1.3:* How many individuals are in the file?

------------------------------------------------------------------------

Now on to the next section: [Performing population assignment with
WGSassign](./02-wgsassign-basics.md)

### Related resources

For a great tutorial on understanding different analyses and the process
of getting genotype likelihood files from ANGSD, check out
[lcwgs-guide-tutorial](https://github.com/nt246/lcwgs-guide-tutorial) on
github, and their paper:
<https://onlinelibrary.wiley.com/doi/10.1111/mec.16077>.
