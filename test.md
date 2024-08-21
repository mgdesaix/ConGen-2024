test
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

### What is a Beagle file?

``` sh
head ./data/tiny.beagle.txt
```

    ## marker   allele1 allele2 SERU2001    SERU2001    SERU2001    SERU2003    SERU2003    SERU2003
    ## scaffold1|size5275185_1104   0   2   0.000000    0.999788    0.000212    0.941168    0.058832    0.000000
    ## scaffold1|size5275185_4870   1   3   0.333333    0.333333    0.333333    0.799889    0.200111    0.000000
    ## scaffold1|size5275185_10574  1   2   0.666649    0.333333    0.000018    0.333333    0.333333    0.333333
    ## scaffold1|size5275185_11007  0   2   0.992246    0.007754    0.000000    0.996108    0.003892    0.000000
    ## scaffold1|size5275185_11167  0   2   0.000058    0.999942    0.000000    0.633983    0.366017    0.000000
