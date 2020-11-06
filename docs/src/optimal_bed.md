# Optimal BED file

![CI](https://github.com/jordiabante/CpelTdm.jl/workflows/CI/badge.svg)
![Docs](https://github.com/jordiabante/CpelTdm.jl/workflows/Docs/badge.svg)
[![MIT license](https://img.shields.io/badge/license-MIT-green.svg)](https://github.com/jordiabante/CpelTdm.jl/blob/master/LICENSE.md)

In order to maximize the usage of RRBS data, it is recommendable to generate
a taylored BED file to pinpoint the regions with data to be analyzed by CpelTdm.
The following command allows to calculate the coverage of the i-th BAM file.

```bash
bedtools genomecov -ibam bam_file_i -bg > cov_i.bed
```

Regions in each `cov_i.bed` can be filtered on the basis of number of overlapping
reads. For instance, the following command allows to retain those regions
overlapping with more than 3 reads:

```bash
cat cov_i.bed | awk '{if($4>3){print $0}}' > cov_i.bed.tmp
mv cov_i.bed.tmp cov_i.bed
```

Next, regions that are within a given distance can be merged. For example,
the following command allows to merge regions that are within 50 bp:

```bash
bedtools merge -d 50 -i cov_i.bed > cov_i.bed.merged
mv cov_i.bed.merged cov_i.bed
```

Finally, one might want to intersect or take the union of the regions obtained
for multiple samples. The following commands create a file containing the union
of the regions across samples:

```bash
cat cov_{1..N}.bed > cov_union.bed.tmp
sort-bed cov_union.bed.tmp > cov_union.bed
bedtools merge -d 50 -i cov_union.bed > cov_union.merged.bed
rm cov_union.bed.tmp cov_union.bed
```
