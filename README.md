# SetariaGWASPipeline
Run MLMM GWAS on large Setaria genotype file

1. Initial filter and parsing of genotype file with vcftools
  * Do a simple inital filter by phenotyped lines, missingness, maf, etc
    * Do additional filtering for lines and MAF on the fly with R tools after splitting into batches
2. Additional filtering to keep only high-quality SNPs
3. LD prune high-quality SNP set for a set of non-redundant SNPs

## Scripts

All paths are relative to the ./src/ folder.

### 1.getSetariaHapmap.bash
* Quick script that:
  1. Downloads and extracts 12.Setaria_598g_8.58M_withRef_imp_phased_maf0.01_FINAL.vcf.bz2 from dropbox into ../data/genotype
  
### 2.setariaIRdataAnalysis.R
* This step is specific to the IR dataset and it generates a file to narrow down the genotype file to just phenotyped lines
 1. Contains code for QC analysis of traits
 2. Does a mixed model analysis to determine influence of different variables on traits.
 3. Calculates BLUEs and BLUPs of the traits based upon models selected by analysis in step 2.
   * BLUEs are calculated from this equation:
    `Trait ~  Genotype:Treatment + (1|Awning)`
   * BLUPs are calculated from this equation:
    `Trait ~  Treatment + (1|Awning) + (1|Genotype)`
 4. New data file is written to:
   * `../data/2.Setaria_IR_2016_datsetset_GWAS.BLUPsandBLUEs.csv`
 5. Trait heritabilities are written to:
   * `../results/2.setariaIR.H2.csv`
 6. Set of lines present in the phenotype file for genotype screening is in:
   * `../data/genotype/keepLines.txt`
 7. Diagnostic plots are written to:
   * `../results/2.setariaIR.BLUEsAndBlupsVsOriginal.pdf`
   
### 3.runVCFhapmapfilter.condor
* Uses vcftools to parse each vcf chromosome file:
  1. Remove minor alleles at 0.1 with --maf 0.1 --max-maf 0.9
  2. Make sure there are only (and at least) 2 alleles --min-alleles 2 --max-alleles 2
  3. Writes SNP file out in 012 formatted matrix with prefix `3.from12.setaria.maf0.1.maxMissing0.1`, rownames (individuals) are in the .012.indv and positions are in .012.pos
  
### 4.convertGenotypeFile.R
* Convert 012 VCF files to a R matrix in format for MLMM
  1. First it filters out SNPs with a large number of heterzygous calls (>25%)
    * This filters 147,824 SNPs and leaves 4,420,660 SNPs remaining on chromosomes 1 to 9
  2. It also calculates MAF and missing numbers, this file was already imputed by Sujan and the vcftools filtered for MAF, so no additional filtering is done
  3. The filtered genotype file and SNP info file is written to `../data/genotype/4.FilteredGenotypeFile.MatrixFormat.noscaffold.hetFilter0.25.maf0.1.rda`
  
### 5.filterGenoforLowQualitySNPs.R
* Filters out bad quality SNPs based on LD to neighboring SNPs
* After this step there are 1,997,780 SNPs remaining
 * The assumption of this step is that neighboring SNPs (those within 2000 base pairs) should be correlated (r<sup>2</sup> > 0.5).
   1. The first step of the filtering is to find a stretch of good SNPs defined as 2 or more SNPs in a row correlated at >0.9 r<sup>2</sup>
   2. If a SNP is within 2000 bp of a good set of SNPs and doesn't have a correlation of r<sup>2</sup>>0.5 with those SNPs it is considered a badly called SNP and removed
   3. If a SNP is not within 2000 bp of another SNP it is still checked for LD to next closest SNP, if <0.5 it is moved to a list of possible missed SNPs
   4. Diagnostic plots of before and after correlations to neighboring SNPs are writen to the `../results` folder.
   5. Output genotype file is written to `../data/genotype/5.filteredSNPs.2kbDistThresh.0.5neighborLD.rda`
    
### 6.performLDfilteringOfSNPs.R
* Filters out highly correlated SNPs to further narrow down genotype file and prevent redundant testing
* After this step there are 1,253,863 SNPs
* Final file has an average distance between SNPs of 314 bp and and average correlation between neighboring SNPs of 0.77.
 1. Tests pairs of SNPs, if correlation between two neighboring SNPs is r<sup>2</sup>>0.975 then only 1 SNP is kept
   * Process is iterated so no two nieghboring SNPs have a r<sup>2</sup> correlation >0.975.
 2. Diagnostic plots are written to `../results/6.FilteredHighSNP.ChromsomeWideNeighboringSNP.LD.subset.pdf`
 3. Output is written to `../data/genotype/6.filteredSNPs.noHighCorSNPs.2kbDistThresh.0.5neighborLD.0.975LDfilter.rda`
 