# ssGSEA3.0

An installable R package for gene-centric **single sample Gene Set Enrichment
Analysis (ssGSEA)** of expression data (e.g. mRNAs, proteins) and site-centric
**PTM Signature Enrichment Analysis (PTM-SEA)** [1] of phosphoproteomics data
using the [**PTM signatures database (PTMsigDB)**](https://proteogenomics.shinyapps.io/ptmsigdb/) [1].

`ssGSEA3.0` replaces `ssGSEA2.0`'s gene-label permutation testing with
feature-level (gene- or site-level) bootstrap resampling for significance.
This fixes a known weight-dependent calibration issue in the original approach
while keeping more detection power than a direct fix to the permutation -- see
`vignette("ssgsea2-vs-ssgsea3-comparison")` for the full story and a
comparison of `ssGSEA2.0`, its fixed variant, and `ssGSEA3.0` on real data.


## Official ssGSEA codebase

`ssGSEA3.0` is an update to `ssGSEA2.0`, developed alongside our [manuscript](https://www.mcponline.org/content/18/3/576.abstract) describing PTM-SEA and [PTMsigDB](https://proteogenomics.shinyapps.io/ptmsigdb/). It is **not the official repository for ssGSEA**. For the official ssGSEA codebase, see [here](https://github.com/GSEA-MSigDB/ssGSEA-gpmodule); for the official GenePattern module, see [here](https://www.genepattern.org/modules/docs/ssGSEAProjection/4).


## Installation

```r
# install.packages("remotes")
remotes::install_github("broadinstitute/ssGSEA3.0", build_vignettes = TRUE)
```


## Quick start

`ssGSEA3.0` covers the full pipeline from a raw GCT file to a summary heatmap:

```r
library(ssGSEA3)

# 1. (optional) preprocess a raw GCT into a gene-centric or site-centric report
# preprocessGCT(...)

# 2. score samples against a gene set database, with bootstrap significance
result <- ssGSEA3(
  input.ds = "my_data.gct",
  output.prefix = "my_run",
  gene.set.databases = "my_gene_sets.gmt"
)

# 3. summarize as a heatmap
ht <- plot_ssgsea_heatmap(result)
ComplexHeatmap::draw(ht)
```

See `vignette("ssgsea2-vs-ssgsea3-comparison")` for three full worked examples
(ordinary ssGSEA on tumor NMF factors against HALLMARK and REACTOME, and
directional PTM-SEA on phosphoproteomics data) using the datasets installed
with the package in `example/`. They show the exact results reported in the
ssGSEA3.0 manuscript; re-running those analyses with the current package gives slightly
different numbers due to bootstrap randomness alone (see the vignette for
why).


## ssGSEA

This is an updated version of the original ssGSEA [2,3] R-implementation. Depending on the input dataset and chosen database (gene sets or PTM signatures), the software performs either ssGSEA or PTM-SEA, respectively. The Molecular Signatures Database ([MSigDB](http://software.broadinstitute.org/gsea/msigdb/)) [4] provides a large collection of curated gene sets. Gene sets are stored as plain text in  [GMT](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GMT:_Gene_Matrix_Transposed_file_format_.28.2A.gmt.29) format. A current version of MSigDB gene set collections can be found in the ```db/msigdb``` subfolder of this repository (not installed with the package). MSigDB gene sets are realeased under [Creative Commons Attribution 4.0 International License](http://software.broadinstitute.org/gsea/msigdb_license_terms.jsp). The license terms can be found in the ```db/msigdb``` folder.

File formats supported are Gene Cluster Text [GCT v1.2](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GCT:_Gene_Cluster_Text_file_format_.28.2A.gct.29) or [GCT v1.3](https://clue.io/connectopedia/gct_format) files. [Morpheus](https://software.broadinstitute.org/morpheus/) provides a convenient way to convert your data tables into GCT format.

For more information about the GSEA method and MSigDB please visit http://software.broadinstitute.org/gsea/.


## PTMsigDB v2.0.0

Please check out our new website for [PTMsigDB](https://proteogenomics.shinyapps.io/ptmsigdb/). We have updated [PTMsigDB](https://proteogenomics.shinyapps.io/ptmsigdb/) to version v2.0.0 in which we provide better and more consistent annotation of each PTM site. We have also inlcuded a *disease* category comprising of signatures associated to certain diseases curated from the table ```Disease-associated_sites``` available at [PhosphoSitePlus (PSP)](https://www.phosphosite.org/homeAction.action) [5].

The PTM signatures database (PTMsigDB) is a database comprised of modification site-specific signatures of perturbations, kinase activities and signaling pathways curated from more than 2,500 publications which provides the foundation to perform PTM-SEA. A unique advantage of PTMsigDB over other pathway databases is the annotation of each PTM site with its reported direction of change upon a specific perturbation or signaling event which is incorporated into the scoring scheme of PTM-SEA. The foundation of PTMsigDB is [PhosphoSitePlus (PSP)](https://www.phosphosite.org/homeAction.action) [5], a comprehensive systems biology resource for PTMs, which provides high-quality curation and annotation of PTMs at the individual residue level. A collection of PTM sites, whose levels are collectively regulated in a curated pathway or upon a perturbation, are defined as a signature set. Signature sets in PTMsigDB can be separated into different categories: 1) **Perturbation signatures** derived from treatment of cells with perturbagens such as small molecules or growth factors; 2) **Signature sets of molecular signaling pathways**; 3) **Kinase-substrate signatures**; and 4) **Disease-associated signature sets**.

To ensure a high degree of compatibility to phosphorylation datasets generated by different software packages and searched against different protein sequence databases, PTMsigDB represents signatures using three different identifiers to represent phosphorylation sites: 1) **PSP site group ID**; 2) **UniProt-centric ID**; 3) **Flanking sequence** (Table 1). While the PSP site group ID provides an unambiguous representation of PTM sites within protein families and across species [5], using this type of identifier restricts the analysis to PTM sites present in PSP. We generally recommend to using the flanking sequence as site identifier, since these are more invariant to updates made to protein sequence databases.

| Database format   | Site accession | Example in PTMsigDB |  Example in dataset | Download
| ----------------- | -------------- | ------------------- | ------------------- | ------------ 
| UniProt-centric   | Uniprot_acc;site-type;direction | Q06609;Y315-p;u | Q06609;Y315-p | [human](https://raw.githubusercontent.com/broadinstitute/ssGSEA2.0/master/db/ptmsigdb/v2.0.0/all/ptm.sig.db.all.uniprot.human.v2.0.0.gmt)<br>[mouse](https://raw.githubusercontent.com/broadinstitute/ssGSEA2.0/master/db/ptmsigdb/v2.0.0/all/ptm.sig.db.all.uniprot.mouse.v2.0.0.gmt)<br>[rat](https://raw.githubusercontent.com/broadinstitute/ssGSEA2.0/master/db/ptmsigdb/v2.0.0/all/ptm.sig.db.all.uniprot.rat.v2.0.0.gmt)
| Flanking sequence | +/-7aa flanking seq-type;direction | ETRICKIYDSPCLPE-p;u | ETRICKIYDSPCLPE-p | [human](https://raw.githubusercontent.com/broadinstitute/ssGSEA2.0/master/db/ptmsigdb/v2.0.0/all/ptm.sig.db.all.flanking.human.v2.0.0.gmt)<br>[mouse](https://raw.githubusercontent.com/broadinstitute/ssGSEA2.0/master/db/ptmsigdb/v2.0.0/all/ptm.sig.db.all.flanking.mouse.v2.0.0.gmt)<br>[rat](https://raw.githubusercontent.com/broadinstitute/ssGSEA2.0/master/db/ptmsigdb/v2.0.0/all/ptm.sig.db.all.flanking.rat.v2.0.0.gmt)
| PSP site group id | site_grp_id-type;direction | 448324-p;u | 448324-p | [human](https://raw.githubusercontent.com/broadinstitute/ssGSEA2.0/master/db/ptmsigdb/v2.0.0/all/ptm.sig.db.all.sitegrpid.human.v2.0.0.gmt)<br>[mouse](https://raw.githubusercontent.com/broadinstitute/ssGSEA2.0/master/db/ptmsigdb/v2.0.0/all/ptm.sig.db.all.sitegrpid.mouse.v2.0.0.gmt)<br>[rat](https://raw.githubusercontent.com/broadinstitute/ssGSEA2.0/master/db/ptmsigdb/v2.0.0/all/ptm.sig.db.all.sitegrpid.rat.v2.0.0.gmt)

Table 1: PTM site representation in PTMsigDB. The direction of change for a PTM site in a signature is indicated by ;u (up-regulation) or ;d (down-regulation). Please note that the annotation of directionality is a feature of PTMsigDB (column: *Example in PTMsigDB*) and must not be included when generating compatible site identifier for a particular dataset (column: *Example in dataset*).


## PTM-SEA

PTM-Signature Enrichment Analysis (PTM-SEA) is a modified version of ssGSEA to perform site-specific signature analysis by scoring PTMsigDB's bi-directional signature-sets. The input to PTM-SEA is a single site-centric data matrix, *m*, stored in [GCT v1.2](https://software.broadinstitute.org/cancer/software/gsea/wiki/index.php/Data_formats#GCT:_Gene_Cluster_Text_file_format_.28.2A.gct.29) or [GCT v1.3](https://clue.io/connectopedia/gct_format)  format and PTM signatures database (PTMsigDB). Each row in *m* represents a single phosphorylation site confidently localized to a specific amino acid residue, with measured abundances across samples specified in columns in *m*. Multiple phosphorylation sites detected on the same peptide have to be converted into separate site-specific entities for every site. While some proteomics software packages, such as  MaxQuant [6], readily produce single site-centric PTM reports, the use of other software packages might require additional preprocessing steps.


## Example datasets

`ssGSEA3.0` installs three worked examples in `example/` (`inst/example/` in
this repository; `system.file("example", package = "ssGSEA3")` once installed),
used by the package vignette:

- **Ordinary ssGSEA (HALLMARK)**: a PANOPLY LUAD tumor NMF factorization
  (`LUAD_matrix_W_combined_signed_gc.gct`) against the MSigDB HALLMARK
  collection (`h.all.v6.2.symbols.gmt`); these results are shown in Figure 1
  of Satpathy et al. [8].
- **Ordinary ssGSEA (REACTOME)**: the same LUAD NMF factorization scored
  against a REACTOME subset (`reactome_subset.gmt`) -- a much larger,
  size/function-heterogeneous database, contrasting with HALLMARK's small,
  curated gene sets.
- **PTM-SEA**: a site-centric phosphoproteomics PI3K-inhibitor perturbation
  dataset, based on Supplemental Table 6 in [1]
  (`PI3K_pert_logP_n2x23936.gct`), against PTMsigDB v2.0.0
  (`ptm.sig.db.all.flanking.human.v2.0.0.gmt`).

Only these example files are installed with the package; the full MSigDB and
PTMsigDB collections stay in this repository's `db/` folder. The MSigDB license
for the two MSigDB files is installed alongside them as `msigdb_license.txt`.


## Legacy standalone scripts

This repository also retains the original standalone R scripts
(`ssgsea-cli.R`, `preprocess_gct-cli.R`) for pipeline integrations (e.g. Terra/
Cromwell WDL) that call them directly rather than installing the R package.
These are not part of the installable `ssGSEA3` package itself. For
command-line use, see ```ssgsea-cli.R --help``` and
```preprocess_gct-cli.R --help```.


## Changes from the original implementation

Original code written by Pablo Tamayo. Adapted with additional modifications by D. R. Mani, Karsten Krug, and Natalie Clark. Adaptions include:

- support of multiple CPU cores (```doParallel``` R-package)
- support of GCT v1.3 format using functions from [cmapR](https://github.com/cmap/cmapR)
- improved handling of missing values
- scoring of directional gene sets (PTMsigDB)
- basic error handling
- improvements in runtime performance
- additional output files like rank plots and parameter files
- replacement of gene-label permutation testing with feature-level (gene- or
  site-level) bootstrap resampling for significance (see `vignette("ssgsea2-vs-ssgsea3-comparison")`)


## License
License Agreement for MSigDB v6.0 and above can be found [here](https://raw.githubusercontent.com/broadinstitute/ssGSEA2.0/master/db/msigdb/license.txt).

## References
1.  Krug, K., Mertins, P., Zhang, B., Hornbeck, P., Raju, R., Ahmad, R., . Szucs, M., Mundt, F., Forestier, D., Jane-Valbuena, J., Keshishian, H., Gillette, M. A., Tamayo, P., Mesirov, J. P., Jaffe, J. D., Carr, S. A., Mani, D. R. (2019). **A curated resource for phosphosite-specific signature analysis**, Molecular & Cellular Proteomics (in Press). http://doi.org/10.1074/mcp.TIR118.000943

1.  Barbie, D. A., Tamayo, P., Boehm, J. S., Kim, S. Y., Susan, E., Dunn, I. F., . Hahn, W. C. (2010). **Systematic RNA interference reveals that oncogenic KRAS- driven cancers require TBK1**, Nature, 462(7269), 108-112. https://doi.org/10.1038/nature08460

1. Abazeed, M. E., Adams, D. J., Hurov, K. E., Tamayo, P., Creighton, C. J., Sonkin, D., et al. (2013).
       **Integrative Radiogenomic Profiling of Squamous Cell Lung Cancer.** Cancer Research, 73(20), 6289-6298.
       http://doi.org/10.1158/0008-5472.CAN-13-1616

1. Subramanian, A., Tamayo, P., Mootha, V. K., Mukherjee, S., Ebert, B. L., Gillette, M. A., et al. (2005).
   **Gene set enrichment analysis: a knowledge-based approach for interpreting genome-wide expression profiles.**
  Proceedings of the National Academy of Sciences of the United States of America, 102(43), 15545-15550. http://doi.org/10.1073/pnas.0506580102

1. Hornbeck, P. V., Zhang, B., Murray, B., Kornhauser, J. M., Latham, V., & Skrzypek, E. (2015). **PhosphoSitePlus, 2014: mutations, PTMs and recalibrations.** Nucleic Acids Research, 43(D1), D512-D520. https://doi.org/10.1093/nar/gku1267
 

1. Cox, J., & Mann, M. (2008). **MaxQuant enables high peptide identification rates, individualized p.p.b.-range mass accuracies and proteome-wide protein quantification.** Nature Biotechnology, 26(12), 1367-1372. https://doi.org/10.1038/nbt.1511


1. Reich, M., Liefeld, T., Gould, J., Lerner, J., Tamayo, P., & Mesirov, J. P. (2006). **GenePattern 2.0.** Nature Genetics, 38(5), 500-501. https://doi.org/10.1038/ng0506-500

1. Satpathy, S., et al. (2025). **Integrative analysis of lung adenocarcinoma across diverse ethnicities and exposures.** Cancer Cell, 43(9), 1731-1757.e17. https://doi.org/10.1016/j.ccell.2025.07.011


***
