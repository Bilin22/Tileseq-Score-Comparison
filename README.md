[//]: # (# Tileseq-Score-Comparison)

[//]: # ()
[//]: # (![newversion]&#40;https://github.com/Bilin22/Tileseq-Score-Comparison/blob/main/SUMO1/mavevis_SUMO1/newver.png&#41;)

[//]: # (![oldversion]&#40;https://github.com/Bilin22/Tileseq-Score-Comparison/blob/main/SUMO1/mavevis_SUMO1/oldver.png&#41;)

[//]: # (> A comparison of the variant effect maps generated by the 2023 &#40;top&#41; and 2019 &#40;bottom&#41; versions of MAVE pipelines for SUMO1 gene.)


**This [BCB330](https://artsci.calendar.utoronto.ca/course/bcb330y1) project aims to re-evaluate the performance of variant effect maps based on different versions of MAVE pipelines
with respect to precision and sensitivity on reliable benchmarks.**

The detailed goals include:
1. Re-process the raw data underlying existing variant effect maps with the latest versions of their respective analysis pipelines.
2. Compile benchmark sets of variants with known pathogenicity from online databases and literature for each map.
3. Compare the predictions made by different versions of variant effect maps using the benchmark sets and use them to infer evidence strength for clinical interpretation.
4. Provide recommendations for optimizing the implementation of MAVEs based on the evaluation result.

**Here's the list of genes we re-processed:**
* Calmodulin 1: [CALM1](https://github.com/Bilin22/Tileseq-Score-Comparison/wiki/CALM1)
* GDP Dissociation Inhibitor 1: [GDI1](https://github.com/Bilin22/Tileseq-Score-Comparison/wiki/GDI1)
* Small Ubiquitin Like Modifier 1: [SUMO1](https://github.com/Bilin22/Tileseq-Score-Comparison/wiki/SUMO1)
* Trans-2,3-Enoyl-CoA Reductase: [TECR](https://github.com/Bilin22/Tileseq-Score-Comparison/wiki/TECR)
* Ubiquitin Conjugating Enzyme E2 I: [UBE2I](https://github.com/Bilin22/Tileseq-Score-Comparison/wiki/UBE2I)
* Methylenetetrahydrofolate reductase: [MTHFR](https://github.com/Bilin22/Tileseq-Score-Comparison/wiki/MTHFR)

Reference: [tileseqMave](https://github.com/rothlab/tileseqMave#joining-variant-counts-and-computing-marginal-frequencies)
