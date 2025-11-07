---
title: 'RustSASA: A Rust Crate for Accelerated Solvent Accessible Surface Area Calculations'
tags:
  - Rust
  - Bioinformatics
  - Proteomics
  - Molecular Dynamics
authors:
  - name: Maxwell J. Campbell
    orcid: 0000-0002-0959-1164
    affiliation: "1"
affiliations:
 - name: University of California, San Francisco, United States
   index: 1
date: 16 July 2025
bibliography: paper.bib

---

# Summary

Solvent accessible surface area (SASA) calculations are fundamental for understanding protein structure, function, and dynamics in computational biology. These calculations quantify the surface area of biomolecules accessible to solvent molecules, providing insights into protein folding, stability, and intermolecular interactions. The Shrake-Rupley algorithm has served as the standard for SASA calculations since 1973, but existing implementations often become computational bottlenecks when analyzing large protein datasets. As proteomics datasets continue to grow with initiatives like AlphaFold producing hundreds of millions of predicted protein structures the need for efficient SASA calculation tools has increased dramatically. RustSASA addresses this challenge by providing a high-performance implementation of the Shrake-Rupley algorithm written in pure Rust, delivering a 7× speed improvement over Freesasa while maintaining calculation accuracy and providing interfaces for multiple programming languages and frameworks (i.e: MDAnalysis).

# Statement of need

Current SASA calculation tools represent a significant computational bottleneck in structural biology workflows, particularly for molecular dynamics simulations and high-throughput analyses. Popular implementations such as those in Biopython and Freesasa, while accurate, become prohibitively slow when processing large protein datasets.

RustSASA addresses this performance gap by leveraging Rust's zero-cost abstractions and memory safety guarantees to create a SASA calculation crate that is significantly faster than Freesasa and Biopython. Benchmarking on representative protein structures demonstrates that RustSASA achieves a 7× improvement over Freesasa, and a 46× performance improvement over Biopython. This performance advantage reduces computational costs for high-throughput structural analyses and makes large-scale comparative studies feasible. Furthermore, RustSASA's multi-language support (Rust and Python), command-line interface, and MDAnalysis package ensure broad accessibility across the computational biology community.

# Results

## Calculation Quality

![Comparing Freesasa and RustSasa on E. coli proteome\label{fig:example}](eval/figures/sasa_chain_comparison_E_coli.pdf){ width=50% }
![Comparing Freesasa and RustSasa on Freesasa comparison dataset\label{fig:example}](eval/figures/sasa_chain_comparison_freesasa_ds.pdf){ width=50% }

To evaluate the accuracy of RustSASA calculations, we compared results to Freesasa [@Mitternacht_2016] on both the predicted E. coli proteome from AlphaFold [@Jumper2021] and the Freesasa evaluation dataset. RustSASA produces SASA values that closely match those from Freesasa, achieving Pearson correlation coefficients > 0.98 on both datasets.

## Performance

![Comparing Freesasa, RustSasa, and Biopython performance on E. coli proteome\label{fig:example}](eval/figures/performance_comparison.pdf){ width=100% }

We evaluated the performance of Freesasa, RustSASA, and Biopython [@biopython] on the predicted E. coli proteome using Hyperfine [@Hyperfine] with three runs and three warmup iterations on a 2024 Apple MacBook Air with an M3 processor and 24GB of unified memory. All methods utilized parallel processing across eight cores. GNU parallel [@Tange2011a] was used to parallelize Freesasa and Biopython, while RustSASA's utilized its internal parallelization.

RustSASA processed the entire proteome in ~8 seconds compared to ~55 seconds for Freesasa and ~368 seconds for Biopython, representing 7× and 46× speed improvements, respectively.

# Conclusion

RustSASA provides a significant advancement in SASA calculation performance while maintaining accuracy, addressing a bottleneck in computational structural biology. The 7× speed improvement over current standards enables previously intractable analyses of large protein datasets and molecular dynamics simulations. By providing interfaces for multiple programming languages alongside a command-line tool and MDAnalysis package, RustSASA ensures broad accessibility across the research community. As structural biology datasets continue to expand, efficient computational tools like RustSASA become essential for advancing our understanding of protein structure and function.

# Acknowledgements

We would like to thank Rodrigo Honorato and Niccolò Bruciaferri for their valuable contributions to this project.

# References
