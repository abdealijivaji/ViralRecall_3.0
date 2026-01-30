---
title: 'Viralrecall - A tool to identify giant viruses integrated into eukaryotic genomes'
tags:
  - Giant virus
  - Python
  - Genomics
authors:
 - name: Abdeali M. Jivaji
   orcid: 0000-0003-2474-9561
   affiliation: 1
   corresponding: true
 - name: Frank O. Aylward
   orcid: 0000-0002-1279-4050
   affiliation: "1, 2"
affiliations:
 - name: Dept. of Biological Sciences, Virginia Tech
   index: 1
 - name: Department of Biological Sciences and Center for Emerging, Zoonotic, and Arthropod-borne Pathogens, Virginia Tech, Blacksburg, Virginia, USA
   index: 2
date: 1 February 2027
bibliography: paper.bib
---

# Summary

Viralrecall 3.0 is a python tool to identify Giant Endogenous Viral Elements (GEVEs) integrated in the genome of eukaryotes. The current version is an update on the original tool by Dr. Aylward and uses the same GVOG HMM database to detect signatures of giant viruses [@aylward:2021]. The key motivation for updating `viralrecall` was to make it more efficient at processing the larger eukaryotic genomes that are being published with the rise in popularity of long-read sequencing.

# Statement of Need

This update was needed since the older version (v2.0) was only accessible as a python script that may not be beginner-friendly and easy to use. It was based on an older version of python that required specific and outdated versions of crucial packages that were not easy to setup/troubleshoot.
With v3.0, we provide an easy to install option via conda as it is widely used by bioinformaticians worldwide. We provide a direct script to install and setup the database automatically reducing errors. We also include an experimental database of key Mirusvirus proteins to identify these mysterious viruses that share similarities with giant viruses.
With the rise of long-read DNA sequencing and the increased abundance of high-quality and complete eukaryotic genomes, this version of the tool has been optimized to be memory-efficient and blazingly fast compared to the v2.0. Specifically in batch mode, where all genome files in a given directory are processed, the tool maximizes the utilization of available cores present in modern high-performance compute clusters to process thousands of genomes in an efficient and quick manner.
With this, we hope to bring inexperienced folks and proficient bioinformaticians alike into the field of giant viruses and identify these mysterious viruses throughout the tree of life [@moniruzzaman:2020].

# State of the Field

Although several tools exist to identify giant viruses based on DNA sequences, these tools are generally broad-ranged and only include identification of giant viruses as a minor feature. Specifically, the tool `genomad` is an excellent resource to identify various phages and mobile genetic elements [@camargo:2024]. Since giant viruses are significantly different from most previously identified viruses, their dataset required special categories for giant virus identification. While `genomad` is an excellent resource, our tool is specifically designed to identify giant viruses, either present in metagenomes as independent contigs, or integrated into the genome of eukaryotes as GEVEs. Thus, we provide a more curated database and tool specifically designed to identify these enigmatic viruses.

# AI usage disclosure

We used Github copilot for basic autocomplete and the chatbot for suggestions regarding some functions, and for writing tests. Generative AI was not used for wrtiting this manuscript.

# Author Contributions

**Abdeali M. Jivaji:** Methodology, Software, Writing - Original Draft. **Frank O. Aylward:** Review, Conceptualization, Resources, Supervision.

# Acknowledgements

We thank the members of the Aylward Lab for their support and comments during the development of this tool.
