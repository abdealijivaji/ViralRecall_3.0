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

Viralrecall 3.0 is a python tool to identify Giant Endogenous Viral Elements (GEVEs) integrated in the genome of eukaryotes. The current version is an update on the original tool by Dr. Aylward and uses the same GVOG HMM database to detect signatures of giant viruses. The key motivation for updating Viralrecall was to make it more efficient at processing the larger euykaryotic genomes that are being published with the rise in popularity of long-read sequencing.

# Statement of Need

This update was needed since the older version (v2.0) was only accessible as a python script that may not be beginner-friendly and easy to use. It was based on an older version of python that required specific and outdated versions of crucial packages that were not easy to setup/troubleshoot.
With v3.0, we provide an easy to install option via conda as it is widely used by bioinformaticians worldwide. The required database has also been refined to eliminate the PFAM database. We also provide a direct script to install and setup the database automatically reducing errors. We also include an experimental database of key Mirusvirus proteins to identify these mysterious viruses that share similarities with giant viruses.
With the rise of long-read DNA sequencing and the increased abundance of high-quality and complete eukaryotic genomes, this version of the tool has been optimized to be memory-efficient and blazingly fast compared to the v2.0. Specifically in batch mode, where all genome files in a given directory are processed, the tool maximizes the utilization of available cores present in high-performance compute clusters to process thousands of genomes in an efficient and quick manner.
With this, we hope to bring inexperienced folks and proficient bioinformaticians alike into the field of giant viruses and identify these mysterious viruses throughout the tree of life.

# References