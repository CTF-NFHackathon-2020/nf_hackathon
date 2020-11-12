# nf_hackathon
Jupyter notebook and accompanying code for the 2020 NF Hackathon

The `processing_functions.R` library contains the functions needed for loading all requisite data into the notebook for analysis, as well as the helper functions used for computing trajectories and deriving marker genes for evaluating tumor samples. Some of the datasets needed to run this notebook are quite large and do not fit into the data/ folder available here. We will make this data available in a shareable form as soon as we can.

The notebook contains all the commands and calls needed to compute the trajectory tree computing for all developing cells as well as the signature heatmap of all MPNST tumor samples that we processed. 

Data Sources : 

[1] Soldatov, Ruslan, et al. "Spatiotemporal structure of cell fate decisions in murine neural crest." Science 364.6444 (2019). https://science.sciencemag.org/content/364/6444/eaas9536.abstract . This was the source of the E9.5 trunk neural crest data. 

[2] Furlan, Alessandro, et al. "Multipotent peripheral glial cells generate neuroendocrine cells of the adrenal medulla." Science 357.6346 (2017).https://science.sciencemag.org/content/357/6346/eaal3753.abstract . This was the source of the E12.5 and E13.5 Schwann cell precursors.

[3] Wolbert, Jolien, et al. "Redefining the heterogeneity of peripheral nerve cells in health and autoimmunity." Proceedings of the National Academy of Sciences 117.17 (2020): 9466-9476. https://www.pnas.org/content/117/17/9466.short . This was the source of the nmSC and mySC (non-myelinating and myelinating Schwann cells). 

[4] Fletcher, Jonathan S., et al. "Cxcr3-expressing leukocytes are necessary for neurofibroma formation in mice." JCI insight 4.3 (2019). https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6413799/ . This was the source of the Scwhann cell single-cell RNA-seq profiles from a mouse plexiform neurofibroma tumor.

[5] https://www.ncbi.nlm.nih.gov/projects/gap/cgi-bin/study.cgi?study_id=phs000792.v1.p1 . Source of the bulk human MPNST tumor samples.

Methods References : 

[6] Aibar, Sara, et al. "SCENIC: single-cell regulatory network inference and clustering." Nature methods 14.11 (2017): 1083-1086. The AUCell package was used for scoring neural crest and Schwann cell signatures in Schwann cells from [4].

[7] Stuart, Tim, et al. "Comprehensive integration of single-cell data." Cell 177.7 (2019): 1888-1902. The Seurat package.

[8] Qiu, Xiaojie, et al. "Reversed graph embedding resolves complex single-cell trajectories." Nature methods 14.10 (2017): 979. Monocle3 package.
