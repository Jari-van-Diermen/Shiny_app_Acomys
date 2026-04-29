# Shiny_app_Acomys
 
## About the R Shiny Web Resource
 
This repository hosts the code for the R Shiny app accompanying the manuscript titled **Decoding mammalian phenotypes through protein domain evolution in African spiny mice**.
 
The purpose of this app is to aid in the exploration of the episodic positive selection results presented in the manuscript. Some useful use-cases include:
 
- **Exploration of the positively selected gene list** (i.e. HyPhy aBSREL results)
    - Exploring the positive selection likelihood-ratio test (LRT) results of the *Acomys cahirinus* (African spiny mice) positively selected gene list.
    - Exploring the functional overrepresentation analysis results (i.e. GO-terms) of the positively selected STRING clusters.
- **Exploration of positively selected codon positions** (i.e. HyPhy MEME results)
    - Exploring the **positive selection LRT results** of codon positions across genes in the positively selected gene list.
    - Exploring **codon-level positive selection signals** in genes from the positively selected gene list.
    - Exploring *A. cahirinus* and *H. sapiens* **amino acid substitutions** and how these relate to codon-level positive selection signals.
    - Exploring the **protein domains** in which codons with positive selection signals reside.

## Web Resource Structure
 
### aBSREL Results
 
Exploration in the web resource is divided between aBSREL and MEME results. The **aBSREL results** represent the estimated positively selected gene list, and allow exploration of the aBSREL LRT results for specific genes of interest. It is also possible to inspect the GO-terms associated with each STRING cluster in the positively selected gene list.
 
### MEME Results
 
The MEME results represent codon-level positive selection signals across the aBSREL gene list, calculated using HyPhy MEME. Several tabs provide different views of the positive selection dataset:
 
- A **MEME LRT results** tab displaying ortholog-specific LRT results and other estimated model parameters, such as synonymous and nonsynonymous substitution rates.
- A **multiple sequence alignment (MSA)** tab for exploring the MSAs used in the aBSREL and MEME evolutionary models.
    - These MSAs are annotated with *MEME site p-values* (a measure of codon-level positive selection) and the *A. cahirinus* empirical Bayes factor (EBF), a rough exploratory measure that indicates the branches where the positive selection signal may have originated. These metrics help identify potentially interesting regions in relation to positive selection.
- An **evolutionary tree viewer** for exploring MEME-estimated amino acid substitutions at specific tree branches and codon positions.
    - Filter options allow selection of codon positions with:
        - Significant MEME LRT p-values.
        - Amino acid substitutions in the *A. cahirinus* or *H. sapiens* branch.
- An **EBF tab** displaying an interactive view of the empirical Bayes factor (EBF) across tree branches and codon positions for selected genes. This visualization uses EBF values to provide insight into the positive selection signals of specific genes — particularly which codons and branches may be driving selection.
- A **UniProt sequence annotations tab** integrating MEME codon-level positive selection results with the UniProtKB/Swiss-Prot database. This interactive visualization places codon-level positive selection signals within specific protein domains, potentially informing which aspects of protein function are affected by amino acid substitutions. The same filtering options as the tree viewer are available to help focus on specific codons of interest.
## Where to Find the Web Resource
 
The web resource can be visited [here](https://jari-van-diermen-shiny-app-acomys.share.connect.posit.cloud/).
 
If you encounter any problems with the web resource, feel free to contact me or open an Issue.