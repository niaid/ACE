##  Exercises: Gene Ontology, Enrichment, Pathways, Networks

####    Gene Ontology 
> The Gene Ontology started as collaboration between the databases for model organisms of mouse (MGI), fruit fly (FlyBase) and baker’s yeast (SGD).  It now receives contributions from over [30 groups](http://geneontology.org/page/go-consortium-contributors-list), which include WormBase, dictyBase, J.Craig Venter, Reactome and others.
    •	GO describes how gene products behave in a cellular context.
    The GO consists of three controlled vocabularies
        o	Molecular function (what does the gene product do?)_ 
        o	Biological process (why does it perform the activity?)
        o	Cellular component (where does it act?) 
Example:   6 – phosphofructokinase

    Molecular Function
    •	GO:0003872, “6-phophofructokinase activity”
    Biological Process
    •	GO:0006096, “glycolysis”
    Cellular Component
    •	GO:0005737, “cytoplasm”


### Excercise #1 - Search Gene Ontology for Gene of interest

> In the exercises below we will explore the protein interactions of genes AOAH and HLA class II genes (HLA-DRB1) in hopes to better understand a "trans" associatiion between them. 
Background info: A paper published in [Nature Genetics, 2012](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3437404/) reported that HLA antigens could form "trans" associations with the expression of AOAH and ARHGAP24 in monocytes (white blood cell) but not in B cells but the mechanism is still unknown.
[Figure 5c](https://www.ncbi.nlm.nih.gov/core/lw/2.0/html/tileshop_pmc/tileshop_pmc_inline.html?title=Click%20on%20image%20to%20zoom&p=PMC3&id=3437404_ukmss-40993-f>0005.jpg) shows that the number of alleles carried in the genome is significantly associated with reduced expression of AOAH and increase in ARHGAP24

##### 1. Go to [Panther.org](http://pantherdb.org) and enter the gene AOAH in the box, then select Organism, Analysis and proceed to submit.  View results: Click on the geneID, then expand section GENE ONTOLOGY DATABASE ANNOTATIONS
---
#### Optional: Explore gene information
- AOAH - [browser view](https://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr7%3A36340943-36975760&hgsid=228955823_ecmjMB4Q2o4Dk8WMKrI0IwNvyANl)
- HLA-DRB1 - [browser view of DRB1](https://genome-euro.ucsc.edu/cgi-bin/hgTracks?db=hg19&lastVirtModeType=default&lastVirtModeExtraState=&virtModeType=default&virtMode=0&nonVirtPosition=&position=chr6_ssto_hap7%3A3889108-3908157&hgsid=228955823_ecmjMB4Q2o4Dk8WMKrI0IwNvyANl)


####    The SNP reference in the paper with ID "rs28366298" is within the binding site for TF YY1
rs28366298 is an A/G/T single-nucleotide variation on human chromosome 6. [View more(https://www.ncbi.nlm.nih.gov/SNP/snp_ref.cgi?rs=28366298&pt=1fS4Q8icZO4m2YHtJFTQs_ZcypnH1xjfuY226v2zAj-fqzJwUP)


### 2. Let's now see if these three genes (HLA-DRB1, YY1 and AOAH) are known to interact using https://string-db.org/)
- First explore STRING using the "More" feature to bring connections between proteins. Export network file.
- Also visit genemania.org, view interactions for the same genes and export network.

### Exercise #2 - Use [DAVID Bioinformatics](https://david.ncifcrf.gov/) to run functional and pathway enrichment by explore the list of genes identified above with STRING.  Which pathways are associated?  
ENSG00000248993
KIR2DL4
LILRB2
CD8A
B2M
TROVE2
HLA-DMA
HLA-DQA2
CD4
PTPN22
AOAH
HLA-DRA
TRIM21
HLA-G
CSK
HLA-DPB1
HLA-DRB1
CD74
HLA-DRB5
HLA-DOA
HLA-DMB
HLA-DQA1

### Exercise #3 GSEA  
1- Explore Diabetes expression dataset using GSEA.  
For reference view [tutorial here](http://software.broadinstitute.org/gsea/doc/desktop_tutorial.jsp)
    - Upload data following demonstration and/or tutorial 
    - Download data from [gsea datasets](http://software.broadinstitute.org/gsea/datasets.jsp)
    
    
### Exercise #4: Cytoscape

1. Import the network data downloaded in previous steps into Cytoscape
- Perform import by using menu File --> Import --> Network

2. (Optional) Follow the [official cytoscape tutorial](http://manual.cytoscape.org/en/stable/Basic_Expression_Analysis_Tutorial.html)

    

