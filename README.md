# miRNA-Analysis
This pipeline analyzes smallRNA-seq data by culminating several well established smallRNA-seq analysis packages into one automated workflow.

![miRNA-Analysis Pipeline Overview](bin/miRNA-Analysis.png)

1. Preprocessing
    1. FastQC: assesses raw and post-trimming and filtering read qualities
    2. Cutadapt: the perform trimming
    3. UniVec filtering: this is from the exceRpt pipeline
2. Known Expression
    1. exceRpt: performs necessary alignments to reference datas in the following order: rRNA...
    2. HTSeq-count: takes the genome alignment bam from exceRpt and obtains the raw counts for known miRNAs (location from miRBase) specifically for DESeq2 analysis
3. Novel miRNA Detection
    1. miRDeep2: the best package currently for detecting novel miRNAs. The subworkflow for this was borrowed from the nf-core community, with some modifications.
    2. mergemiRDeep2: script to combine the detected novel miRNA across all samples based on the predicted hairpin sequence.
4. Differential Expression
    1. DESeq2: will automatically perform differential expression analysis of every condition combination based on the metadata file given. Plots will be created for each combination as well.
5. Functional Analysis
    1. Targets of miRNAs: obtains the gene targets of the differentially expressed miRNAs. The user will need to specify which DESeq2 output or list of miRNAs they are interested int.
    2. Intersect miRNA and RNA-seq: If the user chooses to provide Bulk RNA-seq counts from the same tissue, this process will select only the miRNA targeted genes from the previous process that appear in the Bulk RNA-seq. The user has the option to provide a minimum expression level for each gene in the Bulk RNA-seq. A count with the number of miRNAs targeting the gene is provided in the second column of the output file.
    3. PANTHER: Uses the PANTHER API to obtain gene enrichment of the intersected genes and the number of differentially miRNAs that target it. **To be replaced with clusterProfiler**
    4. Interpreting PANTHER output: takes PANTHER's output files and creates bubble charts and table for the enrichment results.

## Installation and Setup of Pipeline
Download and setup conda (can use anaconda or miniconda) in your Linux machine's directory. You can find instructions on how: [anaconda download and setup](https://docs.anaconda.com/anaconda/install/)
Download the repository to your Linux machine using:
```
git clone https://github.com/utdal/miRNA-Analysis.git
```
Create conda environment for pipeline:
```
conda create --name mirna
conda activate mirna
conda install nextflow
conda install singularity
```

## Running pipeline
I recommend running it this way:
Go into the conf folder and open the user_sample.conf file. Configure parameters as needed.

To run miRNA_Expression portion of pipeline:
```
nextflow run . -profile user_sample,singularity --outdir <output_directory_name>
```
Other optional commands you can add
```
-bg # can add after nextflow in previous command to run the background (does not require terminal in which the command was executed to be open)
-resume # If you stop a run mid way, add this to the end of the command you used to run the pipeline and it will cache all the processes that have already run.
--email # to email when the pipeline is done running.
```


## Sources
This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).
 
> The nf-core framework for community-curated bioinformatics pipelines.
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> Nat Biotechnol. 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://doi.org/10.1038/s41587-020-0439-x).


This pipeline utilizes the following programs:
> Rozowsky, J., Kitchen, R. R., Park, J. J., Galeev, T. R., Diao, J., Warrell, J., Thistlethwaite, W., Subramanian, S. L., Milosavljevic, A., & Gerstein, M. (2019). exceRpt: A Comprehensive Analytic Platform for Extracellular RNA Profiling. Cell Systems, 8(4), 352-357.e3. https://doi.org/10.1016/j.cels.2019.03.004

> miRDeep2 accurately identifies known and hundreds of novel microRNA genes in seven animal clades | Nucleic Acids Research | Oxford Academic. (n.d.). Retrieved February 26, 2025, from https://academic.oup.com/nar/article/40/1/37/1275937?login=true

> Love MI, Huber W, Anders S (2014). “Moderated estimation of fold change and dispersion for RNA-seq data with DESeq2.” Genome Biology, 15, 550. doi:10.1186/s13059-014-0550-8.

> Simon Anders, Paul Theodor Pyl, Wolfgang Huber HTSeq — A Python framework to work with high-throughput sequencing data Bioinformatics (2014), in print, online at doi:10.1093/bioinformatics/btu638

> Thomas, P. D., Ebert, D., Muruganujan, A., Mushayahama, T., Albou, L.-P., & Mi, H. (2022). PANTHER: Making genome-scale phylogenetics accessible to all. Protein Science, 31(1), 8–22. https://doi.org/10.1002/pro.4218

> Li JH, et al.starBase v2.0: decoding miRNA-ceRNA, miRNA-ncRNA and protein-RNA interaction networks from large-scale CLIP-Seq data , Nucleic Acids Res. 2014 Jan;42:D92-7. (A new manuscript has been submitted: Zhou KR, Huang JH, Liu S, Zheng WJ, Liu SR, et al. An Encyclopedic Regulatory and Functional Atlas of RNA Interactomes.)