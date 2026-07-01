# Non-coding Mutations in Multiple Myeloma
Repository for analysis in "Recurrent non-coding mutational hotspots drive chromosomal
rearrangement and disease progression in multiple myeloma" from Blaney, P. et al.

## Primary Analysis and Visualization
The crux of all data analysis and visualization were performed in `R` with the source
code found in `index.qmd` and the rendered HTML at [`_manuscript/index.html`](https://github.com/pblaney/myelomaNoncodingGenome/blob/main/_manuscript/index.html).

## Data Generation

### WGS Data Processing
All matched tumor/normal whole genome sequencing data was preprocessed and analyzed
for somatic variants uniformly using our Nextflow pipeline, the [`mgp1000`](https://github.com/pblaney/mgp1000).

### Mutation Enrichment Identificaiton
Enrichment of SNVs and InDels across the genome was identified using [`fishHook`](https://github.com/mskilab-org/fishHook).
The model was deployed on NYU's HPC UltraViolet due to the large number of samples,
mutations to analyze, and covariates included. The two scripts utilized for this
are `angler.R` and `sushi.R` (kept with in the `scripts/` directory).

Example of execution of scripts for Pan-Myeloma fishHook model:
```bash
# Run the model
sbatch \
  --job-name=goneSnvFishing \
  --time=5:00:00 \
  --mem=96G \
  --cpus-per-task=1 \
  --partition=fn_medium \
  --wrap="singularity exec -B $PWD:/data/ --pwd /data/ fishhook-0.1.simg \
            ./angler.R mgp1000_snv_mutations.rds mgp1000_snv.tuna.rds"

# Compute the P-values
sbatch \
  --job-name=makingSnvSushi \
  --time=3:00:00 \
  --mem=96G \
  --cpus-per-task=1 \
  --partition=fn_short \
  --wrap="singularity exec -B $PWD:/data/ --pwd /data/ fishhook-0.1.simg \
            ./sushi.R mgp1000_snv.tuna.rds mgp1000_snv"
```
>[`fishhook-0.1` Docker container<sup>#</sup>](https://hub.docker.com/layers/patrickblaneynyu/tools/fishhook-0.1/images/sha256-55b84d013a1062938558513deaec1f515be7f96bacb7bcc36a19d477e1c68819)

### Mutational Signature Extraction
De novo extraction of mutational signatures was performed using [`SigProfilerExtractor`](https://github.com/AlexandrovLab/SigProfilerExtractor) with `COSMIC` v3.3 reference signatures. SNVs that fell within the mutationally enriched regions (MERs) identified by `fishHook` were isolated per sample and used as input for extraction.

Code used for extraction from within SigProfiler Docker container:
```python
# BCFtools to extract SNVs within MERs for each patient
# for vcf in `ls -1 *.vcf.gz`; do bcftools view -O v -R MERs.bed ${vcf} | bgzip > merSubset/"${vcf}"; done

# Now run the extraction
sig.sigProfilerExtractor(
  input_type="vcf",output="/data/cortado",input_data="/data/mgp/",
  reference_genome="GRCh38",opportunity_genome="GRCh38",context_type="default",
  exome=False,minimum_signatures=1,maximum_signatures=3,nmf_replicates=100,
  resample=True,batch_size=1,cpu=-1,gpu=False,nmf_init="random",precision="single",
  matrix_normalization="gmm",seeds="random",min_nmf_iterations=10000,
  max_nmf_iterations=1000000,nmf_test_conv=10000,nmf_tolerance=1e-15,
  get_all_signature_matrices=False,stability=0.80,min_stability=0.2,
  combined_stability=1.0,allow_stability_drop=False,cosmic_version=3.3,
  make_decomposition_plots=True,collapse_to_SBS96=True)
``` 
>[`sigprofiler` Docker container](https://hub.docker.com/layers/dceoy/sigprofiler/latest/images/sha256-7ab802ead7d72b0bd5ab5a3a507fdbcdcfc3e42d575d122cc6387dc3c850c26a)

### Chromatin Topography Assessment
Enrichment of SNVs and InDels at mutationally enriched regions (MERs) in the context
of chromatin topography was assessed using [`SigProfilerTopography`](https://github.com/AlexandrovLab/SigProfilerTopography)
and Pan-Myeloma set of histone modification marks H3K27ac, H3K4me1, H3K27me3, H3K4me3,
H4K20me1, CTCF binding sites, DNase I hypersensitive sites and transposase-accessible
chromatin (ATAC) derived from myeloma cell lines and MM patients.

Code used for topography assessment from within SigProfiler Docker container:
```python
import os
os.system('apt-get update && apt-get install -y wget')
os.system('pip install SigProfilerTopography')
from SigProfilerTopography import Topography as topography
from SigProfilerMatrixGenerator import install as genInstall
genInstall.install(genome="GRCh38",offline_files_path="/data/ref/",rsync=False,bash=True)

topography.runAnalyses(
  genome="GRCh38",inputDir="/data/mgp/",outputDir="/data/topography/",
  jobname="PanMyeloma-MER",numofSimulations=5,gender="male",epigenomics=True,
  plot_epigenomics=True,
  epigenomics_files=[
    '/data/ENCFF653YKZ_MM1S_H3K4me1.bed','/data/ENCFF825UPB_MM1S_H3K4me3.bed',
    '/data/ENCFF577YXE_MM1S_H3K27me3.bed','/data/ENCFF121HFO_MM1S_H3K36me3.bed',
    '/data/ENCFF148MKR_MM1S_DNase.bed','/data/ENCFF539QXX_MM1S_CTCF.bed',
    '/data/ENCFF362XCO_MM1S_H4K20me1.bed','/data/ENCFF452IWV_KMS11_H3K4me1.bed',
    '/data/ENCFF400AMJ_KMS11_H3K4me3.bed','/data/ENCFF293AAL_KMS11_H3K27me3.bed',
    '/data/ENCFF838NMF_KMS11_H3K36me3.bed','/data/ENCFF496TKW_KMS11_CTCF.bed',
    '/data/ENCFF683CGI_KMS11_H4K20me1.bed','/data/ENCFF874LYU_KMS11_H3K27ac.bed',
    '/data/RefEpigenome_Myeloma_ATAC.bed','/data/RefEpigenome_Myeloma_H3K4me1.bed',
    '/data/RefEpigenome_Myeloma_H3K4me3.bed','/data/RefEpigenome_Myeloma_H3K27ac.bed',
    '/data/RefEpigenome_Myeloma_H3K27me3.bed','/data/RefEpigenome_Myeloma_H3K36me3.bed'
  ],
  epigenomics_dna_elements=[
    'H3K27ac','H3K4me1','H3K4me3','H3K27me3',
    'H3K36me3','DNase','ATAC','CTCF','H4K20me1'
  ],
  epigenomics_biosamples=['Myeloma','MM1S','KMS11'],
  sbs_probabilities="/data/MER_De_Novo_MutationType_Probabilities.txt",
  exceptions={'SBS1':0.0424,'SBS5':0.8068,'SBS9':0.1139,'SBS84':0.0370},
  average_probability=0.01,num_of_sbs_required=100)
```
>[`sigprofiler` Docker container](https://hub.docker.com/layers/dceoy/sigprofiler/latest/images/sha256-7ab802ead7d72b0bd5ab5a3a507fdbcdcfc3e42d575d122cc6387dc3c850c26a)

### Clustered Mutation Event Detection
Clustered mutation events were found and classified using [`SigProfilerClusters`](https://github.com/AlexandrovLab/SigProfilerClusters) where intramutational distance is compared between SNVs to a sample-dependent threshold and VAF is used to determine the relative timing of the mutations. The sample-dependent IMD threshold is determined by simulations using [`SigProfilerSimulator`](https://github.com/AlexandrovLab/SigProfilerSimulator).

```python
from SigProfilerMatrixGenerator import install as genInstall
from SigProfilerMatrixGenerator.scripts import SigProfilerMatrixGeneratorFunc as matGen
genInstall.install(genome="GRCh38",offline_files_path="/temp/data/ref/",rsync=False,bash=True)

from SigProfilerSimulator import SigProfilerSimulator as sigSim
sigSim.SigProfilerSimulator(
  project="clusters",project_path="/temp/data/mgp/",genome="GRCh38",
  contexts=["96","ID"],simulations=100,chrom_based=True,vcf=True,gender="male")

import os
os.system("pip install seaborn")
from SigProfilerClusters import SigProfilerClusters as hp
hp.analysis(
  project="clusters",genome="GRCh38",contexts="96",simContext=["96"],
  input_path="/temp/data/mgp/",analysis="all",sortSims=False,calculateIMD=True,
  max_cpu=6,subClassify=True,includedVAFs=True,standardVC=False,TCGA=True,
  sanger=False,plotRainfall=False,probability=True)
```
>[`sigprofiler` Docker container](https://hub.docker.com/layers/dceoy/sigprofiler/latest/images/sha256-7ab802ead7d72b0bd5ab5a3a507fdbcdcfc3e42d575d122cc6387dc3c850c26a)


#### <sup>#</sup>Containerized Execution
Docker containers are not prohibited from being used on HPCs so each <sup>#</sup>Docker
container was converted to a Singularity `.simg` container using `singularityware/docker2singularity:latest`
