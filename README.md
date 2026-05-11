Vexomiser — Vector-extended Exomiser
=====================================

> Fork of [Exomiser 15.0.0](https://github.com/exomiser/Exomiser/releases/tag/15.0.0) extended with embedding-based gene–disease association (GDA) scoring.

**Added scoring methods:**

- **INDIGENA** — inductive, similarity/embedding-based GDA scoring
- **Multihop-GDA** — link-prediction / logical query answering for GDA

Both methods are integrated into Exomiser's variant prioritisation pipeline as additional `GeneScorer` implementations.

## Setup

### Requirements

- Java 21 (required to build and run Exomiser)
- [uv](https://github.com/astral-sh/uv) (Python package manager)
- Maven (bundled via `mvnw`)

On Ubuntu/Debian:
```bash
sudo apt install openjdk-21-jdk
export JAVA_HOME=/usr/lib/jvm/java-21-openjdk-amd64   # add to ~/.bashrc
curl -LsSf https://astral.sh/uv/install.sh | sh
```

### Build

```bash
git clone https://github.com/ferzcam/vexomiser.git
cd vexomiser
chmod +x mvnw
./mvnw install -pl phenix-repository        # installs local JARs not on Maven Central
./mvnw package -pl exomiser-cli -am         # builds the CLI (runs tests)
```

To skip tests (faster):
```bash
./mvnw package -pl exomiser-cli -am -DskipTests
```

### Python environment

```bash
uv sync   # creates .venv and installs all dependencies
```

---

## Evaluation data

### Benchmark cases (PAVS)

The benchmark data lives on Ibex and is **not tracked in git**:

```
/ibex/scratch/projects/c2014/exomiser_extension/
```

To sync locally:
```bash
rsync -av zhapacfp@ilogin.ibex.kaust.edu.sa:/ibex/scratch/projects/c2014/exomiser_extension/ data/
```

Structure after sync:
```
data/
├── track1_cases.tsv          # 5,076 phenotype-only cases
├── track2_cases.tsv          # 2,815 phenotype + genotype cases (subset of Track 1)
├── pavs/phenopackets/        # one JSON per case (HPO terms + variant)
├── spiked_vcfs/              # Track 2 input VCFs (GIAB HG001 + causal variant)
└── causal_vcfs/              # causal variant only, one VCF per case
```

**Track 1** — phenotype-only gene ranking. Input: HPO terms. Ground truth: causal gene.

**Track 2** — variant prioritisation. Input: HPO terms + spiked VCF (~94K background variants from GIAB HG001 with the known causal variant inserted). Ground truth: causal variant + gene.

### Exomiser phenotype data

Download separately (~6.2 GB, not in git):
```bash
mkdir -p exomiser-data
wget -P exomiser-data https://data.monarchinitiative.org/exomiser/latest/2406_phenotype.zip
unzip exomiser-data/2406_phenotype.zip -d exomiser-data/
```

---

## Data splits

Train/val/test splits are generated with `eval/generate_splits.py`. The split is:

- **Disease-disjoint**: all cases sharing the same MONDO disease ID go to the same partition. Cases without a MONDO ID (mostly the Saudi cohort) receive a unique synthetic ID so each is treated independently.
- **Stratified by cohort** (Saudi / DDD / Mixed) to preserve population distribution.
- **Ratio**: 8 / 1 / 1

```bash
uv run python eval/generate_splits.py --track 1
```

Output: `data/splits/{train,val,test}.tsv`

| Split | Cases | Track 2 subset |
|-------|-------|----------------|
| Train | 4,062 | 2,261 |
| Val   |   496 |   270 |
| Test  |   518 |   284 |

**Important:** Train on Track 1 splits. The same trained model is evaluated on both Track 1 and Track 2 test splits — no retraining needed between tracks.

---

## Evaluation

Runs Exomiser's built-in phenotype prioritisers (HiPhive, Phive, PhenIX) on all cases and computes per-case gene ranking scores. Results are used as baselines for INDIGENA and Multihop-GDA.

```bash
uv run python eval/exomiser_eval.py \
    --phenotype-data-dir exomiser-data/2406_phenotype \
    --track 1
```

To run a subset of prioritisers (e.g. to resume a partial run):
```bash
PRIORITISERS=phenix uv run python eval/exomiser_eval.py \
    --phenotype-data-dir exomiser-data/2406_phenotype \
    --track 1
```

Output: `data/results/exomiser_{hiphive,phive,phenix}_track{1,2}.tsv`

Each TSV row: `gene_symbol <TAB> case_id <TAB> gene_index <TAB> score_0 ... score_N`
(one score per gene in the eval set, 2,258 genes total)

Metrics computed: MR, MRR, Hits@1, Hits@3, Hits@10, Hits@100, AUC (macro).

### Baseline results (2,258 genes in pool)

**Full dataset — 5,076 cases**

| Prioritiser | MR | MRR | Hits@1 | Hits@3 | Hits@10 | Hits@100 | AUC |
|---|---|---|---|---|---|---|---|
| HiPhive | 500.5 | 0.263 | 0.225 | 0.272 | 0.333 | 0.512 | 0.779 |
| Phive | 826.9 | 0.035 | 0.018 | 0.031 | 0.063 | 0.222 | 0.634 |
| PhenIX | 567.3 | 0.260 | 0.225 | 0.270 | 0.319 | 0.483 | 0.749 |

**Test split — 518 cases**

INDIGENA-* rows use INDIGENA embedding BMA as a drop-in replacement for the
IC-based similarity in each method, with identical gene-phenotype model associations.

| Method | MR | MRR | Hits@1 | Hits@3 | Hits@10 | Hits@100 | AUC |
|---|---|---|---|---|---|---|---|
| HiPhive | 509.7 | 0.269 | 0.228 | 0.280 | 0.344 | 0.492 | 0.775 |
| INDIGENA-HiPhive | 506.2 | 0.268 | 0.226 | 0.278 | 0.347 | 0.510 | 0.777 |
| PhenIX | 552.5 | 0.261 | 0.230 | 0.266 | 0.320 | 0.461 | 0.756 |
| INDIGENA-PhenIX | 503.2 | 0.268 | 0.226 | 0.278 | 0.347 | 0.510 | 0.778 |
| Phive | 847.2 | 0.026 | 0.010 | 0.025 | 0.048 | 0.193 | 0.626 |
| INDIGENA-Phive | 886.7 | 0.017 | 0.010 | 0.017 | 0.027 | 0.058 | 0.608 |
| INDIGENA G4 (HPO only) | 761.2 | 0.116 | 0.079 | 0.118 | 0.191 | 0.375 | 0.664 |
| INDIGENA G4 + MGI | 997.8 | 0.019 | 0.012 | 0.017 | 0.027 | 0.069 | 0.559 |

The gene pool is fixed to all 2,258 genes in both cases so results are directly comparable. Small differences between full and test split are sampling variance.

### Track 2 results (226 genes in pool)

Track 2 combines variant pathogenicity scores from Exomiser's full pipeline (gnomAD
MAF filter + REVEL/MVP pathogenicity) with phenotype similarity.
Each spiked VCF contains ~94K GIAB HG001 background variants plus the known causal variant.
The gene pool is restricted to the 226 unique causal genes in the test split.

To prepare the VCFs (fixes the non-standard REF="." format of the spiked causal
variants so HTSJDK can parse them):
```bash
~/miniforge3/envs/indigena/bin/python eval/fix_spiked_vcfs.py test
```

To run Track 2 (first time — includes ~10 min Exomiser CLI phase):
```bash
nohup ~/miniforge3/envs/indigena/bin/python eval/track2_eval.py \
    --phenotype-data-dir exomiser-data/2406_phenotype \
    --app-props exomiser-data/application.properties \
    --split test \
    --embeddings data/models/indigena_track1_graph4_embeddings.tsv \
    > data/results/track2_eval.log 2>&1 &
```

To add/re-run phenotype scoring without re-running the CLI phase:
```bash
nohup ~/miniforge3/envs/indigena/bin/python eval/track2_eval.py \
    --phenotype-data-dir exomiser-data/2406_phenotype \
    --app-props exomiser-data/application.properties \
    --split test \
    --embeddings data/models/indigena_track1_graph4_embeddings.tsv \
    --skip-cli \
    > data/results/track2_eval.log 2>&1 &
```

**Test split — 284 cases, 226 genes in pool**

Combined score = `EXOMISER_GENE_VARIANT_SCORE × phenotype_score` for all methods
except HiPhive, which uses Exomiser's native `EXOMISER_GENE_COMBINED_SCORE` directly.

| Method | MR | MRR | Hits@1 | Hits@3 | Hits@10 | Hits@100 | AUC |
|---|---|---|---|---|---|---|---|
| HiPhive | 20.1 | 0.580 | 0.415 | 0.669 | 0.877 | 0.905 | 0.914 |
| INDIGENA-HiPhive | 11.9 | 0.468 | 0.271 | 0.479 | 0.926 | 0.951 | 0.951 |
| PhenIX | 26.5 | 0.531 | 0.391 | 0.616 | 0.789 | 0.870 | 0.886 |
| INDIGENA-PhenIX | 24.2 | 0.459 | 0.268 | 0.665 | 0.817 | 0.891 | 0.896 |
| Phive | 34.9 | 0.338 | 0.134 | 0.451 | 0.768 | 0.831 | 0.849 |
| INDIGENA-Phive | 9.7 | 0.405 | 0.144 | 0.465 | 0.940 | 0.972 | 0.961 |

INDIGENA methods consistently improve top-10 and top-100 recall at the cost of
top-1 precision. INDIGENA-Phive achieves the best MR (9.7) and Hits@10 (0.940)
despite Phive being the weakest baseline — the embedding space already encodes
cross-species HP↔MP similarity without explicit cross-species lookup.

---

## INDIGENA training

INDIGENA uses TransD KGE embeddings on the UPheno ontology graph, augmented with
gene–phenotype and disease–phenotype associations derived from the training split.
At inference, each test case is scored against all eval genes using BMA
(Best-Match Average) over phenotype embeddings — no disease entity is required
in the graph at test time (inductive).

### Graph structures

| Graph | Additional edges | Notes |
|-------|-----------------|-------|
| G1 | UPheno OWL2VecStar projection | HP + MP + UPHENO taxonomy |
| G2 | gene → MP/HP | Primary: MGI mouse-knockout phenotypes (human ortholog → MP terms); fallback: PAVS training case HPO terms for genes with no mouse ortholog |
| G3 | disease → HP (training cases) | inductive: test/val diseases excluded |
| G4 | gene ↔ disease (training cases) | supervised association signal |

**Gene–phenotype associations (Graph 2)** use MGI mouse-knockout MP phenotypes
as the primary source, mapped to human genes via the HomoloGene orthology table
(`HOM_MouseHumanSequence.rpt`). Because MP terms are embedded in the UPheno
graph alongside HP terms, BMA comparison between gene (MP-based) and case (HP-based)
phenotype embeddings is valid within the same latent space.

For human genes with no mouse ortholog in MGI, the PAVS training-case HPO terms
are used as a fallback. This keeps self-contained coverage for genes unique to
the human dataset.

This design avoids using OMIM disease–gene associations as a gene-phenotype
source, which would leak the GDA signal we are trying to predict.

### Setup

The training script requires the pre-computed UPheno OWL2VecStar edges from the
INDIGENA repo and the `indigena` conda environment on the workstation:

```bash
# On the workstation (10.74.250.168)
INDIGENA=~/Git/indigena/data

~/miniforge3/envs/indigena/bin/python eval/indigena_train.py \
    --upheno-edges $INDIGENA/upheno_owl2vecstar_edges.tsv \
    --mgi-gene-phenotypes $INDIGENA/gene_phenotypes.csv \
    --hom-file $INDIGENA/HOM_MouseHumanSequence.rpt \
    --graph2 --graph3 --graph4 \
    --track 1 --eval-split test
```

To evaluate a saved checkpoint without retraining:

```bash
~/miniforge3/envs/indigena/bin/python eval/indigena_train.py \
    --upheno-edges $INDIGENA/upheno_owl2vecstar_edges.tsv \
    --mgi-gene-phenotypes $INDIGENA/gene_phenotypes.csv \
    --hom-file $INDIGENA/HOM_MouseHumanSequence.rpt \
    --graph2 --graph3 --graph4 \
    --track 1 --eval-split test --only-eval
```

Output: `data/results/indigena_transd_track1_graph4_mgi_seed0_dim100_bs2048_lr0.001_test.tsv`

---

The Exomiser - A Tool to Annotate and Prioritize Exome Variants
===============================================================

[![GitHub release](https://img.shields.io/github/release/exomiser/Exomiser.svg)](https://github.com/exomiser/Exomiser/releases)
[![GitHub Discussions](https://img.shields.io/github/discussions/exomiser/exomiser?label=data%20release&link=https%3A%2F%2Fgithub.com%2Fexomiser%2FExomiser%2Fdiscussions%2Fcategories%2Fdata-release)](https://github.com/exomiser/Exomiser/discussions/categories/data-release)
[![CircleCI](https://circleci.com/gh/exomiser/Exomiser/tree/development.svg?style=shield)](https://circleci.com/gh/exomiser/Exomiser/tree/development)
[![Codecov](https://img.shields.io/codecov/c/github/exomiser/exomiser)](https://app.codecov.io/github/exomiser/Exomiser)
[![Documentation](https://readthedocs.org/projects/exomiser/badge/?version=latest)](http://exomiser.readthedocs.io/en/latest)
#### Overview:

The Exomiser is a Java program that finds potential disease-causing variants from whole-exome or whole-genome sequencing data.

Starting from a [VCF](https://samtools.github.io/hts-specs/VCFv4.3.pdf) file and a set of phenotypes encoded using the [Human Phenotype Ontology](http://www.human-phenotype-ontology.org) (HPO) it will annotate, filter and prioritise likely causative variants. The program does this based on user-defined criteria such as a variant's predicted pathogenicity, frequency of occurrence in a population and also how closely the given phenotype matches the known phenotype of diseased genes from human and model organism data.

The functional annotation of variants is handled by [Jannovar](https://github.com/charite/jannovar) and uses any of [UCSC](http://genome.ucsc.edu), [RefSeq](https://www.ncbi.nlm.nih.gov/refseq/) or [Ensembl](https://www.ensembl.org/Homo_sapiens/Info/Index) KnownGene transcript definitions and hg19 or hg38 genomic coordinates.

Variants are prioritised according to user-defined criteria on variant frequency, pathogenicity, quality, inheritance pattern, and model organism phenotype data. Predicted pathogenicity data is extracted from the [dbNSFP](http://www.ncbi.nlm.nih.gov/pubmed/21520341) resource. Variant frequency data is taken from the [1000 Genomes](http://www.1000genomes.org/), [ESP](http://evs.gs.washington.edu/EVS), [TOPMed](http://www.uk10k.org/studies/cohorts.html), [UK10K](http://www.uk10k.org/studies/cohorts.html), [ExAC](http://exac.broadinstitute.org) and [gnomAD](http://gnomad.broadinstitute.org/) datasets. Subsets of these frequency and pathogenicity data can be defined to further tune the analysis. Cross-species phenotype comparisons come from our PhenoDigm tool powered by the OWLTools [OWLSim](https://github.com/owlcollab/owltools) algorithm.

The Exomiser was developed by the Computational Biology and Bioinformatics group at the Institute for Medical Genetics and Human Genetics of the Charité - Universitätsmedizin Berlin, the Mouse Informatics Group at the Sanger Institute and other members of the [Monarch initiative](https://monarchinitiative.org).

#### Download and Installation

The prebuilt Exomiser binaries can be obtained from the [releases](https://github.com/exomiser/Exomiser/releases) page and supporting data files can be downloaded from the [Exomiser FTP site](http://data.monarchinitiative.org/exomiser/latest).

It is possible to use the same data sources for multiple versions, in order to avoid having to download the data files for each software point release. We recommend maintaining a dedicated exomiser data directory where you can extract versions of the hg19, hg38 and phenotype data. To do this, edit the ```exomiser.data-directory``` field in the ```application.properties``` file to point to the dedicated data directory. The version for the data releases should also be specified in the ```application.properties``` file:
    
For example, if you have an exomiser installation located at ```/opt/exomiser-cli-11.0.0``` and you have extracted the data files to the directory ```/opt/exomiser-data```. When there is a new data release, you can change the data versions by specifying the version in the ```/opt/exomiser-cli-11.0.0/application.properties``` from
```properties
# root path where data is to be downloaded and worked on
# it is assumed that all the files required by exomiser listed in this properties file
# will be found in the data directory unless specifically overridden here.
exomiser.data-directory=data

# old data versions
exomiser.hg19.data-version=1802
...
exomiser.hg38.data-version=1802
...
exomiser.phenotype.data-version=1802
```
to
```properties
# overridden data-directory containing multiple data versions
exomiser.data-directory=/opt/exomiser-data

# updated data versions
exomiser.hg19.data-version=1805
...
exomiser.hg38.data-version=1805
...
exomiser.phenotype.data-version=1807
```

We strongly recommend using the latest versions of both the application and the data for optimum results.

For further instructions on installing and running, please refer to the [online documentation](https://exomiser.readthedocs.io/en/latest/)..

#### Running it

Please refer to the [manual](https://exomiser.readthedocs.io/en/latest/) for details on how to configure and run the Exomiser.

#### Demo site

There is a limited [demo version](http://exomiser.monarchinitiative.org/exomiser/) of the exomiser hosted by the [Monarch Initiative](https://monarchinitiative.org/). This instance is for teaching purposes only and is limited to small exome analysis.

#### Using The Exomiser in your code

The exomiser can also be used as a library in Spring Java applications. Add the ```exomiser-spring-boot-starter``` library to your pom/gradle build script.

In your configuration class add the ```@EnableExomiser``` annotation
 
 ```java
@EnableExomiser
public class MainConfig {
    
}
```

Or if using Spring boot for your application, the exomiser will be autoconfigured if it is on your classpath.

```java
@SpringBootApplication
public class Application {
    public static void main(String[] args) {
        SpringApplication.run(Application.class, args);
    }
}
```

In your application use the AnalysisBuilder obtained from the Exomiser instance to configure your analysis. Then run the Analysis using the Exomiser class.
Creation of the Exomiser is a complicated process so defer this to Spring and the exomiser-spring-boot-starter. Calling the ```add``` prefixed methods 
will add that analysis step to the analysis in the order that they have been defined in your code.

Example usage:
```
@Autowired
private final Exomiser exomiser;

...
           
    Analysis analysis = exomiser.getAnalysisBuilder()
                .genomeAssembly(GenomeAssembly.HG19)
                .vcfPath(vcfPath)
                .pedPath(pedPath)
                .probandSampleName(probandSampleId)
                .hpoIds(phenotypes)
                .analysisMode(AnalysisMode.PASS_ONLY)
                .modeOfInheritance(EnumSet.of(ModeOfInheritance.AUTOSOMAL_DOMINANT, ModeOfInheritance.AUTOSOMAL_RECESSIVE))
                .frequencySources(FrequencySource.ALL_EXTERNAL_FREQ_SOURCES)
                .pathogenicitySources(EnumSet.of(PathogenicitySource.POLYPHEN, PathogenicitySource.MUTATION_TASTER, PathogenicitySource.SIFT))
                .addPhivePrioritiser()
                .addPriorityScoreFilter(PriorityType.PHIVE_PRIORITY, 0.501f)
                .addQualityFilter(500.0)
                .addRegulatoryFeatureFilter()
                .addFrequencyFilter(0.01f)
                .addPathogenicityFilter(true)
                .addInheritanceFilter()
                .addOmimPrioritiser()
                .build();
                
    AnalysisResults analysisResults = exomiser.run(analysis);
```
 
#### Memory usage

Analysing whole genomes using the ``AnalysisMode.FULL`` will use a lot of RAM (~16GB for 4.5 million variants without any extra variant data being loaded) the standard Java GC will fail to cope well with these.
Using the G1GC should solve this issue. e.g. add ``-XX:+UseG1GC`` to your ``java -jar -Xmx...`` incantation. 

#### Caching

Since 9.0.0 caching uses the standard Spring mechanisms.
 
To enable and configure caching in your Spring application, use the ```@EnableCaching``` annotation on a ```@Configuration``` class, include the required cache implementation jar and add the specific properties to the ```application.properties```.

For example, to use [Caffeine](https://github.com/ben-manes/caffeine) just add the dependency to your pom:

```xml
<dependency>
    <groupId>com.github.ben-manes.caffeine</groupId>
    <artifactId>caffeine</artifactId>
</dependency>
```
and these lines to the ```application.properties```:
```properties
spring.cache.type=caffeine
spring.cache.caffeine.spec=maximumSize=300000
```

#### Recognition

The Exomiser is proud to be recognised by the International Rare Diseases Research Consortium ([IRDiRC](http://www.irdirc.org/)) as an [IRDiRC Recognized Resource](http://www.irdirc.org/research/irdirc-recognized-resources/). This is *'a quality indicator, based on a specific set of criteria, that was created to highlight key resources which, if used more broadly, would accelerate the pace of translating discoveries into clinical applications.'* These resources *'must be of fundamental importance to the international rare diseases research and development community'*.

[![IRDiRC recognised resource](LOGO_IRDIRCRR_RVB.jpg)](http://www.irdirc.org/)
