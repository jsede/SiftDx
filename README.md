# SiftDx
SiftDx: A clinical metagenomics pipeline for sifting through complex sequencing data to detect pathogens with precision

## 🚧 Pipeline Under Construction 🚧

## Installation
Please follow the guide in the ![wiki](https://github.com/jsede/SiftDx/wiki/Installation-Guide)

## Usage
Full pipeline in DNA/RNA mode, just swap out the `--na DNA` in the code below to `--na RNA`

```
nextflow run main.nf -resume --na DNA --r1 PATH_TO_R1 --r2 PATH_TO_R2 --output PATH_TO_OUTPUT -with-trace
```

Full pipeline in RNA mode with ERCC controls, if you use MetaSequins for DNA spike-ins, use the `--sequins` flag instead

```
nextflow run main.nf -resume --na RNA --r1 PATH_TO_R1 --r2 PATH_TO_R2 --output PATH_TO_OUTPUT -with-trace --ercc
```

If you have a negative control, we'd advise you run the pipeline on it first using any of the above commands then run your sample with the negative control using the following.

```
nextflow run main.nf -resume --na RNA --r1 PATH_TO_R1 --r2 PATH_TO_R2 --output PATH_TO_OUTPUT -with-trace --ndata PATH_TO_NEGATIVE
```
Note: the `PATH_TO_NEGATIVE` is the path to the output, rather than the work directory

### TO-DO
- A preprocessing only mode