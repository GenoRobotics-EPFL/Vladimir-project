### Summary

This repo contains code to analyse the sequence of reads from the minION sequencer in real time.
A pipeline will create a consensus for one gene, outputing the quality of the consensus and the result of the identification.
Two pipelines were implemented: the naive pipeline and the bestx pipeline

### Usage

To launch the pipeline, go to the root directory of the repo and execute the pipeline:

```console
foo@bar:~$ python3 <pathToPipelineFile> <nameOfGene>
```

where "<pathToPipelineFile>" is either "src/pipelineBestSequence.py" (for the best-x pipeline) or "src/pipelineNaïve.py (for the naive pipeline), and "<nameOfGene>" is either "matK" or "rbcL" or "psbA-trnH" or "ITS".
It is important to note that the genename is the one that corresponds to the reads given to the pipeline. ALso, that name will be used as name of the databest in the blastn query, so that must match with the db installed locally. If other names are used for the blastn databases, those db must be renamed

The pipeline will take the reads that are in the folder "
