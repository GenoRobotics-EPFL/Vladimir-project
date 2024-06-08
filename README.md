### Summary

This repo contains code to analyse the sequence of reads from the minION sequencer in real time.
A pipeline will create a consensus for one gene, outputing the quality of the consensus and the result of the detection.
Two pipelines were implemented: the naive pipeline and the bestx pipeline

### Usage

To launch the pipeline, go to the root directory of the repo and execute the pipeline:

```console
foo@bar:~$ python3 <pathToPipelineFile> <nameOfGene>
```

where `<pathToPipelineFile>` is either `src/pipelineBestSequence.py` (for the best-x pipeline) or `src/pipelineNaïve.py` (for the naive pipeline), and `<nameOfGene>` is either `matK` or `rbcL` or `psbA-trnH` or `ITS`.
For instance:
```console
foo@bar:~$ python3 src/pipelineBestSequence rbcL
```

The pipelines will take the reads that are in the folder "fastqpass". The pipelines will take each file as one "iteration" of reads coming from the sequencer. In order to populate that folder, the following script can be used:
```console
foo@bar:~$ python3 src/simulateRealTimeOutput.py <pathToFastqFile> <minutesBeforeNextFile>
```
where `<minutesBeforeNextFile>` must be a positive integer that represents the number of minutes the script will sleep before creating the next fastq file. This is done so that we can simulate the sequencer outputting reads at every time interval. If 0 is supplied then the script will create the files directly.
for instance:
```console
foo@bar:~$ python3 src/simulateRealTimeOutput.py allData/Allium_Ursinum_ITS.fastq 0
```
The script can be used in the background to simulate the sequencer creating a fastq file every minute.

### output

The pipelines will output all their results in their respective output file `outputPipelineBest` and `outputPipelineNaive`.
The main output file `result.txt` contains the detection result of each iteration (and the consensus). The folder will also contain graph of the depth covreage for instance.

stdout will show the progression of the pipeline.

### Installation

The first program that the pipelines uses is [Medaka](https://github.com/nanoporetech/medaka). Follow their installation instructions on their repo. The method that worked for me was using the conda channel.
To test that the installation works you should be able to execute
```console
foo@bar:~$ medaka_consensus -h
```

Another program required is [SPOA](https://github.com/rvaser/spoa). I included a binary in the src/ directory, but maybe it won't work on your machine so you can follow their instructions to compile it from source, and replace the binary of this repo (must be placed in the src/ directory)
You should be able to execture
```console
foo@bar:~$ ./src/spoa --version
```

Another program is [Mosdepth](https://github.com/brentp/mosdepth). Same thing as for spoa
You should be able to execture
```console
foo@bar:~$ ./src/mosdepth --version
```

The last program is blastn. Follow the instructions in the Identification repo. You must be able to execute:
```console
foo@bar:~$ blastn -version
```
In addition, you need 4 databases, one for each gene mentioned above. Their name must be exactmu the same as the geneName used as parameter for the pipeline. For instance you should be able to execture:
```console
foo@bar:~$ blastdbcmd -db rbcL -info
```