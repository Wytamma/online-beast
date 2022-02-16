# online-BEAST
[![PyPi](https://img.shields.io/pypi/v/online-beast.svg)](https://pypi.org/project/online-beast/)
[![tests](https://github.com/Wytamma/online-beast/actions/workflows/test.yml/badge.svg)](https://github.com/Wytamma/online-beast/actions/workflows/test.yml)
[![cov](https://codecov.io/gh/Wytamma/online-beast/branch/master/graph/badge.svg)](https://codecov.io/gh/Wytamma/online-beast)

This command line tool can be used to add sequences to an ongoing analysis in BEAST2. This framework is called online Bayesian phylodynamic inference (see [Gill et al., 2020](https://academic.oup.com/mbe/article/37/6/1832/5758268?login=false)).

## Install
Install `online-beast` with pip (requires python -V >= 3.6.2).

```bash
pip install online-beast
```

## Usage 

Give `online-beast` beast the path to a XML file from a previous BEAST2 run (i.e. one that has an associated `.state` file) and a fasta file of sequences to add to the analysis. Sequences in the fasta file must be aligned (i.e. to the sequences in the XML file) and the same length as the other sequences in the XML file. Only new sequences (new descriptors) will be added to the analysis, so new sequences can be append to the fasta file as they are acquired. 

```bash
online-beast data/testGTR.xml data/samples.fasta
```

![](images/output.png)

The new sequences will by added to the XML file and the associated `.state` file (produced automatically by BEAST2).

The analysis can then be resumed (with the additional sequence data) using the BEAST2 resume flag. 

```bash
beast -resume testGTR.xml
```

The online analysis can be visualised in real-time using [Beastiary](https://beastiary.wytamma.com/). The jumps in the trace show where new sequences have been added. 

![](images/beastiary.png)

Date trait data will be automatically parsed. The format of the date trait data (in the fasta descriptor) can be set with the `--date-format` (default `%Y-%m-%d`) and `--delimiter` (default `_`) flags. If there is no date trait in the xml use the `--no-date-trait` flag.

```
online-beast data/ebola.xml data/ebola.fasta --dateformat %d/%m/%Y --date-delimiter _
```

If there is trait data in the XML file you need to specify how to extract it from the fasta descriptor line using the `--trait` flag. The format is `'traitname delimiter group'` e.g. a string separated by spaces. For example to get the `location` trait from `sample_wuhan_2022-04-05` you would use `--trait 'location _ 1'`. The `--trait` flag can be used multiple times to specify multiple traits. 

```bash
online-beast covid.xml data/covid.fasta --trait 'location _ 1'
```

By default the new sequences will be appended to the input XML and state files. Output file names can be specified using the `--output` flag. This will also create a new `.state` file.

```bash
online-beast testGTR.xml samples.fasta --output new_testGTR.xml 
```

If you use the BEAST2 `-statefile` flag to specify the filename of the state (i.e. it is not `xml_filename + .state`). Use the flag `--state-file` to specify the state file path. 

```bash
online-beast testGTR.xml samples.fasta --state-file beast.state 
```

## Explanation

A Markov chain started anywhere near the center of the stationary distribution needs no burn-in ([Geyer 2011](http://www.mcmchandbook.net/)). Online Bayesian phylodynamic inference is akin to transfer learning in the deep learning field. By starting our MCMC with reasonable states (obtained from a previous run) we reduce the amount of optimisation (burn-in) that must be performed to reach convergence. 

Online-beast loosely follows the implementation of [Gill et al., 2020](https://academic.oup.com/mbe/article/37/6/1832/5758268?login=false) for BEAST1. However, most of the implementation of online-beast is handled by the default state system in BEAST2. New sequences are added from the fasta file one at a time. The hamming distance is calculated between the new sequence and all the other sequences in the XML file. The new sequence is grafted onto the tree in the `.state` file, half way along the branch of the closest sequence in the XML file. The new sequence is append to the BEAST XML file. 

## Ebola example

In this example we will make use of a publicly available dataset of sequences from the 2013-2016 *Zaire ebolavirus* outbreak in Sierra Leone. 

In the `data/` folder you'll find a `ebola.xml` file and several fasta files that contain sequences from the outbreak broken up by date. The script below will run an online Bayesian phylodynamic analysis adding new sequences after each run finishes. 


```bash
#!/bin/bash

# Run beast with initial samples
beast data/ebola.xml 
# Update analysis with new samples
online-beast data/ebola.xml data/ebola1.fasta --date-format "%d/%m/%Y" --state-file ebola.xml.state --output ebola.xml
# Resume the analysis
beast -resume ebola.xml 
# Update analysis with new samples
online-beast ebola.xml data/ebola2.fasta --date-format "%d/%m/%Y" --output ebola.xml
# Resume the analysis
beast -resume ebola.xml 
# Update analysis with new samples
online-beast ebola.xml data/ebola3.fasta --date-format "%d/%m/%Y" --output ebola.xml
# Resume the analysis
beast -resume ebola.xml 
```




