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

Give `online-beast` beast the path to a xml file from an previous BEAST run (i.e. one that have been stopped/killed/crashed) and a fasta of sequence to add to the analysis. Sequences in the fasta file must be aligned and the same length as the other sequences in the XML file. Only new sequences (new descriptors) will be added to the analysis, so new sequences can be append to the fasta file as they are acquired. 

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

By default the new sequences will be appended to the input XML and Sate files. Output file names can be specified using the `--output` flag. This will also create a new `.state` file.

```bash
online-beast testGTR.xml samples.fasta --output new_testGTR.xml 
```

If you use the BEAST2 `-statefile` flag to specify the filename of the state (i.e. it is not `xml_filename + .state`). Use the flag `--state-file` to specify the state file path. 

```bash
online-beast testGTR.xml samples.fasta --state-file beast.state 
```

## Explanation

Online-beast loosely follows the implementation of [Gill et al., 2020](https://academic.oup.com/mbe/article/37/6/1832/5758268?login=false) for BEAST1. However, most of the implementation of online-beast is handled by the default state system in BEAST2. New sequences are added from the fasta file one at a time. The pairwise distance is calculated between the new sequence and all the other sequences in the XML file. The new sequence is grafted onto the tree in the `.state` file, half way along the branch of the closest sequence in the XML file. The new sequence is append to the BEAST XML file. 






