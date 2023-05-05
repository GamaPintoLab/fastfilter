# FastFilter

Fast, multi-threaded filtering of FASTQ files with adjustable parameters written in Python.

---

Author: [Gil Poiares-Oliveira](mailto:gpo@ciencias.ulisboa.pt) |
PI: [Margarida Gama-Carvalho](mailto:mhcarvalho@ciencias.ulisboa.pt)

[RNA Systems Biology Lab](https://rnasysbio.rd.ciencias.ulisboa.pt)\
BioISI - Biosystems and Integrative Sciences Institute, Faculty of Sciences, University of Lisbon

---

## Introduction

FastFilter allows you to filter out:
* Homopolymers
* Sequences below a given quality threshold
* Sequences below a given length
* Sequences with `N` or `.` characters

And generate CSV reports to add to your paper.

## Guided Mode

FastFilter can be run in "guided mode" by not specifying any runtime parameters,
i.e. just running:
```sh
fastfilter
```
For this mode to run properly, you have to first assure:
1. Projects are separated by folder and are all contained in the folder
   specified in `PROJECTS_DIR_DEFAULT`, you should edit the script's source code
   to change the string inside the `Path()` constructor to the folder in which
   you save your projects.
2. Your project's folder has a child folder called `cutadapt` which contains
   the sequences that will be filtered.

After assuring these criteria, you can run the script, which will prompt you
which of the folders inside `PROJECTS_DIR_DEFAULT` corresponds to your project,
the script will then create a new folder inside it called `fastfilter` and place
all the output files inside it.

If these criteria aren't met, you can't run the script in Guided Mode, and have
to then specify the input and output folders at runtime, by supplying the `-i`
and `-o` arguments, respectively.

## Options

To list all arguments available, run:
```sh
fastfilter -h
```

`-l` 	Length threshold for sequence

`-s`	Quality score threshold for sequence

`-p` 	Number of repeated adjacent nucleotides from which a sequence is
		considered a homopolymer

`-d`	Dry run (does not export files, primarily for debugging)

`-i` 	Folder which contains the FASTQ sequences for filtering

`-o` 	Folder in which output files will be placed (will be created if it
		doesn't exist)

## Multithreading

This script uses multithreading (running processes in parallel using different
CPU threads) to increase efficiency. The program allocates one CPU thread per
pair of pair-ended FASTQ files.

## Examples

Filter sequences with quality score lower than 50:
```sh
fastfilter -s 50
```

Place output files in folder with path `/foo/bar/fastfilter-output`:
```sh
fastfilter -o "/foo/bar/fastfilter-output"
```
