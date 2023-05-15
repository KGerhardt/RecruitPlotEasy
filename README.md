# RecruitPlotEasy

A tool for interactive Recruitment Plot generation and viewing

# Requirements

- Python >= 3.9
- Pyrodigal >= 2.0
- Plotly
- Numpy 


# Installation

RecruitPlotEasy and its dependencies are installed via pip. Simply run:

```bash
pip install RecruitPlotEasy
```

# Use

RecruitPlotEasy has two interchangable ways to access its behavior:

```bash
recruitploteasy

rpe
```

There are no differences between these methods of calling RecruitPlotEasy - rpe is just a convenient shorthand.

Regardless of which prefix you use, you then have to select an action. The two actions are "build" and "plot". Use looks like this:

```bash
rpe build

rpe plot
```

Each has its own set of options that will print to the screen if you use either action without any other arguments, or add --help. 

## Build

The "build" module controls the creation and modification of a RecruitPlotEasy database. All options in build require a target -d / --database argument to specify a path to a database file to create or modify. Reads are added with the -r / --reads argument, while genomes/MAGs are added with -g / --genome. Only reads are strictly necessary to build and plot with RecruitPlotEasy.

Genomes may be further modified by supplying the --mag and --predict flags. The --mag flag indicates that the file supplied with --genome is a MAG, and that all of the sequences contained within are to be concatenated and treated as part of a single genome. The --predict flag enables protein visualization behavior in the plotting module of RecruitPlotEasy by predicting gene sequences for the current --genome file and storing them in the database.

To add additional reads or genomes, just specify an already existing database with --database and add the reads or genomes with --reads or --genome as desired.

Sequences added with --genome will always overwrite data in a RecruitPlotEasy database, while data added with --reads never will. If the same sequence ID appears in two --genome file additions, the last file added will have its data reflected in the database. This extends to MAG goup and to predicted proteins. Please use caution: there are no warnings for this.

```bash
#example build
rpe build -d my_db.db -r reads.sam -g genomes.fna --mag
```

Planned features:

* Adding multiple reads, genomes at once
* Specifying MAG identity for multiple MAGs at once

## Plot

RecruitPlotEasy's plotting functionality is designed to produce a recruitment plot for every single genome/MAG in every single sample in the database. A directory called "recruitment_plots" will be created in your working directory, which will contain a subdirectory for each sample. Inside each sample directory, you will find a recruitment plot HTML file for each genome/MAG.

Three decisions have to be made at the time when you are creating plots: genome bin width, percent identity bin height, and protein binning.

Genome bin width, controlled by -w / --width determines how many regions a genome will be broken into. RecruitPlotEasy attempts to evenly divide contigs/genomes into segments of --width base pairs. A contig of length 3500 with the default width of 1000 would be divided into 3 bins of approximately 1166 bases each. This behavior is overwritten by the use of -p / --proteins, which determines genome bins at protein boundaries instead.

Percent identity bin height, controlled by -i / --id_step, determines how many percent identity bins will be created. The range always spans 100 to 70 percent identity (inclusive), and the ID bin height determines how many steps will exist along this range. The default 0.5% identity step divides this range into 61 bins at 100, 99.5, 99.0, ... 70.0.

Last, -p / --proteins instructs RecruitPlotEasy to plot proteins as the genome bins. Each protein will be a bin unto itself based on its start and stop index in the genome, and intergenic regions will be created to span any gaps in the genome uncovered by a protein. If proteins overlap, their overlap will be divided in half, and each protein will be given half of the overlapping region. -w / --width is ignored with proteins are being plotted. You cannot use --proteins if you did not add genomes with the --predict flag during database building.

The recruitment plots produced by RecruitPlotEasy are otherwise completely standalone. They are simple (albeit large) HTML files and will open in any web browser without any other kind of installation. They can be passed to other uses who do not have RecruitPlotEasy installed (or even Python) and they will still function identically.

```bash
#example plot
rpe plot -d my_db.db -w 3000
```

Planned features:

* A filter for genomes/MAGs with insufficient coverage to prevent pointless plots.

## Reading a Recruitment Plot

A recruitment plot conists of one main panel and 5 additional panels that either summarize or annotate the data within the main panel. Check out the How-To-Read slide deck to learn all about it!

https://docs.google.com/presentation/d/1XB7J-XDCLsvIUKzjw32PMktL2fXzNQOH61rjSCP309A/edit?usp=sharing
