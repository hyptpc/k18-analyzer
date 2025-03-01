k18-analyzer
============

K1.8(BR) analyzer

## Install and build

By default, usr is empty and is created by copying from example.
User part is removed from git management because it will be changed frequently by each user.
First, copy files from example to usr. And make.

```sh
git clone --branch e72 --recurse-submodules git@github.com:hyptpc/k18-analyzer.git e72
cd e72
cp Makefile.org Makefile
cp -pr example usr
make
```

## Environment variables

Example setting in .bashrc.

```sh
. /sw/packages/root/6.34.38/bin/thisroot.sh
. /sw/packages/geant4/11.2.2/bin/geant4.sh
. /sw/packages/geant4/11.2.2/share/Geant4/geant4make/geant4make.sh
conda activate myenv
export G4WORKDIR=$HOME/work/geant4 # set as you like
export PATH=/group/had/sks/software/unpacker/e70/bin:$PATH
export PATH=$G4WORKDIR/bin/Linux-g++:$PATH
```

Since PATH to unpacker is specified in Makefile, it is unnecessary to be set.

## Anaconda setting

To use Python,
it is necessary to build the Anaconda local environment once using the `conda` command as follows.
Note that it is recommended to use `conda install` instead of `pip install` in the anaconda environment.

```sh
$ conda create -n myenv python=3.9 # myenv is an example name
$ conda activate myenv
$ conda install numpy psutil pyyaml rich
```

Add the following line in .bashrc to activate your environment.

```sh
conda activate myenv
```

If the prompt header of conda is annoying, add the following line in .condarc.

```yaml
changeps1: False
```

## Run analyzer

The data directory at kekcc is as below.

| Period | Path |
| -- | -- |
| E73 in Run91 (2024) | /group/had/knucl/e15/e73_data/run91 |
| E72 in Run92 (2025) | /hsm/had/sks/E72/JPARC2025Feb/rawdata |

The following is run command.

```sh
./bin/RawHist param/conf/analyzer_e73_2024.conf /group/had/knucl/e15/e73_data/run91/run00117.dat.gz tmp.root
```

Ctrl-C will terminate the analysis.

## Use runmanager

__runmanager__ is a script for managing jobs on KEKCC.
Prepare runlist.yml by referring to the example.yml.
Note that the indentation determines the nest depth in yaml.

```yml
#
# RUN LIST (YAML format)
#
# The allowed keys are
#   queue  <- bsub queue (eg. s, l, etc...)
#   unit   <- dividing event unit
#   nporc  <- number of process for merging
#             (must be less than 20)
#   buff   <- intermediate root files will be placed here
#             if nproc is more than 2
#   bin    <- path to executable binary
#   conf   <- path to conf directory/file from the work directory
#   data   <- path to data directory/file
#   root   <- path to the output ROOT output directory/file
#   fig    <- path to figure pdf/png files only used for macro
#
# Be careful of the indent rule
# because YAML format is sensitive to it.
# Some problems will happen
# if there is no default setting declaration.
#

#____________________________________________________
# work directory path
# The following paths must be relative to this path

WORKDIR: ~/k18analyzer/pro

#____________________________________________________
# default setting
# default setting MUST have all items.
# This setting will be inherited unless you explicitly set values
# for the individual cases.

DEFAULT:
  queue: s
  unit:  100000
  nproc: 1
  buff:  /group/had/sks/Users/user/buffer
  bin:   ./bin/Hodoscope
  conf:  ./param/conf/analyzer_default.conf
  data:  ../data
  root:  ../rootfile
  fig:   ../fig/hodo

#____________________________________________________
# Individual settings
# If you want to adapt the default setting to some runs,
# you only have to list keys.

RUN:

  test: # Any keys are OK unless it overlaps.
    bin:  ./bin/Hodoscope
    conf: ./param/conf/analyzer_default.conf
    data: ./test.dat.gz
    root: ./test.root

  3838:
  # The default setting will be adapted.
  # In this case run03838.dat(.gz) will be loaded.

  3800: # ./param/conf/analyzer_03800.conf will be loaded.
    bin:  ./bin/KuramaTracking
    conf: ./param/conf
    data: ./tmp_data

  hodo:
    data: ../data/run03800.dat.gz
```

Run run.py.
Note that the process runs until all jobs have finished.
Ctrl-C kills all jobs and terminates the process.

```sh
./runmanager/run.py runmanager/runlist/foo.yml
```

Run monitor.py on another tty to see the progress of the jobs.
The job status is updated in the "stat" directory, using the same name as the runlist in json format.

```sh
./runmanager/monitor.py runmanager/stat/foo.json
```

## Macro

Macro to list RawHist histograms in single pdf.

```sh
python3 ./macro/rawhist.py runmanager/runlist/foo.yml
```
