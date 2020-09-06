---
stylesheet: https://cdnjs.cloudflare.com/ajax/libs/github-markdown-css/2.10.0/github-markdown.min.css
body_class: markdown-body
pdf_options:
  format: A4
  margin: 24mm 16mm
  displayHeaderFooter: true
  headerTemplate: |-
    <style>
      section {
        margin: 0 auto;
        font-family: system-ui;
        font-size: 11px;
      }
    </style>
    <section>
    </section>
  footerTemplate: |-
    <section>
      <div>
        <span class="pageNumber"></span>
        <!-- /<span class="totalPages"></span> -->
      </div>
    </section>
---

K1.8 analyzer README
====================

<div style="text-align: right;">
2020.09.06
</div><br>

It is assumed to work on KEKCC, CentOS 7.7.1908 (Core).
Use gcc 8.3.0 and root 6.22.02.

## Environment setting

See KEKCC.pdf and set environment variables.
Anaconda/Python environment is needed to use runmanager.

## Unpacker

The unpacker compiled with gcc 8.3.0 is placed in the group directory, /group/had/sks/software/unpacker/unpacker.gcc830. This is updated constantly.

If you want to install in local, install as follows.

```sh
$ git clone ssh://sks@www-online.kek.jp:8022/~/public_html/git/unpacker.git
$ cd unpacker/src
$ cp Makefile.org Makefile
$ make
```

Check if the "unpacker-config" command is available.

```sh
$ unpacker-config --version
2020-01-21
```

## K1.8 analyzer

Install the K1.8 analyzer.

```sh
$ git clone \
ssh://sks@www-online.kek.jp:8022/~/public_html/git/k18-analyzer.git
$ cd k18-analyzer
$ git checkout e40 # choose branch as you like
$ cp Makefile.org Makefile
$ make
```

e.g.) Hodoscope,
Usage: Hodoscope [analyzer config file] [data input stream] [output root file]

```sh
$ ./bin/Hodoscope param/conf/analyzer_2019apr_0.conf \
/group/had/sks/E40/JPARC2019Feb/e40_2019feb/run07334.dat.gz hoge.root
```

### runmanager

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
$ ./runmanager/run.py runmanager/runlist/foo.yml
```

Run monitor.py on another tty to see the progress of the jobs.
The job status is updated in the "stat" directory, using the same name as the runlist in json format.

```sh
$ ./runmanager/monitor.py runmanager/stat/foo.json
```

![runmanager.png](runmanager.png "caption")
