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

Environment setting in KEKCC
============================

<div style="text-align: right;">
Rev. 2024.10.24<br>
Rev. 2023.06.01<br>
2020.09.06
</div><br>

The new KEKCC is Red Hat Enterprise Linux release 9.3 (Plow).
The default compiler is gcc 11.4.1.

## module

The same gcc 8.3.0 as the old system can be used with module,
but is basically unnecessary.

```sh
$ module av

--------------------------- /opt/Modules/modulefiles ---------------------------
gcc/850  gcc/1230  intel/2024  nvhpc/24.1
```

## /sw/packages

You can use the software installed in /sw/packages.
Everything is installed with the default compiler.

| package                |
| :--------------------------- |
| /sw/packages/root/5.34.38    |
| /sw/packages/root/6.32.04    |
| /sw/packages/geant4/9.6.4    |
| ...                          |
| /sw/packages/geant4/11.2.2   |

The ROOT on the group directory is no longer used.

## unpacker

The unpacker is placed in the group directory,
/group/had/sks/software/unpacker/e70
This is updated constantly.

See below for local installation.

```sh
$ git clone ssh://sks@www-online.kek.jp:8022/~/public_html/git/unpacker.git
$ cd unpacker/src
$ cp Makefile.org Makefile
$ make
```

Check if the `unpacker-config` command is available.

```sh
$ unpacker-config --version
2024-07-03
```

## Environment variables

The following is an example of environment setting in .bashrc.

```sh
. /sw/packages/root/6.32.04/bin/thisroot.sh
. /sw/packages/geant4/11.2.2/bin/geant4.sh
. /sw/packages/geant4/11.2.2/share/Geant4/geant4make/geant4make.sh
conda activate myenv
export G4WORKDIR=$HOME/work/geant4 # set as you like
export PATH=/group/had/sks/software/unpacker/e70/bin:$PATH
export PATH=$G4WORKDIR/bin/Linux-g++:$PATH
```

## Anaconda setting

To use Python,
it is recommended to build the Anaconda local environment once using the `conda` command as follows.
Note that it is recommended to use `conda install` instead of `pip install` in the anaconda environment.
Newer versions of python, such as 3.12, are available,
but to run PyROOT under /sw/package, the python version must be default, 3.9.

```sh
$ conda create -n myenv python=3.9 # myenv is an example name
$ conda activate myenv
$ conda install numpy psutil pyyaml rich # and other modules you like
```

Add the following line in .bashrc to activate your environment.

```sh
conda activate myenv
```

If the prompt header of conda is annoying, add the following line in .condarc.

```yaml
changeps1: False
```

## K1.8 analyzer

Install the K1.8 analyzer as follows.

```sh
$ git clone \
ssh://sks@www-online.kek.jp:8022/~/public_html/git/k18-analyzer.git
$ cd k18-analyzer
$ git checkout e70 # choose branch as you like
$ cp Makefile.org Makefile
$ make
```

See README for more information on how to use the analyzer.
