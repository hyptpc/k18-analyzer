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

RAVE INSTALLATION README
====================

<div style="text-align: right;">
2023.01.05
</div><br>

This document is for RAVE installation and linking of it with current K1.8 Analyzer

## Prerequesites
RAVE toolkit requires Boost and CLHEP library.
If not, please download thise prerequesites.	

Boost library is available in KEKCC. 
CLHEP version succeeded in RAVE install : 2.4.6.2

## RAVE INSTALLATION
Unzip rave*.tar.gz into somewhere and run ./boostrap and ./configure.

Please set environment variables for CLHEP. Otherwise configure it with 

```sh
$ ./configure \
	    --with-clhep=/your/path/to/CLHEP 
	    --with-clhep-libpath=/your/path/to/CLHEPLIB
	    --with-clhep-incpath=/your/path/to/CLHEPINC
	    --prefix=/your/path/to/install
	    --disable-java 
	    CXXFLAGS=-std=c++11 
	    LDFLAGS=-L/usr/lib64
```

More details and options are available on ./configure --help

## Linking it with GenFit and K1.8Analyzer

Sample Makefile is summarized in MakeFile.gfRave.
Please modify RAVE_DIR in Makefile to your path to RAVE Installation.

## Environment variables

Path to RAVE and CLHEP is required.
For example, add some lines to bash_profile.

```sh
$ export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:~/path/to/CLHEP/lib
$ export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:~/path/to/RAVE/lib
```
