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

K1.8 analyzer
=============

## Install

```sh
$ git clone -r -b tohoku-counter git@github.com:hyptpc/k18-analyzer.git <dir-name>
$ cd <dir-name>
$ cp Makefile.org Makefile
$ make
```

## Usage

The arguments are `TohokuBench [analyzer config file] [data input stream] [output root file]`.

```sh
$ ./bin/TohokuBench \
param/conf/analyzer_20241029.conf \
/group/had/sks/E72/tohoku-counter-2024dec/data/run00118.dat \
foobar.root
```
