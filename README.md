---
stylesheet: https://cdnjs.cloudflare.com/ajax/libs/github-markdown-css/2.10.0/github-markdown.min.css
body_class: markdown-body
pdf_options:
  format: A4
  margin: 30mm 20mm
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

2020.03.12 shuhei hayakawa

It is assumed to work on KEKCC, Scientific Linux release 6.10 (Carbon).

## Environment variables

Add the following lines to your .bashrc.

```a
module load git/2171
module load gcc/485
module load python/2.7
module load python/3.5
. /sw/packages/root/6.14.06/bin/thisroot.sh
export PATH=$PATH:$HOME/unpacker/bin # Set as appropriate
```

## Unpacker

Install the Unpacker.

```a
$ git clone \
ssh://sks@www-online.kek.jp:8022/~/public_html/git/unpacker.git
$ cd unpacker/src
$ cp Makefile.org Makefile
$ make
```

Check if the "unpacker-config" command is available.

```a
$ unpacker-config --version
2020-01-21
```

## K1.8 analyzer

Install the K1.8 analyzer.

```a
$ git clone \
ssh://sks@www-online.kek.jp:8022/~/public_html/git/k18-analyzer.git
$ cd k18-analyzer
$ git checkout e42
$ cp Makefile.org Makefile
$ make
```

e.g.) Hodoscope,  
Usage: Hodoscope [analyzer config file] [data input stream] [output root file]

```a
$ ./bin/Hodoscope param/conf/analyzer_2019apr_0.conf \
/group/had/sks/E40/JPARC2019Feb/e40_2019feb/run07334.dat.gz hoge.root
```
