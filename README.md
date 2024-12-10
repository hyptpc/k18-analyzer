K1.8 analyzer
=============

## Install

```sh
git clone -r -b tohoku-counter git@github.com:hyptpc/k18-analyzer.git <dir-name>
cd <dir-name>
cp Makefile.org Makefile
make
```

## Usage

The arguments are `TohokuBench [analyzer config file] [data input stream] [output root file]`.

```sh
./bin/TohokuBench param/conf/analyzer_20241029.conf /group/had/sks/E72/tohoku-counter-2024dec/data/run00118.dat foobar.root
```
