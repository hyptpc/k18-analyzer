k18-analyzer
============

K1.8(BR) analyzer

## Quick start

```sh
cd e73_2024
sh script/copy-example-to-user.sh
make
./bin/RawHist param/conf/analyzer_e73_2024.conf /group/had/knucl/e15/e73_data/run91/run00117.dat.gz tmp.root
```

## Use runmanager

```sh
./runmanager/run.py runmanager/runlist/foo.yml
```

```sh
./runmanager/monitor.py runmanager/stat/foo.json
```

## Macro

For RawHist

```sh
python3 ./macro/rawhist.py runmanager/stat/foo.json
```
