#!/bin/bash

work_dir=$(readlink -e $(dirname $0)/..)
example_dir=$work_dir/example
usr_dir=$work_dir/usr
example_files="$example_dir/*.cc"

mkdir -p $usr_dir
for example_file in $example_files
do
  usr_file=`echo $example_file|sed "s/example/usr/g"|sed "s/Example/User/g"`
  # echo $example_file $usr_file
  cp -v $example_file $usr_file
done
