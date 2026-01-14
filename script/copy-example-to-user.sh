#!/bin/bash

work_dir=$(readlink -e $(dirname $0)/..)
example_dir=$work_dir/example
usr_dir=$work_dir/usr
dst_dir=$work_dir/dst
example_files="$example_dir/*.cc"

mkdir -p $usr_dir
for file in $example_dir/User*.cc; do
  [ -e "$file" ] || continue
  cp -avi "$file" "$usr_dir/"
done

mkdir -p $dst_dir
for file in $example_dir/Dst*.cc; do
  [ -e "$file" ] || continue
  cp -avi "$file" "$dst_dir/"
done