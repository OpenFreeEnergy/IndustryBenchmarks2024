#!/bin/zsh

for file in network_setup/transformations/*.json; do
  relpath=${file:14}  # strip off "network_setup/"
  dirpath=${relpath%.*}  # strip off final ".json"
  jobpath="network_setup/${dirpath}.job"
  for repeat in {0..2}; do
    cmd="openfe quickrun $file -o results_${repeat}/$relpath -d results_${repeat}/$dirpath"
    echo -e "#!/usr/bin/env bash\n${cmd}" > $jobpath
    cat $jobpath
  done
done
