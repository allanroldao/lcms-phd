#!/bin/bash

# Convert files to Unix format
tr -d '\r' < ExpMW.txt > ExpMW_unix.txt
tr -d '\r' < PoolMap.txt > PoolMap_unix.txt

# Prepare ExpMW_unix.txt for awk by creating a lookup table
awk '
BEGIN { FS="[ \t]+"; OFS="\t" }
{
  id_to_mw[$1] = $2
}
END {
  for (id in id_to_mw) {
    print id, id_to_mw[id]
  }
}' ExpMW_unix.txt > ExpMW_lookup.txt

# Read PoolMap_unix.txt and create the final dataset
awk '
BEGIN { FS="[ \t]+"; OFS="\t" }
FNR==NR {
  id_to_mw[$1] = $2
  next
}
{
  poolID = $1
  ids[1] = $2
  ids[2] = $3
  ids[3] = $4
  ids[4] = $5
  ids[5] = $6
  ids[6] = $7
  mw_values[1] = id_to_mw[ids[1]]
  mw_values[2] = id_to_mw[ids[2]]
  mw_values[3] = id_to_mw[ids[3]]
  mw_values[4] = id_to_mw[ids[4]]
  mw_values[5] = id_to_mw[ids[5]]
  mw_values[6] = id_to_mw[ids[6]]
  print poolID, mw_values[1], mw_values[2], mw_values[3], mw_values[4], mw_values[5], mw_values[6]
}' ExpMW_lookup.txt PoolMap_unix.txt > FinalDataset.txt

# Clean up temporary files
rm ExpMW_unix.txt PoolMap_unix.txt ExpMW_lookup.txt
