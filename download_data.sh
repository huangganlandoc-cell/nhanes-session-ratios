#!/bin/sh
# Download the public NHANES files used in the analysis (about 540 MB) into $NHANES_DATA_DIR
# (default: ~/data/nhanes_public). Checksums of the files used for the reported results: SHA256SUMS_raw_data.txt
set -e
HERE=$(cd "$(dirname "$0")" && pwd)
D=${NHANES_DATA_DIR:-$HOME/data/nhanes_public}
mkdir -p "$D/nhanes3" "$D/continuous_xpt" "$D/mortality_2019"
N3=https://wwwn.cdc.gov/nchs/data/nhanes3/1a
LM=https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality
for f in lab.dat lab.sas exam.dat exam.sas adult.dat adult.sas; do
  curl -fsSL --retry 3 -o "$D/nhanes3/$f" "$N3/$f"; echo "nhanes3/$f"
done
for f in NHANES_III_MORT_2019_PUBLIC.dat SAS_ReadInProgramAllSurveys.sas; do
  curl -fsSL --retry 3 -o "$D/nhanes3/$f" "$LM/$f"; echo "nhanes3/$f"
done
for c in 1999_2000 2001_2002 2003_2004 2005_2006 2007_2008 2009_2010 2011_2012 2013_2014 2015_2016 2017_2018; do
  f=NHANES_${c}_MORT_2019_PUBLIC.dat
  curl -fsSL --retry 3 -o "$D/mortality_2019/$f" "$LM/$f"; echo "mortality_2019/$f"
done
while read y f; do
  [ -n "$f" ] || continue
  curl -fsSL --retry 3 -o "$D/continuous_xpt/$f.xpt" "https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/$y/DataFiles/$f.xpt"
  echo "continuous_xpt/$f.xpt"
done < "$HERE/data_files_continuous.txt"
echo "Done. To verify: cd \"$D\" && shasum -a 256 -c \"$HERE/SHA256SUMS_raw_data.txt\""
