#!/bin/bash

# this unzips a file that is presumed to contain a single json file of the same name as the zip (without .zip extension)
# this does NOT have error handling for any missing or overwriting of files. 
# it does not take into consideration that other files with the same names may end up getting processed in secondary steps

# temp file is used to support counter
TEMPFILE=/tmp/$$.tmp
echo 0 > $TEMPFILE

# unzip the zipped json file containing the vci records and split into separate files with prefix 'vci'
zipfile=$1
jsonfile=$(echo ${zipfile%.*})
splitfileprefix="vci"
unzip $zipfile $jsonfile
gcsplit --digits=4 --prefix=$splitfileprefix $jsonfile "/^    },$/" "{*}"

# declare and capture the files to process in an array
declare -a files
files=(vci*)
pos=$(( ${#files[*]} - 1 ))
last=${files[$pos]}

# loop through the vci* files, prep and process to dmwg#.json form...
for FILE in "${files[@]}"
do 

  # output the filename being processed
  printf "%s\n" "$FILE"

  # Fetch, increase and store the counter
  counter=$[$(cat $TEMPFILE) + 1]
  echo $counter > $TEMPFILE

  # define the vcixxx.json input filename and dmwgxxx.json output filenames
  input=${FILE}.json
  output="dmwg${counter}.json"

  # The last file split out contains an additional 2 lines at the end 
  # which must be removed before continuing
  if [[ $FILE == $last ]]
  then
   sed -e :a -e '$d;N;2,2ba' -e 'P;D' "$FILE" > tmpfile; mv tmpfile "$FILE"
  fi

  # fix split files by removing the first line and adding the closing paren at the end.
  sed '1d' "$FILE" > tmpfile; mv tmpfile "$input"; echo "    }" >> "$input"
  printf "%s\n" "Transforming $input to $output"

  # transform the vci json to dmwg using project's python script
  python ../VCI2DMWG.py $input $output

done 

# delete the tempfile
unlink $TEMPFILE
