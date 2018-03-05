#!/bin/bash

# this unzips a file that is presumed to contain a single json file of the same name as the zip (without .zip extension)
# this does NOT have error handling for any missing or overwriting of files. 
# it does not take into consideration that other files with the same names may end up getting processed in secondary steps
i
TEMPFILE=/tmp/$$.tmp
echo 0 > $TEMPFILE

zipfile=$1
jsonfile=$(echo ${zipfile%.*})
splitfileprefix="vci"

unzip $zipfile $jsonfile
gcsplit --digits=4 --prefix=$splitfileprefix $jsonfile "/^    },$/" "{*}"

# Loop through product of gcsplit operation above...
for filename in ${splitfileprevix}*; do
    printf "%s\n" "$filename"
    testIt=$(head -3 "$filename" | pcregrep -M '},\n    {\n        \"modeIn')
    if [ -z $testIt ]; then
       echo -n ""
    else
       # Fetch the counter value and increase it
       counter=$[$(cat $TEMPFILE) + 1]

       # Store the new value
       echo $counter > $TEMPFILE

       #define the vcixxx.json input filename and dmwgxxx.json output filenames
       input=${filename}.json
       output="dmwg${counter}.json"
       
       # fix split files by removing the first line and adding the closing paren at the end.
       sed '1d' "$filename" > tmpfile; mv tmpfile "$input"; echo "    }" >> "$input"
       printf "%s\n" "Transforming $input to $output"
       
       # transform the vci json to dmwg using project's python script
       python ../VCI2DMWG.py $input $output 
    fi
done

# delete the tempfile
unlink $TEMPFILE
