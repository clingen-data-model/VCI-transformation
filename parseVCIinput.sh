#!/bin/bash
usage()
{
  echo "Usage: $0 [-o --overwrite] [-d PUBLISH_DATE (yyyy-mm-dd)] [-g PUBLISH-GROUP]"
  exit 2
}

exists()
{
  for file in "$@"; do
    if [[ -f $file ]]; then
      # file found, return true
      return 0
    fi
  done
  return -1
}

#########################
# Main script starts here

unset PUBLISH_DATE PUBLISH_GROUP
OVERWRITE=false

options=':od:g:h'
while getopts $options option
do
    case $option in
        o  ) OVERWRITE=true;;
        d  ) PUBLISH_DATE=$OPTARG;;
        g  ) PUBLISH_GROUP=$OPTARG;;
        h  ) usage; exit;;
        \? ) echo "Unknown option: -$OPTARG" >&2; exit 1;;
        :  ) echo "Missing option argument for -$OPTARG" >&2; exit 1;;
        *  ) echo "Unimplemented option: -$OPTARG" >&2; exit 1;;
    esac
done

shift $(($OPTIND - 1))

# make sure the publish_date and publish_group argument are passed
if [ -z "$PUBLISH_DATE" -o -z "$PUBLISH_GROUP" ]
then
    echo "Both -d publish_date and -g publish_group are required."
    exit 1
fi

PUBLISH_GROUP=$(echo $PUBLISH_GROUP | awk '{print toupper($0)}')

outputgzip="$PUBLISH_GROUP-CGSEPIO-$PUBLISH_DATE.gzip"
splitfileprefix="vci-"
vcimask="$splitfileprefix*"
cgsepiomask="cg-sepio-*"

if exists $vcimask -o exists $cgsepiomask -o exists $outputgzip
then
  if ! $OVERWRITE; then
    echo "FOUND EXISTING FILES - No Overwrite option provided."
    exit 1
  else
    echo "FOUND EXISTING FILES and REMOVING - Overwrite option provided"
    if exists $vcimask; then rm vci-*; fi
    if exists $cgsepiomask; then rm cg-sepio-*; fi
    if exists $outputgzip; then rm $outputgzip; fi
  fi
fi

# this unzips a file that is presumed to contain a single json file of the same name as the zip (without .zip extension)

# unzip the zipped json file containing the vci records and split into separate files with prefix 'vci'
jsonfile="$PUBLISH_GROUP-INTERPRETATIONS-$PUBLISH_DATE.json"
zipfile=$jsonfile.zip

if [ ! -f $jsonfile ]; then
  if [ ! -f $zipfile ]; then
    unzip $zipfile $jsonfile
    if [ $? -ne 0 ]; then
      echo "Error: Unzip of $zipfile failed!"
      exit 1
    fi
  else
    echo "Error: VCI json or zip files ($jsonfile, $zipfile) not found."
    exit 1
  fi
fi

gcsplit --quiet --digits=4 --prefix=$splitfileprefix $jsonfile "/^    },$/" "{*}"
if [ $? -ne 0 ]; then
  echo "Error: gcsplit of $jsonfile failed!"
  exit 1
fi

# declare and capture the files to process in an array
declare -a files
files=($vcimask)
pos=$(( ${#files[*]} - 1 ))
last=${files[$pos]}
# use regex to pull the matching suffix of the vci* files to build the output filenames that match.
regex="$splitfileprefix([0-9]+)"

# loop through the vci* files, prep and process to dmwg#.json form...
for FILE in "${files[@]}"
do
  if [[ $FILE =~ $regex ]]
  then
    suffix="${BASH_REMATCH[1]}"
  else
    echo "Error: $FILE doesn't match the expected vci-*.json format!"
    exit 1
  fi

  # The last file split out contains an additional 2 lines at the end
  # which must be removed before continuing
  if [[ $FILE == $last ]]
  then
    sed -e :a -e '$d;N;2,2ba' -e 'P;D' "$FILE" > tmpfile; mv tmpfile "$FILE"
    if [ $? -ne 0 ]; then
      echo "Error: sed failed removing the additional 2 lines from the last file $FILE."
      exit 1
    fi
  fi

  # define the vci-xxx.json input filename and cg-sepio-acmg-xxx.json output filenames
  input=${FILE}.json

  # fix split files by removing the first line and adding the closing paren at the end.
  sed '1d' "$FILE" > tmpfile; mv tmpfile "$input"; echo "    }" >> "$input"
  if [ $? -ne 0 ]; then
    echo "Error: sed failed fixing vci file first line and adding closing paren at end of $FILE."
    exit 1
  fi
  rm $FILE

  output="cg-sepio-${suffix}.json"
  printf "%s\n" "Transforming $input to $output"

  # transform the vci json to dmwg using project's python script
  python3 ../VCI2cgsepio.py $input $output -p $PUBLISH_DATE
  if [ $? -ne 0 ]; then
    echo "Error: VCI to ClinGen-SEPIO-ACMG format failed for $input."
    exit 1
  fi

  # remove vci-*.json input file
  rm $input

done

# gzip files and remove vci-* files
tar -czf $outputgzip cg-sepio*
rm $cgsepiomask
