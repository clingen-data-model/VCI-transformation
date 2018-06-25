# VCI-transformation
Scripts for transforming VCI JSON-LD into DMWG Interpretation JSON-LD

This project consists of 

1) A series of python classes representing classes in the DMWG interpretation data model
2) A script for dynamically generating these python classes 
3) Static components for handling serialization
4) A script for converting VCI JSON files into DMWG JSON files using the  classes above.
5) Unit tests
6) A folder of sample input files with an archive of historical real data file sets.
7) A script for parsing and transforming a single multi-record input json file.

## Setup Python Environment
VCI2DMWG has a dependency on the clingen_interpretation library, which is used to serialize DMWG-style
JSON files.  The github repository can be found [here](https://github.com/clingen-data-model/interpretation_json).
Follow the instructions on that page to install the library, or just install this and all other dependencies using
`pip install -r requirements.txt`.

### alternative option with virtualenv
Some may prefer to work in a python virtualenv like this:

```
virtualenv venv
. venv/bin/activate
pip install -r requirements.txt
```

(If you do this, you will need to activate the virtualenv in whichever shell you want to use these scripts)

## Execute Script 
### Standard - Single One-Record Input File
To convert a VCI json file (input.json) into a DMWG JSON file (output.json) use the script VCI2DMG.py:
```python VCI2DMWG.py input.json output.json```

Sample input and output files are found in the `test_data` directory.

### Batch - Single Multi-Record Input File
When given a single zip file containing a single json file containing multiple VCI json records
(see test_data/archive/DMWG-PAH-INTERPRETATIONS.json.zip for example), you can automate the
unzipping, parsing and transformations of each record with a single bash script called 
parseVCIinput.sh found in the root directory.

To parse and transform a single mulit-record input zip file, run
```./parseVCIinput.sh <zipfilename> &> <outputlogfilename>```
where 
```<zipfilename>``` is a reference to the zip file containing a single json file with one ormore VCI records lists, and
```<outputlogfilename>``` is a reference to the output (errors included) from process.  

#### saving results to archive folder (for ClinGen use only)
The output of the VCI transformation run and its logged output should be passed to the clinvar-submitter process.

The resulting dmwg*.json files should be zipped and passed on.
The input VCI zip file should be stored in the test_data/archive folder and checked into github for future 
users to reference as needed. 

## Generating Dependent Classes
To generate the classes, run
```python generate_interpretation_library.py <flattened>```
where ```<flattened>``` is the directory containing the flattened data model.  

This will create several files:
  * interpretation_generated.py: Python classes with getters and setters for properties
  * coding_generated.py: Python classes generated for Coding and CodeableConcept
  * interpretation_constants.py: property and other names
  * ValueSets/*: JSON files containing the defined value sets

In addition there are several static files:
  * node.py: base class for the generated python classes
  * interpretation_extras.py: Some utility functions for creating DMWG style interpretations and serializing them
  * CodingFactory: decorators to handle codings and codeable concepts.

## Running Unit Tests
The tests are to be run using `pytest`.
