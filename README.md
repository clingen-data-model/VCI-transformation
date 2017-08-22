# VCI-transformation
Scripts for transforming VCI JSON-LD into DMWG Interpretation JSON-LD

This project consists of 

1) A series of python classes representing classes in the DMWG interpretation data model
2) A script for dynamically generating these python classes 
3) Static components for handling serialization
4) A script for converting VCI JSON files into DMWG JSON files using the  classes above.
5) Unit tests

VCI2DMWG has a dependency on the clingen_interpretation library, which is used to serialize DMWG-style
JSON files.  The github repository can be found [here](https://github.com/clingen-data-model/interpretation_json).
Follow the instructions on that page to install the library.

To convert a VCI json file (input.json) into a DMWG JSON file (output.json) use the script VCI2DMG.py:
```python VCI2DMWG.py input.json output.json```

Sample input and output files are found in the `test_data` directory.

To generate the classes, run
```python generate_interpretation_library.py <flattened>```
where <flattened> is the directory containing the flattened data model.  

This will create several files:
  * interpretation_generated.py: Python classes with getters and setters for properties
  * coding_generated.py: Python classes generated for Coding and CodeableConcept
  * interpretation_constants.py: property and other names
  * ValueSets/*: JSON files containing the defined value sets

In addition there are several static files:
  * node.py: base class for the generated python classes
  * interpretation_extras.py: Some utility functions for creating DMWG style interpretations and serializing them
  * CodingFactory: decorators to handle codings and codeable concepts.



The tests are to be run using `pytest`.
