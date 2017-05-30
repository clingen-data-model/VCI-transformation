# VCI-transformation
Scripts for transforming VCI JSON-LD into DMWG Interpretation JSON-LD

This project consists of 

1) A series of python classes representing classes in the DMWG interpretation data model
2) A script for dynamically generating these python classes 
3) Static components for handling serialization
4) A script for converting VCI JSON files into DMWG JSON files using the  classes above.
5) Unit tests



To generate the classes, run
```python generate_interpretation_library.py <flattened>```
where <flattened> is the directory containing the flattened data model.  

This will create several files:
  * interpretation_generated.py: Python classes with getters and setters for properties
  * coding_generated.py: 
