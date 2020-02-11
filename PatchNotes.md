# Patch Notes for simPLEX

## 7th June
- Added ability category float parameter.
  - Runs significantly slower than a continuous float.
  - Should be renamed to discrete, rather then category float.
  - There are duplicates of the rate categories currently.
- Changed the form of the lua call to make new parameters in the model file.
- When empty path (aka "") is specified for output file, this will result in no data being saved for this output file.
- Found a bug regarding substituion counts.
- Currently doesn't run in release mode.
- Added basic logging funtionality.
- Tracks substitution locations.
