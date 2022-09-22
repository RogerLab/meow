#### v0.5.2 2022-09-22
- Changed output file name to match desired ouput when given zero high partition classes or zero low partition classes
  (formerly would add ".high.nex" or ".low.nex" to the end of the files)

#### v0.5.1 2022-09-04
- Added option to include invariant class in output nexus file
- Changed remove invariant flag from -I to -ri
- Added flags:
	- Invariant class (-I)
	- Sort sequence file (-sort)
	- Add epsilon to H-Clust (-e)
- Added "total classes" count to output log
- Various bug fixes
- Updated tests to reflect changes

#### v0.5.0 2022-07-08 
- Initial commits