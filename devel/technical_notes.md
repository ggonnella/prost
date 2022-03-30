# Snakemake
- in previous versions of snakemake the `ancient()` function did
   not work properly here (i.e. on files created by other rules);
   in this case, update snakemake (version 5.32.2 works)
- when updating snakemake, conda messed up with libffi; the
  workaround is to go to the directory lib under the active env location
  (see with conda info, in my case it is `~/tools/miniconda/3/lib`) and do
  `ln -s libffi.so.7 libffi.so.6`

