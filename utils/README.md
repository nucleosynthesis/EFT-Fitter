#  A bunch of helpful scripts to plot whats going on 

## `checkinputs.py`

check your inputs with this script

`python checkinputs.py module.input`

  -  if the input is a data (eg likelhood scan or best fits) it will plot those
  -  if the input is a theory uncertainty file, it will plot the uncertainties and correlations

## `pickle2text.py`

cript to covert piclked results into a simple txt file (can then be used with `overlay_scans.py`

`python pickle2text.py results.pkl parameter outfilename`

## `root2text.py`

A script to take a ROOT (higgsCombineX.root) file of a likelihood scan and turn into plain text

`python root2text.py input.root outname [brach_names]`

  - By default, all branches will be plotted, select only parameters you want using branch_names options - brn1 brnn2 brn3 ...


## `overlay_scans.py`

Very quick tool to plot likelihood scans from txt files - i.e those produced from above (1D only)

`python overlay_scans.py outname input1.txt [input2.txt input3.txt....]`

  -  set outname=`show` to just show the results rather than save as pdf/png.
  -  overlay as many scans  as you like. The script will plot whichever parameter it finds though

