# Generation tidy spectra table for *basepeakfit*

## 1.) extracting spectrum_table with msaccess.exe (proteowizard)

+ Visit <http://proteowizard.sourceforge.net/>, download proteowizard and install.
+ run msconvert on RAW file

```
msconvert.exe 210421OFc1_PM38_025.raw -x "spectrum_table" (-o x:\workingdir)
```

+ A spectrum table file `210421OFc1_PM38_025.raw.spectrum_table.txt` is exported

## 2.) tidy exported spectrum table file

The exported files needs needs following modifications to be loaded into R:

1.) Remove first line (contains filename), e.g. windows command line command: 

```
more +1 210421OFc1_PM38_025.raw.spectrum_table.txt > 210421OFc1_PM38_025.raw.txt
```

2.) Rename hash of column ID (first column) and convert to tab separated values file, e.g. with powershell:


```
$tbl = Get-Content .\210421OFc1_PM38_025.raw.txt
$tbl = $tbl -replace("#", "NULL")
$tbl = $tbl -replace(" {1,}", "`t")
Set-Content .\210421OFc1_PM38_025.tsv $tbl
```

The exported file `210421OFc1_PM38_025.tsv` can be used as input file basepeakfit.R
