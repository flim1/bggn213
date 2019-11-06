Class 11 Structural Bioinformatics
================

# What is in the PDB database

Download PDB statistics summary sorted by experimental method and
molecular type

``` r
data <- read.csv("Data Export Summary.csv")
```

Q. Determine the percentage of structures solved by X-Ray and Electron
Microscopy

``` r
sum(data$Total)
```

    ## [1] 157530

``` r
ans <- data$Total/sum(data$Total)*100
names(ans) <- data$Experimental.Method
ans
```

    ##               X-Ray                 NMR Electron Microscopy 
    ##         89.06176601          8.13495842          2.51444169 
    ##               Other        Multi Method 
    ##          0.19234432          0.09648956

``` r
round(ans,2)
```

    ##               X-Ray                 NMR Electron Microscopy 
    ##               89.06                8.13                2.51 
    ##               Other        Multi Method 
    ##                0.19                0.10

Q. What proportion of structures are protein

``` r
round(sum(data$Proteins)/sum(data$Total)*100, 2)
```

    ## [1] 92.71

Q2: Type HIV in the PDB website search box on the home page and
determine how many HIV-1 protease structures are in the current PDB?
3048

# WOrking with biomolecular data in R

``` r
library(bio3d)
pdb <- read.pdb("1hsg")
```

    ##   Note: Accessing on-line PDB file

``` r
pdb
```

    ## 
    ##  Call:  read.pdb(file = "1hsg")
    ## 
    ##    Total Models#: 1
    ##      Total Atoms#: 1686,  XYZs#: 5058  Chains#: 2  (values: A B)
    ## 
    ##      Protein Atoms#: 1514  (residues/Calpha atoms#: 198)
    ##      Nucleic acid Atoms#: 0  (residues/phosphate atoms#: 0)
    ## 
    ##      Non-protein/nucleic Atoms#: 172  (residues: 128)
    ##      Non-protein/nucleic resid values: [ HOH (127), MK1 (1) ]
    ## 
    ##    Protein sequence:
    ##       PQITLWQRPLVTIKIGGQLKEALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYD
    ##       QILIEICGHKAIGTVLVGPTPVNIIGRNLLTQIGCTLNFPQITLWQRPLVTIKIGGQLKE
    ##       ALLDTGADDTVLEEMSLPGRWKPKMIGGIGGFIKVRQYDQILIEICGHKAIGTVLVGPTP
    ##       VNIIGRNLLTQIGCTLNF
    ## 
    ## + attr: atom, xyz, seqres, helix, sheet,
    ##         calpha, remark, call

``` r
pdb$atom[1,"resid"]
```

    ## [1] "PRO"

``` r
aa321(pdb$atom[pdb$calpha,"resid"])
```

    ##   [1] "P" "Q" "I" "T" "L" "W" "Q" "R" "P" "L" "V" "T" "I" "K" "I" "G" "G"
    ##  [18] "Q" "L" "K" "E" "A" "L" "L" "D" "T" "G" "A" "D" "D" "T" "V" "L" "E"
    ##  [35] "E" "M" "S" "L" "P" "G" "R" "W" "K" "P" "K" "M" "I" "G" "G" "I" "G"
    ##  [52] "G" "F" "I" "K" "V" "R" "Q" "Y" "D" "Q" "I" "L" "I" "E" "I" "C" "G"
    ##  [69] "H" "K" "A" "I" "G" "T" "V" "L" "V" "G" "P" "T" "P" "V" "N" "I" "I"
    ##  [86] "G" "R" "N" "L" "L" "T" "Q" "I" "G" "C" "T" "L" "N" "F" "P" "Q" "I"
    ## [103] "T" "L" "W" "Q" "R" "P" "L" "V" "T" "I" "K" "I" "G" "G" "Q" "L" "K"
    ## [120] "E" "A" "L" "L" "D" "T" "G" "A" "D" "D" "T" "V" "L" "E" "E" "M" "S"
    ## [137] "L" "P" "G" "R" "W" "K" "P" "K" "M" "I" "G" "G" "I" "G" "G" "F" "I"
    ## [154] "K" "V" "R" "Q" "Y" "D" "Q" "I" "L" "I" "E" "I" "C" "G" "H" "K" "A"
    ## [171] "I" "G" "T" "V" "L" "V" "G" "P" "T" "P" "V" "N" "I" "I" "G" "R" "N"
    ## [188] "L" "L" "T" "Q" "I" "G" "C" "T" "L" "N" "F"

``` r
# Select all C-alpha atoms (return their indices)
ca.inds <- atom.select(pdb, "calpha")
ca.inds
```

    ## 
    ##  Call:  atom.select.pdb(pdb = pdb, string = "calpha")
    ## 
    ##    Atom Indices#: 198  ($atom)
    ##    XYZ  Indices#: 594  ($xyz)
    ## 
    ## + attr: atom, xyz, call

``` r
# Print details of the first few selected atoms
head(pdb$atom[ca.inds$atom,])
```

    ##    type eleno elety  alt resid chain resno insert      x      y     z o
    ## 2  ATOM     2    CA <NA>   PRO     A     1   <NA> 30.307 38.663 5.319 1
    ## 9  ATOM     9    CA <NA>   GLN     A     2   <NA> 30.158 36.492 2.199 1
    ## 18 ATOM    18    CA <NA>   ILE     A     3   <NA> 29.123 33.098 3.397 1
    ## 26 ATOM    26    CA <NA>   THR     A     4   <NA> 29.774 30.143 1.062 1
    ## 33 ATOM    33    CA <NA>   LEU     A     5   <NA> 27.644 27.003 1.144 1
    ## 41 ATOM    41    CA <NA>   TRP     A     6   <NA> 30.177 24.150 1.279 1
    ##        b segid elesy charge
    ## 2  40.62  <NA>     C   <NA>
    ## 9  41.30  <NA>     C   <NA>
    ## 18 34.13  <NA>     C   <NA>
    ## 26 30.14  <NA>     C   <NA>
    ## 33 30.12  <NA>     C   <NA>
    ## 41 30.82  <NA>     C   <NA>

``` r
# Selected xyz coordinates
head(pdb$xyz[, ca.inds$xyz])
```

    ## [1] 30.307 38.663  5.319 30.158 36.492  2.199

Q8. Use the Bio3D write.pdb() function to write out a protein only PDB
file for viewing in VMD. Also write out a second separate PDB file for
the ligand with residue name MK1

First select a protein and write out a file “1hsg\_protein.pdb”

``` r
prot <- atom.select(pdb, "protein", value=TRUE)
write.pdb(prot, file="1hsg_protein.pdb")
```

Select ligand and write out a file

``` r
lig <- atom.select(pdb, "ligand", value=TRUE)
write.pdb(lig, file="1hsg_ligand.pdb")
```

``` r
# The 'devtools' package allows us to install development versions
#install.packages("devtools")
# Install the bio3d.view package from bitbucket
#devtools::install_bitbucket("Grantlab/bio3d-view")
```

``` r
# Load the package
library("bio3d.view")
# view the 3D structure
view(pdb, "overview", col="sse")
```

    ## Computing connectivity from coordinates...

``` r
# Load the package
pdb <- read.pdb("1hel")
```

    ##   Note: Accessing on-line PDB file

``` r
# Normal mode analysis calculation
modes <- nma(pdb)
```

    ##  Building Hessian...     Done in 0.04 seconds.
    ##  Diagonalizing Hessian...    Done in 0.17 seconds.

``` r
m7 <- mktrj(modes,
 mode=7,
 file="mode_7.pdb")
view(m7, col=vec2color( rmsf(m7) ))
```

    ## Potential all C-alpha atom structure(s) detected: Using calpha.connectivity()
