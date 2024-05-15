# MetaboRI

The retention time is a crucial piece of information for metabolite identification. However, the retention time of the same metabolite can vary under different chromatographic conditions and parameters, leading to discrepancies between experimental and library retention times. ***MetaboRI*** is a retention index-based program designed to address these issues.

The R package can be installed and used through the following methods:

```R
devtools::install_github("dongdongdong7/MetaboRI")
library(MetaboRI)
```

## Endogenous retention time corrects metabolites

Endogenous retention time correction metabolite have three types: Carnitine, Dicarboxylic acid, and Hydroxy fatty acid. These metabolites have several advantages: (1) wide coverage of retention time; (2) endogenous biological metabolites; (3) retention time is related to chain length, making it easier to predict.

The following format (```lc_tibble```) is used to store endogenous retention time-corrected metabolites: 

### Carnitine

```R
> data("standard_car", )
> standard_car
# A tibble: 16 × 11
   id          name                      rt    lc    ri exactmass formula   smiles             inchi inchikey predicted
   <chr>       <chr>                  <dbl> <dbl> <dbl>     <dbl> <chr>     <chr>              <chr> <chr>    <lgl>    
 1 HMDB0000062 L-Carnitine             0.52     0     0      162. C7H15NO3  C[N+](C)(C)C[C@H]… InCh… PHIQHXF… FALSE    
 2 HMDB0000201 L-Acetylcarnitine       0.92     2   200      204. C9H17NO4  CC(=O)O[C@H](CC(O… InCh… RDHQFKQ… FALSE    
 3 HMDB0000824 Propionylcarnitine      2.24     3   300      218. C10H19NO4 CCC(=O)O[C@H](CC(… InCh… UFAHZIU… FALSE    
 4 HMDB0002013 Butyrylcarnitine        3.08     4   400      231. C11H21NO4 CCCC(=O)O[C@H](CC… InCh… QWYFHHG… FALSE    
 5 HMDB0013128 Valerylcarnitine        3.92     5   500      246. C12H23NO4 CCCCC(=O)O[C@H](C… InCh… VSNFQQX… FALSE    
 6 HMDB0000756 Hexanoylcarnitine       4.56     6   600      260. C13H25NO4 CCCCCC(=O)O[C@H](… InCh… VVPRQWT… FALSE    
 7 HMDB0000791 Octanoylcarnitine       5.52     8   800      288. C15H29NO4 CCCCCCCC(=O)O[C@H… InCh… CXTATJF… FALSE    
 8 HMDB0000651 Decanoylcarnitine       6.29    10  1000      315. C17H33NO4 CCCCCCCCCC(=O)O[C… InCh… LZOSYCM… FALSE    
 9 HMDB0002250 Dodecanoylcarnitine     6.91    12  1200      344. C19H37NO4 CCCCCCCCCCCC(=O)O… InCh… FUJLYHJ… FALSE    
10 HMDB0005066 Tetradecanoylcarnitine  7.48    14  1400      372. C21H41NO4 CCCCCCCCCCCCCC(=O… InCh… PSHXNVG… FALSE    
11 HMDB0000222 Palmitoylcarnitine      7.97    16  1600      400. C23H45NO4 CCCCCCCCCCCCCCCC(… InCh… XOMRRQX… TRUE     
12 HMDB0000848 Stearoylcarnitine       8.42    18  1800      428. C25H49NO4 CCCCCCCCCCCCCCCCC… InCh… FNPHNLN… FALSE    
13 HMDB0006460 Arachidyl carnitine     8.8     20  2000      455. C27H53NO4 CCCCCCCCCCCCCCCCC… InCh… SVJLJQB… FALSE    
14 HMDB0062468 Docosanoylcarnitine     9.12    22  2200      483. C29H57NO4 CCCCCCCCCCCCCCCCC… InCh… IUMXSSO… FALSE    
15 HMDB0240665 Lignoceroylcarnitine    9.58    24  2400      511. C31H61NO4 CCCCCCCCCCCCCCCCC… InCh… YDUFZFU… FALSE    
16 HMDB0006347 Hexacosanoyl carnitine  9.68    26  2600      539. C33H65NO4 CCCCCCCCCCCCCCCCC… InCh… KOCKWDD… TRUE     
```

### Dicarboxylic acid

```R
> data("standard_dic")
> standard_dic
# A tibble: 25 × 11
   id          name                  rt    lc    ri exactmass formula  smiles                inchi   inchikey predicted
   <chr>       <chr>              <dbl> <dbl> <dbl>     <dbl> <chr>    <chr>                 <chr>   <chr>    <lgl>    
 1 HMDB0002329 Oxalic acid         0.54     2   200      90.0 C2H2O4   OC(=O)C(O)=O          InChI=… MUBZPKH… FALSE    
 2 HMDB0000691 Malonic acid        0.73     3   300     104.  C3H4O4   OC(=O)CC(O)=O         InChI=… OFOBLEO… FALSE    
 3 HMDB0000254 Succinic acid       1.4      4   400     118.  C4H6O4   OC(=O)CCC(O)=O        InChI=… KDYFGRW… FALSE    
 4 HMDB0000661 Glutaric acid       2.39     5   500     132.  C5H8O4   OC(=O)CCCC(O)=O       InChI=… JFCQEDH… FALSE    
 5 HMDB0000448 Adipic acid         3.15     6   600     146.  C6H10O4  OC(=O)CCCCC(O)=O      InChI=… WNLRTRB… FALSE    
 6 HMDB0000857 Pimelic acid        3.91     7   700     160.  C7H12O4  OC(=O)CCCCCC(O)=O     InChI=… WLJVNTC… FALSE    
 7 HMDB0000893 Suberic acid        4.52     8   800     174.  C8H14O4  OC(=O)CCCCCCC(O)=O    InChI=… TYFQFVW… FALSE    
 8 HMDB0000784 Azelaic acid        5.05     9   900     188.  C9H16O4  OC(=O)CCCCCCCC(O)=O   InChI=… BDJRBEY… FALSE    
 9 HMDB0000792 Sebacic acid        5.53    10  1000     202.  C10H18O4 OC(=O)CCCCCCCCC(O)=O  InChI=… CXMXRPH… FALSE    
10 HMDB0000888 Undecanedioic acid  5.98    11  1100     216.  C11H20O4 OC(=O)CCCCCCCCCC(O)=O InChI=… LWBHHRR… FALSE    
# ℹ 15 more rows
# ℹ Use `print(n = ...)` to see more rows
```

### Hydroxy fatty acid

```R
> data("standard_hyd")
> standard_hyd
# A tibble: 21 × 11
   id          name                        rt    lc    ri exactmass formula  smiles            inchi inchikey predicted
   <chr>       <chr>                    <dbl> <dbl> <dbl>     <dbl> <chr>    <chr>             <chr> <chr>    <lgl>    
 1 HMDB0000115 Glycolic acid             0.63     2   200      76.0 C2H4O3   OCC(O)=O          InCh… AEMRFAO… FALSE    
 2 HMDB0000190 L-Lactic acid             0.9      3   300      90.0 C3H6O3   C[C@H](O)C(O)=O   InCh… JVTAAEK… FALSE    
 3 CID11266    2-Hydroxybutyric acid     1.82     4   400     104.  C4H8O3   CCC(C(=O)O)O      InCh… AFENDNX… FALSE    
 4 HMDB0001863 2-Hydroxyvaleric acid     3.19     5   500     118.  C5H10O3  CCCC(O)C(O)=O     InCh… JRHWHSJ… FALSE    
 5 HMDB0001624 2-Hydroxycaproic acid     4.32     6   600     132.  C6H12O3  CCCCC(O)C(O)=O    InCh… NYHNVHG… FALSE    
 6 CID2750949  2-hydroxyheptanoic acid   5.16     7   700     146.  C7H14O3  CCCCCC(C(=O)O)O   InCh… RGMMREB… FALSE    
 7 HMDB0000711 2-Hydroxyoctanoic acid    5.83     8   800     160.  C8H16O3  CCCCCCC(O)C(O)=O  InCh… JKRDADV… FALSE    
 8 CID5282897  2-Hydroxynonanoic acid    6.43     9   900     174.  C9H18O3  CCCCCCCC(C(=O)O)O InCh… BTJFTHO… FALSE    
 9 HMDB0242148 2-Hydroxycapric acid      6.96    10  1000     188.  C10H20O3 CCCCCCCCC(O)C(O)… InCh… GHPVDCP… FALSE    
10 HMDB0059736 2-hydroxyundecanoic acid  7.43    11  1100     202.  C11H22O3 CCCCCCCCCC(O)C(O… InCh… MNRBGFK… FALSE    
# ℹ 11 more rows
# ℹ Use `print(n = ...)` to see more rows
```

## Input

***MetaboRI*** requires the user to prepare the following ```tibble``` as input. The correction metabolite carnitine is used as an example:

### experimental_data

Experimental data, including feature id and corresponding retention time.

```R
> data("experimental_data")
> experimental_data
# A tibble: 4,629 × 2
   feature              rt
   <chr>             <dbl>
 1 0.05_97.9687m/z  0.0501
 2 0.14_136.0734m/z 0.142 
 3 0.21_391.8773n   0.214 
 4 0.43_308.2184m/z 0.431 
 5 0.45_234.9509m/z 0.451 
 6 0.45_188.0227m/z 0.451 
 7 0.45_161.0119m/z 0.451 
 8 0.45_200.9741m/z 0.451 
 9 0.45_170.0658m/z 0.451 
10 0.45_146.1651m/z 0.451 
# ℹ 4,619 more rows
# ℹ Use `print(n = ...)` to see more rows
```

### experimental_car

The retention time correction metabolite under experimental chromatographic conditions, using carnitine as an example. The retention time correction metabolite with unknown retention time can be set to NA and the predicted value set to TRUE. 

```R
> data("experimental_car")
> experimental_car
# A tibble: 16 × 11
   id          name                      rt    lc    ri exactmass formula   smiles             inchi inchikey predicted
   <chr>       <chr>                  <dbl> <int> <dbl>     <dbl> <chr>     <chr>              <chr> <chr>    <lgl>    
 1 HMDB0000062 L-Carnitine             0.55     0     0      162. C7H15NO3  C[N+](C)(C)C[C@H]… InCh… PHIQHXF… FALSE    
 2 HMDB0000201 L-Acetylcarnitine       0.91     2   200      204. C9H17NO4  CC(=O)O[C@H](CC(O… InCh… RDHQFKQ… FALSE    
 3 HMDB0000824 Propionylcarnitine      2.03     3   300      218. C10H19NO4 CCC(=O)O[C@H](CC(… InCh… UFAHZIU… FALSE    
 4 HMDB0002013 Butyrylcarnitine        3.03     4   400      231. C11H21NO4 CCCC(=O)O[C@H](CC… InCh… QWYFHHG… FALSE    
 5 HMDB0013128 Valerylcarnitine        3.9      5   500      246. C12H23NO4 CCCCC(=O)O[C@H](C… InCh… VSNFQQX… FALSE    
 6 HMDB0000756 Hexanoylcarnitine       4.63     6   600      260. C13H25NO4 CCCCCC(=O)O[C@H](… InCh… VVPRQWT… FALSE    
 7 HMDB0000791 Octanoylcarnitine       5.65     8   800      288. C15H29NO4 CCCCCCCC(=O)O[C@H… InCh… CXTATJF… FALSE    
 8 HMDB0000651 Decanoylcarnitine       6.48    10  1000      315. C17H33NO4 CCCCCCCCCC(=O)O[C… InCh… LZOSYCM… FALSE    
 9 HMDB0002250 Dodecanoylcarnitine     7.17    12  1200      344. C19H37NO4 CCCCCCCCCCCC(=O)O… InCh… FUJLYHJ… FALSE    
10 HMDB0005066 Tetradecanoylcarnitine  7.76    14  1400      372. C21H41NO4 CCCCCCCCCCCCCC(=O… InCh… PSHXNVG… FALSE    
11 HMDB0000222 Palmitoylcarnitine      8.28    16  1600      400. C23H45NO4 CCCCCCCCCCCCCCCC(… InCh… XOMRRQX… FALSE    
12 HMDB0000848 Stearoylcarnitine       8.74    18  1800      428. C25H49NO4 CCCCCCCCCCCCCCCCC… InCh… FNPHNLN… FALSE    
13 HMDB0006460 Arachidyl carnitine    NA       20  2000      455. C27H53NO4 CCCCCCCCCCCCCCCCC… InCh… SVJLJQB… TRUE     
14 HMDB0062468 Docosanoylcarnitine    NA       22  2200      483. C29H57NO4 CCCCCCCCCCCCCCCCC… InCh… IUMXSSO… TRUE     
15 HMDB0240665 Lignoceroylcarnitine   NA       24  2400      511. C31H61NO4 CCCCCCCCCCCCCCCCC… InCh… YDUFZFU… TRUE     
16 HMDB0006347 Hexacosanoyl carnitine NA       26  2600      539. C33H65NO4 CCCCCCCCCCCCCCCCC… InCh… KOCKWDD… TRUE
```

## Workflow

### Load data

```R
data("experimental_data")
data("experimental_car")
data("standard_car")
```

### Retention time prediction

```R
ResList <- lc2rt(lc_tibble = experimental_car, targetRt = max(experimental_data$rt))
experimental_car <- ResList$lc_tibble
ResList$p
```

<img src=".\assets\image-20240515153642876.png" alt="image-20240515153642876" style="zoom:50%;" />

```R
ResList <- lc2rt(lc_tibble = standard_car, targetRt = max(experimental_data$rt))
standard_car <- ResList$lc_tibble
ResList$p
```

<img src=".\assets\image-20240515153730135.png" alt="image-20240515153730135" style="zoom:50%;" />

### Set dead time

It is recommended that a lower limit be set for the range of retention time corrections, which is typically the value with the minimum retention time of experimental metabolite feature or the dead time. 

```R
experimental_car <- deadTime_add(lc_tibble = experimental_car, deadTime = min(experimental_data$rt))
standard_car <- deadTime_add(lc_tibble = standard_car, deadTime = min(experimental_data$rt))
```

### Calculate retention index

```R
experimental_data$ri <- sapply(experimental_data$rt, function(x) {
  rt2ri(x, experimental_lc = experimental_car)
})
```

```R
> experimental_data
# A tibble: 4,629 × 3
   feature              rt     ri
   <chr>             <dbl>  <dbl>
 1 0.05_97.9687m/z  0.0501 -100  
 2 0.14_136.0734m/z 0.142   -81.7
 3 0.21_391.8773n   0.214   -67.2
 4 0.43_308.2184m/z 0.431   -23.7
 5 0.45_234.9509m/z 0.451   -19.9
 6 0.45_188.0227m/z 0.451   -19.9
 7 0.45_161.0119m/z 0.451   -19.9
 8 0.45_200.9741m/z 0.451   -19.9
 9 0.45_170.0658m/z 0.451   -19.9
10 0.45_146.1651m/z 0.451   -19.9
# ℹ 4,619 more rows
# ℹ Use `print(n = ...)` to see more rows
```

### Retention time correction

```R
# Make the experimental lc_tibble and standard lc_tibble have the same number of lc.
ResList <- compare_lc_tibble(experimental_lc = experimental_car, standard_lc = standard_car)
experimental_car <- ResList$experimental_lc
standard_car <- ResList$standard_lc
# Calculate new retention time
experimental_data <- orignRT2newRT(experimental_data = experimental_data,
                                   experimental_lc = experimental_car,
                                   standard_lc = standard_car, thread = 1)
```

```R
> experimental_data 
# A tibble: 4,629 × 4
   feature          orign_rt     ri new_rt
   <chr>               <dbl>  <dbl>  <dbl>
 1 0.05_97.9687m/z    0.0501 -100   0.0501
 2 0.14_136.0734m/z   0.142   -81.7 0.136 
 3 0.21_391.8773n     0.214   -67.2 0.204 
 4 0.43_308.2184m/z   0.431   -23.7 0.409 
 5 0.45_234.9509m/z   0.451   -19.9 0.427 
 6 0.45_188.0227m/z   0.451   -19.9 0.427 
 7 0.45_161.0119m/z   0.451   -19.9 0.427 
 8 0.45_200.9741m/z   0.451   -19.9 0.427 
 9 0.45_170.0658m/z   0.451   -19.9 0.427 
10 0.45_146.1651m/z   0.451   -19.9 0.427 
# ℹ 4,619 more rows
# ℹ Use `print(n = ...)` to see more rows
```



