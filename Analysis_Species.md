# 10-year rotation study Arkansas
Alejandro Rojas  
January 30, 2015  



Read the data into R for multivariate analysis using a text file that is tab-delimited and converted into presence/absence data:

|                | contRice | Corn_Rice | Rice_Corn | Rice_Corn_Soy | Rice_Soy | Rice_Soy_Corn | Rice.Wheat_Rice.Wheat | Rice.Wheat_Soy.Wheat | Soy_Rice | Soy.Wheat_Rice.Wheat |
|:---------------|:--------:|:---------:|:---------:|:-------------:|:--------:|:-------------:|:---------------------:|:--------------------:|:--------:|:--------------------:|
|Py_irregulare   |    6     |     9     |     8     |      11       |    4     |       3       |           4           |          16          |    10    |          8           |
|Py_paroecandrum |    6     |     4     |     5     |       6       |    6     |       4       |           2           |          15          |    5     |          13          |
|Py_sp.          |    2     |     0     |     3     |       2       |    3     |       0       |           0           |          3           |    3     |          0           |
|Py_spinosum     |    8     |    14     |    14     |      11       |    5     |       9       |           4           |          14          |    5     |          17          |
|Py_sylvaticum   |    7     |     1     |     5     |       4       |    10    |       6       |           1           |          4           |    3     |          2           |



|                | contRice | Corn_Rice | Rice_Corn | Rice_Corn_Soy | Rice_Soy | Rice_Soy_Corn | Rice.Wheat_Rice.Wheat | Rice.Wheat_Soy.Wheat | Soy_Rice | Soy.Wheat_Rice.Wheat |
|:---------------|:--------:|:---------:|:---------:|:-------------:|:--------:|:-------------:|:---------------------:|:--------------------:|:--------:|:--------------------:|
|Py_irregulare   |    1     |     1     |     1     |       1       |    1     |       1       |           1           |          1           |    1     |          1           |
|Py_paroecandrum |    1     |     1     |     1     |       1       |    1     |       1       |           1           |          1           |    1     |          1           |
|Py_sp.          |    1     |     0     |     1     |       1       |    1     |       0       |           0           |          1           |    1     |          0           |
|Py_spinosum     |    1     |     1     |     1     |       1       |    1     |       1       |           1           |          1           |    1     |          1           |
|Py_sylvaticum   |    1     |     1     |     1     |       1       |    1     |       1       |           1           |          1           |    1     |          1           |

Relative abundance was established for each of the samples:

|                | contRice  | Corn_Rice | Rice_Corn | Rice_Corn_Soy | Rice_Soy  | Rice_Soy_Corn | Rice.Wheat_Rice.Wheat | Rice.Wheat_Soy.Wheat | Soy_Rice  | Soy.Wheat_Rice.Wheat |
|:---------------|:---------:|:---------:|:---------:|:-------------:|:---------:|:-------------:|:---------------------:|:--------------------:|:---------:|:--------------------:|
|Py_irregulare   | 0.2068966 | 0.3214286 | 0.2285714 |   0.3235294   | 0.1428571 |   0.1363636   |       0.3636364       |      0.3076923       | 0.3846154 |        0.200         |
|Py_paroecandrum | 0.2068966 | 0.1428571 | 0.1428571 |   0.1764706   | 0.2142857 |   0.1818182   |       0.1818182       |      0.2884615       | 0.1923077 |        0.325         |
|Py_sp.          | 0.0689655 | 0.0000000 | 0.0857143 |   0.0588235   | 0.1071429 |   0.0000000   |       0.0000000       |      0.0576923       | 0.1153846 |        0.000         |
|Py_spinosum     | 0.2758621 | 0.5000000 | 0.4000000 |   0.3235294   | 0.1785714 |   0.4090909   |       0.3636364       |      0.2692308       | 0.1923077 |        0.425         |
|Py_sylvaticum   | 0.2413793 | 0.0357143 | 0.1428571 |   0.1176471   | 0.3571429 |   0.2727273   |       0.0909091       |      0.0769231       | 0.1153846 |        0.050         |

Calculate distances using bray-curtis for abundance data and Jaccard for presence/absence data:

```r
samplePA.dist = vegdist(x = t(dataPA), method="jaccard")
sampleREL.dist = vegdist(t(dataREL), method = "bray")
```

Calculation of hierarchical clustering for 10-year rotation data based on presence/absence data and relative abundance data:

```r
samplePA.cl = hclust(samplePA.dist)
plot(samplePA.cl, main = "Hierarchical Clustering based on presence/absence", 
     xlab="Method Complete", sub="")
```

![](Analysis_Species_files/figure-html/unnamed-chunk-5-1.png) 

```r
sampleREL.cl = hclust(sampleREL.dist)
plot(sampleREL.cl, main = "Hierarchical clustering based on relative abundance",
     xlab = "Method Complete", sub = "")
```

![](Analysis_Species_files/figure-html/unnamed-chunk-5-2.png) 

The heatmap revealed some of the differences observed among samples based on relative abundance, and some of the groupings also observed in hierarchical clustering.

```r
par(mar=c(5,7,1,1))
hcol <- brewer.pal(n = 9, name = "BuGn")
heatmap.2(dataREL,scale = "none", col=hcol, key=TRUE, trace="none", 
          density.info = "none", keysize = 0.7,  margins=c(15,15),
          cexRow = 0.9, cexCol = 0.9)
```

![](Analysis_Species_files/figure-html/unnamed-chunk-6-1.png) 

Calculation of principal coordinate analysis, the presence/absence did not show any difference among the samples on the hierachical clustering, therefore it will not show any differences in the analysis below.  However, relative abundance has an stronger effect, thus it is easier to observe differences among samples.

```r
samplePA.pcoa = cmdscale(samplePA.dist)
sampleREL.pcoa = cmdscale(sampleREL.dist)
```

Principal coordinate analysis for 10-year rotation system, using presence/absence and relative abundance:

```r
samplePA.fr <- as.data.frame(samplePA.pcoa)
ggplot(samplePA.fr, aes(x=samplePA.fr[,1], y=samplePA.fr[,2])) + geom_point() + geom_text(aes(label=row.names(samplePA.fr), size=3, vjust=-0.8), position=position_jitter(width=0.05), show_guide = FALSE) + labs(title="Principal Coordinate Analysis of 10-year rotation data - Presence/Absence", x="PCoA Axis 1", y="PcoA Axis 2")
```

![](Analysis_Species_files/figure-html/unnamed-chunk-8-1.png) 

```r
Rotation <- factor(c(1,2,2,3,2,3,4,4,2,4))
sampleREL.fr <- as.data.frame(sampleREL.pcoa)
ggplot(sampleREL.fr, aes(x=sampleREL.fr[,1], y=sampleREL.fr[,2])) + geom_point(aes(fill=Rotation), colour="black", pch=21, size=3) + geom_text(aes(label=row.names(sampleREL.fr)), size=3, vjust=-1.5, show_guide = FALSE) + labs(title="Principal Coordinate Analysis of 10-year rotation data - Relative Abundance", x="PCoA Axis 1", y="PcoA Axis 2")
```

![](Analysis_Species_files/figure-html/unnamed-chunk-8-2.png) 
