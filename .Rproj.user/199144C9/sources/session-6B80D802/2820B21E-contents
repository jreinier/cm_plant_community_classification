Plant Community Classification
================
John Reinier
7/21/2022

## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax
for authoring HTML, PDF, and MS Word documents. For more details on
using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that
includes both content as well as the output of any embedded R code
chunks within the document. You can embed an R code chunk like this:

## Query data

``` r
vibi.combined.new <- dbGetQuery(conn, "SELECT * FROM nr_misc.plot_data_for_veg_comm_analysis_newest_sample;")
```

## Reshape into site X species matrix

``` r
vibi.combined.matrix.new <- cast(vibi.combined.new, plot_id ~ acronym, value='abundance', fun="mean")
```

## Change NULLS to 0

``` r
vibi.combined.matrix.new[is.na(vibi.combined.matrix.new)] <- 0
```

## Apply Hellinger standardization and create euclidean distance matrix

``` r
vibi.combined.new.hell <- decostand(vibi.combined.matrix.new, "hel", margin = 1)
vibi.combined.new.hell.dist <- vegdist(vibi.combined.new.hell, "euc")
```

<div class="warning"
style="padding:0.1em; background-color:#E9D8FD; color:#69337A">

<span>
<p style="margin-top:1em; text-align:center">
<b>The Hellinger Transformation</b>
</p>
<p style="margin-left:1em;">
This (above) is called the Hellinger transformation. The Euclidian
distance function applied to Hellinger-transformed data produces a
Hellinger distance matrix. Useful before PCA and RDA and K-means
partitioning.
</p>

</span>

</div>

## k-means clustering with 7 groups

``` r
vibi.combined.new.kmeans7 <- kmeans(vibi.combined.new.hell, centers=7, nstart=100)
```

``` r
vibi.combined.new.kmeans7.indval.out = indval(vibi.combined.matrix.new[,c(-1,-647)], vibi.combined.new.kmeans7$cluster, numitr=10000)

vibi.new.kmeans7.gr <- vibi.combined.new.kmeans7.indval.out$maxcls[vibi.combined.new.kmeans7.indval.out$pval <= 0.05]
vibi.new.kmeans7.iv <- vibi.combined.new.kmeans7.indval.out$indcls[vibi.combined.new.kmeans7.indval.out$pval <= 0.05]
vibi.new.kmeans7.pv <- vibi.combined.new.kmeans7.indval.out$pval[vibi.combined.new.kmeans7.indval.out$pval <= 0.05]
vibi.new.kmeans7.fr <- apply(vibi.combined.matrix.new[,c(-1,-647)] > 0, 2, sum)[vibi.combined.new.kmeans7.indval.out$pval <= 0.05]
vibi.new.kmeans7.fidg <- data.frame(group=vibi.new.kmeans7.gr, indval=vibi.new.kmeans7.iv, pvalue=vibi.new.kmeans7.pv, freq=vibi.new.kmeans7.fr)
vibi.new.kmeans7.fidg <- vibi.new.kmeans7.fidg[order(vibi.new.kmeans7.fidg$group, -vibi.new.kmeans7.fidg$indval),]
```

## Write IndVal results to csv \*\*NOT RUN

``` r
write.csv(vibi.new.kmeans7.fidg, file="vibi_new_indval_kmeans7_results_07222022.csv")
```

## Combine k-means clusters (7 groups) with site X spp matrix

``` r
kmeans7_cluster <- vibi.combined.new.kmeans7$cluster
vibi.new.kmeans7.cluster <- cbind(vibi.combined.matrix.new, kmeans7_cluster)
```

## Write matrix w/ group assignment to csv \*\*NOT RUN

``` r
write.csv(vibi.new.kmeans7.cluster, file="vibi_new_cluster_kmeans7_results_04182022.csv")
```

## Try optsil function on k-means result from above

``` r
vibi.combined.new.kmeans7.optsil <- optsil(kmeans7_cluster, vibi.combined.new.hell.dist)
summary(vibi.combined.new.kmeans7.optsil)
```

    ## Number of clusters  =  7 
    ## 
    ##   1   2   3   4   5   6   7 
    ##  58 157 131  41  88 123  45 
    ## 
    ## call     =  optsil.default(x = kmeans7_cluster, dist = vibi.combined.new.hell.dist) 
    ## created  =  Wed Aug 17 10:24:25 2022

## Extract optsil results and combine with matrix

## Write optsil results to csv \*\*NOT RUN

``` r
write.csv(vibi.new.kmeans7.optsil.cluster, file="vibi_new_cluster_kmeans7_optsil_results_04182022.csv")
```

``` r
vibi.new.kmeans7.optsil.indval.out = indval(vibi.combined.matrix.new[,c(-1,-647)], vibi.combined.new.kmeans7.optsil$clustering, numitr=10000)

vibi.new.kmeans7.optsil.gr <- vibi.new.kmeans7.optsil.indval.out$maxcls[vibi.new.kmeans7.optsil.indval.out$pval <= 0.05]
vibi.new.kmeans7.optsil.iv <- vibi.new.kmeans7.optsil.indval.out$indcls[vibi.new.kmeans7.optsil.indval.out$pval <= 0.05]
vibi.new.kmeans7.optsil.pv <- vibi.new.kmeans7.optsil.indval.out$pval[vibi.new.kmeans7.optsil.indval.out$pval <= 0.05]
vibi.new.kmeans7.optsil.fr <- apply(vibi.combined.matrix.new[,c(-1,-647)] > 0, 2, sum)[vibi.new.kmeans7.optsil.indval.out$pval <= 0.05]
vibi.new.kmeans7.optsil.fidg <- data.frame(group=vibi.new.kmeans7.optsil.gr, indval=vibi.new.kmeans7.optsil.iv, pvalue=vibi.new.kmeans7.optsil.pv, freq=vibi.new.kmeans7.optsil.fr)
vibi.new.kmeans7.optsil.fidg <- vibi.new.kmeans7.optsil.fidg[order(vibi.new.kmeans7.optsil.fidg$group, -vibi.new.kmeans7.optsil.fidg$indval),]
```

## Write Indval results from optsil groups to csv \*\*NOT RUN

``` r
write.csv(vibi.new.kmeans7.optsil.fidg, file="vibi_new_indval_kmeans7_optsil_results_04182022.csv")
```

## Run analyses above just on “Ruderal Wet-Mesic Thicket…” type to try and break out subtypes

``` r
vibi.new.thicket <- dbGetQuery(conn, "SELECT * FROM nr_misc.plot_data_for_veg_comm_analysis_newest_sample
WHERE plot_id IN (SELECT plot_id FROM nr_misc.plot_data_for_upland_vibi_development WHERE group_description = 'Ruderal Wet-Mesic Shrubland and Thicket');")
```

## Reshape into site X species matrix

``` r
vibi.matrix.new.thicket <- cast(vibi.new.thicket, plot_id ~ acronym, value='abundance', fun="mean")
```

## Change NAs to 0

``` r
vibi.matrix.new.thicket[is.na(vibi.matrix.new.thicket)] <- 0
```

``` r
##Hellinger standardization and euclidean distance matrix
vibi.matrix.new.thicket.hell <- decostand(vibi.matrix.new.thicket, "hel", margin = 1)
vibi.matrix.new.thicket.hell.dist <- vegdist(vibi.matrix.new.thicket.hell, "euc")
```

## K-means clustering

``` r
vibi.new.thicket.hell.dist.kmeans2 <- kmeans(vibi.matrix.new.thicket.hell, centers=2, nstart=100)
```

## Run indicator species analysis on k-means clusters

``` r
vibi.new.thicket.kmeans2.indval.out = indval(vibi.matrix.new.thicket[,c(-1)], vibi.new.thicket.hell.dist.kmeans2$cluster, numitr=10000)

vibi.new.thicket.kmeans2.gr <- vibi.new.thicket.kmeans2.indval.out$maxcls[vibi.new.thicket.kmeans2.indval.out$pval <= 0.05]
vibi.new.thicket.kmeans2.iv <- vibi.new.thicket.kmeans2.indval.out$indcls[vibi.new.thicket.kmeans2.indval.out$pval <= 0.05]
vibi.new.thicket.kmeans2.pv <- vibi.new.thicket.kmeans2.indval.out$pval[vibi.new.thicket.kmeans2.indval.out$pval <= 0.05]
vibi.new.thicket.kmeans2.fr <- apply(vibi.matrix.new.thicket[,c(-1)] > 0, 2, sum)[vibi.new.thicket.kmeans2.indval.out$pval <= 0.05]
vibi.new.thicket.kmeans2.fidg <- data.frame(group=vibi.new.thicket.kmeans2.gr, indval=vibi.new.thicket.kmeans2.iv, pvalue=vibi.new.thicket.kmeans2.pv, freq=vibi.new.thicket.kmeans2.fr)
vibi.new.thicket.kmeans2.fidg <- vibi.new.thicket.kmeans2.fidg[order(vibi.new.thicket.kmeans2.fidg$group, -vibi.new.thicket.kmeans2.fidg$indval),]
```

## Write IndVal results from K-means to csv \*\*NOT RUN

``` r
write.csv(vibi.new.thicket.kmeans2.fidg, file="vibi_new_thicket_indval_kmeans2_results_04192022.csv")
```

## Combine K-means clusters with plot matrix

``` r
write.csv(vibi.new.thicket.kmeans2.cluster, file="vibi_new_thicket_kmeans2_results_04192022.csv")
```

## Run optsil function on K-means result from above

``` r
vibi.new.thicket.kmeans2.optsil <- optsil(kmeans2_thicket_cluster, vibi.matrix.new.thicket.hell.dist)
```

## Combine cluster results from meadow analysis with spp matrix

``` r
kmeans2_thicket_optsil_cluster <- vibi.new.thicket.kmeans2.optsil$clustering 
vibi.new.thicket.kmeans2.optsil.cluster <- cbind(vibi.matrix.new.thicket, kmeans2_thicket_optsil_cluster)
```

``` r
write.csv(vibi.new.thicket.kmeans2.optsil.cluster, file="vibi_new_thicket_kmeans2_optsil_results_04192022.csv")
```

## Run indicator species analysis on K-means + optsil results

``` r
vibi.new.thicket.kmeans2.optsil.indval.out = indval(vibi.matrix.new.thicket[,c(-1)], vibi.new.thicket.kmeans2.optsil$clustering, numitr=10000)

vibi.new.thicket.kmeans2.optsil.gr <- vibi.new.thicket.kmeans2.optsil.indval.out$maxcls[vibi.new.thicket.kmeans2.optsil.indval.out$pval <= 0.05]
vibi.new.thicket.kmeans2.optsil.iv <- vibi.new.thicket.kmeans2.optsil.indval.out$indcls[vibi.new.thicket.kmeans2.optsil.indval.out$pval <= 0.05]
vibi.new.thicket.kmeans2.optsil.pv <- vibi.new.thicket.kmeans2.optsil.indval.out$pval[vibi.new.thicket.kmeans2.optsil.indval.out$pval <= 0.05]
vibi.new.thicket.kmeans2.optsil.fr <- apply(vibi.matrix.new.thicket[,c(-1)] > 0, 2, sum)[vibi.new.thicket.kmeans2.optsil.indval.out$pval <= 0.05]
vibi.new.thicket.kmeans2.optsil.fidg <- data.frame(group=vibi.new.thicket.kmeans2.optsil.gr, indval=vibi.new.thicket.kmeans2.optsil.iv, pvalue=vibi.new.thicket.kmeans2.optsil.pv, freq=vibi.new.thicket.kmeans2.optsil.fr)
vibi.new.thicket.kmeans2.optsil.fidg <- vibi.new.thicket.kmeans2.optsil.fidg[order(vibi.new.thicket.kmeans2.optsil.fidg$group, -vibi.new.thicket.kmeans2.optsil.fidg$indval),]
```

``` r
write.csv(vibi.new.thicket.kmeans2.optsil.fidg, file="vibi_new_thicket_indval_kmeans2_optsil_results_04192022.csv")
```

Note that the `echo = FALSE` parameter was added to the code chunk to
prevent printing of the R code that generated the plot.
