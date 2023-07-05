# demo

``` r
library(HCCPS)

demo_data = data.frame(sample1 = runif(4,1,10),sample2 = runif(4,1,10),sample3 = runif(4,1,10))
rownames(demo_data) = HCCPS:::getID(c("SQSTM1","BIRC5","HSPB8","NCF2"))

score(demo_data)
```

    ##             sample1    sample2     sample3
    ## score1  4.078288178 2.11790009 3.032441887
    ## score2  1.929340123 0.33282570 1.172263634
    ## score3  0.000000000 0.00000000 0.000000000
    ## score4  0.000000000 0.00000000 0.000000000
    ## score5  0.000000000 0.00000000 0.000000000
    ## score6  0.000000000 0.00000000 0.000000000
    ## score7  0.000000000 0.00000000 0.000000000
    ## score8  0.000000000 0.00000000 0.000000000
    ## score9  0.000000000 0.00000000 0.000000000
    ## score10 0.000000000 0.00000000 0.000000000
    ## score11 0.000000000 0.00000000 0.000000000
    ## score12 0.000000000 0.00000000 0.000000000
    ## score13 0.000000000 0.00000000 0.000000000
    ## score14 0.000000000 0.00000000 0.000000000
    ## score15 0.000000000 0.00000000 0.000000000
    ## score16 0.000000000 0.00000000 0.000000000
    ## score17 0.000000000 0.00000000 0.000000000
    ## score18 0.000000000 0.00000000 0.000000000
    ## score19 0.000000000 0.00000000 0.000000000
    ## score20 0.000000000 0.00000000 0.000000000
    ## score21 0.000000000 0.00000000 0.000000000
    ## score22 0.007067867 0.01057132 0.003580759
    ## score23 0.000000000 0.00000000 0.000000000
    ## score24 0.185229104 0.03195342 0.112544875
    ## score25 1.253020410 0.21615546 0.761332977
    ## score26 0.000000000 0.00000000 0.000000000
    ## score27 0.000000000 0.00000000 0.000000000
    ## score28 0.000000000 0.00000000 0.000000000
    ## score29 0.000000000 0.00000000 0.000000000
    ## score30 1.909105010 0.32933499 1.159968815
    ## score31 0.000000000 0.00000000 0.000000000
    ## score32 0.000000000 0.00000000 0.000000000
    ## score33 0.000000000 0.00000000 0.000000000
    ## score34 0.000000000 0.00000000 0.000000000
    ## score35 0.000000000 0.00000000 0.000000000
    ## score36 0.000000000 0.00000000 0.000000000
    ## score37 0.000000000 0.00000000 0.000000000
    ## score38 0.000000000 0.00000000 0.000000000
    ## score39 0.000000000 0.00000000 0.000000000
    ## score40 0.000000000 0.00000000 0.000000000
    ## score41 0.000000000 0.00000000 0.000000000
    ## score42 0.000000000 0.00000000 0.000000000
    ## score43 0.000000000 0.00000000 0.000000000
    ## score44 1.000000000 1.00000000 1.000000000
    ## score45 1.000000000 1.00000000 1.000000000
    ## score46 1.790158687 1.10567043 1.424480291
    ## score47 1.000000000 1.00000000 1.000000000
    ## score48 1.256822170 1.10226633 1.347673917

``` r
score(demo_data[, 1, drop=FALSE])
```

    ##             sample1
    ## score1  4.078288178
    ## score2  1.929340123
    ## score3  0.000000000
    ## score4  0.000000000
    ## score5  0.000000000
    ## score6  0.000000000
    ## score7  0.000000000
    ## score8  0.000000000
    ## score9  0.000000000
    ## score10 0.000000000
    ## score11 0.000000000
    ## score12 0.000000000
    ## score13 0.000000000
    ## score14 0.000000000
    ## score15 0.000000000
    ## score16 0.000000000
    ## score17 0.000000000
    ## score18 0.000000000
    ## score19 0.000000000
    ## score20 0.000000000
    ## score21 0.000000000
    ## score22 0.007067867
    ## score23 0.000000000
    ## score24 0.185229104
    ## score25 1.253020410
    ## score26 0.000000000
    ## score27 0.000000000
    ## score28 0.000000000
    ## score29 0.000000000
    ## score30 1.909105010
    ## score31 0.000000000
    ## score32 0.000000000
    ## score33 0.000000000
    ## score34 0.000000000
    ## score35 0.000000000
    ## score36 0.000000000
    ## score37 0.000000000
    ## score38 0.000000000
    ## score39 0.000000000
    ## score40 0.000000000
    ## score41 0.000000000
    ## score42 0.000000000
    ## score43 0.000000000
    ## score44 1.000000000
    ## score45 1.000000000
    ## score46 1.790158687
    ## score47 1.000000000
    ## score48 1.256822170