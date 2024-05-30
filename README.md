# dynetR
R implementation of the DyNet network rewiring algorithm as described in [Goenawan et al., 2016](https://academic.oup.com/bioinformatics/article/32/17/2713/2450724).

The main function, `dynetR` takes a named list of weighted edge lists or adjacency matrices to calculate rewiring values as described in the paper.

To install the package use
`devtools::install_github('olbeimarton/dynetR')`

Example

`m1 <- matrix(rnorm(36), nrow = 6,  dimnames = list(c("A", "B","C","D","E","F"),
                                                    c("A", "B","C","D","E","F")))`
                                                    
`m2 <- matrix(rnorm(36), nrow = 6,  dimnames = list(c("A", "B","C","D","E","F"),
                                                    c("A", "B","C","D","E","F")))`
                                                    
 `m3 <- matrix(rnorm(36), nrow = 6,  dimnames = list(c("A", "B","C","D","E","F"),
                                                    c("A", "B","C","D","E","F")))`
                                                    
 `m4 <- matrix(rnorm(36), nrow = 6,  dimnames = list(c("A", "B","C","D","E","F"),
                                                    c("A", "B","C","D","E","F")))`

`m5 <- matrix(rnorm(9), nrow = 3,  dimnames = list(c("A", "B","C"), c("A", "B","C")))`
                                                    
 `mList <- list('a' = m1, 'b' = m2, 'c' = m3, 'd' = m4, 'e' = m5)`
 
 `output<-dynetR(mList)`

The main function, `dynetR`, returns a list of items. The first item in the list is a dataframe with the node rewiring values: 
```
> output
[[1]]
  name rewiring degree degree_corrected_rewiring
1    A 11.66336      7                  1.666195
2    B 17.72024      7                  2.531463
3    C 12.43599      7                  1.776569
4    D 18.47138      7                  2.638768
5    E 13.48869      7                  1.926956
6    F 13.05306      5                  2.610611
```
The second element is a simple plot visualising the union network and the rewiring values.

![output plot](example_plot.png)


A typical follow-up analysis usually involves a more in-depth study of the most rewired nodes, to understand why they are rewired. 

The `small_multiples_plot` function returns a small multiples plot, focussing on the interactions of a single node across the compared network states. The function requires two parameters, the input adjacency matrices and the name of the focus node, `F` in this example.

`small_multiples_plot(mList, 'F')`

![output small multiples plot](small_multiples_example.png)
