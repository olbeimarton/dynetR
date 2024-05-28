# dynetR
R implementation of the DyNet network rewiring algorithm as described in [Goenawan et al., 2016](https://academic.oup.com/bioinformatics/article/32/17/2713/2450724).

The main function, `dynetR` takes a list of weighted edge lists or adjacency matrices to calculate rewiring values as described in the paper.

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
                                                    
 `mList <- list(m1, m2, m3, m4)`
 
 `dynetR(mList)`

The main function, `dynetR`, returns a list. The first item in the list is a dataframe with the node rewiring values, and the second element is a simple plot visualising the union network and the rewiring values.
