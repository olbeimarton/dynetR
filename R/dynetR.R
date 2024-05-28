#' @title Helper function to set rownames equal to colnames for a single matrix
#'
#' @description
#' Helps to set matching row and colum nnames
#' @param mat A single input matrix
#'

.set_rownames_to_colnames <- function(mat) {
  if (!is.null(colnames(mat))) {
    rownames(mat) <- colnames(mat)
  }
  return(mat)
}

#' @title Helper function to format weighted edge list dataframe
#'
#' @description Formats edge list to adjacency matrix
#'
#' @param el An input weighted edge list with from, to, and weight columns
#'
#'
.edgelist_to_adjacency <- function(el) {
  g1 <- graph.data.frame(el)
  gdf1 <- get.adjacency(g1, attr = "weight", sparse = F)
}

#' @title Helper function to format input data, matrix or edge list df
#' @description Formats input list data depending on input type
#' @param input_list A list of matrices or edge lists
#'
#'
format_indata <- function(input_list) {
  # checking if elements are matrices
  if (all(sapply(input_list, is.matrix))) {
    # do stuff
    formatted_matrix_list <- lapply(input_list, .set_rownames_to_colnames)
  } else if (all(sapply(input_list, is.data.frame))) {
    formatted_matrix_list <- lapply(input_list, .edgelist_to_adjacency)
  } else {
    stop("Input must be a matrix or an edge list data frame.")
  }
  return(formatted_matrix_list)
}


#' @title Helper function to expand adjacency matrices to include all unique nodes
#' @description Expands formatted adjacency matrices to matching sizes
#' @param adj_matrices An input list of adjacency matrices
#'
#'
#'
.expand_adjacency_matrices <- function(adj_matrices) {
  # Identify all unique nodes across all adjacency matrices
  all_nodes <- unique(unlist(lapply(adj_matrices, rownames)))

  # Function to expand a single adjacency matrix
  expand_matrix <- function(matrix, all_nodes) {
    current_nodes <- rownames(matrix)
    missing_nodes <- setdiff(all_nodes, current_nodes)

    # Create an expanded matrix with all nodes
    expanded_matrix <- Matrix(0,
      nrow = length(all_nodes), ncol = length(all_nodes),
      dimnames = list(all_nodes, all_nodes), sparse = F
    )

    # Fill in the existing values
    expanded_matrix[current_nodes, current_nodes] <- matrix[current_nodes, current_nodes]

    expanded_matrix <- as.matrix(expanded_matrix)

    return(expanded_matrix)
  }

  # Apply the expand_matrix function to each adjacency matrix in the list
  expanded_matrices <- lapply(adj_matrices, expand_matrix, all_nodes = all_nodes)

  return(expanded_matrices)
}


#' @title Helper function to create a union of multiple adjacency matrices
#' @description Builds a union adjacency matrix from a list of adjacency matrices
#' @param adj_matrices An input list of adjacency matrices
#'

.union_adjacency_matrices <- function(inlist) {
  # List of input adjacency matrices
  matrices <- inlist

  # Get the union of all nodes across all matrices
  all_nodes <- unique(unlist(lapply(matrices, rownames)))

  # Create an empty adjacency matrix for the union
  union_matrix <- matrix(0, nrow = length(all_nodes), ncol = length(all_nodes))
  rownames(union_matrix) <- all_nodes
  colnames(union_matrix) <- all_nodes

  # Populate the union matrix with edges from each input matrix
  for (adj in matrices) {
    nodes <- rownames(adj)
    for (i in seq_along(nodes)) {
      for (j in seq_along(nodes)) {
        if (adj[i, j] > 0) {
          union_matrix[nodes[i], nodes[j]] <- 1
        }
      }
    }
  }

  return(union_matrix)
}



#' Helper function to summarise matrices
.sumMatrices <- function(matrices) {
  # Initialize a variable to hold the result
  result <- matrix(0, nrow = nrow(matrices[[1]]), ncol = ncol(matrices[[1]]))
  # Initialize a counter matrix with the same dimensions as the result matrix
  counter <- matrix(0, nrow = nrow(result), ncol = ncol(result))
  # Loop through the list of matrices and sum them
  for (matrix in matrices) {
    result <- result + matrix
    counter <- counter + (matrix != 0)
  }
  # Divide the result matrix by the counter matrix element-wise
  result <- result / counter
  result[!is.finite(result)] <- 0
  # Return the result
  return(result)
}

#' Helper function to get the square values of the values in the list of matrices
.squareMatrices <- function(matrices) {
  # Initialize a list to hold the squared matrices
  squared_matrices <- list()

  # Loop through the list of matrices and square each one
  for (matrix in matrices) {
    squared_matrices[[length(squared_matrices) + 1]] <- matrix^2
  }

  # Return the list of squared matrices
  return(squared_matrices)
}

#' Helper function to calculate the centroid distances
.centroidDistance <- function(matrices) {
  result <- list()
  for (matrixx in matrices) {
    result[[length(result) + 1]] <- sqrt(rowSums(matrixx))
  }
  return(result)
}


#' @title Calculate the DyNet rewiring values for multiple graphs (adjacency matrices)
#'
#' @description Calculates the Dn rewiring values for the nodes included in the networks
#' as described in `Goenawan et al., 2016.`
#' @import igraph
#' @import readr
#' @import ggplot2
#' @import ggraph
#' @import tidygraph
#' @import dplyr
#' @import tibble
#' @import Matrix
#' @param matrix_list List of adjacency matrices corresponding to the input networks.
#'
#' @export

dynetR <- function(matrix_list) {
  # format input data
  matrix_list <- format_indata(matrix_list)

  # Match matrix sizes
  expanded_matrices <- .expand_adjacency_matrices(matrix_list)

  # Calculate non-zero means of matrices
  nonZeroMean <- .sumMatrices(expanded_matrices)

  # Calculate standardised weights (divided by non-zero means)
  standardised <- lapply(expanded_matrices, "/", nonZeroMean)
  standardisedNoNan <- lapply(standardised, function(x) replace(x, is.nan(x), 0))
  # Calculate centroids of all node states by summing
  standardSums <- Reduce("+", standardisedNoNan)

  # Calculate euclidean distance from node to centroid
  centroid <- standardSums / length(standardisedNoNan)

  # Subtract centroid values from standardised matrix
  standardMinusCentroid <- lapply(standardisedNoNan, "-", centroid)
  sqr <- .squareMatrices(standardMinusCentroid)
  centDist <- .centroidDistance(sqr)
  centDistBound <- do.call(rbind, centDist)
  centFinalDist <- colSums(centDistBound)
  rewiring <- centFinalDist / (length(centDist) - 1)

  # calculate degree and correct for it
  unimatrix <- .union_adjacency_matrices(matrix_list)
  unimatrix_dataframe <- graph_from_adjacency_matrix(unimatrix, mode = "undirected", diag = T) |>
    degree(mode = "all") |>
    as.data.frame() |>
    rownames_to_column()
  # note that in the case of undirected graphs, an edge that starts and ends in the same node
  # increases the corresponding degree value by 2 (i.e. it is counted twice)
  unimatrix_dataframe <- unimatrix_dataframe |>
    rename("name" = "rowname", "degree" = `degree(graph_from_adjacency_matrix(unimatrix, mode = "undirected", diag = T), mode = "all")`) |>
    tibble()

  rewiring <- rewiring |>
    as.data.frame() |>
    tibble::rownames_to_column() |>
    rename("name" = "rowname")
  output_dataframe <- left_join(rewiring, unimatrix_dataframe, by = "name") |>
    mutate(degree_corrected_rewiring = rewiring / degree)

  # plot

  unimatrix_graph <- graph_from_adjacency_matrix(unimatrix, mode = "undirected", diag = T) |>
    as_tbl_graph() |>
    activate(nodes) |>
    left_join(output_dataframe) |>
      ggraph(layout = "nicely") +
      geom_edge_link() +
      geom_node_point(aes(fill = rewiring, size = degree), shape = 21) +
      scale_fill_distiller(palette = "Reds", direction = 1) +
      geom_node_text(aes(label = name, size = degree), repel = T, fontface = "bold", max.overlaps = 15) +
      scale_size_continuous(range = c(6, 10)) +
      theme_graph()

  output <- list(output_dataframe, unimatrix_graph)
  return(output)
}

#' Helper function to pretty print matrices
.mat_op_print <- function(..., width = 0) {
  # get arguments
  args <- list(...)
  chars <- sapply(args, is.character)

  # auxilliary function to create character of n spaces
  spaces <- function(n) paste(rep(" ", n), collapse = "")

  # convert vectors to row matrix
  vecs <- sapply(args, is.vector)
  args[vecs & !chars] <- lapply(args[vecs & !chars], function(v) matrix(v, ncol = 1))

  # convert all non-characters to character with format
  args[!chars] <- lapply(args[!chars], format, width = width)

  # print names as the first line, if present
  arg_names <- names(args)
  if (!is.null(arg_names)) {
    get_title <- function(x, name) {
      if (is.matrix(x)) {
        paste0(name, spaces(sum(nchar(x[1, ])) + ncol(x) - 1 - nchar(name)))
      } else {
        spaces(nchar(x))
      }
    }
    cat(mapply(get_title, args, arg_names), "\n")
  }

  # auxiliary function to create the lines
  get_line <- function(x, n) {
    if (is.matrix(x)) {
      if (nrow(x) < n) {
        spaces(sum(nchar(x[1, ])) + ncol(x) - 1)
      } else {
        paste(x[n, ], collapse = " ")
      }
    } else if (n == 1) {
      x
    } else {
      spaces(nchar(x))
    }
  }

  # print as many lines as needed for the matrix with most rows
  N <- max(sapply(args[!chars], nrow))
  for (n in 1:N) {
    cat(sapply(args, get_line, n), "\n")
  }
}

#' Helper function that prints the results of the subsequent steps
.dynet_verbose <- function(matrix_list) {
  # format input data
  matrix_list <- format_indata(matrix_list)

  # Match matrix sizes
  expanded_matrices <- .expand_adjacency_matrices(matrix_list)

  # get union matrix

  matrix_union <- .union_adjacency_matrices(expanded_matrices)

  # Calculate non-zero means of matrices
  nonZeroMean <- .sumMatrices(expanded_matrices)

  # Calculate standardised weights (divided by non-zero means)
  standardised <- lapply(expanded_matrices, "/", nonZeroMean)
  standardisedNoNan <- lapply(standardised, function(x) replace(x, is.nan(x), 0))

  # Calculate centroids of all node states by summing and subsequent division
  standardSums <- Reduce("+", standardisedNoNan)
  centroid <- standardSums / length(standardisedNoNan)

  # Subtract centroid values from standardised matrix
  standardMinusCentroid <- lapply(standardisedNoNan, "-", centroid)
  sqr <- .squareMatrices(standardMinusCentroid)
  centDist <- .centroidDistance(sqr)
  centDistBound <- do.call(rbind, centDist)
  # print(centDistBound)
  centFinalDist <- colSums(centDistBound)
  print(length(centFinalDist))
  rewiring <- centFinalDist / (length(centDist) - 1)
  print("Input matrices:")
  .mat_op_print(A = matrix_list[[1]], " , ", B = matrix_list[[2]])
  print("Non-zero mean:")
  print(nonZeroMean)
  print("Standardised weights (divided by non-zero means:")
  .mat_op_print(A = standardisedNoNan[[1]], " , ", B = standardisedNoNan[[2]])
  print("Centroid:")
  print(centroid)
  print("Standardised weight - centroid")
  .mat_op_print(A = standardMinusCentroid[[1]], " , ", B = standardMinusCentroid[[2]])
  print("Square of standard minus centroid")
  .mat_op_print(A = sqr[[1]], " , ", B = sqr[[2]])
  print("Square root of standard minus centroid row sums")
  print(centDistBound)
  print("Final rewiring value")
  print(rewiring)
  print(matrix_union)
}
