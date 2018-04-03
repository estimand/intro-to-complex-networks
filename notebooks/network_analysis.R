#' ---
#' title: "Network analysis using R"
#' author: "Gianluca Campanella"
#' date: ""
#' output:
#'   pdf_document:
#'     latex_engine: xelatex
#'     toc: true
#'     highlight: pygments
#' documentclass: article
#' classoption: oneside
#' geometry: margin=2cm
#' fontsize: 12pt
#' lang: british
#' mainfont: "Lucida Sans OT"
#' sansfont: "Lucida Sans OT"
#' monofont: "Lucida Console DK"
#' urlcolor: blue
#' ---

#+ include=FALSE
rm(list=ls())

set.seed(42)

library(knitr)
opts_chunk$set(dev="png", dpi=300)

#' # Setup
#'
#' Let's start by loading the packages we'll be using later.
#'
#' You can install any missing packages using `install.packages` or the RStudio
#' GUI.

#+ message=FALSE
library("corpcor")
library("igraph")
library("poweRlaw")
library("tidyverse")

#' # Brief introduction to `igraph`
#'
#' [`igraph`](http://igraph.org/) is an efficient network analysis toolkit
#' implemented in C, with a convenient R interface.
#'
#' We'll briefly review some of the main functions of `igraph`, before delving
#' into an '--omics' application.
#' For a longer introduction to `igraph`, check out Katherine Ognyanova's
#' [*Network Analysis and Visualization with R and `igraph`*](http://kateto.net/netscix2016).
#'
#' Graphs can be specified using the function `graph`.
#' For example, let's create an undirected graph with three vertices and three
#' edges.

g1 <- graph(edges = c(1, 2, 2, 3, 3, 1), directed = FALSE)

#' Graphs can also be specified using `graph_from_literal`.

g2 <- graph_from_literal(A-B-C-D-E-F, A-C-E, B-F:G:H)

#' **Small** graphs can be plotted using `plot`.

#+ fig.width=5, fig.height=5
plot(g2)

#' Vertices and edges can be accessed using the functions `V` and `E`,
#' respectively.

V(g2)

E(g2)

#' Alternatively, the adjacency matrix representation is also available.

g2[]

#' # Correlation networks
#'
#' We'll use data from the study *Autoantibody profiling of sera from
#' individuals with Systemic Lupus Erythematosus (SLE) by protein microarray*
#' ([E-MTAB-5900](https://www.ebi.ac.uk/arrayexpress/experiments/E-MTAB-5900/)
#' on [ArrayExpress](https://www.ebi.ac.uk/arrayexpress/)).
#'
#' The original data can be found in the `input_data` folder.
#' The dataset we'll be using was generated using the script `build_dataset.R`,
#' and is stored in the `datasets` folder.
#'
#' Let's start by loading the data.

sle <- read_rds("datasets/sle_proteomics.rds")

#' Next, we perform some basic exploratory analyses.

sle %>%
    select(-starts_with("P")) %>%
    head()

sle %>%
    group_by(country, ethnicity) %>%
    summarize(
        n = n(),
        n_cases = sum(case),
        cases_pct = mean(case)
    )

#' We now focus on the protein-protein correlation network (irrespective of
#' disease status), which we estimate using the shrinkage estimator implemented
#' in the `corpcor` package.

cor_all <- sle %>%
           select(starts_with("P")) %>%
           as.matrix() %>%
           cor.shrink()

#' We have a few options to turn this correlation matrix into a graph:
#'
#' * Use the correlation coefficients directly as weights
#' * Use the absolute value of the correlation coefficients as weights
#' * Dichotomise with respect to an arbitrary (absolute) correlation threshold
#'
#' We'll illustrate the last method using an arbitrary threshold based on the
#' quantiles of the absolute correlation coefficient distribution.

abs_cors <- abs(cor_all[upper.tri(cor_all)])

#+ fig.width=5, fig.height=5
ggplot(tibble(cor = abs_cors), mapping = aes(x = cor)) +
    geom_density() +
    theme_bw()

adj_matrix <- ifelse(abs(cor_all) > quantile(abs_cors, 0.95), 1, 0)

g3 <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")

#' The graph we've obtained is not necessarily connected.
#'
#' Let's define a function `largest_component` to extract the largest connected
#' component of the graph, and apply it to our object.

largest_component <- function(graph, ...) {
    cs <- components(graph, ...)
    induced_subgraph(graph, which(cs$membership == which.max(cs$csize)))
}

lc3 <- largest_component(g3)

#' We'll also define a function to plot the degree distribution (in log-log
#' scale), and a function to fit a discrete power law and a discrete log-normal
#' distribution to the degrees (using the package `poweRlaw`).

plot_degree_distribution <- function(graph) {
    ggplot(tibble(degree = degree(graph)), mapping = aes(x = degree)) +
        geom_point(stat = "count") +
        scale_x_continuous("Degree", trans = "log10") +
        scale_y_continuous("Frequency", trans = "log10") +
        theme_bw()
}

fit_degree_distribution <- function(graph) {
    degrees <- degree(graph)

    pl <- displ$new(degrees)
    pl$setXmin(estimate_xmin(pl))

    lnorm <- dislnorm$new(degrees)
    lnorm$setXmin(pl$getXmin())
    lnorm$setPars(estimate_pars(lnorm))

    plot(pl, pch = 20, cex = 0.5)
    lines(pl, col = "red", lwd = 2)
    lines(lnorm, col = "blue", lwd = 2)

    compare_distributions(pl, lnorm)$p_two_sided
}

#' Let's apply these functions to `lc3`.

#+ fig.width=5, fig.height=5
plot_degree_distribution(lc3)

#+ fig.width=5, fig.height=5
fit_degree_distribution(lc3)

#' It appears that the discrete log-normal distribution provides a better
#' fit in this case.
#'
#' # Differential correlation networks
#'
#' In case-control studies we're often interested in identifying correlates of
#' disease status.
#' We can extend this idea to networks by investigating differences between
#' the protein-protein correlation networks observed under each condition.

cor_cases <- sle %>%
             filter(case) %>%
             select(starts_with("P")) %>%
             as.matrix() %>%
             pcor.shrink()

cor_controls <- sle %>%
                filter(!case) %>%
                select(starts_with("P")) %>%
                as.matrix() %>%
                pcor.shrink()

cor_diff <- cor_cases - cor_controls

#' We'll again dichotomise with respect to an arbitrary threshold based on the
#' quantiles of the absolute correlation difference distribution.

abs_diffs <- abs(cor_diff[upper.tri(cor_diff)])

#+ fig.width=5, fig.height=5
ggplot(tibble(cor = abs_diffs), mapping = aes(x = cor)) +
    geom_density() +
    theme_bw()

adj_matrix <- ifelse(abs(cor_diff) > quantile(abs_diffs, 0.95), 1, 0)

g4 <- graph_from_adjacency_matrix(adj_matrix, mode = "undirected")

#' We'll also retrieve the largest component, and plot its degree distribution.

lc4 <- largest_component(g4)

#+ fig.width=5, fig.height=5
plot_degree_distribution(lc4)

#' Finally, we'll fit a discrete power law and a discrete log-normal
#' distribution to the degrees.

#+ fig.width=5, fig.height=5
fit_degree_distribution(lc4)

#' It appears that both distributions provide similar fits in this case.
#'
#' We'll now compute different centrality measures.

centralities <- tibble(
    id = names(V(lc4)),
    degree = degree(lc4, normalized = TRUE),
    closeness = closeness(lc4, normalized = TRUE),
    betweenness = betweenness(lc4, normalized = TRUE),
    eigen = eigen_centrality(lc4)$vector
)

#' First of all, note how different measures of centrality are correlated.

cor(centralities %>% select(-id), method = "spearman")

#' Let's investigate some of the more 'central' vertices.

centralities %>%
    arrange(desc(degree))

#' Alternatively, we can compare ranks.

centrality_ranks <- centralities %>%
                    mutate_if(is.numeric, funs(min_rank(desc(.))))

centrality_ranks %>%
    arrange(degree)

#' Let's focus on proteins that appear amongst the top 3 for any centrality
#' measure.

centrality_ranks %>%
    filter_if(is.numeric, any_vars(. <= 3))

#' Consulting the
#' [array design file](https://www.ebi.ac.uk/arrayexpress/files/A-MTAB-620/A-MTAB-620.adf.txt)
#' on ArrayExpress, we arrive at the following table:
#'
#' | `id`        | UniProt                                          | Gene name |
#' |-------------|--------------------------------------------------|-----------|
#' | P000131_1   | [Q13526](https://www.uniprot.org/uniprot/Q13526) | *PIN1*    |
#' | P000165_1   | [Q92934](https://www.uniprot.org/uniprot/Q92934) | *BAD*     |
#' | P000201_1   | [O75081](https://www.uniprot.org/uniprot/O75081) | *CBFA2T3* |
#' | P000213_1   | [Q9UHB7](https://www.uniprot.org/uniprot/Q9UHB7) | *AFF4*    |
#' | P001280_1   | [P54252](https://www.uniprot.org/uniprot/P54252) | *ATXN3*   |
#' | P002239.2_1 | [Q9H999](https://www.uniprot.org/uniprot/Q9H999) | *PANK3*   |
#' | P003029_1   | [Q5JRK9](https://www.uniprot.org/uniprot/Q5JRK9) | *PAGE2B*  |
#' | P003218_1   | [Q06330](https://www.uniprot.org/uniprot/Q06330) | *RBPJ*    |
#'
#' Whilst interpretation of these results is beyond the scope of this workshop,
#' we note that:
#'
#' * *PIN1* has been found to be abnormally activated in SLE by
#'   [Wei et al.](https://doi.org/10.1002/art.39741)
#' * *PANK3* codes for isoform 3 of pantothenate kinase, the first enzyme in the
#'   CoA biosynthetic pathway.
#'   CoA deficiency has been hypothesised to be involved in the pathogenesis of
#'   SLE by [Leung](https://doi.org/10.1016/j.mehy.2004.01.017).
#' * *RBPJ* is part of the Notch signalling pathway, which has been explored as
#'   a potential target for SLE treatment by
#'   [Teachey et al.](https://doi.org/10.1182/blood-2007-05-087353)

