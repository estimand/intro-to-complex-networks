
set.seed(42)

library(igraph)

igraph.options(list(
    vertex.size = 1,
    vertex.label = NA
))

G <- barabasi.game(1000, directed = FALSE)
G$layout <- layout.fruchterman.reingold(G)

sp <- get.shortest.paths(G, 29, 481)$vpath[[1]]

e.cols <- rep("#263238", ecount(G))
e.widths <- rep(1, ecount(G))
for (i in 2:length(sp)) {
    idx <- get.edge.ids(G, sp[(i-1):i])
    e.cols[idx] <- "#FF5722"
    e.widths[idx] <- 2
}

pdf("shortest_path.pdf")
par(mar=rep(0, 4))
plot(G, vertex.color = "#263238", edge.color = e.cols, edge.width = e.widths)
dev.off()

