
set.seed(42)

library(igraph)

n <- 500
p <- 2 / n

igraph.options(list(
    vertex.size = 1,
    vertex.color = "#263238",
    vertex.label = NA,
    edge.color = "#263238"
))

pdf("models.pdf", width = 14, height = 7)
par(mar = c(1, 1, 2, 1), mfrow = c(1, 2))

G1 <- erdos.renyi.game(n, p)
G1$layout <- layout.fruchterman.reingold

G2 <- barabasi.game(n, directed = FALSE)
G2$layout <- layout.fruchterman.reingold

plot(G1, main = "Random")
plot(G2, main = "Prefential attachment")

dev.off()

