
set.seed(42)

library(igraph)

igraph.options(list(
    vertex.size = 1,
    vertex.color = "#263238",
    vertex.label = NA,
    edge.color = "#263238"
))

G <- barabasi.game(1000, directed = FALSE)
G$layout <- layout.fruchterman.reingold(G)

pdf("degree_dist.pdf", width = 14, height = 7)
par(mfrow = c(1, 2))

par(mar = rep(0, 4))
plot(G)

par(mar = c(5, 5, 1, 1))
x <- seq(1, 40, 0.01)
dd <- degree.distribution(G)
plot(x + 1, x**-3, type = "l", lty = 3, lwd = 2, col = "#FF5722",
     xlim = range(x), ylim = c(0, 1), axes=FALSE,
     xlab = "Degree", ylab = "Probability")
points(2:length(dd), dd[-1], pch = 20)
axis(1, at = c(1, seq(5, max(x), 5)))
axis(2, at = seq(0, 1, 0.1))

dev.off()

