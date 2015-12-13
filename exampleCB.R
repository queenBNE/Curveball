# Erdos Renyi tests
# Create an Erdos Renyi random network
require('igraph')
source('curveball.versions.R')
g <- erdos.renyi.game(20, 0.1, directed = T)
I <- get.adjlist(g, mode = "out")

# Randomize simple directed network
sample <- curveball.randomise.dir(I, 1000, F, T)
h <- graph.adjlist(sample[[1000]])
prod(degree(g, mode = "in") == degree(h, mode = "in"))
prod(degree(g, mode = "out") == degree(h, mode = "out"))
is.simple(h)

# Randomize directed network 
sample <- curveball.randomise.dir(I, 1000, T, T)
h <- graph.adjlist(sample[[1000]])
prod(degree(g, mode = "in") == degree(h, mode = "in"))
prod(degree(g, mode = "out") == degree(h, mode = "out"))
is.simple(h)

# Create undirected Erdos Renyi random network 
g <- erdos.renyi.game(20, 0.2)
I <- get.adjlist(g)
sample <- curveball.randomise.undir(I, 1000, T)
h <- graph.adjlist(sample[[1000]], mode = "all")
prod(degree(g, mode = "in") == degree(h, mode = "in"))
prod(degree(g, mode = "out") == degree(h, mode = "out"))
is.simple(h)


