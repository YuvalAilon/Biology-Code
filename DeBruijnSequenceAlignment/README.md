# De-Bruijn Sequence Assembly
De-Bruijn Assembly works by splitting every DNA fragment into smaller overlapping fragments, and finding similarities with a graph
## k-mers
A k-mer is a substring of size k contained within a biological sequence. For example all 4-mers of `ACCGCTC` are `ACCG`, `CCGC`, `CGCT`, and `GCTC`

When finding similarities k-mer size has tradeoffs. shorter k-mers allow for stringing together of sequences with less overlap, but may lead to ambiguities or incorrect matching
## Building the Graph
Every sequence is first divided into its k-mers. For example `ACCGCTC` with `k = 4` is split into `ACCG`, `CCGC`, `CGCT`, and `GCTC`. Then, every one of the resulting k-mers is split into it's k-1-mers, which are connected by an edge on the graph. For example `ACCG` would result in an `ACC` node connected to `CCG` node. The entire sequence would produce this graph:
```
ACC -> CCG -> CGC -> GCT -> CTC
```
Since edges of this graph must be unique, several DNA fragments can build one long graph together with a traversable Eulerian Path. For example `GCTCTATA` could "latch on" to the `GCT -> CTC` node producing the following larger graph:
```
ACC -> CCG -> CGC -> GCT -> CTC -> TCT -> CTA -> TAT -> ATA
```
## Reconstructing a Sequence
To reconstruct a sequnce from a De-Bruijn graph, first a valid Eulerian Path must be found, and then the nodes reconstructed using a sliding window. With the example path above,
a valid Eulerian Path is trivial, and the sequence is reconstructed in this fashion
```
   1 2 3 4 5 6 7 8 9 A B
1: A C C
2: | C C G
3: |   C G C
4: |   | G C T
5: |   |   C T C
6: |   |   | T C T
7: |   |   |   C T A
8: |   |   |   | T A T
9: |   |   |   |   A T A
   |   |   |   |   |   |
   A C C G C T C T A T A
```

From the two pieces `ACCGCTC` and `GCTCTATA`

## Limitations
Sometimes a graph contains multiple valid Eulerian Paths, for example:
!["A:](/DeBruijnSequenceAlignment/De-Bruijn-Graph-1.png)

can refer to either `AAATAAAGCAAA` or `AAAGCAAATAA`. One way to avoid this is to use larger k values, since there are more possibilites for nodes, and less chance for an overlap.

A limitation of my implementation is that it assumes 100% perfect data. More complex sequence aligners may weigh the edges based on how often each k-mer appears and prune edges with sufficiently small weight. I hope to implement this in the future. 