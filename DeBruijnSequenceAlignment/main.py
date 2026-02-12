# Returns a list of all kmers of a specified length of sequence
def kmers_of(sequence, k):
    if len(sequence) < k:
        raise Exception("k must be less than or equal to the length of the sequence")

    kmers = []
    for i in range(len(sequence) - (k - 1)):
        kmers.append(sequence[i: i + k])

    return kmers


# Returns a list of length 2 with format [left k-1mer, right k-1mer]
def l_r_k_minus_1_mers(kmer):
    return [
        kmer[0: len(kmer) - 1],
        kmer[1: len(kmer)]
    ]


def reconstruct_graph_from_node(graph, node):
    # Since currently bubbles are filtered out, this function is quite simple
    # as I continue to add complexity, this function will become less trivial
    if not node in graph:
        return [node]
    return [node] + reconstruct_graph_from_node(graph, next(iter(graph[node]["edges"])))

# Takes a list of k-1mer nodes and reconstructs them into the final DNA sequence
def dna_sequence_from_nodes(nodes):
    sequence = nodes[0]
    reconstruction = sequence + "\n"
    for i in range(1, len(nodes)):
        sequence += nodes[i][-1]
        reconstruction += " " * i + nodes[i] + "\n"
    return {"sequence": sequence, "reconstruction": reconstruction}


def construct_de_bruijn_graph(kmers):
    de_bruijn_graph = {}
    for kmer in kmers:
        k_minus_1_mers = l_r_k_minus_1_mers(kmer)
        if k_minus_1_mers[0] not in de_bruijn_graph:
            de_bruijn_graph[k_minus_1_mers[0]] = {"edges": {k_minus_1_mers[1]}, "in_degree": 0}
        else:
            de_bruijn_graph[k_minus_1_mers[0]]["edges"].add(k_minus_1_mers[1])
    for key in de_bruijn_graph:
        for node in de_bruijn_graph[key]["edges"]:
            if node in de_bruijn_graph:
                de_bruijn_graph[node]["in_degree"] += 1
                if de_bruijn_graph[node]["in_degree"] > 1:
                    raise Exception(
                        "Bubble on node " + node + ", reassembly stopped. There may be no valid reconstruction, or k-value may be too small")
    return de_bruijn_graph


def de_bruijn_alignment(fragments, k):
    length_k_mers = lambda frags: kmers_of(frags, k)
    kmers_from_fragments = map(length_k_mers, fragments)
    kmers = []
    for list in kmers_from_fragments:
        kmers += list

    de_bruijn_graph = construct_de_bruijn_graph(kmers)

    starting_node = ""
    for node in de_bruijn_graph:
        if de_bruijn_graph[node]["in_degree"] - len(de_bruijn_graph[node]["edges"]) == -1:
            if starting_node == "":
                starting_node = node
            else:
                raise Exception("More than one node with out degree > in degree, no valid Eulerian Path")

    for key in de_bruijn_graph:
        print(key + " : " + str(de_bruijn_graph.get(key)))
    reconstructed_sequence = dna_sequence_from_nodes(reconstruct_graph_from_node(de_bruijn_graph, starting_node))
    print(reconstructed_sequence["reconstruction"])
    print(reconstructed_sequence["sequence"])
    print("^^^ Reconstructed Sequence ^^^")


de_bruijn_alignment(["CGTGTAGCATC", "GCACGCTCAA", "AGCATCCCACATCGCACGCT", "CACGCTCAATCGA"], 6)

"""

For the above example
Reconstructed Strand:
          CGTGTAGCATCCCACATCGCACGCTCAATCGA
From Segments:
    1.    CGTGTAGCATC       
    2.                      GCACGCTCAA
    3.         AGCATCCCACATCGCACGCT
    4.                       CACGCTCAATCGA

"""
