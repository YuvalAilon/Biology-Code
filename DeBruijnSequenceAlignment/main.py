def kmers_of(sequence, k):
    if len(sequence) < k:
        raise Exception("k must be less than or equal to the length of the sequence")

    kmers = []
    for i in range(len(sequence) - (k - 1)):
        kmers.append(sequence[i: i + k])

    return kmers


def l_r_k_minus_1_mers(kmer):
    return [
        kmer[0: len(kmer) - 1],
        kmer[1: len(kmer)]
    ]


def eulerian_path(graph):
    starting_node = ""
    for node in graph:
        if graph[node]["in_degree"] - len(graph[node]["edges"]) == -1:
            starting_node = node
    print(starting_node)


def de_bruijn_alignment(fragments, k):
    length_k_mers = lambda frags: kmers_of(frags, k)
    kmers_from_fragments = map(length_k_mers, fragments)
    kmers = []
    for list in kmers_from_fragments:
        kmers += list
    print(kmers)

    de_bruijn_graph = construct_de_bruijn_graph(kmers)

    eulerian_path(de_bruijn_graph)

    for key in de_bruijn_graph:
        print(key + " : " + str(de_bruijn_graph.get(key)))


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
    return de_bruijn_graph


de_bruijn_alignment(["TCCAGG", "TCACGTCCA"], 3)
