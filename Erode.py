import gudhi as gd
import networkx as nx
import sys
import numpy as np
from subprocess import call


tetra_to_node = dict()
node_to_tetra = dict()
filtration = dict()
points = dict()


def readOFF(filename):
    stree = gd.SimplexTree()
    f = open(filename, 'r')
    header = f.readline().strip('\n')
    if header != 'OFF':
        print('Invalid OFF file')
    num = f.readline().strip('\n').split(' ')
    num_v = int(num[0])
    num_f = int(num[1])
    for ii in range(num_v):
        v = f.readline().strip('\n').split(" ")
        ver = list()
        # print(v)
        for c in v:
            if c != '':
                ver.append(float(c))
        points[ii] = (ver[0], ver[1], ver[2])
        # point_to_ind[(ver[0], ver[1], ver[2])] = ii
    for ii in range(num_f):
        line = f.readline().strip('\n').split(" ")
        # line = line[:-1]
        face_ind = list()
        for jj in range(1, int(float(line[0])) + 1):
            face_ind.append(int(float(line[jj])))
        stree.insert(face_ind)
        key = tuple(sorted(face_ind))
        tetra_to_node[key] = ii
        node_to_tetra[ii] = key
    f.close()
    print('Read ', num_v, 'vertices, ', num_f, ' Faces')
    return stree, num_f


def createGraph(stree: gd.SimplexTree):
    G = nx.Graph()
    for f in stree.get_skeleton(1):
        if len(f[0]) != 2:
            continue
        adj = stree.get_cofaces(f[0], 1)
        if len(adj) == 1:
            continue
        edge_list = list()
        for i in range(len(adj)):
            u = tetra_to_node[tuple(adj[i][0])]
            for j in range(i+1, len(adj)):
                v = tetra_to_node[tuple(adj[j][0])]
                edge_list.append((u, v))
        G.add_edges_from(edge_list)
    return G


def read_filtration(filein):
    f = open(filein, 'r')
    dim = int(f.readline().strip('\n'))
    each_dim = 1
    for i in range(dim):
        each_dim = each_dim * int(f.readline().strip('\n'))
    for ii in range(each_dim):
        filtration[ii] = float(f.readline().strip('\n'))
    print('Read ', len(filtration), ' Values')


def get_filtration(f: list):
    val = max([filtration[vind] for vind in f])
    return val


if __name__ == '__main__':
    data = sys.argv[1]
    delta = 10
    if len(sys.argv) > 2:
        delta = float(sys.argv[2])
    else:
        print('Setting delta to default 10')
    st, numf = readOFF(data + '/' + data + '_pruned_jacobi.off')
    read_filtration(data + '/' + data + '.txt')
    print('Creating dual graph to find connected components ...')
    graph = createGraph(st)
    # groups = nx.connected_components(graph)
    # print(numf, len(graph.nodes))
    # print('Found ', len(list(groups)), ' connected components ...')
    f_to_export = list()
    print('Eroding with threshold ', delta)
    count = 0
    # print([len(c) for c in nx.connected_components(graph)])
    for g in nx.connected_components(graph):
        # print(list(g))
        min_val = min([get_filtration(node_to_tetra[n]) for n in g])
        print(count, ' : ', min_val)
        for n in g:
            facet = node_to_tetra[n]
            if get_filtration(facet) > (min_val + delta):
                f_to_export.append(facet)
        count += 1
    if len(f_to_export) > 0:
        np.savetxt(data + '/' + data + '_pruned_eroded_jacobi.txt', f_to_export, fmt='%d')
        call(['python3', 'WriteVTP.py', data, str(delta)])
    else:
        print('No face to write')