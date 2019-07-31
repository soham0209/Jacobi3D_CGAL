import gudhi as gd
import networkx as nx
import sys
import numpy as np

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


def boundary(f: list):
    bd_edges = list()
    for i in range(len(f)):
        bd_edges.append(sorted(f[:i] + f[i+1:]))
    return bd_edges


def getnbrs(f: list, stree: gd.SimplexTree):
    nbrs = list()
    bd_edges = boundary(f)
    for e in bd_edges:
        adj = stree.get_cofaces(e, 1)
        for tri in adj:
            if tri[0] != f:
                nbrs.append(tri[0])
    return nbrs


def read_filtration(filein):
    f = open(filein, 'r')
    dim = int(f.readline().strip('\n'))
    each_dim = 1
    for i in range(dim):
        each_dim = each_dim * int(f.readline().strip('\n'))
    for ii in range(each_dim):
        filtration[ii] = abs(float(f.readline().strip('\n')))
    print('Read ', len(filtration), ' Values')


def get_filtration(f: list):
    val = min([filtration[vind] for vind in f])
    return val


if __name__ == '__main__':
    data = sys.argv[1]
    pers_file = data + '/' + data + '.txt'
    delta = 10
    argc = len(sys.argv)
    if argc > 2:
        delta = float(sys.argv[2])
    else:
        print('Setting delta to default 10')
    if argc == 4:
        resample_flag = sys.argv[3]
        if resample_flag == '-r':
            pers_file = data + '/' + data + '_resampled.txt'
    st, numf = readOFF(data + '/' + data + '_pruned_walked.off')
    read_filtration(pers_file)
    faces = list()
    for f in st.get_skeleton(2):
        if len(f[0]) != 3:
            continue
        faces.append(f[0])
    sorted(faces, key = get_filtration)
    min_val = get_filtration(faces[0])
    print(min_val)
    rem_faces = list()
    for f in faces:
        if get_filtration(f) > (min_val + delta):
            rem_faces.append(f)
    print(len(rem_faces))
    np.savetxt(data + '/' + data + '_final.txt', rem_faces)
