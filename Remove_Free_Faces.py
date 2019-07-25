import gudhi as gd
import sys
import queue
import operator
import numpy as np
from subprocess import call


point_to_ind = dict()
dim = [-1, -1, -1]
rem_faces = list()

def readOFF(filename):
    stree = gd.SimplexTree()
    f = open(filename, 'r')
    header = f.readline().strip('\n')
    if header != 'OFF':
        print('Invalid OFF file')
    num = f.readline().strip('\n').split(' ')
    num_v = int(num[0])
    num_f = int(num[1])
    points = dict()
    for ii in range(num_v):
        v = f.readline().strip('\n').split(" ")
        ver = list()
        # print(v)
        for c in v:
            if c != '':
                ver.append(float(c))
        points[ii] = (ver[0], ver[1], ver[2])
        point_to_ind[(ver[0], ver[1], ver[2])] = ii
    for ii in range(num_f):
        line = f.readline().strip('\n').split(" ")
        # line = line[:-1]
        face_ind = list()
        for jj in range(1, int(float(line[0])) + 1):
            face_ind.append(int(float(line[jj])))
        stree.insert(face_ind)
    f.close()
    print('Read ', num_v, 'vertices, ', num_f, ' Faces')
    dim[0] = max(point_to_ind.items(), key=operator.itemgetter(1))[0][0]
    dim[1] = dim[0]
    dim[2] = dim[0]
    return stree, points, num_f


def WriteOFF(vertices, stree: gd.SimplexTree, tag: dict, outfile, num_faces):
    f = open(outfile, 'w')
    f.write('OFF\n')
    f.write(str(len(vertices)) + ' ' + str(num_faces) + ' ' + '0\n')
    count = 0
    for i in vertices:
        v = list(vertices[i])
        line = str(v[0]) + ' ' + str(v[1]) + ' ' + str(v[2]) + '\n'
        f.write(line)
    for t in stree.get_skeleton(2):
        if len(t[0]) != 3:
            continue
        if tuple(t[0]) in tag:
            continue
        face = t[0]
        rem_faces.append(face)
        line = str(len(face)) + ' ' + str(face[0]) + ' ' + str(face[1]) + ' ' + str(face[2]) + '\n'
        f.write(line)
        count = count + 1
    print('Wrote ' + str(count) + ' faces to ' + outfile)


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


def check_Free(f: list, stree: gd.SimplexTree, tag, ver_dict:dict, dim = dim):
    bd_edges = boundary(f)
    for e in bd_edges:
        edge = [ver_dict[u] for u in e]
        if is_boundary_edge(edge, dim):
            return False
        adj = stree.get_cofaces(e, 1)
        count = 0
        for t in adj:
            if tuple(t[0]) in tag:
                count = count + 1
        if (len(adj) - count) == 1:
            return True
    return False



def Prune(stree: gd.SimplexTree, ver_dict: dict):
    Q = queue.Queue()
    tag = dict()
    for e in stree.get_skeleton(1):
        if len(e[0]) != 2:
            continue
        edge = [ver_dict[u] for u in e[0]]
        adj = stree.get_cofaces(e[0], 1)
        if len(adj) == 1 and not is_boundary_edge(edge, dim):
            Q.put(adj[0][0])
            # tag[tuple(adj[0][0])] = True
    while not Q.empty():
        f = Q.get()
        if tuple(f) in tag:
            continue
        tag[tuple(f)] = True
        for t in getnbrs(f, stree):
            if check_Free(t, stree, tag, ver_dict):
                Q.put(t)
    return tag


def is_boundary_edge(ed: list, dim):
    i = 0
    # print(ed)
    same_ind = set()
    while i < 3:
        if ed[0][i] == ed[1][i]:
            same_ind.add(i)
        i = i + 1
    for i in same_ind:
        if ed[0][i] == 0 or ed[0][i] == dim[i]:
            return True
    return False


if __name__ == '__main__':
    data = sys.argv[1]
    cpx, pts, nF = readOFF(data + '/' + data + '_jacobi.off')
    tags = Prune(cpx, pts)
    print('Total faces: ', nF)
    print('Faces Deleted: ', len(tags))
    faces_to_write = nF - len(tags)
    print('Faces to Write: ', faces_to_write)
    WriteOFF(pts, cpx, tags, data + '/' + data + '_pruned_jacobi.off', faces_to_write)
    np.savetxt(data + '/' + data + '_pruned_jacobi.txt', rem_faces)
    call(['python3', 'WriteVTP.py', data])
