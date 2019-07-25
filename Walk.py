import gudhi as gd
import sys
import networkx as nx
import queue
import numpy as np


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
        # point_to_ind[(ver[0], ver[1], ver[2])] = ii
    for ii in range(num_f):
        line = f.readline().strip('\n').split(" ")
        # line = line[:-1]
        face_ind = list()
        for jj in range(1, int(float(line[0])) + 1):
            face_ind.append(int(float(line[jj])))
        stree.insert(face_ind)
    f.close()
    print('Read ', num_v, 'vertices, ', num_f, ' Faces')
    return stree, points, num_f


data = sys.argv[1]
st, ver, _ = readOFF(data+'/'+data+'_triangulated.off')
pruned_st, _, _ = readOFF(data + '/' + data + '_pruned_jacobi.off')
tetra_ind = dict()
ind_to_face = dict()
count = 0
G = nx.Graph()
status = dict()
pruned_faces = dict()
for f in pruned_st.get_skeleton(2):
    if len(f[0]) != 3:
        continue
    pruned_faces[tuple(sorted(f[0]))] = True


def getstart():
    counter = 1
    face = list()
    for t in st.get_skeleton(3):
        if len(t[0]) != 4:
            continue
        if counter == 1:
            face = sorted(t[0])
            counter = counter + 1
        keyed = tuple(sorted(t[0]))
        status[keyed] = False
    return face


def getbdry(t: list):
    boundary = list()
    for i in range(len(t)):
        boundary.append(t[:i]+t[i+1:])
    return boundary


def getnbrs(t):
    nbrs = list()
    bd = getbdry(t)
    for face in bd:
        adj = st.get_cofaces(face, 1)
        for n in adj:
            if t != n[0]:
                nbrs.append(n[0])
    return nbrs


def getFace(t, n):
    t1 = set(t)
    t2 = set(n)
    face = t1.intersection(t2)
    return list(face)


f_to_print = set()
start_v = getstart()
# status[tuple(sorted(start_v))] = True
# Q = queue.Queue()
# Q.put(start_v)
stack = list()
stack.append(start_v)
while len(stack) != 0:
    v = stack.pop()
    key = tuple(sorted(v))
    if not status[key]:
        status[key] = True
        for w in getnbrs(v):
            f = tuple(sorted(getFace(v, w)))
            if f in pruned_faces:
                f_to_print.add(f)
                continue
            stack.append(w)
# mesh = trimesh.load(data + '/' + data + '_pruned_jacobi.off')
# print(len(f_to_print))
# new_mesh = trimesh.Trimesh(mesh.vertices, list(f_to_print))
# new_mesh.export(data + '/' + data + '_walked.off')
print(len(f_to_print))
np.savetxt(data + '/' + data + '_walked.txt', list(f_to_print), fmt='%d')
