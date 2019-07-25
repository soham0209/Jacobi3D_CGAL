import vtk
import numpy as np
import sys
import os


def writeVtk(Vert, JacobiFaces, filename):
    count = 0
    vpoints = vtk.vtkPoints()
    for ii in range(len(Vert)):
        vpoints.InsertNextPoint(Vert[ii][0], Vert[ii][1], Vert[ii][2])
    vpfaces = vtk.vtkCellArray()
    for e in JacobiFaces:
        count = count + 1
        t = e
        polygon = vtk.vtkPolygon()
        polygon.GetPointIds().SetNumberOfIds(3)
        polygon.GetPointIds().SetId(0, int(t[0]))
        polygon.GetPointIds().SetId(1, int(t[1]))
        polygon.GetPointIds().SetId(2, int(t[2]))
        vpfaces.InsertNextCell(polygon)
    writeTofile(vpoints, vpfaces, filename, count= len(JacobiFaces))


def writeTofile(Points, Faces, filename, count):
    polydata = vtk.vtkPolyData()
    polydata.SetPoints(Points)
    polydata.SetPolys(Faces)
    # polydata.GetPointData().SetScalars(value)
    gw = vtk.vtkXMLPolyDataWriter()
    # dirname = cwd + '/' + dataName + dirName
    # if not os.path.exists(dirname):
    #  os.mkdir(dirname)
    # print('Writing Unstable Manifolds')
    gw.SetFileName(filename + '.vtp')
    gw.SetInputData(polydata)
    gw.Write()
    print('Wrote ' , count, 'faces to ' + filename + '.vtp')


def writeOFF(Vert, JacobiFaces, filename):
    print('Writing to ' + filename+'.off'+'...')
    f = open(filename+'.off', 'w')
    f.write('OFF\n')
    f.write(str(len(Vert))+' '+str(len(JacobiFaces))+' '+'0\n')
    for ii in range(len(Vert)):
        line = str(int(Vert[ii][0]))+' '+str(int(Vert[ii][1]))+' '+str(int(Vert[ii][2]))+'\n'
        f.write(line)
    for ii in range(len(JacobiFaces)):
        line = str(len(JacobiFaces[ii]))+' '+ str(int(JacobiFaces[ii][0])) + ' '+ str(int(JacobiFaces[ii][1]))+' '+str(int(JacobiFaces[ii][2]))+'\n'
        f.write(line)
    f.close()

    
if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Usage: python <dataName>')
    datafile = sys.argv[1]
    delta = ''
    if len(sys.argv) > 2:
        delta = str(float(sys.argv[2]))
    # Simplex_tree, Vertices, num_face = readOFF(datafile+'/'+datafile+'_triangulated.off')
    # readperseus(datafile+'/'+datafile+'.txt')
    # j_faces = computejacobi(Simplex_tree)
    # print(len(j_faces))
    # writeVtk(j_faces, datafile+'/'+datafile+'_jacobi')
    resampled = False
    vertfile = datafile+'/'+datafile + '_vert.txt'
    if resampled:
        vertfile = datafile + '/' + datafile + '_resampled_vert.txt'
    facefile = datafile+'/'+datafile + '_jacobi.txt'
    pruned_facefile = datafile+'/'+datafile+'_pruned_jacobi.txt'
    vert = np.loadtxt(vertfile)
    p_j_facefile = datafile+'/'+datafile+'_pruned_eroded_jacobi.txt'
    p_w_facefile = datafile+'/'+datafile+'_pruned_walked.txt'
    final_facefile = datafile+'/'+datafile+'_final.txt'
    if os.path.exists(facefile):
        j_faces = np.loadtxt(facefile)
        writeVtk(vert, j_faces, datafile+'/'+datafile+'_jacobi')
        writeOFF(vert, j_faces, datafile+'/'+datafile+'_jacobi')
    if os.path.exists(pruned_facefile):
        jp_faces = np.loadtxt(pruned_facefile)
        writeOFF(vert, jp_faces, datafile+'/'+datafile+'_pruned_jacobi')
        writeVtk(vert, jp_faces, datafile+'/'+datafile+'_pruned_jacobi')
    if os.path.exists(p_j_facefile):
        p_faces = np.loadtxt(p_j_facefile)
        writeOFF(vert, p_faces, datafile+'/'+datafile+'_pruned_eroded_jacobi_' + delta)
        writeVtk(vert, p_faces, datafile+'/'+datafile+'_pruned_eroded_jacobi_' + delta)
    if os.path.exists(p_w_facefile):
        j_faces = np.loadtxt(p_w_facefile)
        writeVtk(vert, j_faces, datafile+'/'+datafile+'_pruned_walked')
        writeOFF(vert, j_faces, datafile+'/'+datafile+'_pruned_walked')
    if os.path.exists(final_facefile):
        final_faces = np.loadtxt(final_facefile)
        writeVtk(vert, final_faces, datafile + '/' + datafile + '_final')
        writeOFF(vert, final_faces, datafile + '/' + datafile + '_final')
