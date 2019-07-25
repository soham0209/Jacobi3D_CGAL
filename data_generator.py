import math
import numpy as np
import os


def gaussian(X, var=1):
    v = np.sum(np.multiply(X, X)/var,axis=1)
    v = np.exp(-v)
    val= np.sum(v)
    return val

def createfin(points,centre,r,r_sqr,dirctn=1):
    normal = np.cross(np.array(points[0])-np.array(centre),np.array(points[1])-np.array(centre))
    #normal = np.array([1, 1, 1])
    d = - normal.dot(np.array(centre))
    newpts = list()
    (xstart, xstop) = (int(centre[0]),int(centre[0] + r + 3))
    (ystart, ystop) = (int(centre[1]), int(centre[1] + r + 3))
    (zstart, zstop) = (int(centre[2]), int(centre[2] + r + 3))
    for z in range(zstart,zstop):
        for y in range(ystart,ystop):
            for x in range(xstart,xstop):
                if intersectsplane(x, y, z,normal,d,centre,r_sqr):
                    newpts.append((x, y, z))
    # normal = np.cross(normal,np.array(points[0])-np.array(points[1]))
    # last_pt = newpts[-1]
    # d = -normal.dot(newpts[-1])
    # (xstart,xstop) = (int(last_pt[0]-6),int(last_pt[0]+11))
    # (ystart, ystop) = (int(last_pt[1]-5), int(last_pt[1]+11))
    # (zstart, zstop) = (int(last_pt[2]), int(last_pt[2]+11))
    # for z in range(zstart,zstop):
    #     for y in range(ystart,ystop):
    #         for x in range(xstart,xstop):
    #             if intersectsplane(x, y, z,normal,d,centre,r_sqr):
    #                 newpts.append((x,y,z))
    return  newpts


def checknear(max_point:list, v):
    for p in max_point:
        if (p[0]-v[0])**2 + (p[1]-v[1])**2 + (p[2]-v[2])**2 < 4:
            return True
    return False


def writedata3D(nx, ny, nz, pts, vertexfile, databinfile):
    values = list()

    # artificial_max = max_pos
    # pts = np.array(artificial_max + list(pts))
    for z in range(nz):
        for y in range(ny):
            for x in range(nx):
                # if (x, y, z) in max_pos:
                #      val = 300
                #
                v = pts - np.array([x, y, z])
                val = gaussian(v, 0.2) * 10.0
                # filebin.write(val)
                values.append(val)
                s = '%.4f' % val
                databinfile.write(s + '\n')
                vertexfile.write(str(x) + ' ' + str(y) + ' ' + str(z) + '\n')
    return values

def createSphere(incl_slice,azimuth_slice,r,centre):
    azm = (2 * math.pi) / azimuth_slice
    incl = math.pi / incl_slice
    mean = set()
    for k in range(0, incl_slice + 1):
        for p in range(0, azimuth_slice):
            x = centre[0] + r * np.sin(k * incl) * np.cos(p * azm)
            y = centre[1] + r * np.sin(k * incl) * np.sin(p * azm)
            z = centre[2] + r * np.cos(k * incl)
            mean.add((int(x), int(y), int(z)))
    return mean


def intersectSphere(nx,ny,nz,centre,r_sqr):
    mean = list()
    for z in range(int(centre[2]) + 1,nz+1):
        for y in range(1,ny+1):
            for x in range(1,nx+1):
                if doesIntersect(x,y,z,r_sqr,centre):
                    mean.append((x,y,z))

    return mean

def intersectsplane(x,y,z,norm,d,centre,r_sqr):
    corners = list()
    is_outside = 0
    outside_sphere = 0
    corners.append((x, y, z))
    corners.append((x - 1, y, z))
    corners.append((x, y - 1, z))
    corners.append((x - 1, y - 1, z))
    corners.append((x, y, z - 1))
    corners.append((x - 1, y, z - 1))
    corners.append((x, y - 1, z - 1))
    corners.append((x - 1, y - 1, z - 1))
    for coor in corners:
        dist = norm.dot(np.array(coor))+ d
        dist1 = (coor[0]-centre[0])**2 + (coor[1]-centre[1])**2 + (coor[2]-centre[2])**2 - r_sqr
        if dist > 0:
            is_outside = is_outside + 1
        else:
            is_outside = is_outside - 1
        if dist1 > 0:
            outside_sphere = outside_sphere + 1
        else:
            outside_sphere = outside_sphere - 1


    if is_outside == 8 or is_outside == -8:
        return False
    elif outside_sphere > 0:
        return True
    else:
        return False

def doesIntersect(x, y, z,r_sqr,centre):
    corners = list()
    is_outside = 0
    corners.append((x, y, z))
    corners.append((x - 1, y, z))
    corners.append((x, y - 1, z))
    corners.append((x - 1, y - 1, z))
    corners.append((x, y, z - 1))
    corners.append((x - 1, y, z - 1))
    corners.append((x, y - 1, z - 1))
    corners.append((x - 1, y - 1, z-1))
    for coor in corners:
        dist = (coor[0]-centre[0])**2 + (coor[1]-centre[1])**2 + (coor[2]-centre[2])**2 - r_sqr
        if dist > 0:
            is_outside = is_outside + 1
        else:
            is_outside = is_outside - 1
    if is_outside == 8 or is_outside ==-8:
        return False
    else:
        return True


def createcircle(slice , r):
    step = (2 * math.pi) / slice
    mean = set()
    for theta in range(0,slice):
        x = centre[0] + r * np.cos(theta * step)
        y = centre[1] + r * np.sin(theta * step)
        mean.add((int(x), int(y)))
    return mean

def writedata2D(nx,ny,pts,vertexfile,databinfile):
    values = list()
    # artificial_max = [(45, 43), (50, 32), (48, 18), (37, 12)]
    # pts = pts + artificial_max
    for y in range(ny):
        for x in range(nx):
            # if (x, y) in artificial_max:
            #     val = -2
            # else:
            v = pts - np.array((x, y))
            val = gaussian(v, 4) * 10.0
            # filebin.write(-val)
            values.append(val)
            s = '%.4f' % val
            databinfile.write(s + '\n')
            vertexfile.write(str(x) + ' ' + str(y)  + ' ' + str(0)+'\n')
    return values

def applyRot(pts,X, Y, Z,centre):
    rotated_pts = list()
    for p in pts:
        T_minus = np.matrix([[1,0,0,-centre[0]],[0,1,0,-centre[1]],[0,0,1,-centre[2]],[0,0,0,1]],dtype=np.float32)
        Rx = np.matrix([[1,0,0,0],[0,math.cos(X),-math.sin(X),0],[0,math.sin(X),math.cos(X),0],[0,0,0,1]],dtype=np.float32)
        Ry = np.matrix(
            [[math.cos(Y), 0, -math.sin(Y), 0], [0, 1, 0, 0], [math.sin(Y), 0,  math.cos(Y), 0], [0, 0, 0, 1]],
            dtype=np.float32)
        Rz = np.matrix(
            [[math.cos(Z),-math.sin(Z), 0, 0], [math.sin(Z), math.cos(Z), 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]],
            dtype=np.float32)
        T = np.matrix([[1,0,0,centre[0]],[0,1,0,centre[1]],[0,0,1,centre[2]],[0,0,0,1]],dtype=np.float32)
        v = np.matrix([p[0],p[1],p[2],1],dtype=np.float32)
        v = v.transpose()
        f = np.matmul(T_minus,v)
        s_x = np.matmul(Rx,f)
        s_y = np.matmul(Ry,s_x)
        s_z = np.matmul(Rz, s_y)
        newpt = np.matmul(T, s_z)
        newpt = np.array(newpt.transpose()).reshape(4)
        rotated_pts.append((newpt[0],newpt[1],newpt[2]))
    return rotated_pts

def createMax(CornerPoint:list):
    corners = list()
    for x,y,z in CornerPoint:
        corners.append((x, y, z))
        corners.append((x - 1, y, z))
        corners.append((x, y - 1, z))
        corners.append((x - 1, y - 1, z))
        corners.append((x, y, z - 1))
        corners.append((x - 1, y, z - 1))
        corners.append((x, y - 1, z - 1))
        corners.append((x - 1, y - 1, z - 1))
    return corners


def createHair(cen:list, r, deg):
    hairpt = list()
    x1 = cen[0] + 9*math.cos(math.radians(deg))
    y1 = cen[1] + 9*math.sin(math.radians(deg)) - 1
    for x in range(30,51):
        for y in range(5, 51):
            if (x - cen[0]) ** 2 + (y - cen[1]) ** 2 >= r ** 2:
                if intersectLine(cen, x1, y1, x, y):
                    hairpt.append((x, y))
    return hairpt


def intersectLine(cen, x1, y1, x, y):
    corners = list()
    is_left = 0
    corners.append((x, y))
    corners.append((x + 1, y))
    corners.append((x + 1, y + 1))
    corners.append((x, y + 1))
    for co_or in corners:
        dist = ((x1 - cen[0]) * (co_or[1] - cen[1])) - ((y1 - cen[1])* (co_or[0]-cen[0]))
        if dist < 0:
            is_left = is_left + 1
        else:
            is_left = is_left - 1
    if is_left == 4 or is_left == -4:
        return False
    else:
        return True



cwd = os.getcwd()
a = 30
nx = a
ny = a
nz = a
dim=[nx, ny, nz]
centre = [a/2, a/2, a/2]
dim2 = [nx, ny]
centre2 = [nx/2, ny/2]

#dataName='bone'
dataName = 'halfsphere_' + str(nx) + '_' + str(ny) + '_' + str(nz)
# dataName = 'seeded_circle_'+str(nx)+'_'+str(ny)+'_only_fin'
dirname = cwd + '/' + dataName



if not os.path.exists(dirname):
    os.mkdir(dirname)

binfilename = dirname + '/' + dataName + '.raw'
datafilename = dirname + '/' + dataName +'.txt'
vertfilename = dirname + '/' + dataName +'_vert.txt'

radius_sqr = 16
r = math.ceil(math.sqrt(radius_sqr))
eps = 1e-06
#
pts = intersectSphere(nx,ny,nz,centre,radius_sqr)
finpts = [(18, 15, 15), (17, 17, 16)]
# max_pos = [(20, 31, 42), (28, 48, 40), (15, 12, 41), (18, 10, 23)]
# pts = max_pos
#pts = list()
#fin = createfin(finpts,centre,r,radius_sqr,dirctn=1)
#pts = pts + fin # + applyRot(fin,math.radians(100),math.radians(5),math.radians(5),centre)

# createfin(pts,centre,r,dirctn=-1)
#pts = np.array(list(createSphere(16,16,r,centre)))
# pts = np.array(list(createcircle(16,r)) + createHair(centre2, r, 30) + createHair(centre2, r, -45))
# pts = np.array(createHair(centre2, r, 30) + createHair(centre2, r, -45))

datafile = open(datafilename,'w')
vertfile = open(vertfilename,'w')
datafile.write(str(len(dim))+'\n')
for n in range(len(dim)):
    datafile.write(str(dim[n])+'\n')
values = writedata3D(nx,ny,nz,np.array(pts),vertfile,datafile)
# for n in range(len(dim2)):
  #  datafile.write(str(dim2[n])+'\n')
# values = writedata2D(nx,ny,pts,vertfile,datafile)


# values = writedata2D(nx,ny,pts,vertfile,datafile)
vertfile.close()
datafile.close()

filebin=open(binfilename,'wb')
dt= np.dtype('d')
np.array(values,dt).tofile(filebin,format="%.4f")
filebin.close()
