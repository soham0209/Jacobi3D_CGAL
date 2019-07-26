# -*- coding: utf-8 -*-
"""
Created on Tue Jul 23 18:42:15 2019

@author: mustafa
"""
import random
import numpy as np
import sys
from sklearn.neural_network import MLPRegressor
from scipy.interpolate import RegularGridInterpolator


def resample(f,sampled_points,number_of_desired_points):
    
    #number_of_desired_points=max(number_of_desired_points, len(sampled_points))
    
    index = int(random.random() * number_of_desired_points)
    beta = 0.0
    
    P=[p for p in  sampled_points]
    
    W=[f[p] for p in sampled_points]
    #print(W)
    W=W/np.sum(W)
    
   
    
    mw = max(W)
    
    newpoints=[]
    for i in range(number_of_desired_points):
        
        beta += random.random() * 2.0 * mw
        while beta > W[index]:
            beta -= W[index]
            index = (index + 1) % number_of_desired_points
            
        newpoints.append(P[index])
    #print([f[p ] for p in  newpoints] )    
    return newpoints



def resample_cont(f, sampled_points, number_of_desired_points, simulated_annealing_fac=0.15):
    
    
    index = int(random.random() * number_of_desired_points)
    beta = 0.0
    
    P=[p for p in  sampled_points]
    
    W=[f(p[0],p[1],p[2] ) for p in sampled_points]
    #print(W)
    #print(P)
    W=W/np.sum(W)
    
   
    
    mw = max(W)
    
    newpoints=[]
    for i in range(number_of_desired_points):
        
        beta += random.random() * 2.0 * mw
        while beta > W[index]:
            beta -= W[index]
            index = (index + 1) % number_of_desired_points
            
        newpoints.append(P[index])
        
    if simulated_annealing_fac==0:
        return newpoints
    else:
            
        alpha=min( int(number_of_desired_points*simulated_annealing_fac)+1,number_of_desired_points)
        
        beta= min(int(number_of_desired_points*(1-simulated_annealing_fac))+1,number_of_desired_points )  
        
        
        
        original_sampled=random.sample(sampled_points,alpha )
        
        sample_new_points = random.sample(newpoints, beta)
        
        out = random.sample(original_sampled+sample_new_points,number_of_desired_points)
    
        return out


def generate_grid(lower_corner,higher_corner):
    grid=[]
    for k in range(lower_corner[0],higher_corner[0]):
        for j in range(lower_corner[1],higher_corner[1]):
            for i in range(lower_corner[2],higher_corner[2]):
                grid.append((i, j, k))
    return grid            


def random_point_around_p(point,R,r):
    
    x=random.uniform(0.1, 2)
    y=random.uniform(0.1, 2)
    z=random.uniform(0.1, 2)
    P=np.array([x,y,z])
    
    scaler=random.uniform(R, r)
    
    P=P/np.linalg.norm(P)
    
    
    return P*scaler+np.array(point)

def generate_random_cloud_around_p(point,R,r,pointcld_size=30):
    out=[]
    for _ in range(pointcld_size):
        p=random_point_around_p(point,R,r)
        out.append(p)
    return out    
    

def sample_density(lower_corner,higher_corner,f, f_threshold=0.01,number_of_sampless=10,number_of_iterations=100):
    # there are 26 neighbors for (i,j,k)
    neighbors=[(1,0,0),(0,1,0),(0,0,1),
               
               (-1,0,0),(0,-1,0),(0,0,-1),
               
               (-1,-1,0),(0,-1,-1),(-1,0,-1),
               
               (1,1,0),(1,1,0),(1,0,1),
               
               (-1,1,0),(-1,1,0),(-1,0,1),
               
               (1,-1,0),(1,-1,0),(1,0,-1),
               
               (-1,1,1),(1,-1,1),(1,1,-1),
               
               (-1,1,-1),(-1,-1,1),(1,-1,-1),
               
               (-1,-1,-1),(1,1,1)
               ]
    
    # the grid must be the domain of f--that is  same as f.keys()
    
    #grid=generate_grid(lower_corner,higher_corner)
    
    grid=f.keys()
    #assert( set(grid)==set(f.keys()))
    samples=random.sample(grid,number_of_sampless)
    
    
    
    for _ in range(number_of_iterations):
        newsamples=[]
        for sample in samples:
            
            for n in neighbors:
                
                current=np.array(n)+np.array(sample)

                if current[0]>=lower_corner[0] and current[0]<higher_corner[0]:
                    if current[1]>=lower_corner[1] and current[1]<higher_corner[1]:
                        if current[2]>=lower_corner[2] and current[2]<higher_corner[2]:
                            if f[ tuple( current) ]>=f_threshold: # only consider these points larger than a certain value
                            
                            
                                newsamples.append( tuple(current ))
                            
        sumsamples=set(samples).union(set(newsamples))  
        
        samples=resample(f,sumsamples,number_of_sampless)
                 
    return samples
    


def sample_density_cont(lower_corner,higher_corner,f,R=1.5,r=0.005,local_point_cld_size=8, f_threshold=0.1,number_of_sampless=100,number_of_iterations=20,simulated_annealing_fac=0.15):
    
    # the grid must be the domain of f--that is  same as f.keys()
    
    grid=generate_grid(lower_corner,higher_corner)
    

    samples=random.sample(grid,min(number_of_sampless,len(grid)-1))
    
    
    newsamples=[] 
    for _ in range(number_of_iterations):
        
        for sample in samples:
            
            neighbors=generate_random_cloud_around_p(sample,R,r,local_point_cld_size)
            
            
            for n in neighbors:
                
                current=np.array(n)

                if current[0]>=lower_corner[0] and current[0]<higher_corner[0]:
                    if current[1]>=lower_corner[1] and current[1]<higher_corner[1]:
                        if current[2]>=lower_corner[2] and current[2]<higher_corner[2]:
                            
                            if f( current[0],current[1],current[2]  )>=f_threshold: # only consider these points larger than a certain value
                            
                            
                                newsamples.append( tuple(current ))
                            
        sumsamples=set(samples).union(set(newsamples))  

        #print(number_of_sampless)
        samples=resample_cont(f, sumsamples, number_of_sampless, simulated_annealing_fac)
                 
    return samples


def generate_scalar_function_from_dictionary(f_dic):
    from scipy.interpolate import Rbf
    
    X=np.array([ list(x) for x in f_dic.keys()])
    
    y=np.array([ f_dic[x] for x in  f_dic ])
    
    rbfi = Rbf(X[:,0], X[:,1], X[:,2], y)

    
    return rbfi


def read_pers_file(filename, f_val: dict, grid):
    pers_file = open(filename, 'r')
    total_dim = int(pers_file.readline().strip('\n'))
    total_ver = 1
    for i in range(total_dim):
        total_ver = total_ver * int(pers_file.readline().strip('\n'))
    for i in range(total_ver):
        val = float(pers_file.readline().strip('\n'))
        f_val[grid[i]] = val
    pers_file.close()
    print('Read ', total_ver, ' values')


def write_resampled_file(f_dict, f,  filename, num_v, newpts):
    fout = open(resam_pers_file, 'w')
    fout.write('3\n')
    fout.write(str(num_v) + '\n')
    fout.write('1\n')
    fout.write('1\n')
    count = 0
    for x in f_dict:
        line = str(f_dict[x]) + '\n'
        fout.write(line)
        count = count + 1
    for p in newpts:
        val = f(p[0] ,p[1], p[2])
        line = str(val) + '\n'
        fout.write(line)
        count = count + 1
    assert count == num_v
    print('Wrote ', num_v, ' values')
    fout.close()


if __name__ == "__main__":
    
    f_dic = dict()
    if len(sys.argv) < 3:
        print('Usage <dataName> <gridDimension>')
    data = sys.argv[1]
    dim = int(sys.argv[2])
    grid=generate_grid(lower_corner=(0,0,0),higher_corner=(dim,dim,dim))
    pers_file_path = data + "/" + data + ".txt"
    vert_file_path = data + "/" + data + "_resampled_vert.txt"
    resam_pers_file = data + "/" + data + "_resampled.txt"
    read_pers_file(pers_file_path, f_dic, grid)
    f = generate_scalar_function_from_dictionary(f_dic)
    out = sample_density_cont((0,0,0),(20,20,20), f, simulated_annealing_fac = 0.3)
    finalpts = grid + list(out)
    np.savetxt(vert_file_path, finalpts, fmt = '%4f')
    num_ver = len(finalpts)
    print('Wrote ', num_ver, ' vertices')
    write_resampled_file(f_dic, f, resam_pers_file, num_ver, out)

