#!/afs/ipp/.cs/python_modules/amd64_sles11/python27/python/2.7.10/bin/python
import Illustris_python as il
import matplotlib.pyplot as plt
import numpy as np

nsnap = np.array([68,85,103,120])

for nn in range(0,nsnap.size) :
    path1 = '/virgo/simulations/Illustris/Illustris-1/'    


    file1 = '../Data/Illustris_data/Illustris1/Relaxed/m200_halos_hydr_rel_'+(str(nsnap[nn]))+'.txt'
    id1 = np.loadtxt(file1, dtype = 'int',unpack=True, usecols=[0])
    mvir1 = np.loadtxt(file1, unpack=True, usecols=[1])
    x1 = np.loadtxt(file1, unpack=True, usecols=[4])
    y1 = np.loadtxt(file1, unpack=True, usecols=[5])
    z1 = np.loadtxt(file1, unpack=True, usecols=[6])
    r1 = np.loadtxt(file1, unpack=True, usecols=[13])

    hbox = 37500.
    box = 75000.
    dim = mvir1.size

    
    for i in range(0,dim) :
        print i,nsnap[nn]
        bin = np.log10(2*r1[i])/100
        
        rhodm = []
        rhost = []
        rhoga = []
        xrho = []
        cc = 0
        halo = il.snapshot.loadHalo(path1,nsnap[nn],id1[i],'dm')
        x = halo['Coordinates'][:,0]
        y = halo['Coordinates'][:,1]
        z = halo['Coordinates'][:,2]
        if ((x.max()-x.min())> hbox) :
            if(x1[i] <= hbox):
                x1[i] = x1[i] + box
            for j in range(0,x.size) :
                if(x[j] <= hbox):
                    x[j] = x[j] +box
        if ((y.max()-y.min())> hbox) :
            if(y1[i] <= hbox):
                y1[i] = y1[i] + box
            for j in range(0,y.size) :
                if(y[j] <= hbox):
                    y[j] = y[j] +box
        if ((z.max()-z.min())> hbox) :
            if(z1[i] <= hbox):
                z1[i] = z1[i] + box
            for j in range(0,z.size) :
                if(z[j] <= hbox):
                    z[j] = z[j] +box
        dx = x- x1[i]
        dy = y- y1[i]
        dz = z- z1[i]
        dist = np.sqrt(dx**2 + dy**2 + dz**2)
        dist = np.sort(dist)
        for j in range(0,100) :
            count = 0
            rad = 10**(bin*(j+1))
            volume = (4./3.)*np.pi*rad**3
            if (cc < dist.size) :
                for k in range(cc,dist.size):
                    if (dist[k]< rad):
                        count = count+1
                rho = (cc+count)*6.3*1.e+6/volume
                rhodm.append(rho)
                xrho.append(rad)
                cc = cc + count
            if (cc == dist.size and j<99) :
                xrho.append(rad)
                rhodm.append(rho)
        
        rhodm = np.array(rhodm)
        xrho = np.array(xrho)
        
        print 'dm done'
        
        
        cc=0
        buf = 0
        halo = il.snapshot.loadHalo(path1,nsnap[nn],id1[i],'gas')
        x = halo['Coordinates'][:,0]
        y = halo['Coordinates'][:,1]
        z = halo['Coordinates'][:,2]
        ms = halo['Masses']*1.e+10
        if ((x.max()-x.min())> hbox) :
            for j in range(0,x.size) :
                if(x[j] <= hbox) :
                    x[j] = x[j] + box
        if ((y.max()-y.min())> hbox) :
            for j in range(0,y.size) :
                if(y[j] <= hbox):
                    y[j] = y[j] + box
        if ((z.max()-z.min())> hbox) :
            for j in range(0,z.size) :
                if(z[j] <= hbox):
                    z[j] = z[j] + box

        dx = x- x1[i]
        dy = y- y1[i]
        dz = z- z1[i]
        dist = np.sqrt(dx**2 + dy**2 + dz**2)
        indmass = dist.argsort()
        dist = np.sort(dist)
        for j in range(0,100) :
            count = 0
            mass = 0.
            rad = 10**(bin*(j+1))
            volume = (4./3.)*np.pi*rad**3
            if (cc < dist.size) :
                for k in range(cc,dist.size):
                    if (dist[k]< rad):
                        count = count+1
                        mass = mass + ms[indmass[k]]                        
                rho = (buf+mass)/volume
                rhoga.append(rho)
                cc = cc + count
                buf = buf +mass
            if (cc == dist.size and j<99) :
                rhoga.append(rho)
        
        rhoga = np.array(rhoga)
        
        
        print 'gas done'
        
        
        
        cc=0
        buf = 0
        halo = il.snapshot.loadHalo(path1,nsnap[nn],id1[i],'stars')
        print 'reading ok'
        x = halo['Coordinates'][:,0]
        y = halo['Coordinates'][:,1]
        z = halo['Coordinates'][:,2]
        ms = halo['Masses']*1.e+10
        print x
        print ms
        if ((x.max()-x.min())> hbox) :
            for j in range(0,x.size) :
                if(x[j] <= hbox) :
                    x[j] = x[j] + box
        if ((y.max()-y.min())> hbox) :
            for j in range(0,y.size) :
                if(y[j] <= hbox):
                    y[j] = y[j] + box
        if ((z.max()-z.min())> hbox) :
            for j in range(0,z.size) :
                if(z[j] <= hbox):
                    z[j] = z[j] + box
        dx = x- x1[i]
        dy = y- y1[i]
        dz = z- z1[i]
        dist = np.sqrt(dx**2 + dy**2 + dz**2)
        indmass = dist.argsort()
        dist = np.sort(dist)
        for j in range(0,100) :
                #print 'j=',j
            count = 0
            mass = 0.
            rad = 10**(bin*(j+1))
            volume = (4./3.)*np.pi*rad**3
            if (cc < dist.size) :
                for k in range(cc,dist.size):
                    #print 'cc=',cc
                    if (dist[k]< rad):
                        count = count+1
                        mass = mass + ms[indmass[k]]
                rho = (buf+mass)/volume
                rhost.append(rho)
                cc = cc + count
                buf = buf + mass
            if (cc == dist.size and j<99) :
                rhost.append(rho)

        rhost = np.array(rhost)
            
        print 'stars done'
               #print rhohalo
        filename='profile_hydr_rel_'+str(id1[i])+'_'+str(nsnap[nn])+'.txt'
        np.savetxt(filename,np.c_[xrho,rhodm,rhoga,rhost],fmt='%f   %1.4e   %1.4e   %1.4e')
        

