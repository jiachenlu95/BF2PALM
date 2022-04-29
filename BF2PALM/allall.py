import numpy as np
from scipy.io import savemat
import netCDF4 as nc
import scipy.interpolate as spi

### Calculate important vertical profiles from 3D field
def c_prof(u): # Calculate profile using intrinsic average 
    tmp = np.zeros(np.shape(u)[2])
    for i in range(np.shape(u)[2]):
        tmp[i] = np.nanmean(u[:,:,i])
    return(tmp)
def c_profA(u):
    tmp = np.zeros(np.shape(u)[2])
    for i in range(np.shape(u)[2]):
        tmp[i] = np.nanmean(abs(u[:,:,i]))
    return(tmp)


def c_std(u): # Calculate profile using intrinsic average 
    tmp = np.zeros(np.shape(u)[2])
    for i in range(np.shape(u)[2]):
        tmp[i] = np.nanstd(u[:,:,i])
    return(tmp)

def c_disp(u,w): # Calculate dispersive flux
    tmp = np.zeros(np.shape(u))
    for i in range(np.shape(u)[2]):
        tmp[:,:,i] = (u[:,:,i] - c_prof(u)[i])*(w[:,:,i] - c_prof(w)[i])
    return(tmp)

def c_disp1(u): # Dispersive parameter for a single variable.
    tmp = np.zeros(np.shape(u))
    for i in range(np.shape(u)[2]):
        tmp[:,:,i] = u[:,:,i] - c_prof(u)[i]
    return(tmp)


def c_tf(uw,u,w): # Calculate turbulent flux
    return(uw-u*w)
    
    
def c_dpdx(p_ens,dx): # Calculate dp/dx
    disp_ens = np.zeros(np.shape(p_ens))
    tmp = np.zeros(np.shape(p_ens))
    for i in range(np.shape(p_ens)[2]):
        disp_ens[:,:,i] = p_ens[:,:,i] - np.nanmean(p_ens[:,:,i])
    for i in range(1,np.shape(p_ens)[0]):
        for j in range(np.shape(p_ens)[1]):
            for k in range(np.shape(p_ens)[2]):
                tmp[i,j,k] = (disp_ens[i,j,k] - disp_ens[i-1,j,k])/dx
    return(tmp)

def c_dpdy(p_ens,dy): # Calculate dp/dy
    disp_ens = np.zeros(np.shape(p_ens))
    tmp = np.zeros(np.shape(p_ens))
    for i in range(np.shape(p_ens)[2]):
        disp_ens[:,:,i] = p_ens[:,:,i] - np.nanmean(p_ens[:,:,i])
    for i in range(np.shape(p_ens)[0]):
        for j in range(1,np.shape(p_ens)[1]):
            for k in range(np.shape(p_ens)[2]):
                tmp[i,j,k] = (disp_ens[i,j,k] - disp_ens[i,j-1,k])/dy
    return(tmp)


def c_dpdz(p_ens,z): # Calculate dp/dz
    #dz = np.insert(np.diff(z*16, n=1, axis=0),0,0.5)
    dz = 0.5
    disp_ens = np.zeros(np.shape(p_ens))
    tmp = np.zeros(np.shape(p_ens))
    for i in range(np.shape(p_ens)[2]):
        disp_ens[:,:,i] = p_ens[:,:,i] - np.nanmean(p_ens[:,:,i])
    for i in range(np.shape(p_ens)[0]): # x-dir
        for j in range(np.shape(p_ens)[1]):  # y-dir
            for k in range(1,np.shape(p_ens)[2]):  # z-dir
                tmp[i,j,k] = (disp_ens[i,j,k] - disp_ens[i,j,k-1])/dz
    return(tmp)



def c_tri(uvw_ens,u_ens,v_ens,w_ens,vw_ens,uw_ens,uv_ens): # Triple correlation
    tmp = np.zeros(np.shape(uvw_ens))
    tmp = uvw_ens - u_ens*v_ens*w_ens
    - c_tf(vw_ens,v_ens,w_ens)*u_ens
    - c_tf(uw_ens,u_ens,w_ens)*v_ens
    - c_tf(uv_ens,u_ens,v_ens)*w_ens
    return(tmp)


def c_dz(prof,z): # treat dz as constant, not valid for layers near the top.
    dz = np.insert(np.diff(z*16, n=1, axis=0),0,0.5)
    tmp = np.zeros(np.size(prof))
    for i in range(1,np.size(prof)):
        tmp[i] = (prof[i]-prof[i-1])/dz[i]
    return(tmp)


def c_tkeb():
    

    # shear production
    z = zw_3d

    phi = np.ones(z.size) # probability to encounter a fluid particle within the constant elevation z plane.
    phi[:32] = 1 - lp # Below the canopy we have probability of 1-lp
    dphi = np.gradient(phi,z.ravel()*16)



    shearP = np.zeros([3,z.shape[0]])
    shearP[0,:] = -c_prof(c_tf(uw,u,w))*np.gradient(c_prof(u),z.ravel()*16)
    shearP[1,:] = -c_prof(c_tf(vw,v,w))*np.gradient(c_prof(v),z.ravel()*16)
    shearP[2,:] = -c_prof(c_tf(ww,w,w))*np.gradient(c_prof(w),z.ravel()*16)
    # axis=0 -- x
    # axis=1 -- y
    # axis=2 -- z

    shearPV = np.zeros([3,z.shape[0]])
    shearPV[0,:] = -c_prof(c_tf(uw,u,w))*c_prof(u)/(phi)*dphi
    shearPV[1,:] = -c_prof(c_tf(vw,v,w))*c_prof(v)/(phi)*dphi
    shearPV[2,:] = -c_prof(c_tf(ww,w,w))*c_prof(w)/(phi)*dphi


    # wake production
    wakeP = np.zeros([9,z.shape[0]])
    # --j=x
    wakeP[0,:] = -c_prof(c_disp1(uu - u*u)*np.gradient(c_disp1(u),0.5,axis=0))
    wakeP[1,:] = -c_prof(c_disp1(uv - u*v)*np.gradient(c_disp1(v),0.5,axis=0))
    wakeP[2,:] = -c_prof(c_disp1(uw - u*w)*np.gradient(c_disp1(w),0.5,axis=0))
    # --j=y
    wakeP[3,:] = -c_prof(c_disp1(uv - v*u)*np.gradient(c_disp1(u),0.5,axis=1))
    wakeP[4,:] = -c_prof(c_disp1(vv - v*v)*np.gradient(c_disp1(v),0.5,axis=1))
    wakeP[5,:] = -c_prof(c_disp1(vw - v*w)*np.gradient(c_disp1(w),0.5,axis=1))
    # --j=z
    wakeP[6,:] = -c_prof(c_disp1(uw - w*u)*np.gradient(c_disp1(u),z.ravel()*16,axis=2))
    wakeP[7,:] = -c_prof(c_disp1(vw - w*v)*np.gradient(c_disp1(v),z.ravel()*16,axis=2))
    wakeP[8,:] = -c_prof(c_disp1(ww - w*w)*np.gradient(c_disp1(w),z.ravel()*16,axis=2))
    # form induced shear stress


    shearA = np.zeros([3,z.shape[0]])
    shearA[0,:] = -c_prof(c_tf(uw,u,w))*c_prof(np.gradient(c_disp1(u),z.ravel()*16,axis=2))
    shearA[1,:] = -c_prof(c_tf(vw,v,w))*c_prof(np.gradient(c_disp1(v),z.ravel()*16,axis=2))
    shearA[2,:] = -c_prof(c_tf(ww,w,w))*c_prof(np.gradient(c_disp1(w),z.ravel()*16,axis=2))
    # turbulent transport uuw,vvw,www
    


    turbTV = np.zeros([3,z.shape[0]])
    turbTV[0,:] = -1/2*c_prof(c_tri(uuw,u,u,w,uw,uw,uu))/(phi)*dphi
    turbTV[1,:] = -1/2*c_prof(c_tri(vvw,v,v,w,vw,vw,vv))/(phi)*dphi
    turbTV[2,:] = -1/2*c_prof(c_tri(www,w,w,w,ww,ww,ww))/(phi)*dphi

    turbT = np.zeros([3,z.shape[0]])
    turbT[0,:] = -1/2*np.gradient(c_prof(c_tri(uuw,u,u,w,uw,uw,uu)),z.ravel()*16)
    turbT[1,:] = -1/2*np.gradient(c_prof(c_tri(vvw,v,v,w,vw,vw,vv)),z.ravel()*16)
    turbT[2,:] = -1/2*np.gradient(c_prof(c_tri(www,w,w,w,ww,ww,ww)),z.ravel()*16) 
    


    
    # dispersive transport

    disTV = np.zeros([3,z.shape[0]])

    disTV[0,:] = -1/2*c_prof(c_disp1(w) * c_disp1(uu - u*u))/(phi)*dphi
    disTV[1,:] = -1/2*c_prof(c_disp1(w) * c_disp1(vv - v*v))/(phi)*dphi
    disTV[2,:] = -1/2*c_prof(c_disp1(w) * c_disp1(ww - w*w))/(phi)*dphi

    disT = np.zeros([3,z.shape[0]])
    disT[0,:] = -1/2*np.gradient(c_prof(c_disp1(w) * c_disp1(uu - u*u)),z.ravel()*16)
    disT[1,:] = -1/2*np.gradient(c_prof(c_disp1(w) * c_disp1(vv - v*v)),z.ravel()*16)
    disT[2,:] = -1/2*np.gradient(c_prof(c_disp1(w) * c_disp1(ww - w*w)),z.ravel()*16)



    # pressure transport

    presTV =  -c_prof(c_tf(wp,w,p))/(phi)*dphi
    presT = -np.gradient(c_prof(c_tf(wp,w,p)),z.ravel()*16)

    # Residue for sgs part

    sgsT = 0 - (np.nansum(turbT,axis=0)+np.nansum(disT,axis=0)+np.nansum(disTV,axis=0)+np.nansum(turbTV,axis=0)+presT+presTV)
    sgsD = 0 - (np.nansum(shearP,axis=0)+np.nansum(shearPV,axis=0)+np.nansum(wakeP,axis=0)+np.nansum(shearA,axis=0))
               
    return(shearP,shearA,wakeP,turbT,disT,presT,sgsT,sgsD,presTV,turbTV,shearPV,disTV)

# Everything interpolate to zu_3d !!

def interpZZ(u,zold,znew):
    uNew = np.zeros([u.shape[0],u.shape[1],u.shape[2]])
    for i in range(u.shape[0]):
        for j in range(u.shape[1]):
            f = spi.interp1d(zold,u[i,j,:],axis=0, fill_value="extrapolate")
            #f = interpolate.interp1d()
            uNew[i,j,:] = f(znew) 
    return(uNew)

def readDataNo(var):
    tmp = np.mean(data.variables[var][:].data,axis=0)
    tmp[tmp<=-1000] = np.nan
    tmp = np.moveaxis(tmp,0,2)
    tmp = np.moveaxis(tmp,1,0)[:,:,:]
    return(tmp)

def readDataI(var,zw_3d,zu_3d):
    tmp = np.mean(data.variables[var][:].data,axis=0)
    tmp[tmp<=-1000] = np.nan
    tmp = np.moveaxis(tmp,0,2)
    tmp = np.moveaxis(tmp,1,0)[:,:,:]
    interpZZ(tmp,zw_3d,zu_3d)
    return(tmp)

def pressureDef(topo):    
    def pressureDefPre(topo): # Count how much pb and pf was sampled from topo
        pfNN = np.zeros(topo.shape[1])
        pbNN = np.zeros(topo.shape[1])

        for j in range(topo.shape[1]):
            for i in range(topo.shape[0]):
                try:
                    if np.isnan(topo[i,j]) and ~np.isnan(topo[i-1,j]): # fontface
                        pfNN[j] +=1
                        #print('cao')
                    if np.isnan(topo[i,j]) and ~np.isnan(topo[i+1,j]): # back face
                        pbNN[j] +=1
                except:
                    0  
        return(pfNN,pbNN)

    pfNN,pbNN = pressureDefPre(topo)
    
    pf = []
    pb = []
    pfO = []
    pbO = []
    distO = []
    pfN = 0
    pbN = 0
    dist = []
    o1 = 0; o2 = 0; o3 = 0
    for j in range(topo.shape[1]):
        
        if ~np.isnan(topo[0,j]) and ~np.isnan(topo[-1,j]): # no B grid in the first and the end - normal row
            for i in range(topo.shape[0]):
                try:
                    if np.isnan(topo[i,j]) and ~np.isnan(topo[i-1,j]): # fontface
                        pf.append(topo[i-1,j])
                        pfN += 1
                        itmp = i
                    if np.isnan(topo[i,j]) and ~np.isnan(topo[i+1,j]): # back face
                        pb.append(topo[i+1,j])
                        pbN += 1
                        dist.append(i-itmp+1) # index of the paired pf and pb
                except:
                    0
            
        if ~np.isnan(topo[0,j]) and np.isnan(topo[-1,j]): # no B grid in the first but the end
            count = 0
            #print('Outlier 1 found at' + str(j) + ' th row')
            o1+=1
            for i in range(topo.shape[0]):

                try:
                    if np.isnan(topo[i,j]) and ~np.isnan(topo[i-1,j]): # fontface normal except for the last
                        if  count == pfNN[j]-1: # if that's the last frontal face 
                            
                            pf.append(topo[i-1,j])
                            pb.append(topo[0,j])
                            pfN += 1
                            pbN += 1
                            dist.append(topo.shape[0]-i) # should be the length of the last continued building grids
                            #print(topo.shape[0]-i)
                            
                        else:
                            pf.append(topo[i-1,j])
                            pfN += 1
                            itmp = i
                            count += 1              
                    if np.isnan(topo[i,j]) and ~np.isnan(topo[i+1,j]): # back face normal
                        pb.append(topo[i+1,j])
                        pbN += 1
                        dist.append(i-itmp+1) # index of the paired pf and pb
                        
                except:
                    1
                          
        if np.isnan(topo[0,j]) and ~np.isnan(topo[-1,j]): # no B grid in the end but the first
            
            first = True
            #print('Outlier 2 found at' + str(j) + ' th row')
            o2+=1
            
            for i in range(topo.shape[0]):
                try:
                    if np.isnan(topo[i,j]) and ~np.isnan(topo[i-1,j]): # fontface normal
                        #print('sampling pf'+str(i)+str(j))
                        pf.append(topo[i-1,j])
                        pfN += 1
                        itmp = i
                        
                    if np.isnan(topo[i,j]) and ~np.isnan(topo[i+1,j]): # back face normal except for the first
                        if first: 
                            #print('sampling pb'+str(i)+str(j))
                            pb.append(topo[i+1,j])
                            pbN += 1
                            dist.append(i) # grid count to the end - distance fixed
                            first = False
                        else:
                            pb.append(topo[i+1,j])
                            #print('sampling pb'+str(i)+str(j))
                            pbN += 1
                            dist.append(i-itmp+1) # index of the paired pf and pb     
                except:
                    2

        if np.isnan(topo[0,j]) and np.isnan(topo[-1,j]): # B grid in the end and the first, record them in a seperate array
            count = 0
            first = True
            o3+=1
            #print('Outlier 3 found at' + str(j) + ' th row')
            
            for i in range(topo.shape[0]):
                
                try:
                    if np.isnan(topo[i,j]) and ~np.isnan(topo[i-1,j]): # fontface normal except for the last
                        
                        if  count == pfNN[j]-1: # if that's the last frontal face
                            pfO.append(topo[i-1,j])
                            upSize = topo.shape[0]-i
                            
                        else:
                            pf.append(topo[i-1,j])
                            pfN += 1
                            itmp = i
                            count += 1
                except:
                    3
                            
            for i in range(topo.shape[0]):
                
                try:
                    if np.isnan(topo[i,j]) and ~np.isnan(topo[i+1,j]): # back face normal except for the first
                        
                        if first:
                            #print('qppen')
                            pbO.append(topo[i+1,j])
                            lowSize = i
                            distO.append(lowSize+upSize+1) # index of the paired pf and pb

                            first = False
                            
                        else:
                            pb.append(topo[i+1,j])
                            pbN += 1
                            dist.append(i-itmp+1) # index of the paired pf and pb
                        
                except:
                    3
                    
            
    pf.extend(pfO)
    pb.extend(pbO)
    dist.extend(distO)
                    
    return(np.array(pf),np.array(pb),np.array(dist))


# for idealized cases


names = globals()
utau = 0.21105008

for num in range(2,3):
    dirName = names['dir' + str(num)]
    dataName = names['dataName' + str(num)]
    for i in range(np.size(dataName)):          


        data_dir = dirName+dataName[i]+'/OUTPUT/'
        fileName = dataName[i]+'_av_3d.nc'
        print(fileName)

        dxdy = 0.5 # for icealized cases
        data = nc.Dataset(data_dir+fileName)
        h = 16
        zu_3d = data.variables['zu_3d'][:].data/h
        zw_3d = data.variables['zw_3d'][:].data/h


        # read data
        e = readDataNo('e')
        s = readDataNo('s')
        p = readDataNo('p')
        u = readDataNo('u')
        v = readDataNo('v')
        wp = readDataNo('wpp')
        
        
        w = readDataI('w',zw_3d,zu_3d)

        
        uw = readDataNo('uw')
        vw = readDataNo('vw')
        uv = readDataNo('uv')

        uu = readDataNo('uu')
        vv = readDataNo('vv')
        ww = readDataNo('ww')
        
        wc = readDataNo('wc')

        uuw = readDataNo('uuw')
        vvw = readDataNo('vvw')
        www = readDataNo('www')

        
        Utmp = np.nanmean(c_prof(u)[1:33]**2*abs(c_prof(u)[1:33])/c_prof(u)[1:33])
        pTrue = p+2/3*e

        dp = np.zeros(33)
        for k in range(1,33):
            pf,pb,dist = pressureDef(pTrue[:,:,k]) # front face of building
            #pf,pb,pfN[k],pbN[k],dist = pressureDefS(pTrue[:,:,k]) # front face of building
            
            volP = dxdy*np.array(dist)*0.0003828125 # # dxdy * dpdx * distance between these two pairs
            f = pf.flatten()
            b = pb.flatten()
            dp[k] = np.nanmean(f-b+volP)
            #cdeq[k] = np.mean(f-b+volP)/Utmp[k]
        cdeq = np.nanmean(dp)/Utmp




        U = np.sqrt(u**2+v**2) # spatially averaged velocity


        # domain shape

        xp = u.shape[2]
        yp = u.shape[1]
        zp = u.shape[0]

        # Momentum
        Tuv = c_tf(uv,u,v)
        Tuw = c_tf(uw,u,w)
        Tvw = c_tf(vw,v,w)

        Twc = c_tf(wc,w,s)
        Dwc = c_disp(w,s)

        Duv = c_disp(u,v)
        Duw = c_disp(u,w)
        Dvw = c_disp(v,w)


        speed = c_prof(U)
        vel = np.sqrt(c_prof(u)**2+c_prof(v)**2)  # wind speed

        stdS = c_std(U)

        TKE = 1/2*(c_tf(uu,u,u)+c_tf(vv,v,v)+c_tf(ww,w,w))+e

        dp = c_prof(np.gradient(p,0.5,axis=0))*16/utau**2


        delta = 0.5
        l = 0.5
        sgsDD = -c_prof((0.19+0.74*l/delta)*e**(3/2)/l)
        cm = 0.1
        
        l = 1.8*zu_3d
        
        #K_m = np.zeros(np.shape(e))
        
        K_m = l*cm*np.sqrt(e)
        
        t = 0
        for i in range(xp):
            for j in range(yp):
                if np.isnan(u[i,j,5]):
                    t+=1
        lp = t/xp/yp 
        
        [shearP,shearA,wakeP,turbT,disT,presT,sgsT,sgsD,presTV,turbTV,shearPV,disTV]=c_tkeb()



        mdic = {
            'name': fileName,
            ## TKE budget
            'shearProduction': shearP,
            'formA': shearA,
            'wakeProduction': wakeP,
            'turbTransport': turbT,
            'dispTransport': disT,
            'presTransport': presT,
            'sgsT': sgsT,
            'presTV':presTV,
            'turbTV':turbTV,
            'shearPV':shearPV,
            'disTV':disTV,
            'sgsDD':sgsDD,

            ## Profiles

            'u': c_prof(u),
            'v': c_profA(v), # abs value of crossstream
            'w': c_prof(w),
            's': c_prof(s),
            # momentum fluxes
            'Tuv':c_prof(Tuv),
            'Tuw':c_prof(Tuw),
            'Tvw':c_prof(Tvw),
            'Duv':c_prof(Duv),
            'Duw':c_prof(Duw),
            'Dvw':c_prof(Dvw),

            'Twc':c_prof(Twc),
            'Dwc':c_prof(Dwc),


            'U':c_prof(U), # wind speed
            'vel':vel, # spatially averaged velocity
            'stdS':stdS, # std of 
            'TKE':c_prof(TKE),

            # ped height value
            'U125': U[:,:,3],
            'u125': u[:,:,3],
            'v125': v[:,:,3],
            'TKE125': TKE[:,:,3],
            's125': s[:,:,3],
            'Tuw125': Tuw[:,:,3],
            'Tuv125': Tuv[:,:,3],
            'Tvw125': Tvw[:,:,3],
            'Duw125': Duw[:,:,3],
            'Duv125': Duv[:,:,3],
            'Dvw125': Dvw[:,:,3],

    # ped heigh No 2
            'U175': U[:,:,4],
            'u175': u[:,:,4],
            'v175': v[:,:,4],
            'TKE175': TKE[:,:,4],
            's175': s[:,:,4],
            'Tuw175': Tuw[:,:,4],
            'Tuv175': Tuv[:,:,4],
            'Tvw175': Tvw[:,:,4],
            'Duw175': Duw[:,:,4],
            'Duv175': Duv[:,:,4],
            'Dvw175': Dvw[:,:,4],

    # Canopy top
            'Uct': U[:,:,31],
            'uct': u[:,:,31],
            'vct': v[:,:,31],
            'TKEct': TKE[:,:,31],
            'sct': s[:,:,31],
            'Tuwct': Tuw[:,:,31],
            'Tuvct': Tuv[:,:,31],
            'Tvwct': Tvw[:,:,31],   
            'Duwct': Duw[:,:,31],
            'Duvct': Duv[:,:,31],
            'Dvwct': Dvw[:,:,31],

            # height levels

            'zu':zu_3d,
            'zw':zw_3d,

            # dp along with z-axis

            'dp':dp,
            'cdeq':cdeq,
        }

        # not good for shear production
        savemat(fileName+'all.mat',mdic)

