from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib.pyplot as plt
import numpy as np


from . import evalgeo


def print_domain(domain):
    fig, ax = plt.subplots(figsize=(6,6))
    im = ax.imshow(domain)
    cb = plt.colorbar(im)
    cb.set_label('Estimated building height [m]')
    plt.tight_layout()
    plt.savefig('domain.pdf')



def print_entropy(x,y):


    xtick_font = {
        "family": "DejaVu Sans",
        "size": 10,
        "weight": "bold",
        "alpha": 1.0,
        "zorder": 3,
    }

    color="#003366"
    edgecolor="k"
    linewidth=0.5
    alpha=0.7

    #dir_topo = '/Users/jiachenlu/forSimulationOnly/highRes/'

    fig, ax = plt.subplots(figsize=(6,6), subplot_kw={"projection": "polar"})

    ax.set_theta_zero_location("N") # 
    ax.set_theta_direction("clockwise")
    ax.set_ylim(top=y.max())

    # configure the y-ticks and remove their labels
    ax.set_yticks(np.linspace(0, y.max(), 5))
    ax.set_yticklabels(labels="")

    # configure the x-ticks and their labels
    xticklabels = ["N", "", "E", "", "S", "", "W", ""]
    ax.set_xticks(ax.get_xticks())
    ax.set_xticklabels(labels=xticklabels, fontdict=xtick_font)
    ax.tick_params(axis="x", which="major", pad=-2)

    ax.bar(
        x,
        height=y,
        width=0.1,
        align="center",
        bottom=0,
        zorder=2,
        color=color,
        edgecolor=edgecolor,
        linewidth=linewidth,
        alpha=alpha,
    )
    ax.set_thetamin(-90)
    ax.set_thetamax(90)
    plt.tight_layout()
    plt.savefig('polarO.png',dpi=300)
    
    dist = 2 # Edge evenly distributed to just 2 orientations
    hMin = 1/dist*np.log(1/dist)*dist

    dist = 18 # Edge evenly distributed to all 18 orientations
    hMax = 1/dist*np.log(1/dist)*dist

    H_w = np.nansum(y/sum(y)*1*np.log(y/sum(y)*1))
    phi = 1-(hMax-H_w)/(hMax-hMin)
    
    return(phi)

def print_results(phi,domain,area,areaH,bldH):
    print('$\phi=$'+str(phi))
    print(domain.shape)
    print('Rasterization error = '+str((evalgeo.cal_lp(domain)-np.array(area).sum()/domain.size)/evalgeo.cal_lp(domain)*100)[:4]+'%')
    print('{}%  real building height recognized, the rest is set to {} m'.format(str(np.array(areaH).sum()/np.array(area).sum()*100)[:4],bldH))
    return()



def showDiagram(domain,angle,weighted,phi,angleRotate,area,name):
        # Triming for proper size and fit the computational request for parallelization
    # 300,330
    dir_topo = ''
    fig = plt.figure(figsize=(12, 6))


    ax1 = fig.add_subplot(1,3,1)
    ax2 = fig.add_subplot(1,3,2)

    ax1.contourf(domain)
    ax1.axis('equal')


    lp = evalgeo.cal_lp(domain)
    s = (1-lp)*0.001

    ax1.set_title('Computaional domain ' +'$\\lambda_p$ = '+str(lp)[0:6]+'\n'+ name + ' Domain size = '
                  + str(domain.shape))

    np.savetxt(dir_topo+name+'_topo',domain,fmt='%d')



    ax2.contour(domain,linewidths=0.1,colors='r')    
    ax2.axis('equal')
    ax2.set_title('Top sink of scalar = 1e-7*'+str(s)[5:9])


    ax31= fig.add_subplot(133, polar=True)



    xtick_font = {
        "family": "DejaVu Sans",
        "size": 15,
        "weight": "bold",
        "alpha": 1.0,
        "zorder": 3,
    }

    color="#003366"
    edgecolor="k"
    linewidth=0.5
    alpha=0.7


    x = np.unique(angle)-angleRotate
    

    x = np.deg2rad(x)
    y = np.array(weighted)
    
    x = np.concatenate((x,x+np.pi),axis=0)
    y = np.concatenate((y,y),axis=0)


    ax31.set_theta_zero_location("N")
    ax31.set_theta_direction("clockwise")
    ax31.set_ylim(top=y.max())

    # configure the y-ticks and remove their labels
    ax31.set_yticks(np.linspace(0, y.max(), 5))
    ax31.set_yticklabels(labels="")

    
    
    # configure the x-ticks and their labels
    xticklabels = ["N", "", "E", "", "S", "", "W", ""]
    ax31.set_xticks(ax31.get_xticks())
    ax31.set_xticklabels(labels=xticklabels, fontdict=xtick_font)
    ax31.tick_params(axis="x", which="major", pad=-2)

    ax31.bar(
        x,
        height=y,
        width=0.1,
        align="center",
        bottom=0,
        zorder=2,
        color=color,
        edgecolor=edgecolor,
        linewidth=linewidth,
        alpha=alpha,
    )
    
    #ax31.set_thetamin(-90)
    #ax31.set_thetamax(90)
    
    
    gamma,prof = evalgeo.cal_alignness(domain)
    

    ax32= fig.add_subplot(8,3,24)
    title = 'Orientations of building edges \n Edge entropy $\phi$  = '+ str(phi)[:6]
    title += '\n Rotation angle  = '+str(angleRotate)+'$^{\circ}$'
    title += '\n Average building size $A_0$  = '+str(area.mean())[:6]+'$m^2$'
    ax31.set_title(title)


    ax32.plot(prof,c='k')
    ax32.set_xlim(0,domain.shape[0])
    #ax32.frame('off')
    ax32.set_xlabel('Alignedness $\gamma$ = '+str(gamma)[0:6])


    plt.tight_layout()
    plt.savefig(dir_topo + name + '.png',dpi=300)
