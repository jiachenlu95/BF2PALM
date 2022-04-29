import numpy as np
import matplotlib.pyplot as plt
import BF2PALM as bm
names = globals()

def showDiagram(nn,angle,weighted,phi,area,angleRotate,dirName,name):
        # Triming for proper size and fit the computational request for parallelization
    # 300,330
    fig = plt.figure(figsize=(12, 6))


    ax1= fig.add_subplot(1,3,1)
    ax2= fig.add_subplot(1,3,2)
    ax1.contourf(nn)
    ax1.axis('equal')


    t = 0
    for x in range(nn.shape[0]):
        for y in range(nn.shape[1]):
            if nn[x,y] == 0:
                t+=1



    lp = 1-t/nn.shape[0]/nn.shape[1]
    s = (1-lp)*0.001
    ax1.set_title('Computaional domain ' +'$\\lambda_p$ = '+str(lp)[0:6]+'\n'+ name + ' Domain size = '
                  + str(nn.shape))

    np.savetxt(dirName+name+'_topo',nn,fmt='%d')



    ax2.contour(nn,linewidths=0.1,colors='r')    
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
    
    gamma = bm.cal_alignness(nn)
    

    ax32= fig.add_subplot(8,3,24)
    title = 'Orientations of building edges \n Edge entropy $\phi$  = '+ str(phi)[:6]
    title += '\n Alignnedness $\\gamma$  = '+str(gamma)
    title += '\n Average building size $A_0$  = '+str(area.mean())[:6]+'$m^2$'
    ax31.set_title(title)

    prof,var = bm.cal_variance(nn)
    ax32.plot(prof,c='k')
    ax32.set_xlim(0,nn.shape[0])
    #ax32.frame('off')
    ax32.set_xlabel('Variance = '+str(var)[0:6])

    plt.tight_layout()
    plt.savefig(dirName+name+'.png',dpi=300)
