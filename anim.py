import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation


  

t=[]
theta=[]
with open(file="./data.dat") as f:
    # read l(n)
    f.readline()
    l = np.array(list(map(float, f.readline().split())))
    # read delta_t
    f.readline()
    delta_t=float(f.readline())
    max_frame_num=int(30/delta_t)

    for i in range(max_frame_num):
        lis = np.array(list(map(float, f.readline().split())))
        if len(lis)==0:
            break

        t.append(lis[0])
        theta.append(lis[1:])

interval=0.2
def update(frame):
    ax.cla()  # ax をクリア
    ax.set_xlabel("x[m]")
    ax.set_ylabel("y[m]")
    ax.set_xlim(-lim_abs,lim_abs)
    ax.set_ylim(-lim_abs,lim_abs)

    xdat=[]
    ydat=[]
    xdat.append(0)
    ydat.append(0)
    use_frame=frame*int(0.1/delta_t)
    # print(use_frame)
    ax.set_title(f"t={round(delta_t*use_frame,3):.2f}s")

    for i in range(len(theta[0])):
        xdat.append(l[i]*np.sin(theta[use_frame][i])+xdat[-1])
        ydat.append(-l[i]*np.cos(theta[use_frame][i])+ydat[-1])

    xdat=np.array(xdat)
    ydat=np.array(ydat)
    ax.plot(xdat,ydat,marker=".")

  
g=9.8
lim_abs=sum(l)*(1+1/10)

fig=plt.figure()
ax=fig.subplots()
ax.set_aspect("equal")

anim = FuncAnimation(fig, update, frames=range(min(max_frame_num,len(theta)//int(0.1/delta_t))), interval=10)
anim.save(f"anim{len(theta[0])}.gif")