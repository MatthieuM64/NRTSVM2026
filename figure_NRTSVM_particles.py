#! /usr/bin/python
# -*- coding:utf8 -*-

import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import scipy.integrate as Int
import scipy.optimize as opt
import os
import sys
import time
import multiprocessing

fontsize=20
plt.rc("font",size=fontsize)
plt.rc('font', family='serif')
plt.rc('text', usetex=True)
#matplotlib.use('TkAgg')

color=['#ff0000','#ff6600','#00ff00','#006600','#00ffff','#0000ff','#cc66ff','k']
fmt='os>*^hv'

clock=time.time()

#Physical parameters
v0=0.0015
eta=0.02
rho0=64
JAB=0.25
JBA=-0.25
LX=12
LY=12
init=3
ran=0

#Time intervals
DT=1
tmax=70000
t0=tmax-2000

#Movie creation
NCPU=4
multi=True
movie=False

#Read parameters in command line
for arg in sys.argv[1:]:
	if "-v0=" in arg:
		v0=float(arg[4:])
	elif "-eta=" in arg:
		eta=float(arg[5:])
	elif "-rho0=" in arg:
		rho0=float(arg[6:])
	elif "-JAB=" in arg:
		JAB=float(arg[5:])
	elif "-JBA=" in arg:
		JBA=float(arg[5:])
	elif "-LX=" in arg:
		LX=int(arg[4:])
	elif "-LY=" in arg:
		LY=int(arg[4:])
	elif "-init=" in arg:
		init=int(arg[6:])
	elif "-ran=" in arg:
		ran=int(arg[5:])
	elif "-t0=" in arg:
		t0=int(arg[4:])
	elif "-tmax=" in arg:
		tmax=int(arg[6:])
	elif "-NCPU=" in arg:
		NCPU=int(arg[6:])
	elif "-movie"==arg:
		movie=True
	else:
		print("Bad Argument: ",arg)
		sys.exit(1)
		
if NCPU==1:
	multi=False
elif NCPU<1:
	print("Bad value of NCPU: ",NCPU)
	sys.exit(1)

arrlen=0.15
lw=0.001
dpi=320

#Create one snapshot (one frame of the movie)
def Snapshot(i):
	if not os.path.isfile('snapshots/figure_TSVM_v0=%.8g_eta=%.8g_rho0=%.8g_JAB=%.8g_JBA=%.8g_LX=%d_LY=%d_init=%d_ran=%d_%d.png'%(v0,eta,rho0,JAB,JBA,LX,LY,init,ran,i)):
		t=t0+i*DT
		
		fig=plt.figure(figsize=(6,6))
		gs=matplotlib.gridspec.GridSpec(1,1,width_ratios=[1],height_ratios=[1],left=0.085,right=0.965,bottom=0.06,top=0.94,hspace=0.1,wspace=0.1)

		ax=plt.subplot(gs[0,0])
		
		data=np.fromfile('data_TSVM_particles/TSVM_particles_v0=%.8g_eta=%.8g_rho0=%.8g_JAB=%.8g_JBA=%.8g_LX=%d_LY=%d_init=%d_ran=%d_t=%d.bin'%(v0,eta,rho0,JAB,JBA,LX,LY,init,ran,t),dtype=np.float32).reshape(-1, 3)
		X,Y,THETA=data.T
		Npart=len(data)
		
		ax.quiver(X,Y,arrlen*np.cos(THETA),arrlen*np.sin(THETA),np.arange(Npart)%2,cmap=colors.ListedColormap(['r', 'b']),
		           angles='xy', scale_units='xy', scale=1, pivot='tail', width=lw, headwidth=5, headlength=5, rasterized=True)
		
		plt.axis('equal')
		plt.xlim([0,LX])
		plt.ylim([0,LY])
		plt.xticks([0,0.25*LX,0.5*LX,0.75*LX,LX])
		plt.yticks([0,0.25*LY,0.5*LY,0.75*LY,LY])
		ax.xaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('$%.8g$'))
		ax.yaxis.set_major_formatter(matplotlib.ticker.FormatStrFormatter('$%.8g$'))
		
		plt.text(0.01*LX,0.01*LY,'$t=%d$'%(t),ha='left',va='bottom',fontsize=20,bbox=dict(facecolor="white", alpha=0.5, edgecolor="none"),zorder=1)
		plt.text(0.5*LX,1.03*LY,'$v_0=%.8g$, $\eta=%.8g$, $\\rho_0=%.8g$, $\\mu=%.8g$'%(v0,eta,rho0,JAB),ha='center',va='center',fontsize=20)
		
		plt.savefig('snapshots/figure_TSVM_v0=%.8g_eta=%.8g_rho0=%.8g_JAB=%.8g_JBA=%.8g_LX=%d_LY=%d_init=%d_ran=%d_%d.png'%(v0,eta,rho0,JAB,JBA,LX,LY,init,ran,i),dpi=dpi)
		plt.close()
		
		print('-snap=%d/%d -t=%d -tcpu=%d'%(i+1,Nsnap,t,time.time()-clock))
		del fig,data

#Find the datafile and create the corresponding snapshots
os.system('mkdir -p snapshots')

i=0
t=t0
ARG=[]
while os.path.isfile('data_TSVM_particles/TSVM_particles_v0=%.8g_eta=%.8g_rho0=%.8g_JAB=%.8g_JBA=%.8g_LX=%d_LY=%d_init=%d_ran=%d_t=%d.bin'%(v0,eta,rho0,JAB,JBA,LX,LY,init,ran,t)) and t<=tmax:
	ARG.append(i)
	i+=1
	t+=DT

i=0
SNAP=[]
while os.path.isfile('snapshots/figure_TSVM_v0=%.8g_eta=%.8g_rho0=%.8g_JAB=%.8g_JBA=%.8g_LX=%d_LY=%d_init=%d_ran=%d_%d.png'%(v0,eta,rho0,JAB,JBA,LX,LY,init,ran,i)) and i*DT<=tmax:
	SNAP.append(i)
	i+=1
	
Nsnap=len(ARG)
print('%d Snapshots (%d already created)'%(Nsnap,len(SNAP)))

if multi:
	pool=multiprocessing.Pool(NCPU)
	results=pool.imap_unordered(Snapshot,ARG)
	pool.close()
	pool.join()
else:
	for i in ARG:
		Snapshot(i)

if movie:	
	os.system('mkdir -p movies')
	os.system('ffmpeg -v quiet -stats -y -r 25/1 -i snapshots/figure_TSVM_v0=%.8g_eta=%.8g_rho0=%.8g_JAB=%.8g_JBA=%.8g_LX=%d_LY=%d_init=%d_ran=%d_%%01d.png -c:v h264 -r 25 -s %dx%d -crf 20 movies/movie_TSVM_v0=%.8g_eta=%.8g_rho0=%.8g_JAB=%.8g_JBA=%.8g_LX=%d_LY=%d_init=%d_ran=%d.mp4'%(v0,eta,rho0,JAB,JBA,LX,LY,init,ran,6*dpi,6*dpi,v0,eta,rho0,JAB,JBA,LX,LY,init,ran))
	if Nsnap>(tmax-t0)/DT:
		os.system('rm snapshots/figure_TSVM_v0=%.8g_eta=%.8g_rho0=%.8g_JAB=%.8g_JBA=%.8g_LX=%d_LY=%d_init=%d_ran=%d_*.png'%(v0,eta,rho0,JAB,JBA,LX,LY,init,ran))
	
print('OK - time=%d sec'%(time.time()-clock))
