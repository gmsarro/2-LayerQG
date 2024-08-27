#Written by Joonsuk M. Kang
#Modified by Giorgio M. Sarro
def lwatend(lwa,dt):
	import numpy as np
	# time-tendency of LWA
	# calculated using leapfrog except at t=0 and -1.
	lwatend=np.zeros(np.shape(lwa))
	lwatend[1:-1,:,:]=(lwa[2:,:,:]-lwa[:-2,:,:])/(2*dt)
	lwatend[0,:,:]=(lwa[1,:,:]-lwa[0,:,:])/(dt)
	lwatend[-1,:,:]=(lwa[-1,:,:]-lwa[-2,:,:])/(dt)
	return lwatend


# 2. Advection by Uref
def urefadv(lwa,uref,dx,filt=True):
	import numpy as np
	# advection of LWA by uref
	# second-order finite differencing on periodic boundary
	# filt: if True, apply 1-2-1 filter in time
	lwagradx=(np.roll(lwa,-1,axis=2)-np.roll(lwa,1,axis=2))/(2*dx)
	out=-lwagradx*uref[:,:,np.newaxis]
	if filt == True:	out[1:-1,:]=out[:-2,:]*0.25+out[1:-1,:]*0.5+out[2:,:]*0.25
	return out


# 3. Advection by ue
def ueadv(q,qref,u,uref,dx,dy,filt=True):
	import numpy as np
	# advection of qe by ue d dx int^{eta}_{0} ue qe dy'
	# second-order finite differencing on periodic boundary
	# filt: if True, apply 1-2-1 filter in time
	tn,yn,xn=np.shape(q)
	Iuq=np.zeros(np.shape(q))
	for t in range(tn):
		for y1 in range(yn):
			q_e=q[t,:,:]-qref[t,y1]
			u_e=u[t,:,:]-uref[t,y1]
			for x in range(xn):
				for y2 in range(yn):
					if y2 <  y1 and q_e[y2,x] > 0:
						Iuq[t,y1,x]+=u_e[y2,x]*q_e[y2,x]*dy
					if y2 >= y1 and q_e[y2,x] <= 0:
						Iuq[t,y1,x]+=u_e[y2,x]*q_e[y2,x]*(-dy)
	out=(np.roll(Iuq,-1,axis=2)-np.roll(Iuq,1,axis=2))/(2*dx)
	if filt == True:	out[1:-1,:]=out[:-2,:]*0.25+out[1:-1,:]*0.5+out[2:,:]*0.25
	return -out
	

# 4 zonal eddy flux convergence
def eddyflux_x(ue,ve,dx,filt=True):
	import numpy as np
	# zonal eddy flux convergence, -1/2 d/dx (v^2 - u^2)
	# second-order finite differencing on periodic boundary
	# filt: if True, apply 1-2-1 filter in time
	v2_u2=0.5*(ve[:,:,:]**2-ue[:,:,:]**2)
	out=-(np.roll(v2_u2,-1,axis=2)-np.roll(v2_u2,1,axis=2))/(2*dx)
	if filt == True:	out[1:-1,:]=out[:-2,:]*0.25+out[1:-1,:]*0.5+out[2:,:]*0.25
	return out

# 5 meridional eddy flux convergence
def eddyflux_y(ue,ve,dy,filt=True):
	import numpy as np
	# meridional eddy flux convergence,  d/dy (uv)
	# second-order finite differencing, uv = 0 added at north and south
	# filt: if True, apply 1-2-1 filter in time
	uv=np.pad(ue[:,:,:]*ve[:,:,:],((0,0),(1,1),(0,0)),mode='constant',constant_values=0)
	out=(uv[:,2:,:]-uv[:,:-2,:])/(2*dy)
	if filt == True:	out[1:-1,:]=out[:-2,:]*0.25+out[1:-1,:]*0.5+out[2:,:]*0.25
	return out


# 6 vertical eddy flux convergence (heat flux)
def eddyflux_z(ve,te,Ld,filt=True):
	import numpy as np
	# heat flux convergence,  vT/Ld**2
	# filt: if True, apply 1-2-1 filter in time
	out=ve*te/(Ld**2)
	if filt == True:	out[1:-1,:]=out[:-2,:]*0.25+out[1:-1,:]*0.5+out[2:,:]*0.25
	return out


# 7 all EP flux
def eddyflux(ve,qe,filt=True):
	import numpy as np
	# EP flux
	out=-ve*qe
	if filt == True:	out[1:-1,:]=out[:-2,:]*0.25+out[1:-1,:]*0.5+out[2:,:]*0.25
	return out

# 3. Effect of LH
def LH(p,q,qref,L,dx,dy,filt=True):
	import numpy as np
	# effect of LH
	# second-order finite differencing on periodic boundary
	# filt: if True, apply 1-2-1 filter in time
	tn,yn,xn=np.shape(p)
	out=np.zeros(np.shape(p))
	for t in range(tn):
		for y1 in range(yn):
			q_e=q[t,:,:]-qref[t,y1]
			for x in range(xn):
				for y2 in range(yn):
					if y2 <  y1 and q_e[y2,x] > 0:
						out[t,y1,x]+=L*p[t,y2,x]*dy
					if y2 >= y1 and q_e[y2,x] <= 0:
						out[t,y1,x]+=L*p[t,y2,x]*(-dy)
	if filt == True:	out[1:-1,:]=out[:-2,:]*0.25+out[1:-1,:]*0.5+out[2:,:]*0.25
	return -out 

