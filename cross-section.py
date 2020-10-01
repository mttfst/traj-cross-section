#!/usr/bin/env python
# coding: utf-8

# # Traj cross-section - Plotting of a cross-section along a trajectory
# ## Used packages

# In[8]:


import numpy as np
from netCDF4 import Dataset
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import datetime
from matplotlib.colors import LinearSegmentedColormap
from mpl_toolkits.axes_grid1.inset_locator import inset_axes, zoomed_inset_axes


# ## Functions
# Simple bi-linear interpolation. 

# In[9]:


def rot2grid(rlon,rlat,rlon_min,rlat_min,dlon,dlat):
    x=(rlon-rlon_min)/dlon
    y=(rlat-rlat_min)/dlat
    return x,y

def interpol(x,y,var):
    x0=int(x)
    x1=x0+1
    y0=int(y)
    y1=y0+1

    w0=(x-x0)*(y-y0)
    w1=(x1-x)*(y-y0)
    w2=(x1-x)*(y1-y)
    w3=(x-x0)*(y1-y)
    vel = w0*var[y1,x1]+w1*var[y1,x0]+w2*var[y0,x0]+w3*var[y0,x1]
    return vel


# ## read atmospheric data from NetCDF file

# In[10]:


f=Dataset('l.nc','r')

rglon = f.variables['rlon'][:]
rglat = f.variables['rlat'][:]
gdlon = round(rglon[1]-rglon[0],4)
gdlat = round(rglat[1]-rglat[0],4)

times = f.variables['time'][:]
time_unit=f.variables['time'].units
surfpres = f.variables['PS'][:]
surftemp = f.variables['T_2M'][:]
temp = f.variables['T'][:]
pres = f.variables['P'][:]
qv = f.variables['QV'][:]
td_2m = f.variables['TD_2M'][:]
tke = f.variables['TKE'][:]
hhl = f.variables['HHL'][0,:,0,0]

f.close()



# ## prepare the output variable(s)
# 
# In this example, I want to plot the virtual potential temperature together with the TKE.
# The virtual potential temperature needs to be calculated from temperature, pressure, and specific moisture.
# At the surface, I add an additional level with 2 m altitude. In the COSMO-model the levels are counted from top to bottom so add the extra level as level 51. 
# 
# 
# In the COSMO-model, the temperature is located on the so-called full levels, in contrast to that the TKE is located on the half levels. This means we have to interpolate the TKE from the half levels to the full levels.  
# 

# In[11]:


# init zlev and add one extra level for the 2m values
zlevs=np.zeros(temp.shape[1]+1)
zlevs[temp.shape[1]]=2.

# fill zlev with values
for i in range(len(hhl)-1):
    zlevs[i]=hhl[i+1]+(hhl[i]-hhl[i+1])/2. - hhl[len(hhl)-1]

    
# def theta   and tke
theta=np.zeros((temp.shape[0],temp.shape[1]+1,temp.shape[2],temp.shape[3]))
tke_fl=np.zeros(theta.shape)

# calc theta 
for i in range(len(zlevs)-1):
    theta[:,i,:,:]=(temp[:,i,:,:]*(1+0.61*qv[:,i,:,:]))*(pres[:,len(zlevs)-2,:,:]/pres[:,i,:,:])**(287./1005.)
    #theta[:,i,:,:]=(temp[:,i,:,:])*(pres[:,len(zlevs)-2,:,:]/pres[:,i,:,:])**(287./1005.)
    tke_fl[:,i,:,:]=(tke[:,i,:,:]+tke[:,i+1,:,:])/2.



# calc  2m values
e_2m=6.1078*np.exp(17.1*(td_2m-273.15)/(235+(td_2m-273.15)))
qv_2m=0.66*e_2m/surfpres
theta[:,len(zlevs)-1,:,:]=surftemp*(1+0.61*qv_2m)
tke_fl[:,len(zlevs)-1,:,:]=(tke_fl[:,len(zlevs)-2,:,:]+tke[:,len(zlevs)-1,:,:])/2.





# ## Read in the trajectory data 

# In[12]:


# load trajectory data
mean_traj = np.loadtxt('mean_traj.txt')
stda_traj = np.loadtxt('stda_traj.txt')

# save the column with the dates
traj_timestemps = [str(int(mean_traj[i,0])) for i in range(len(mean_traj))]


# delete the dates from the rest of the data
mean_traj=np.delete(mean_traj,0,axis=1)
stda_traj=np.delete(stda_traj,0,axis=1)

#for l in zlevs:
#    print(int(l))


# ## convert time

# In[13]:



ref_year  = int(time_unit[14:18])
ref_month = int(time_unit[19:21])
ref_day   = int(time_unit[22:24])
ref_h     = int(time_unit[25:27])
ref_min   = int(time_unit[28:30])
ref_s     = int(time_unit[31:33])



refdate = datetime.datetime(ref_year,ref_month,ref_day,ref_h,ref_min,ref_s)
reftime = datetime.time(ref_h,ref_min,ref_s)

mpl_dates = []
timestemps=[]
dates=[]
for t in times:
    time_diff = datetime.timedelta(seconds=int(t))
    date = refdate + time_diff  
    mpl_dates.append(mpl.dates.date2num(date))
    dates.append(date)
    timestemps.append(date.strftime("%Y%m%d%H%M%S"))
    

traj_dates = []
traj_mpl_dates = []
for t in traj_timestemps:
    date = datetime.datetime.strptime(t,'%Y%m%d%H%M%S')
    traj_dates.append(date)
    traj_mpl_dates.append(mpl.dates.date2num(date))

# i=0
# for d in mpl_dates:
#     print(i,mpl.dates.num2date(d),timestemps[i])
#     i+=1


# ## define the time range for the plot

# In[14]:


start_date = '20170531080000'
end_date   = '20170531230000'


# ## extend the trajectory
# 
# I want to start the cross-section plot bevor the trajectory starts, therefore I need to extend the trajectory at the start and at the end. 
# 

# In[15]:


# I store the data in individual lists and use the .append function because it is handy. 

# trajectory
trax  = [] 
tray  = [] 
traz  = [] 
trazr = [] 

trazr_std1 = [] 
trazr_std2 = [] 

# standard deviation
stdx  = [] 
stdy  = [] 
stdz  = [] 
stdzr = []

sdate = datetime.datetime.strptime(start_date,'%Y%m%d%H%M%S')
edate = datetime.datetime.strptime(end_date,'%Y%m%d%H%M%S')

print (sdate, edate)

# extend the begin of the tray
while True:
    if traj_dates[0] == sdate:
        break
    else:
        traj_dates = [traj_dates[0] - datetime.timedelta(seconds=30)] + traj_dates
        
        trax.append( mean_traj[0,0])
        tray.append( mean_traj[0,1])
        traz.append( None) # here i add None because i dont want to have this line in the plot later. 
        trazr.append(None)
        trazr_std1.append(None)
        trazr_std2.append(None)
        
        stdx.append( stda_traj[0,0])
        stdy.append( stda_traj[0,1])
        stdz.append( None) # here i add None because i dont want to have this line in the plot later. 
        stdzr.append(None)
        
# add the original traj values        
for m,s in zip(mean_traj,stda_traj):
    trax.append( m[0])
    tray.append( m[1])
    traz.append( m[2]) # here i add None because i dont want to have this line in the plot later. 
    trazr.append(m[3])
    trazr_std1.append(m[3]+s[3])
    trazr_std2.append(m[3]-s[3])
    
    stdx.append( s[0])
    stdy.append( s[1])
    stdz.append( s[2]) # here i add None because i dont want to have this line in the plot later. 
    stdzr.append(s[3])

    
# extend the end of the traj
while True:
    if traj_dates[-1] == edate:
        break
    else:
        traj_dates.append(traj_dates[-1] + datetime.timedelta(seconds=30)) 
        
        trax.append( mean_traj[-1,0])
        tray.append( mean_traj[-1,1])
        traz.append( None) # here i add None because i dont want to have this line in the plot later. 
        trazr.append(None)
        trazr_std1.append(None)
        trazr_std2.append(None)
        
        stdx.append( stda_traj[-1,0])
        stdy.append( stda_traj[-1,1])
        stdz.append( None) # here i add None because i dont want to have this line in the plot later. 
        stdzr.append(None)


    
traj_mpl_dates = []    
for d in traj_dates:
    traj_mpl_dates.append(mpl.dates.date2num(d))
    

    
trax  = np.asarray(trax)
tray  = np.asarray(tray)
traz  = np.asarray(traz)
trazr = np.asarray(trazr)

stdx  = np.asarray(stdx)
stdy  = np.asarray(stdy)
stdz  = np.asarray(stdz)
stdzr = np.asarray(stdzr)


# ## interpolate fields to the trajectory

# In[16]:


# count number of time steps for the cross section
i=0
for t in timestemps:
    if t == start_date:
        first_time_step=i
    if t == end_date:
        last_time_step=i
    i+=1


# init the array that be plotted
plotvar1 = np.zeros((23,last_time_step-first_time_step+1))
plotvar2 = np.zeros(plotvar1.shape)

# loop of time steps in the plot
c=0  
for t in range(first_time_step,last_time_step+1):
    
    # loop to find traj time that matches the polt time
    for i in range(len(traj_dates)):
        if traj_dates[i] == dates[t]:
            
            # transform rot coord. to grid coord.
            xg,yg=rot2grid(trax[i],tray[i],rglon[0],rglat[0],gdlon,gdlat)
            
            # interpolate to trajectory position
            for level in range(plotvar1.shape[0]):
                inter=interpol(xg, yg, theta[t,28+level,:,:])
                plotvar1[level,c]=inter
                
                inter=interpol(xg, yg, tke_fl[t,28+level,:,:])
                plotvar2[level,c]=inter
                
            
            c+=1
            break


# ### Emission time 
# The trajectory I use here is the mean trajectory of a huge set of particle trajectories. I want to add the emission time of these particles in the plot. I simply hardcoded this. 

# In[17]:


# an array that contains values on the emission times
emission_time=np.zeros(len(traj_dates))
emission_time[:]=None

emission_time[100:210]=0
emission_time[440:540]=0

#for i in range(len(traj_dates)):
#    print(i,traj_dates[i])


# ## plot

# In[18]:


# original

# create figure
fig = plt.figure(figsize=(11,5))

ax = fig.add_subplot() # main plot 

# grid for the plot. I only use the bottom-most levels [28:]
X, Y = np.meshgrid(dates[first_time_step:last_time_step+1], zlevs[28:])

# color maps
cmap = LinearSegmentedColormap.from_list(name='mycmap',colors =['darkblue','royalblue','lightskyblue','#c6dbef','#fee6ce','#fdd0a2','#fdae6b','#fd8d3c', '#f16913','#d94801'],N=64)
cmap2 = LinearSegmentedColormap.from_list(name='mycmap',colors =['black','darkgray'],N=64)

# The emission time as plot_date
ax.plot_date(traj_mpl_dates,emission_time,fmt='-',linewidth=4,color='green',label='emission time')

# filled contour
con  = ax.contourf(X,Y,plotvar1[:,:], levels=np.arange(290,306),extend='both',cmap=cmap)

# non-filled contour
con2 = ax.contour(X,Y,plotvar2,cmap=cmap2,levels=np.arange(0,7),extend='both')

# mean trajectory and standard deviation
ax.plot_date(traj_mpl_dates,trazr,fmt='-',linewidth=1.5,color='crimson',label='mean trajectory')
ax.plot_date(traj_mpl_dates,trazr_std1,fmt='-',linewidth=1,color='palevioletred',label='standard deviation')
ax.plot_date(traj_mpl_dates,trazr_std2,fmt='-',linewidth=1,color='palevioletred')


# axis settings
ax.set_ylim(0,3000)
ax.set_ylabel('Altitude [m]',size='x-large')
ax.xaxis.set_major_formatter(mpl.dates.DateFormatter('%H:%M'))
ax.tick_params(labelsize='x-large')


# color bar
cb = fig.colorbar(con,ax=ax)
cb.ax.tick_params(labelsize='x-large')
cb.set_label('Virtual potential temperature [K]',size='x-large',labelpad=15)



#legend
custom_lines = [mpl.lines.Line2D([0], [0], color='crimson', lw=3),
                mpl.lines.Line2D([0], [0], color='palevioletred', lw=3),
                mpl.lines.Line2D([0], [0], color='green', lw=3),
                mpl.lines.Line2D([0], [0], color='black', lw=3),
                mpl.lines.Line2D([0], [0], color='darkgray', lw=3)]
ax.legend(custom_lines, ['mean trajectory', 'standard deviation', 'emission time','TKE > 0 J kg$^{-1}$','TKE > 5 J kg$^{-1}$'], loc='upper left')




plt.tight_layout()

# output
plt.savefig('cross.png' ,dpi=500 )


# ### Nested plot
# The same plot as before but with a zoom of the bottom-most layer in the upper left corner

# In[19]:


# nested

# create figure
fig = plt.figure(figsize=(11,5))

ax = fig.add_subplot() # main plot 

# grid for the plot. I only use the bottom-most levels [28:]
X, Y = np.meshgrid(dates[first_time_step:last_time_step+1], zlevs[28:])

# color maps
cmap = LinearSegmentedColormap.from_list(name='mycmap',colors =['darkblue','royalblue','lightskyblue','#c6dbef','#fee6ce','#fdd0a2','#fdae6b','#fd8d3c', '#f16913','#d94801'],N=64)
cmap2 = LinearSegmentedColormap.from_list(name='mycmap',colors =['black','darkgray'],N=64)

# The emission time as plot_date
ax.plot_date(traj_mpl_dates,emission_time,fmt='-',linewidth=4,color='green',label='emission time')

# filled contour
con  = ax.contourf(X,Y,plotvar1[:,:], levels=np.arange(290,306),extend='both',cmap=cmap)

# non-filled contour
con2 = ax.contour(X,Y,plotvar2,cmap=cmap2,levels=np.arange(0,7),extend='both')

# mean trajectory and standard deviation
ax.plot_date(traj_mpl_dates,trazr,fmt='-',linewidth=1.5,color='crimson',label='mean trajectory')
ax.plot_date(traj_mpl_dates,trazr_std1,fmt='-',linewidth=1,color='palevioletred',label='standard deviation')
ax.plot_date(traj_mpl_dates,trazr_std2,fmt='-',linewidth=1,color='palevioletred')


# axis settings
ax.set_ylim(0,3000)
ax.set_ylabel('Altitude [m]',size='x-large')
ax.xaxis.set_major_formatter(mpl.dates.DateFormatter('%H:%M'))
ax.tick_params(labelsize='x-large')

# color bar
cb = fig.colorbar(con,ax=ax)
cb.ax.tick_params(labelsize='x-large')
cb.set_label('Virtual potential temperature [K]',size='x-large',labelpad=15)

#legend
custom_lines = [mpl.lines.Line2D([0], [0], color='crimson', lw=3),
                mpl.lines.Line2D([0], [0], color='palevioletred', lw=3),
                mpl.lines.Line2D([0], [0], color='green', lw=3),
                mpl.lines.Line2D([0], [0], color='black', lw=3),
                mpl.lines.Line2D([0], [0], color='darkgray', lw=3)]
ax.legend(custom_lines, ['mean trajectory', 'standard deviation', 'emission time','TKE > 0 J kg$^{-1}$','TKE > 5 J kg$^{-1}$'], loc='upper right')


# nested plot
axins = inset_axes(ax,1.8,1.2,loc='upper left')

# filled contour only, otherwise to messy
con  = axins.contourf(X,Y,plotvar1[:,:], levels=np.arange(290,306),extend='both',cmap=cmap)

# axis settings
axins.set_xlim(dates[first_time_step],dates[first_time_step+16])
axins.set_ylim(2,60)
axins.xaxis.set_major_formatter(mpl.dates.DateFormatter('%H:%M'))
axins.yaxis.set_tick_params(right=True,labelright=True,left=False,labelleft=False)
axins.set_yticks([10,20,30,40,50])
axins.tick_params(axis='x', labelrotation=-30)
axins.set_xticks(axins.get_xticks()[::2])
plt.xticks(ha='left')


plt.tight_layout()

# output
plt.savefig('cross.png' ,dpi=500 )


# In[ ]:




