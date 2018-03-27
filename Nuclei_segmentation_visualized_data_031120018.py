import numpy as np
import skimage.io as sk
from os import path
import os
import matplotlib.pyplot as plt
from scipy import ndimage

####Reading the data into an array (using dictionary)###

drt = '.../Kaggle 2018/stage1_train'  #Insert your own path to the folder here

subdirs = os.listdir(drt)
subdirs=subdirs[1:len(subdirs)] #cut off DS_Store
bigdict={}
names=[]
for s in subdirs:
        name=s[0:4]  ##only want first 4 characters of name in dictionary
        names.append(name)
        bigdict[name] = [[],[]]  #Set the 4-character name equal to two empty arrays, the first for the image array, the second for the masks
        im_dir = path.join(drt, s, 'images') #path to images folder
        mask_dir=path.join(drt, s, 'masks') #path to masks folder
        for i in os.listdir(im_dir):
            os.chdir(im_dir)
            print i
            bigdict[name][0].append(sk.imread(i)) #Append image array
        for i in os.listdir(mask_dir):
            os.chdir(mask_dir)
            bigdict[name][1].append(sk.imread(i)) #Append mask arrays
                      
##################       

### Create dictionary of images ###

imgdict={}
for n in names:
    img=bigdict[name][0][0] 
    imgdict[name]=[]
    imgdict[name].append(img) 
    shap=img.shape
    imgdict[name].append(shap)
    if len(shap)==3:
        imgdict[name].append(True)
    else:
        imgdict[name].append(False)        
                              
### Number of nuclei hist###
nucnum=[]

for n in names:
    nucnum.append(len(bigdict[n][1]))

f,ax=plt.subplots()    
ax.hist(nucnum,bins=25)
ax.set_xlabel('Number of Nuclei in Pic')
ax.set_ylabel('N')
ax.set_xlim([0,400])
ax.set_ylim([0,225])
plt.show()

### Im size hist###
sidelength=[]
for n in names:
    sidelength.append(len(bigdict[n][1][0]))
    
f,ax=plt.subplots()    
ax.hist(sidelength,bins=25)
ax.set_xlabel('Sidelength in Number of Pix')
ax.set_ylabel('N')
#ax.set_xlim([0,400])
#ax.set_ylim([0,225])
plt.show()

### Density hist ###
f,ax=plt.subplots()

ax.hist(np.log10([np.float(x) for x in nucnum]/np.power(sidelength,2)),bins=20)
ax.set_xlabel('Log(Density) (# of Nuclei/$pix^2$)')
ax.set_ylabel('N')
#ax.set_xlim([0,400])
#ax.set_ylim([0,225])
plt.show()

### Max min plot ###
maxs=[]
mins=[]

for n in names:
    maxs.append(max(np.concatenate(bigdict[n][0][0][:,:,0])))
    mins.append(min(np.concatenate(bigdict[n][0][0][:,:,0])))
f,ax=plt.subplots()
ax.scatter(mins,maxs,s=5)
ax.set_xlabel('Min Pix Value')
ax.set_ylabel('Max Pix Value')
plt.show()

### Std hist ###
stds=[np.std(np.concatenate(bigdict[n][0][0][:,:,0])) for n in names]

f,ax=plt.subplots()    
ax.hist(stds,bins=25)
ax.set_xlabel('Standard Deviation R')
ax.set_ylabel('N')
#ax.set_xlim([0,400])
#ax.set_ylim([0,225])
plt.show()

### Find sizes of nuclei (in number of pixels) ###

sizes=[]
for n in names:
    masks=bigdict[n][1]
    for m in masks:
        sizes.append(len(np.nonzero(np.concatenate(m))[0]))
        
### This code is for getting rid of the holes. Beware it only gets rid of holes of size one pixel ### 

### First need to fix inverted masks ###

for n in names:
    masks=bigdict[n][1]
    for m in masks:
        if len(np.nonzero(m)[0])>len(np.concatenate(m)):
            nzval=m[np.nonzero(m)[0][0]][np.nonzero(m)[1][0]] ##Need to know what the original nonzero value is (usually 255)
            #Go through each pixel and replace with opposite of it's original value
            for x in range(len(m)):
                for y in range(len(m[0])):
                    if masks[x][y]!=0:
                        masks[x][y]=0
                    else:
                        masks[x][y]=nzval

### Now it's time to patch the holes ###
               
for n in names:
    masks=bigdict[n][1]
    for m in range(len(masks)):
        print m
        nz_rows=np.nonzero(masks[m])[0]  # e.g (nz_row[i],nz_cols[i]) will be the coordinates of a nonzero pixel
        nz_cols=np.nonzero(masks[m])[1]
        uniqrows=np.unique(nz_rows) #Consolidate to unique rows (These are the rows that we know contain at least one nonzero pixel)
        for u in uniqrows:
            rowinds=np.where(nz_rows==u)[0] 
            col_corr=nz_cols[rowinds]  #These are the isolated column indices of nonzero pixels in the unique row
            coldiffs=np.ediff1d(col_corr) #Find the difference between each consecutive element
            if len(np.where(coldiffs==2)[0])>0: #Jump of 2 in indices indicates pixels to left and right are nonzero
                for i in np.where(coldiffs==2)[0]: #Following condition makes sure pixel is not on the edge and that the ones above and below it are nonzero (i.e it is a hole)
                    if (u!=0 and u!=len(masks[m])-1 and col_corr[i]!=0 and col_corr[i]!=len(masks[m][0])-1 and masks[m][u-1][col_corr[i]+1]!=0 and masks[m][u+1][col_corr[i]+1]!=0):
                        masks[m][u][col_corr[i]+1]=masks[m][u-1][col_corr[i]+1] #Value at hole (zero) now value of pixel directly below (random choice) 
                        
                
##############

### Sobel filtering example###  
img100=bigdict[names[100]][0][0]
plt.imshow(img100)
plt.show()
plt.figure()
plt.imshow(ndimage.sobel(img100))
plt.show()                                                         
        
