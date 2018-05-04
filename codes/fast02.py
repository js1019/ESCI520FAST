#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 21 13:05:36 2018

@author: jpcorn
python 3.6
"""

import obspy
from scipy.signal import spectrogram

from skimage.transform import resize
import numpy as np
import pywt
from numpy import linalg as LA
from scipy import stats, sparse
import matplotlib.pyplot as plt
from collections import Counter

data_downsample = 500
Ni = 800 # number of indices to keep
specwindow = 50
speclag = 2
specimage = 100
specimagelag = 20
#rownew, colnew = np.mgrid[0:32:1, 0:200:200/128]# for interpolation; not used
# pre-allocate arrays
binaryfinger = np.zeros((64,64))
haarspec = np.zeros((32,64))
MSH = np.array([],dtype='int8').reshape(0,500)
hs= np.array([],dtype='int8').reshape(0,64)

# computes fingerprint by replacing 1 with 10, -1 with 01, and 0 with 00
def fingerprint(data,bfp):
    bfp = np.zeros((64,64),dtype='int8')
    for k in range(0,64):
        for i in range(0,32):
            if data[i,k] == 1:
                bfp[i*2:i*2+2,k] = np.array([1,0])
            elif data[i,k] == -1:
                bfp[i*2:i*2+2,k] = np.array([0,1])
    return bfp
        
def haarwavelet(sp,H,binary):
    # H and binary are matrices of zeros
    # compute haar wavelet according to https://pywavelets.readthedocs.io/en/latest/
    (w1,(w2,w3,w4)) = pywt.dwt2(data=sp, wavelet='haar')
    H[0:16,0:32] = w1
    H[0:16,32:64]=w2
    H[16:32,0:32]=w3
    H[16:32,32:64]=w4
    
    # divide by norm
    Hhat=H/LA.norm(H,2,0)
    # Zscores of each element
    Z=stats.zscore(Hhat, axis=1, ddof=1)
    flat_indices = np.argpartition(-1*np.abs(Z.ravel()), Ni-1)[:Ni] # 800 largest values
    row_indices, col_indices = np.unravel_index(flat_indices, Z.shape) #get their indices
    # keep signs only
    binary[row_indices,col_indices] = np.sign(Z[row_indices,col_indices]) 
    return binary


# binfreq bins and averages frequencies from N elements down to 32 elements
# resamples spectrogram down to 32 frequencies
def binfreq(arr):
    newspec = np.zeros((32,arr.shape[1]))
    # if not divisible by 32, last bin will have a smaller interval than the others
    if arr.shape[0]%32 is not 0:
        bins=arr.shape[0]-arr.shape[0]%32
        for i in range(0,arr.shape[1]):
            newspec[0:32,i] = np.mean(arr[0:bins,i].reshape(-1, arr.shape[0]//32), axis=1) #mean of the other bins
            newspec[31,i] = np.mean([newspec[31,i],arr[bins::,i]]) #mean of the last bin
    else: # arr is divisible by 32
        for i in range(0,arr.shape[1]):
            newspec[:,i] = np.mean(arr[:,i].reshape(-1, arr.shape[0]//32), axis=1)
    return newspec #returns downsmapled spectrogram


def hashish(F,h1):
    # input fingerprint and all 500 permutation matrices
    A = sparse.find(F) #find indices of 1's
    p = np.zeros((500,),dtype='int8')
    for i in range(0,500):
        p[i] = np.min(h1[A[0]+i*64,A[1]]) # find first 1 in the permuted matrix
    return p # returns 500-element long minhash for each image

# computes matching minhash , 
def lsh(data):
    u,uind = np.unique(data,axis=0,return_inverse=True)
    # returns array of length(# of images), each element references signature
    # matching images will reference the same signature
    # order matters, but should it???
    return uind 

#sac datafiles
datafile = "./ZA20.SHZ.SAC"
#datafile2="./ZA20.SHZ.SAC"

# read in the sac files: http://docs.obspy.org/tutorial/
st = obspy.read(datafile,format='sac') #+ obspy.read(datafile2)

#filter the data
kargs = {'freqmin':10,'freqmax':200,'corners':4,'zerophase':True}
st.filter('bandpass',**kargs)
st.resample(data_downsample)
print('filter complete')

#compute spectrogram according to 
# https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.spectrogram.html
f, t, Sxx = spectrogram(st[0].data, fs=data_downsample,window=('hamming'),nperseg=specwindow,noverlap=speclag,nfft=512,scaling='spectrum',mode='magnitude')
print('spectrogram complete')

# downsample spectrogram
spec = binfreq(Sxx)
print('binning complete')

# number of windows to expect
Nfp = int(np.floor((spec.shape[1]-(specimage-specimagelag))/specimagelag))

# create permutation matrices
# better as columns? # should it be 0:4096?
# how does the paper achieve 0:255?
for r in range(0,500):
    hs = np.vstack((hs,np.random.randint(0,4096,size=(64,64))))

# how to choose 64 columns of image from 100? randomly? This provides that
#s=np.random.randint(0,100,size=(64,))
    
for i in range(1,Nfp+1):
    if i > -1: #for testing purposes
        window10s200  = spec[:,(i-1)*specimagelag:(i-1)*specimagelag+specimage]
        #window10s = window10s200[:,s] #downsample choosing random columns
        window10s = resize(window10s200 , (32,64), order=3,mode='constant') #downsample by interpolation
        #window10s = window10s200[:,np.arange(0,specimage,specimage/64,'int')] #downsample by equispaced indices
        haarz = haarwavelet(window10s,haarspec*0,haarspec*0)
        fngr = fingerprint(haarz,binaryfinger*0)
        p=hashish(fngr,hs)
        MSH=np.vstack((MSH,p))
        
print('hash table complete') 

buckets = []

# create buckets of matching bins (length 4 signatures at a time - paper uses 5)
for j in range(0,125):
    indx = lsh(MSH[:,j*4:j*4+4])

    if np.size(indx) > 0:
        for k in range(0,max(indx)+1):
            temp = np.where(indx==k)
            if np.size(temp)==2:
                buckets.append(tuple(temp[0])) #create pairs of matching signatures
            elif np.size(temp)>2:
                for n in range(0,np.size(temp)):
                    for m in range(n+1,np.size(temp)):
                        buckets.append((temp[0][n],temp[0][m]))

# count the number of matching pairs
db = Counter(elem for elem in buckets)

#view the 10 most common pairs
q=db.most_common(10)
y=list(map(list, buckets))
flat_list = [item for sublist in y for item in sublist]
db2 = Counter(flat_list)


#view the 14 most common windows
q2 = db2.most_common(14)













# plot data, spectrogram, and downsampled spectrogram    
plt.figure(1)
plt.subplot(311)
plt.plot(np.arange(0,st[0].stats.npts/data_downsample,1/data_downsample),st[0].data)
plt.subplot(312)
plt.pcolormesh(t,f,Sxx)
plt.subplot(313)
#plt.colorbar(orientation='horizontal')
plt.pcolormesh(t,np.linspace(0,32,32),spec)

# plot last binary fingerprint image  
fig=plt.figure(2)
plt.pcolormesh(fngr,figure=fig,cmap='binary_r')
plt.colorbar()
plt.show()

#plot last -1/1/0 fingerprint
fig=plt.figure(3)
plt.pcolormesh(haarz,figure=fig,cmap='binary_r')
plt.colorbar()

plt.figure(4,figsize=(8,11))
for i in range(0,len(q)):
    if i > -1: #for testing purposes
        window10s200  = spec[:,q[i][0][0]*specimagelag:q[i][0][0]*specimagelag+specimage]
        
        plt.subplot(len(q),2,i*2+1)
        plt.gca().set_title(str(q[i]))
        plt.pcolormesh(window10s200 )
        plt.tight_layout(h_pad=.5)
        window10s200  = spec[:,q[i][0][1]*specimagelag:q[i][0][1]*specimagelag+specimage]
        plt.subplot(len(q),2,i*2+2)
        plt.pcolormesh(window10s200 )
        plt.tight_layout(h_pad=.5)
plt.figure(5,figsize=(8,15))
for i in range(0,len(q)):
    if i > -1: #for testing purposes
        S0 = (q[i][0][0])*specimagelag*specwindow-(q[i][0][0]-1)*speclag
        S1 = ((q[i][0][0])*specimagelag+specimage)*specwindow-((q[i][0][0]-1)*specimagelag+specimage)*speclag
        plt.subplot(len(q),2,i*2+1)
        plt.gca().set_title(str(q[i]))
        plt.plot(np.linspace(S0/data_downsample,S1/data_downsample,S1-S0),st[0].data[S0:S1],lw=.2)
        S0 = (q[i][0][1])*specimagelag*specwindow-(q[i][0][0]-1)*speclag
        S1 = ((q[i][0][1])*specimagelag+specimage)*specwindow-((q[i][0][1]-1)*specimagelag+specimage)*speclag
        plt.subplot(len(q),2,i*2+2)
        plt.plot(np.linspace(S0/data_downsample,S1/data_downsample,S1-S0),st[0].data[S0:S1],lw=.2)
        plt.tight_layout(h_pad=.5)


plt.figure(7,figsize=(8,15))
for i in range(0,len(q2)):
    if i > -1: #for testing purposes
        S0 = (q2[i][0])*specimagelag*specwindow-(q2[i][0]-1)*speclag
        S1 = ((q2[i][0])*specimagelag+specimage)*specwindow-((q2[i][0]-1)*specimagelag+specimage)*speclag
        plt.subplot(len(q2),2,i+1)
        plt.gca().set_title(str(q2[i]))
        plt.plot(np.linspace(S0/data_downsample,S1/data_downsample,S1-S0),st[0].data[S0:S1],lw=.2)
       # plt.subplots_adjust(top=np.ceil(i/2)*2+1,bottom=np.ceil(i/2))
        plt.tight_layout(h_pad=.5)
plt.show()