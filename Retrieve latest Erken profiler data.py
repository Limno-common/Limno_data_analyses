# -*- coding: utf-8 -*-
"""
Created on Tue Aug 23 11:33:00 2022

@author: Shuqi Lin
"""
#%% Import libraries required
import urllib.request 
import os
import datetime as dt
from xml.etree.ElementTree import parse
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import time

#%% Define the functions 
def get_observation_url(url_source,url_tool,time):
    url_format= 'xml'
    url_time=str(time)+'T00:00:00'
    url=url_source+url_tool+'&format='+url_format+'&mode=since-time&p1='+url_time+'&p2='
    return url

    time=[]
    for i in root[1]:
        time.append(i.attrib['time'])
    time=np.array(time)
    SWR=[]
    AirP=[]
    AirT=[]
    Hum=[]
    WS=[]
    Prec=[]
    for i in range(len(root[1])):    
        SWR.append(float(eliminate_nan(root[1][i][0].text)))
        AirT.append(float(eliminate_nan(root[1][i][3].text)))
        Hum.append(float(eliminate_nan(root[1][i][4].text)))
        WS.append(float(eliminate_nan(root[1][i][9].text)))
        Prec.append(float(eliminate_nan(root[1][i][17].text)))
        AirP.append(float(eliminate_nan(root[1][i][18].text)))
    return pd.DataFrame({'u10':WS,'AirT':AirT,'SWR':SWR,'AirP':AirP,'Prec':Prec,'Hum':Hum},index=time)

def eliminate_nan(x):
    if x=='"NAN"':
        return np.nan
    else:
        return x

def get_temp_profile(root):
    ## Get time series
    time=[]
    for i in root[1]:
        time.append(i.attrib['time'])
    time=np.array(time)

    ## Get vertical position
    depth=[]
    for i in root[0][1]:
        d=i.attrib['name']
        if len(d.split('WTemp')[0])!=0:
            pos='RefTemp'
        else: 
            pos=d.split('_')[0].split('WTemp')[1]+'.'+d.split('_')[1]
            pos=float(pos)
        depth.append(pos)

    temp=pd.DataFrame(index=time,columns=depth)

    ## Get temperature profile
    for i in range(len(root[1])):
        profile=[]
        for j in root[1][i]:
            profile.append(j.text)
        profile=np.array(profile)
        temp.iloc[i]=profile
    return temp

def get_YSI_var(root):
    ## Get time series
    time=[]
    for i in root[1]:
        time.append(i.attrib['time'])
    time=np.array(time)
    ## Get temp, O2, Chl
    temp=[]
    Chl=[]
    for i in range(len(root[1])):
        temp.append(float(root[1][i][0].text))
        Chl.append(float(root[1][i][9].text))
    return pd.DataFrame({'Time':time,'Temp':temp,'Chl':Chl})

def get_DO(root):
    time=[]
    for i in root[1]:
        time.append(i.attrib['time'])
    time=np.array(time)
    Sur_DO=[]
    Bot_DO=[]
    for i in range(len(root[1])):
        Sur_DO.append(float(root[1][i][17].text))
        Bot_DO.append(float(root[1][i][18].text))
    return pd.DataFrame({'Time':time,'Sur_DO':Sur_DO,
                         'Bot_DO':Bot_DO})
    

#%% Extract daily temperature observations from Erken logger 
now=dt.datetime.now()
start_time=now.date() - dt.timedelta(days=3) 
url_source='http://130.238.87.115:8080/vilandatalogger/?command=DataQuery&uri=Server:'
url_tool='Erken_Underwater_Temp_Sys.Output60min'
url=get_observation_url(url_source,url_tool,start_time)
tree=parse(urllib.request.urlopen(url))
root = tree.getroot()
temp=get_temp_profile(root)
temp=temp.astype('float64')
temp[0.5]=np.mean(temp[['RefTemp','RefTemp']],axis=1)
temp=temp[np.sort(temp.columns[2:])]
temp.plot(y=0.5,figsize=(16,4)) # Only the recent 3 months data are available
temp.to_csv('Erken_temperature_logger.obs',sep='\t')

#%% Format temp data into GOTM readible file and store it in the forecast folder
temp=temp.reset_index().rename(columns={'index':'Time'})
temp['Time']=pd.to_datetime(temp['Time'])
temp['Date']=temp['Time'].apply(lambda d:d.date())
temp=temp.groupby('Date').mean().reset_index()
#%%
temp=pd.melt(temp,id_vars =['Date'],var_name ='Depth', value_name ='temp')
temp=temp.sort_values(by=['Date','Depth'])
temp['Depth']=temp['Depth']*-1
temp['temp']=np.round(temp['temp'],2)
midnight_time=dt.datetime.min.time()
temp['Date']=temp['Date'].apply(lambda d:dt.datetime.combine(d,midnight_time))

#%% Get Surface O2 and temp from YSI sonde
url_tool='YSI_ErkenWiFi.SondeHourly'
start_time=now.date() - dt.timedelta(days=100) 
url=get_observation_url(url_source,url_tool,start_time)
tree=parse(urllib.request.urlopen(url))
root = tree.getroot()
YSI=get_YSI_var(root)
YSI['Time']=pd.to_datetime(YSI['Time'])
YSI.to_csv('YSI data_todate.csv',index=False,sep='\t')

#%% Get Bottom DO from YSI sonde
url_tool='FloatWind_ErkenWiFi.Output60min'
start_time=now.date() - dt.timedelta(days=200) 
url=get_observation_url(url_source,url_tool,start_time)
tree=parse(urllib.request.urlopen(url))
root = tree.getroot()
DO=get_DO(root)
DO['Time']=pd.to_datetime(DO['Time'])
DO.to_csv('YSI DO_todate.csv',index=False,sep='\t')
