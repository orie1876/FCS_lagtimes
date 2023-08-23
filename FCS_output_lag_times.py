import numpy as np
import csv
import multipletau               #Install this packaged
import matplotlib.pyplot as plt
import scipy.optimize as opt
import pandas as pd
import math


########################################### Run all scripts ############################################

k_green=2.6783865980914547 
w_green=2.7584195881055327e-07 
k_red=5.037835549184945 
w_red=1.907238974194528e-07 

number_of_files=5
filename="FCS"

pathlist=[]

rootpath='/Users/Mathew/Documents/Current analysis/20230809_PDZ_FCS/Run_1_PDZ2/'

pathlist.append(r"/Users/Mathew/Documents/Current analysis/20230809_PDZ_FCS/Run_1_PDZ2/Peptide_only_100nM/")
pathlist.append(r"/Users/Mathew/Documents/Current analysis/20230809_PDZ_FCS/Run_1_PDZ2/20uM/")

def runall(path):
    global green
    global red
    

   
                                                 
    
    green = []                                                                       # This is where the red and green data will be saved. Makes an array. Data exists as two columns of 10 us bursts (green and red)
    red = []               
    
    toloadpath=path
    toloadpath1=path+filename
    print(toloadpath1)
    with open(toloadpath1) as csvDataFile:                                                # Opens the file as a CSV
       csvReader = csv.reader(csvDataFile,delimiter='\t')                           # Assigns the loaded CSV file to csvReader. 
       for row in csvReader:
           green.append(row[0])                                                     # For every row in in csvReader, the values are apended to green and red.         
           red.append(row[1])
     
    for i in range(number_of_files-1):
        j=i+2
        if j<10 :
            toload=toloadpath+filename+"_0"+str(j)
            print(toload)
        else:
            toload=toloadpath+filename+"_"+str(j)
            print(toload)
       
        with open(toload) as csvDataFile:                                                # Opens the file as a CSV
            csvReader = csv.reader(csvDataFile,delimiter='\t')                           # Assigns the loaded CSV file to csvReader. 
            for row in csvReader:
                green.append(row[0])                                                     # For every row in in csvReader, the values are apended to green and red.         
                red.append(row[1])
    
        
    autocorrelate()                                                     # Autocorrelated the sample - comment out if you just want to look at the dye.
                                                                           

########################################### Autocorrelate ############################################
    
def autocorrelate():
   global c
   global new_c
   global new_d
  


   
        
   x=np.array(green,dtype=float)                                                    # Convert the csv columns to float values - Why does it not know?
   c=multipletau.autocorrelate(x,m=16, normalize=True,deltat=1e-5)                  # Correlate the data using the multipletau python. The deltat is the bin width. 
   new_c=np.delete(c, (0), axis=0)                                                  # This deletes the first row which is at t = 0. 
   
       
   
   x=np.array(red,dtype=float)                                                    # Convert the csv columns to float values - Why does it not know?
   d=multipletau.autocorrelate(x,m=16, normalize=True,deltat=1e-5)                  # Correlate the data using the multipletau python. The deltat is the bin width. 
   new_d=np.delete(d, (0), axis=0) 





def fungreen(x,n,td):
    k=k_green                  # This value is the beam waist measured from the dye only.
    return (1/n)*((1+x/td)**(-1))*((1+x/(k**2*td))**(-0.5))

########################################### Fit with unknown diffusion coefficient, but known beam waist ############################################
def fitgreen():
    xdata=new_c[:, 0]  
    ydata=new_c[:,1]
    guesses=np.array([20,6e-5])
    (n_, td_), _ = opt.curve_fit(fungreen, xdata, ydata,guesses)
    params= opt.curve_fit(fungreen, xdata, ydata,guesses)
    
   
    
    y_fit = fungreen(xdata, n_, td_)
    
    
   
       # plotting
    fig = plt.figure(figsize=(10, 8))
    fig.canvas.set_window_title('FCS Curve')

    # autocorrelation
    ax1 = fig.add_subplot(211)
    ax1.plot(xdata,ydata, "--k",
         color="grey", label="correlate (numpy)")
    ax1.set_xlabel("Time lag (s)")
    ax1.set_ylabel("Autocorrelation")
    ax1.set_xscale('log')
    ax1.set_xlim(1e-5,10)
    ax1.set_ylim(0,max(new_c[:,1]))
    ax1.plot(xdata, y_fit, '-',color='green')
    
    
    print(("Green_N = %r \r") %params[0][0])
    print(("Green_td = %r \r") %params[0][1])

    
    Diff=(w_green)**2/(4*params[0][1])
    
    print(("Green_D = %r \r") %Diff)
    
    Rh=1.381e-23*298/(6*3.141*8.9e-4*Diff)

    print(("Green_r = %r \r") %Rh)
    
    conf_volume=(math.pi**(3/2)*k_green*w_green**3)*1000 
    concentration=(params[0][0]/6.022e23)/conf_volume
    print(("Green_conc = %r \r") %concentration)
    
    return params[0][1]

def funred(x,n,td):
    k=k_red                  # This value is the beam waist measured from the dye only.
    return (1/n)*((1+x/td)**(-1))*((1+x/(k**2*td))**(-0.5))

########################################### Fit with unknown diffusion coefficient, but known beam waist ############################################
def fitred():
    xdata=new_d[:, 0]  
    ydata=new_d[:,1]
    guesses=np.array([25,5e-5])
    (n_, td_), _ = opt.curve_fit(funred, xdata, ydata,guesses)
    params= opt.curve_fit(funred, xdata, ydata,guesses)
    
   
    
    y_fit = funred(xdata, n_, td_)
    
    
   
       # plotting
    fig = plt.figure(figsize=(10, 8))
    fig.canvas.set_window_title('FCS Curve')

    # autocorrelation
    ax1 = fig.add_subplot(211)
    ax1.plot(xdata,ydata, "--k",
         color="grey", label="correlate (numpy)")
    ax1.set_xlabel("Time lag (s)")
    ax1.set_ylabel("Autocorrelation")
    ax1.set_xscale('log')
    ax1.set_xlim(1e-5,10)
    ax1.set_ylim(0,max(new_d[:,1]))
    ax1.plot(xdata, y_fit, '-',color='red')
    
    
    print(("Red_N = %r \r") %params[0][0])
    print(("Red_td = %r \r") %params[0][1])

    
    Diff=(w_green)**2/(4*params[0][1])
    
    print(("Red_D = %r \r") %Diff)

    Rh=1.381e-23*298/(6*3.141*8.9e-4*Diff)

    print(("Red_r = %r \r") %Rh)
    conf_volume=(math.pi**(3/2)*k_red*w_red**3)*1000 
    concentration=(params[0][0]/6.022e23)/conf_volume
    print(("Red_conc = %r \r") %concentration)
    
    return params[0][1]


Output_all = pd.DataFrame(columns=['Path','Td_green','Td_red'])

for path in pathlist:

    runall(path)
    green_lag=fitgreen()
    red_lag=fitred()
    

    
    
    Output_all= Output_all.append({'Path':path,'Td_green':green_lag,'Td_red':red_lag},ignore_index=True)
    
    Output_all.to_csv(rootpath + 'Results.csv', sep = '\t')
    # def rebin(a,newshape):
    #     '''Rebin an array to a new shape.
    #     '''
    #     assert len(a) == len(newshape)

    #     slices = [ slice(0,old, float(old)/new) for old,new in zip(a.shape,newshape) ]
    #     coordinates = mgrid[slices]
    #     indices = coordinates.astype('i')   #choose the biggest smaller integer index
    #     return a[tuple(indices)]
    
    # new_width=len(green)/1000
    
    
    
    # green_plot=rebin(green,new_width)
    
    
            

