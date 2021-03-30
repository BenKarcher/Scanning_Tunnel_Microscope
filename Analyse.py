#sorry about lack of comments was hacked together in a hurry will try to upload a better version once I add this functionality to my math library

import matplotlib
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from scipy.integrate import quad

path="/home/ben/Documents/usbstick/good pics/"
fileExtension=".txt"

re=np.real
im=np.imag
ab=np.absolute

def read2d(fName):
    data=[]
    with open (path+fName+fileExtension, "r") as file:
        for line in file:
            line=line.split("\t")
            line=[float(s) for s in line]
            data.append(line)
    return np.array(data)

def heatmap2d(arr: np.ndarray):
    lower=np.percentile(arr,2)
    upper=np.percentile(arr,98)
    plt.imshow(arr, cmap='plasma',vmin=lower,vmax=upper)
    plt.colorbar(label="volt")
    plt.show()
    
def plot(arr: np.ndarray,length=1,zLabel="",percentile=2,title="",frequency=False):
    lower=np.percentile(arr,percentile)
    upper=np.percentile(arr,100-percentile)
    plt.imshow(arr, cmap='plasma',vmin=lower,vmax=upper,extent=([-1,1,-1,1] if frequency else [0,length,0,length]))
    plt.colorbar(label=zLabel)
    if title != "":
        plt.title(title)
    plt.show()
    
def surface(func,shape,xmin=-1,xmax=1,ymin=-1,ymax=1,args=[]):
    result=np.zeros(shape)
    for i in range(shape[0]):
        x=xmin+(xmax-xmin)*(i/shape[0])
        for j in range(shape[1]):
            y=ymin+(ymax-ymin)*(j/shape[1])
            result[i,j]=func(x,y, *args)
    return result
def linex(shape):
    result=np.ones(shape)
    a=np.floor(shape[0]/2)
    b=np.ceil(shape[0]/2)
    for x in range(shape[1]):  
        result[int(a),x]=0;
        result[int(b),x]=0;
    c=np.floor(shape[1]/2)
    d=np.ceil(shape[1]/2)
    result[int(a),int(c)]=1;
    result[int(a),int(d)]=1;
    result[int(b),int(c)]=1;
    result[int(b),int(d)]=1;
    return result
def liney(shape):
    result=np.ones(shape)
    a=np.floor(shape[1]/2)
    b=np.ceil(shape[1]/2)
    for y in range(shape[0]):  
        result[y,int(a)]=0;
        result[y,int(b)]=0;
    c=np.floor(shape[0]/2)
    d=np.ceil(shape[0]/2)
    result[int(a),int(c)]=1;
    result[int(a),int(d)]=1;
    result[int(b),int(c)]=1;
    result[int(b),int(d)]=1;
    return result
def simpleDenoise(arr: np.ndarray):
    for x in range(1,arr.shape[1]-1):
        for y in range(1,arr.shape[0]-1):
            avg=(arr[x-1,y-1]+arr[x-1,y]+arr[x-1,y+1]+arr[x,y-1]+arr[x,y+1]+arr[x+1,y-1]+arr[x+1,y]+arr[x+1,y+1])/8
            if ab(avg-arr[x,y])>0.05*ab(avg):
                arr[x,y]=avg
def integrate(arr: np.ndarray):
    result=np.zeros(arr.shape)
    a=0
    for y in range(arr.shape[0]):
        result[y,0]=arr[y,0]
        a+=result[y,0]
    for x in range(1,arr.shape[1]):
        for y in range(arr.shape[0]):
            result[y,x]=result[y,x-1]+arr[y,x]
    b=0
    for y in range(arr.shape[0]):
        b+=result[y,arr.shape[1]-1]
    angle=(a-b)/arr.shape[0]/arr.shape[1]
    for x in range(arr.shape[1]):
        for y in range(arr.shape[0]):
            result[y,x]+=angle*x
    return result
    
        
def gauss(x,y,w=160):
    return np.exp(-w*(x**2+y**2)**4)

def fft(arr: np.ndarray):
    return np.fft.fftshift(np.fft.fft2(arr))
def ifft(arr: np.ndarray):
    return np.fft.ifft2(np.fft.ifftshift(arr))

hbarC=197.3269804#eV nm
electronMass=511000#eV
work=4
constant=2*np.sqrt(2*electronMass*work)/hbarC
print(constant)
def k(I,D):
    return I*np.exp(constant*D)


Current=1e9*read2d("HDCurrent")
Height=1e9*read2d("HDHeight")
CurrentBack=1e9*read2d("HDCurrentBack")
HeightBack=1e9*read2d("HDHeightBack")

Graphit=1e9*read2d("graphit")
GraphitBack=read2d("graphitBack")


plot(Height,length=200,zLabel="[nm]",title="Gold kugeln Z-position [nm]")

denoised=Height.copy()
simpleDenoise(denoised)

plot(denoised,length=200,zLabel="[nm]",title="Reduziertes rauschen Z-position [nm]")

Fourier=fft(Height)
Fourier2=fft(denoised)

plot(ab(Fourier),frequency=True,zLabel="[nm]",title="FFT der Z-position [Nyquist frequenz]")
plot(ab(Fourier2),frequency=True,zLabel="[nm]",title="FFT der Z-position ohne Rauschen [Nyquist frequenz]")

plot(HeightBack,length=200,zLabel="[nm]",title="Z-position Rück­weg[nm]")
plot((Height+HeightBack)/2,length=200,zLabel="[nm]",title="Z-position beide wege Gemittelt[nm]")

plot(Current,length=200,zLabel="[nA]",title="Gold kugeln Strom [nm]")
inter=integrate(Current-np.average(Current)/2)
plot(inter,length=200,zLabel="[nm*nA]",title="Strom Integriert [nm]")
Combined=np.sqrt(ab(Current*Height))
plot(Combined,length=200,zLabel="[sqrt(nA*nm)]",title="Gold kugeln Kombiniert [nm]")

Fourier=Fourier*linex(Combined.shape)
plot(ab(Fourier),frequency=True,zLabel="[nm]",title="FFT der Z-position nach Achsenfilter [Nyquist frequenz]")
Fixed=ifft(Fourier)
plot(re(Fixed),length=200,zLabel="[nm]",title="Z-position nach Achsenfilter [nm]")

Combined=np.sqrt(ab(Current*Fixed))
plot(Combined,length=200,zLabel="[sqrt(nA*nm)]",title="Kombiniert nach Achsenfilter [nm]")

Fourier=fft(Combined)
plot(ab(Fourier),frequency=True,zLabel="[nm]",title="FFT Kombiniert nach Achsenfilter [Nyquist frequenz]")

Blur=surface(gauss, Combined.shape)
plot(Blur,frequency=True,zLabel="",title="tiefpass Filter [Nyquist frequenz]")

Fourier2=Fourier*Blur
plot(ab(Fourier2),frequency=True,zLabel="[nm]",title="FFT Kombiniert nach Achsen+Tiefpass-filter [Nyquist frequenz]")
plot(re(ifft(Fourier2)),length=200,zLabel="[sqrt(nA*nm)]",title="Kombiniert nach Achsen+Tiefpass-filter [nm]")



plot(Graphit,length=1.03,zLabel="[nA]",title="Graphit Strom [nm]")
Fourier=fft(Graphit)
plot(ab(Fourier),frequency=True,zLabel="[nA]",title="FFT von Graphit Strom [Nyquist frequenz]")
Fourier=Fourier*liney(Graphit.shape)*surface(gauss,Graphit.shape)
plot(ab(Fourier),frequency=True,zLabel="[nA]",title="FFT Strom nach Achsen+Tiefpass-filter [Nyquist frequenz]")
plot(re(ifft(Fourier)),length=1.03,zLabel="[nA]",title="Graphit Strom nach Achsen+Tiefpass-filter [nm]")
