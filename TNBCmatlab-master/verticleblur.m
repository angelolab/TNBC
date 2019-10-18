function obj=verticleblur(ifile)

%Distance between true and shadow image in pixels
ts=15;

%Estimated probe size in um
ps=1;
 
%FOV size in um
fs=800
 
%image size in pixels
is=2048;
 
%probesize in pixels
psp=round(ps/fs*is);
 
%Weight for substraction
wt=-0.3;
 
%Generate Kernel for convolution
filterKernel=zeros(2*ts+1+psp*2,2*ts+1+psp*2);
gw=fspecial('gaussian',psp*2,psp)*wt;
[gr,gc]=size(gw) 

 
%Offset in pixels, this is to approximate the angle of shadowing.  If 90
%degrees these values are zero.  Values here are for a 75 degree angle 
%(15 degrees relative to the normal).

Angle=15 
r_offset=round(sind(Angle)*ts);
c_offset=round(cosd(Angle)*ts);
 
filterKernel(1+r_offset:gr+r_offset,(ts+psp+1+c_offset-gc/2):(ts+psp+gc/2+c_offset))=gw;
 
filterKernel(ts+psp+1,ts+psp+1)=1;
  
obj=imfilter(ifile,filterKernel);
 
 
 
end