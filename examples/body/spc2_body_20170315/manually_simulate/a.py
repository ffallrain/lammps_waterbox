#!/usr/bin/python
import sys,os


ox,oy,oz = 9.800,9.050,13.590
hx,hy,hz = 10.340,8.760,14.390
Hx,Hy,Hz = 10.100,9.950,13.300
cx = ox*16+hx+Hx
cx = cx/18.
cy = oy*16+hy+Hy
cy = cy/18.
cz = oz*16+hz+Hz
cz = cz/18.

fx = -5.8419090
fy = -4.5012881
fz = 2.2438824

print cx
print cy
print cz
#shift_x = 0.5 * a* t**2
#shift_x = 0.5 * (fx/mass)* timestep **2
shift_x = 0.5 * (fx/18)* 1 **2 * 4.184 * 10**2 * 1e-4
shift_y = 0.5 * (fy/18)* 1 **2 * 4.184 * 10**2 * 1e-4
shift_z = 0.5 * (fz/18)* 1 **2 * 4.184 * 10**2 * 1e-4

print "%f"%shift_x
print "%f"%shift_y
print "%f"%shift_z

ox,oy,oz = 9.80097,9.04055,13.5868 
hx,hy,hz = 10.3235,8.77227,14.4058 
Hx,Hy,Hz = 9.97989,9.99433,13.3809 

cx2= ox*16+hx+Hx
cx2= cx2/18.
cy2= oy*16+hy+Hy
cy2= cy2/18.
cz2= oz*16+hz+Hz
cz2= cz2/18.

print ">> Step1: %16.8f %16.8f %16.8f "%(cx,cy,cz)
print ">> Shift: %16.8f %16.8f %16.8f "%(shift_x,shift_y,shift_z)
print ">> Step2: %16.8f %16.8f %16.8f "%(cx2,cy2,cz2)
