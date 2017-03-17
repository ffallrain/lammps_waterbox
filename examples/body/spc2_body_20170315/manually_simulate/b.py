#!/usr/bin/python
import sys,os
import numpy as np 
import math

t = 1

inertia_1 = 1.58, 1.05, 1.22,-0.08,-0.32, 0.51 
ix,iy,iz = inertia_1[:3]
inertia_2 = 1.12, 1.12, 1.62, 0.49, 0.08, 0.35

torque =  -5.8419090,  -4.5012881,   2.2438824
tx,ty,tz = torque

ox = 0.5 * tx/ix * t**2 * 4.182* 1e-4
oy = 0.5 * ty/iy * t**2 * 4.182* 1e-4
oz = 0.5 * tz/iz * t**2 * 4.182* 1e-4

print ox,oy,oz
sys.exit()

matrix_x = np.array( ( (1,0,0)                , (0,cos(ox), sin(ox))  , (0,-sin(ox),cos(ox))  ) )
matrix_y = np.array( ( (cos(oy), 0 , -sin(oy)), (0,1,0)               , (sin(oy),0,cos(oy) )  ) )
matrix_z = np.array( ( cos(oz),sin(oz),0)     , (-sin(oz) , cos(oz),0), (0,0,1)                 )


ox,oy,oz = 9.800,9.050,13.590
hx,hy,hz = 10.340,8.760,14.390
Hx,Hy,Hz = 10.100,9.950,13.300

cx = ox*16+hx+Hx
cx = cx/18.
cy = oy*16+hy+Hy
cy = cy/18.
cz = oz*16+hz+Hz
cz = cz/18.

dx = ox-cx
dy = oy-cy
dz = oz-cz


