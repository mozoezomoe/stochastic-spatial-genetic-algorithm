# -*- coding: utf-8 -*-
"""
Created on Tue Dec 26 23:40:52 2017

@author: mexicore
"""

from Module import *


print "Welcome to Stochastic-Spatial-Genetic-Algorithm"
print "First, I need the rates of the genetic expression (here we'll recommend you some standar values)"

K1=float(raw_input("k1: Rate of binding of the transcription factor to the promoter region (2.55): "))
K_1=float(raw_input("k_-1: Unbinding rate of the promoter-transcription factor pair (150.5): "))
K2=float(raw_input("k2: Rate of DNA transcription (600): "))
K3=float(raw_input("k3: Rate of mRNA translation (50): "))
K4=float(raw_input("k4: mRNA degradation rate (10): "))
K5=float(raw_input("k5: Protein degradation rate (1.05): "))
K6=float(raw_input("k6: Basal transcription rate (15): "))


print "\n1Second, I need some values from the initial state"

TF = float(raw_input("TF: Amount of transcription factor: "))
Pm = float(raw_input("Pm: Free couple promoter-factor ratio (1):"))
Pm_A = float(raw_input("Pm_A Binded couple promoter-factor ratio (0): "))
m = float(raw_input("m: Amount of mRNA produced (0): "))
p = float(raw_input("p: Amount of protein produced (0) "))
t = float(raw_input("t: Initial time in hours (0) "))
tmax = float(raw_input("tmax: Total time in simulation (50) "))

print "\n Our virtual cell size is x[-30,30] and y[-10,10]"
print "\n Knowing that... in which area want to calculate the TF expression? (arbitrary magnitudes)"

xmin = float(raw_input("xmin: "))
xmax = float(raw_input("xmax: "))
ymin = float(raw_input("ymin: "))
ymax = float(raw_input("ymax: "))

print "Wait please, we're printing your results!"

TFr=main(xmin,xmax,ymin,ymax)
(trange,prange) = gill(K1, K_1,K2, K3, K4, K5, K6, TFr[0], Pm, Pm_A, m, p, t, t, xmin,xmax,ymin,ymax)
plot(trange,prange, label='Area []')
legend()