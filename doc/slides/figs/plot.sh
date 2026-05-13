#!/bin/bash

plot_meshtal.py 104 flux-c.pdf runtpe.h5 -v z=49 -thr 1e13:k,2e13:y,2.5e13:orange,3e13:r
plot_meshtal.py 104 flux-2D.pdf runtpe.h5 -v z=49 -l "n Flux $\mathrm{\left[\frac{n}{cm^2 \cdot MJ}\right]}$" -c turbo
plot_meshtal.py 104 flux-1D.pdf runtpe.h5 -v y=0,z=49 --grid -ol "n Flux $\mathrm{\left[\frac{n}{cm^2 \cdot MJ}\right]}$" -lw 0.5

compute_pf.py ACRR.out 134 11,71,45,311 -m PF.pdf -p 4.171 -no 2 -ms 0.9 -fs 1.1
