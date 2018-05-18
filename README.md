# Network_processes

## Overview 

`Network_processes` is a small library that implements analytical techniques to study dynamical processes on networks. Epidemic spreading processes on complex networks is an important area of current research. The dynamical properties of a disease can be captured using heterogeneous mean field theory (HMF), while the static, final-state properties are best portrayed through the use of generating functions (GFs). This library enables the study of an SIR process using both HMF and GFs. The HMF can be integrated either through RK4 using the `HMF` class, or through stochastic integration using the `STO` class.

![alt text](https://raw.githubusercontent.com/PeterStAndrews/Network_processes/master/example/three_classes.eps)

## Heterogeneous mean field theory

HMF is based on the degree block approximation that partitions nodes according to their degree as well as their state. For a given network, a dynamical process is simulated by a system of ODEs at each degree k in the degree distribution. The integration of each system can be performed either through RK4 or stochastic integration. To access an RK4 integrator for the SIR process, use the `HMF` class. The class `STO` performs Gillespie stochastic simulation for the SIR process. 

## Generating functions

Generating functions can be used to track probability distributions in an orderly manner. In network science, they can be used to find the outbreak size of a disease (amongst many other applications) by taking advantage of a mapping between epidemiology and percolation theory. To predict the final size of an epidemic using GFs for an SIR process, we use the `GFs` class. 

## Author & license 
Copyright (c) 2017-2018, Peter Mann 
Licensed under the GNU General Public Licence v.2.0 <https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html>.
