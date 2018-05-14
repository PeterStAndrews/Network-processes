# Network-processes

## Overview 

`Network-processes` is a small library that implements analytical techniques to study dynamical processes on networks. Epidemic spreading processes on complex networks is an important area of current research. The dynamical properties of a disease can be captured using heterogeneous mean field theory (HMF), while the static, final-state properties are best portrayed through the use of generating functions (GFs). This library enables the study of an SIR process using both HMF and GFs. 

## Heterogeneous mean field theory

HMF is based on the degree block approximation that partitions nodes according to their degree as well as their state. For a given network, a dynamical process is simulated by a system of ODEs at each degree k in the degree distribution. To access this for the SIR process, use the `HMF` class. 

## Generating functions

Generating functions can be used to track probability distributions in an orderly manner. In network science, they can be used to find the outbreak size of a disease (amongst many other applications) by taking advantage of a mapping between epidemiology and percolation theory. To predict the final size of an epidemic using GFs for an SIR process, we use the `GFs` class. 

## Author & license 
Copyright (c) 2017-2018, Peter Mann 
Licensed under the GNU General Public Licence v.2.0 <https://www.gnu.org/licenses/old-licenses/gpl-2.0.en.html>.
