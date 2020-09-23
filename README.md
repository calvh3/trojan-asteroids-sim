# trojan-asteroids-sim

## Intro

This is a Python program to examine the stability of the L4/5 lagrange points. These points despite being situated at maxima in the effective potential, are stable due to the Coriolis force, which deflects particles around the maxima rather than away from them.  

Jupiter has a large collection of asteroids named trojans and greeks (after the Trojan wars) which librate around its L4 and L5 points respectively. Although Jupiter has by far the most Trojans other planets such as Mars also have asteroids trapped in their lagrange points.

## Overview
The file solver.py numerically solves the reduced three body problem in a rotating reference frame (solver.py).  
The other files all produce an example plot (described at the start of each). The program can be used to map out the stability of the lagrange points in both position and velocity space, as well as investigating the effect of planet mass and orbit radius.
