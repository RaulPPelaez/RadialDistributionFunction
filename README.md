# Raul P. Pelaez 2017 Radial Distribution Function  

## NAME   
rdf -  Computes the Radial Distribution Function (RDF) of a group of positions in a file, averages it for all snapshots in the file.  
rdf can compute in 4 modes:  
	* GPU  
	  - Nbody  
	  - Neighbour list  
	* CPU  
	  - Nbody  
	  - Neighbour list  
	  
rdf will choose a between GPU/CPU according to the number of particles (unless specified with -device) and will choose NBody/Neighbour list according to the number of particles and the factor L/rcut.  
## COMPILE WITH  

```
$ make
```
You may have to change the Makefile to adequate it to the CUDA target architechture, currently set to -arch=sm_35

Use:
```
$ make test 
```

To compile and run several test using random numbers, the resulting rdf will be compared between CPU and GPU implementations. Which should be numerically identical.

## SYNOPSYS  

rdf [OPTIONS]... [FILE]...  

## DESCRIPTION  
   Compute the Radial Distribution Function.  
   
   With no FILE, or when file is - or not specified, reads from standard input.  

   Required options:  

   -N
       Number of particles, all snapshots must have the same number of particles  

   -L, -Lx [lx] -Ly [ly]  -Lz[lz]  
       Box size, positions will be folded to a box of size Lx Ly Lz. -L will make Lx= Ly = Lz = L  

   -rcut  
       Maximum distance in the rdf, distances greater than rcut will be ignored.  
   
   -nbins  
       Number of bins in the position pair histogram (binSize = rcut/nbins)  

   -Nsnapshots   
       Number of snapshots in the file, a snapshot must be separated from the next with a single line  

   -dim [=3D]  
       Dimensionality of the input positions. Affects how the histogram is normalized to compute the rdf.  
      Can be 3D, 2D or q2D (treat as 3D, but normalize as 2D)
	  
   -device [=auto]  
       Switch between GPU/CPU implementations of the algorithm. By default rdf chooses the best according to N  
	
   -outputDecimals [=5]  
	   Number of decimals in the output file, set through cout<<setprecision()  
   -fixBIAS
	   This will weight the distance of a pair in a bin according to the position inside the bin (instead of weighting all distances as 1).  

   
   
 ## FILE FORMAT   
   The file must have at least "dim" columns (the rest will be ignored) and each snapshot (including the first)  
   must be preceded by a line (no matter the content as long as it is a single line). See example.  


## EXAMPLES:  

```
---pos.dat----
#
1 2 3 4 5
6 7 8 9
10 11 12 13 14
#
15 16 17 18
19 20 21
22 23 24
------

$ cat pos.dat | rdf -N 3 -Nsnapshots 2 -L 25 -nbins 10 -rcut 2 > rdf.dat
```
rdf will take the file as 2 snapshots with the positions of 3 particles in 3D.

