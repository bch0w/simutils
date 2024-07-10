# ISOSURFACE README

Short instruction set to make an Isosurface figure as in Chow et al. 2022b 
Figure 6A using the ParaView GUI. 

1) File -> Open the figure and `apply` 
2) Click `Contour` button (looks like a half sphere, next to calculator)
3) Set 'Isosurface' to desired value
4) Apply Filter -> Alphabetical -> `Python Calculator` to contour
5) In calculator 'Expression': inputs[0].Points[:,2] * 1E-3
	This will set the value of each point to its Z (depth) value
6) Set Array Name to 'Z\_values'
7) Set Coloring to 'Z\_values'
8) Scroll down to 'Transforming' and set third column Scale to something like 4
	This will exagerrate the Z dimension 
