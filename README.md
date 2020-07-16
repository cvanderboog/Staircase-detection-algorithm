# Staircase-detection-algorithm
A detection algorithm (python) to detect oceanic thermohaline staircases from a ice tethered profilers and Argo floats.

____________________________________________
To run the algorithm: 
from staircase_detection_algorithm import *

p  = 'pressure profiles'	# pressure profiles with resolution of 1 dbar
ct = 'temperature profiles'  # conservative temperature profile (deg. C)
sa = 'salinity profiles'   #absolute salinity profile (g / kg)

c1 = 0.0005 # mixed layer gradient with gradient detection method
c2 = 0.005  # density difference of mixed layer
c3 = 200    # averaging window to obtain background profiles
c4 = 30     # maximum gradient layer height

gl, ml, masks = get_mixed_layers(p,ct,sa,c1,c2,c3,c4)


____________________________________________
To make your own global dataset:
%run main.py

____________________________________________
The algorithm computes the following variables:

n 			    : number of profiles in the file
pressure		: pressure levels of the profiles (max. = 2000 dbar)
mixed_layers	: size of array to save the properties of each mixed layer (max. 200 mixed layers per profile)
gradient_layers	: size of array to save the properties of each gradient layer (max. 200 gradient layers per profile)  

General variables in each file:
prof			: number of profile, counted per float
FloatID 		: identification number of each float
juld			: Julian date of observation. (Days since 1950-1-1)
lon				: longitude of observation
lat				: latitude of observation
ct				: conservative temperature  (size: ‘pressure’, ’n’)
sa				: absolute salinity profile (size: ‘pressure’, ’n’)

Mixed layer variables: 
ml_T			: average conservative temperature (deg. C)
ml_S			: average absolute salinity (g/kg)
ml_r			: average density anomaly (sigma, reference pressure 1000) 
ml_p			: average pressure 
ml_h			: vertical extent (dbar)
ml_R			: average density ratio
ml_Tu			: average Turner angle 
ml_dT			: conservative temperature variations in mixed layer (deg. C)
ml_dS			: absolute salinity variations in mixed layer (g/kg)
ml_dr			: density anomaly variations in mixed layer (sigma, reference pressure 1000) 


mask_ml_sf_layer	: A mask of all mixed layers (size: ‘mixed_layers’, ’n’) that have been identified as 
				the ‘salt-finger regime’. The value of each sequence indicates the number of 
				mixed layers in that particular sequence. This mask is applicable to all mixed-
				layer variables.

mask_ml_sf		: A mask of all mixed layers (size: ‘pressure’, ’n’) that have been identified as 
				the ‘salt-finger regime’. The value of each sequence indicates the number of 
				mixed layers in that particular sequence. This mask is applicable to the 	
				 variables ‘ct’ and ‘sa’.

mask_ml_dc_layer	: A mask of all mixed layers (size: ‘mixed_layers’, ’n’) that have been identified as 
				the ‘diffusive-convection regime’. The value of each sequence indicates the 
				number of mixed layers in that particular sequence. This mask is applicable to 
				all mixed-layer variables.

mask_ml_dc		: A mask of all mixed layers (size: ‘pressure’, ’n’) that have been identified as 
				the ‘diffusive-convection regime’. The value of each sequence indicates the 
				number of mixed layers in that particular sequence. This mask is applicable to 
				the variables ‘ct’ and ‘sa’.

Gradient layer variables:
gl_dT			: conservative temperature difference across gradient layer (deg. C)
gl_dS			: absolute salinity difference across gradient layer (g/kg)
gl_r			: density anomaly difference across gradient layer (sigma, reference pressure 
				  1000) 
gl_h			: vertical extent (dbar)
gl_R			: average density ratio 
gl_Tu			: average Turner angle 
gl_dTdz			: average vertical temperature gradient within gradient layer
gl_dSdz			: average vertical salinity gradient within gradient layer


mask_gl_sf_layer: A mask of all gradient layers (size: ‘gradient_layers’, ’n’) that have been 
			 	identified as the ‘salt-finger regime’. The value of each sequence indicates the 
				number of mixed layers in that particular sequence. This mask is applicable to 
				all gradient-layer variables.

mask_gl_sf		: A mask of all gradient layers (size: ‘pressure’, ’n’) that have been identified as 
				the ‘salt-finger regime’. The value of each sequence indicates the number of 
				mixed layers in that particular sequence. This mask is applicable to the 				  variables ‘ct’ and ‘sa’.

mask_gl_dc_layer: A mask of all gradient layers (size: ‘gradient_layers’, ’n’) that have been 
				identified as the ‘diffusive-convection regime’. The value of each sequence 
				indicates the number of gradient layers in that particular sequence. This mask 
				is applicable to all gradient-layer variables.

mask_gl_dc		: A mask of all gradient layers (size: ‘pressure’, ’n’) that have been identified as 
				the ‘diffusive-convection regime’. The value of each sequence indicates the 				  
				number of gradient layers in that particular sequence. This mask is applicable 
				to the variables ‘ct’ and ‘sa’.

Sequences of gradient layers are saved to perform quality controls:
mask_qc			: A mask which specifies which requirements are met. This is a binary count, 
				where 
 				  +2 		temperature variation in GL < temperature variation ML
				  +4 		salinity variation in GL < salinity variation ML
				  +8 		density variation in GL < density variation ML
				  +16 	gradient layer height < mixed layer height
				  +32 	maximum gradient layer height is less than specified
				  +64 	no inversions in gradient layer
				If all requirements are met, this adds up to 127. In this case, the gradient layer 
				is considered as a potential gradient layer of a staircase. These gradient layers 
				are used to identify sequences. 

mask_seq_sf_count	: the total number of layers per sequence in the salt-finger regime
mask_seq_sf_layer	: the layer number within sequence in the salt-finger regime (increase with 
				depth)
mask_seq_sf_num	: the number of the sequence per profile (increase with depth) in the salt-finger 
				regime

mask_seq_dc_count	: the total number of layers per sequence in the diffusive-convective regime
mask_seq_dc_layer	: the layer number within sequence in the diffusive-convective regime (increase 
				  with depth)
mask_seq_dc_num	: the number of the sequence per profile (increase with depth) in the diffusive-
				  convective regime


Acknowledgements:
The Ice-Tethered Profiler data were collected and made available by the Ice-Tethered Profiler Program (Krishfield et al., 2008; Toole et al., 2011) based at Woods Hole Oceanogaphic Institution (http://www.whoi.edu/itp). The Argo data were collected and made freely available by the International Argo Program and the national programs that contribute to it (http://www.argo.ucsd.edu, http://argo.jcommops.org). The Argo Program is part of the Global Ocean Observing System (Argo, 2000). 
