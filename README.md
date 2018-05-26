# proj13_FF501_SDU
This is a duplicate of README.txt
   			  _____ _____  _    _ 
			 / ____|  __ \| |  | |
			| (___ | |  | | |  | |
			 \___ \| |  | | |  | |	
			 ____) | |__| | |__| |
			|_____/|_____/ \____/  
          
		     University of Southern Denmark
	_____________________________________________________

	  Particle mass estimation using jackknife resampling

				by krlor17
		vibra17, heroe17, jakal17, ishac16
	_____________________________________________________


Acknowledgements:
This project was made possible through the guidance of 

Benjamin Jäger

John Bulava

As we'd otherwise not know what jackknife resampling was in the first place.


--- Purpose of script ---

The matlab script DataJK.m naively estimates the masses of the 11
particles in the data/ directory.

The data are values of the time dependent correlator from a QCD sim.

The script employs Jackknife resampling to estimate statistical error
- without the bother of considering error propagation.


--- Running the script ---
A MATLAB installation is required. Matlab R2017b or newer is recommended.

If possible, running the script in an environment not using a x.org server
makes the process quicker (e.g. SSH connection to a ubuntu server system),
 as MATLAB will then select SOFTWARE OPENGL rendering. (see bottom)
In a x server environment (or equivalent) MATLAB will run JVM by default, 
putting quite a load on low-performance systems.

It is not possible to run the script using the -nojvm option for MATLAB.


-- Linux --
               
To run the script for all particles on linux, simply execute run.sh
$	./run.sh

Make sure that the script has the needed permissions i.e.
$	chmod 770 run.sh



-- Windows --

The easiest approach is to call DataJK(datafile, name) from the MATLAB IDE.
The $name gives the folder and file name under the generated results. I.e. 

>> DataJK pion_ud_correl_all_data_ascii.dat pion

Makes the script save the results for the pion under results/pion/
i.e. the final results in .tex format will be saved as 

" results/pion/pion_final_table.tex "



-- Mac --

Whatever you feel like. 


-- Tips for running script via. ssh --

