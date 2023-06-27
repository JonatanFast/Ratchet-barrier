# Ratchet-barrier
Analysis for paper: Geometric symmetry breaking and nonlinearity can increase thermoelectric power
Preprint manuscript available at: https://doi.org/10.48550/arXiv.2304.01616

# There are 4 main scripts (IPYNB format, open in VS code) + 1 normal python file containing all my functions:
1. 	"Calculate_transmission.IPYNB" is used to create an array containg the transmission probability of a 
	ramp shaped barrier based on the schrödinger equation (Derives it in the manuscript SI). 
	Typically you do not need to run this script as I have calculated and saved the transmission probability 
	for a number of different settings, found as JSON files in "Ratchet Barrier/Model/Calculated transmission"
2.	"Ratchet_no_heating.IPYNB" calculates a large 2D array of the SSE between experiment without heating and
	L-B model for a series of A and µ, used to extract the range of A and µ that can give a good fit. 
	File from 1 is used in this script.
	You typically don't need to run this script, I have saved the results as JSON file in "Ratchet Barrier/Model/Fit A and mu to dark curve"
3. 	"Fit_DeltaT" does the fitting between experiment and model to determine Delta T as a function of heating voltage. 
	It uses files from 1 and 2. you typically don't need to run this file as I have saved the results in 
	"Ratchet Barrier/Model/Fit A and mu to dark curve" as "Fit_T_L92nm_Utop_340meV_data.JSON".
4. 	"Ratchet_analysis" Is most likely the only script you may be using. It requires the result file from 3, 
	as well as from 1 if you want to compare experiment and L-B model (choose file with 40 mV range). 
5. 	"ratchet functions" contains many different functions used in the scripts, you hopefully won't be doing anything
	to alter them but need to know that in the way it is currently being imported, the file must be in the same folder as the script you are running


