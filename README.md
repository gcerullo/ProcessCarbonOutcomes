# ProcessCarbonOutcomes
Process all carbon outcomes


#Rscripts 

#1. CalculatePlantationCarbonVals 
Fits models based on SSB inventory data to determine the carbon per year in 12- and 6-year albizia and eucalyptus plantations, respectively. 

#2. CalculateALlHabCarbonVals 
#This code calculates carbon (above and belowground) of each habitat type through time 

#Code notes:
Plantation ACD is calculate is from SSB inventory data 
Primary, 1L and R is from Philipson and uses edited code originally written by 
Philipson 2020; Science. I adjust Philipson's code to allow estimates for 60 years (they
estimate out to 30/35 years; I assume that slope and intercept of their models stay the same 
and that 1L and R values plateu once they reach primary forest values)
I estimate ACD in twice-logged forest 
This code also:
1. Adds different harvest delays 
2 Incorporates belowground carbon and necromass