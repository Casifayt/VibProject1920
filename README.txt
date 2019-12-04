-------------------------------------------------------
--------------- THEORY OF VIBRATION -------------------
- Analysis of the dynamic behaviour of a truss bridge -
----------------- ULiège 2019 - 2020-------------------
-------------------- FAYT Casimir ---------------------
-------------------------------------------------------

Dear Mrs Prijot, and other assistants whose names I 
don't know (sorry)

Two main files are available : 
- ThVib_Group40.m
- ThVib_Group40b.m

The first one is the one I used along the project. I 
selected which method to run by commenting or decommenting 
the corresponding function calls.

The main parameters are directly chosen in the script at lines :
	- Nbr of elements : 			line 39
	- Nbr of modes : 			line 66
	- Duration of simulation :		line 93
	- Time step of simulation :		line 94
	- Type of external force :		line 101
	- Nbr of DOFs for model reduction :	line 136

The second one does not need ot be modified directly in 
the script, I wrote an input function each time a parameter 
had to be chosen (apart from time properties attention, 
still at the lines given above)

---------------------------------------------------------

I implemented four different external shaker signals :
	- A null one (no force)
	- A step one (step at first quarter of total duration)
	- A dirac one (peak at first quarter of total duration)
	- A chirped one (the signal asked in the assignment)

And the possibility to choose for thenumber of degrees of freedom 
retained in the model reductions (between four nodes, eight or twelve)

Best regards,
FAYT Casimir
Casimir.Fayt@student.uliege.be