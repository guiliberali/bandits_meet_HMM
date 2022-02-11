# Replication_Morphing_HMM

This repository contains the replication code for the paper 'Morphing for Consumers Dynamics: Bandits meet HMM' by Gui Liberali and Alina Ferecatu from the Rotterdam School of Management at Erasmus University. 

The repository is structured as follows:
- 'Configuration.R' defines global variables necessary for other files, and loads + installs the necessary packages.
- 'Functions.R' contains all helper functions needed for subsequent analysis performed in the repository. 

Anyone aiming to replicate the results from the paper should first run 'Configuration.R' and 'Functions.R'. In the paper, two studies are performed - each has a separate folder, which is structured as follows. 

<p>
    <img src="Repository Overview.PNG" width="800" height="600" />
</p>

To replicate a figure or table, in general the following steps should be followed:

**Step 1.** Check the table or figure to be replicated in the '4. Figures & Tables' map.
So for instance, for Table 7 of study 2, this file would be in study2/4. Figures & Tables/Table7.R. 

**Step 2.** In case it is necessary, run the files from '3. Simulation files' required for the respective figure or table.

To give an example: for Table 7 of study 2, we would have to change our working directory to '3. Simulation code' and then run: 

* Simulation Results for One Morphing Opportunity.R
* Simulation Results for Two Morphing Opportunities.R
* Simulation Results for Random Morphing.R

After these files have been run, there will be several variables/objects/parameters defined in the environment. These can be used to create the respective table. 

**Note** the HMM estimates for study 1 can only be replicated with a working stan and C++ distribution on your computer. 

**Step 3.** Run the files from '4. Figures & Tables' for the respective figure or table, using the generated objects. Do not forget to change the working directory back to '4. Figures & Tables'. 

It is worth noting that some parts of the code will take a bit of time to run: for instance, Gmatrix.out, as well as the files in '3. Simulation code' for both studies. 
