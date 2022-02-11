
--------------------------------------------------------------------------------------------------------------------------------

If you are interested specifically in the detailed documentation and code needed to reproduce the tables and figures in the paper 'Morphing for Consumers Dynamics: Bandits meet HMM' by Gui Liberali and Alina Ferecatu (2022), please refer to the replication files in this link: https://pubsonline.informs.org/doi/suppl/10.1287/mksc.2021.1346. 

--------------------------------------------------------------------------------------------------------------------------------

This GitHub repository is made public to promote the use of multi-armed bandits (MABs) in marketing and is constantly evolving. 
It contains code used in Liberali and Ferecatu (2022), and more. 
The code and data uploaded to the journal website (see link above) can only be used to replicate the findingds of the published article per regulations of the publisher. 
The code and data found in this GitHub repository, on the other hand, can be freely used and adapted by you in your own projects or future research.  

Before using the code in this repository, you should first run 'Configuration.R' and 'Functions.R'. 
- 'Configuration.R' defines global variables necessary for other files, and loads + installs the necessary packages.
- 'Functions.R' contains all helper functions needed for subsequent analysis performed in the repository. 

The data from the two studies reported in Liberali and Ferecatu (2022) are made available here. Each has a separate folder, structured as follows. 

<p>
    <img src="Repository Overview.PNG" width="800" height="600" />
</p>

Main steps:

Step 1. Find the table or figure you are interested in the '4. Figures & Tables' map.
For instance, for Table 7 of study 2, this file would be in study2/4. Figures & Tables/Table7.R. 

Step 2. In case it is necessary, run the files from '3. Simulation files' required for the respective figure or table.
For example: for Table 7 of study 2, we would have to change our working directory to '3. Simulation code' and then run: 

* Simulation Results for One Morphing Opportunity.R
* Simulation Results for Two Morphing Opportunities.R
* Simulation Results for Random Morphing.R

After these files have been run, there will be several variables/objects/parameters defined in the environment. These can be used to create the respective table. 

Note that for the HMM estimates for study 1 you need a working Stan distribution and a C++ distribution on your computer. 

Step 3. Run the files from '4. Figures & Tables' for the respective figure or table, using the generated objects. Do not forget to change the working directory back to '4. Figures & Tables'. 

Some parts of the code will take a bit of time to run: for instance, Gmatrix.out, as well as the files in '3. Simulation code' for both studies. 
