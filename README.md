**Summary**

The code in this repository contains the scripts to run the analyses described in the manuscript " Integrating passive acoustic data from a towed hydrophone array with visual line transect data to estimate abundeca and surface availability: A case study with sperm whales (*Physeter macrocephalus*)" by Sigourney et al. 2022 (in review at Methdods in Ecology and Evolution). Code is provided to recreate the simulations and run the models on the simulation code and code to re-run the analysis of the sperm whale data presented in the text. There are 7 scripts in total with 5 scripts that are used in the simulation analysis and 2 scripts to run the sperm whale analysis

**Simulation code**

*Run_simulations.R* – This script includes the parameter settings for the simulations including the number of whales (N.true), parameters to simulate detection and parameters to choose the size of the window (i.e., truncation distances). Window is currently set to 6 KM for both the perpendicular (x) distance and forward (y) distance and the detection parameters are set so detection is approximately 1 at zero distance decreasing to \<0.05 at 6 KM.

*Simulate_Data.R*- This script can be used to simulate datasets where the output is similar to the raw input for the real model once the acoustic data have been analyzed and summarized by an acoustician. At the start of each simulation each whale is randomly assigned a perpendicular distance (perp.dist). They also are randomly assigned to a specific point in their dive cycle. All whales start at a forward distance of 6 KMs. Detection for visual data is assume to follow a half-normal detection process. Detection for the acoustic data for testing out the CMR component of the methods is based on a logistic model. For the Hybrid method detection for the acoustic data for the DS component of the analysis is also assumed to follow the same half-normal detection process as the visual data.

*Organize_Simulation_Data.R* – This script takes reformats the data to be in the correct format to be analyzed by the CMR methods. Specifically, it takes the simulated observations from below and assigns them to forward distance bins and creates a capture history matrix.

*CMR-DS Models.R* - This script runs the CMR-DS analysis described in the text

*Hybrid Model.R* – This script run both the Hybrid method and for comparison also includes the DS-DS method.

**Sperm whale analysis code**

*Run_Hybrid_Model_with_sperm_whale_data.R* – This script reads in the acoustic data and visual data collected during the 2013 AMAPPS cruises by the NEFSC. The truncation distance was chosen to be 7.6 KM (7600 meters). Data are organized so they can be analyzed using the Hybrid Method described in the text.

*Hybrid_model_for_sperm_whales.R* – This script analyses the sperm whale data using the Hybrid Method. To analyze both the visual data and acoustic data a hazard rate model is used (in lieu of a half-normal detection function). Also, because the visual line transects used a double observer approach mark-recapture distance sampling (MRDS) is use to analyze the visual data to correct for detectability on the trackline, p(0). Finally, an estimate of availability bias is included in the visual analysis to compare the corrected abundance estimates to the estimates from the Hybrid Method.

*Disclaimer*

This repository is a scientific product and is not official communication of the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All NOAA GitHub project code is provided on an ‘as is’ basis and the user assumes responsibility for its use. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.
