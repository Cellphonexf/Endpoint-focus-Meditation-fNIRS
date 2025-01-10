##Endpoint-Focus Meditation Regulates Strategic Present-Orientation: An Exploratory fNIRS Study##

Welcome to the GitHub repository for our research project. This repository contains the code and data used in our study.

##Data and Code Usage##

##Behavioral Data Analysis & Visualization##

Data: behavioral_rawdata.xlsx
Script: BehavioralDataAnalysis.R & BehavioralDataViso.R

- Run the script segment by segment to follow each step of our analysis.
- Generated graphs: pic_timing.pdf, pic_logk.pdf, pic_money.pdf, Long-logk.pdf, Short-logk.pdf
- Final figure: Fig. 2 & Fig. 3 (integrated using Adobe Illustrator)

##Hemodynamic Data Analysis & Visualization##

Data: Folder fNIRS_preprocessed_data (please contact the corresponding author for use)
Script: data_analysis_fNIRS.m

- Run in MATLAB (version R2017b) with Homer3 toolboxes (https://doi.org/10.1364/AO.48.00D280).
- Data info and exclusions: fNIRS_participant_list.xlsx
- Preprocessing summary: ProcStreamFunctionsSummary (in the folder)
- Generated graphs: HbO overall contrast Pair1.pdf, HbO overall contrast Pair2.pdf, FC endpoint.pdf, FC mindfulness.pdf (saved in Folder fNIRS_pics)
- Final figure: Fig. 4 & Fig. 5 (integrated using Adobe Illustrator)

##Overall Brain Activation Visualization##

Data for plotting: fNIRS_visualization.xlsx
Script: visualization_fnirs.m

-Toolboxes: EasyTopo (please refer to: 'EasyTopo 2.0 User Manual.pdf')
-Generated data files: input_contrastFront_data.mat
-Visualized maps: Meditation Contrast Fig.eps (saved in Folder fNIRS_pics)
-Final figure: Fig. 4 (integrated using Adobe Illustrator)

##Overall Functional Connectivity Visualization##

Data for plotting: FC_table.xlsx
Script: data_analysis_fNIRS.m

-Toolboxes: BrainNet Viewer (please refer to: 'BrainNet_viewer.pdf')
-Generated data files: FC_table.xlsx (sheets 'endpoint' & 'mindfulness')
-Generated edge files: FC_endpoint.edge, FC_mindfulness.edge (saved in Folder fNIRS_pics)
-Generated node files: nodes_endpoint.node, nodes_mindfulness.node (saved in Folder fNIRS_pics)
-Visualized maps: Brain_FC_endpoint.eps, Brain_FC_mindfulness (saved in Folder fNIRS_pics)
-Final figure: Fig. 5 (integrated using Adobe Illustrator)

##Activation Dynamic Visualization##

Data for plotting: analyzed data from Script data_analysis_fNIRS.m
Script: data_analysis_fNIRS.m

-Generated figure (all valid channels): saved in Folder Epoch_figs 
-Generated figure (selected channels): saved in Folder Epoch_figs_selected
-Final figure: Fig. 6 (integrated using Adobe Illustrator)

 ##Behavioral-Neural Correlational Analysis & Visualization##

Data: correlation_rawdata.xlsx
Script: CorrelationAnalysis.R

- Run the script segment by segment to follow each step of our analysis.
- Generated graphs: MTG-k_endpoint.pdf, MTG-k_mindfulness.pdf
- Final figure: Fig. 7 (integrated using Adobe Illustrator)

##Spatial Registration Results##

Files: SpatialRegistration_BA.txt, SpatialRegistration_MNI.txt
Toolboxes: SPM8 and NIRS-SPM (MATLAB version R2013b)
Reference: https://doi.org/10.1016/j.neuroimage.2008.08.036

##Edited on January 10, 2025