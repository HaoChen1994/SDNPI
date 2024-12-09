Algorithms  
-------
SDNPI: A simultaneous differential network and pathway identification analysis


Maintainer
-------
Hao Chen   <hao.chen@sdu.edu.cn>
Chencheng Ma   <machencheng@mail.sdu.edu.cn>
Xiang Zhan   <zhanx@seu.edu.cn>


Publication
-------
Chen H, He Y, Ji J, Ma C, Zhan X (2024). Simultaneous differential network and pathway identification for analyzing effects of brain functional connectivities on the pathway from apoe gene to alzheimer's disease. Submitted to the Annals of Applied Statistics.


Abstract
-------
The Apolipoprotein E (ApoE) gene is a well-known genetic factor that increases the risk of Alzheimer's disease (AD). However, the underlying mechanism for how APOE influences AD has yet to be fully elucidated. Prevailing evidence suggests that ApoE gene is associated with alterations in brain functional connectivities, which play an important role in AD's onset and progression. This paper investigates the intermediate role of brain functional connectivities or network as high-dimensional mediators on the causal pathway from ApoE gene to AD. To address analysis with high-dimensional mediators and a binary outcome variable (i.e., AD status), we propose a Simultaneous Differential Network and Pathway Identification (SDNPI) procedure, which first estimates individual-specific brain functional connectivities and then conducts penalization-based mediation analysis for binary outcomes. A novel MM-ADMM computing algorithm which combines majorization-minimization (MM) and alternating direction method of multipliers (ADMM) is proposed to solve the optimization problem involved in SDNPI. After demonstrating SDNPI's effectiveness using simulation studies, we use it to investigate how brain functional connectivities mediate the effect of ApoE gene on Alzheimer's disease, and find that several critical nodes along with their regions of interest in the brain to shed new light on the biological and genetic interpretations of Alzheimer's disease.

KEYWORDS: Alzheimer's disease, Brain connectivities, Differential network, High-dimensional mediators, Pathway Lasso.


Usage
-------
1. EstMediator.R: This R script contains functions designed to estimate the mediators (individual brain functional connectivity) associated with the first step of SDNPI.

2. PathLasso.R: This R script contains functions designed to solve the optimization problem associated with the second step of SDNPI.

2. Example.R: This script provides a practical example of applying the SDNPI methodology, including a demonstration on identifying differential networks and pathways using a simulated dataset tailored to mimic conditions relevant to brain functional connectivity studies.