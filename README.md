# velocity-curvature-power-law-protocol

**Biological kinematics: a detailed review of the velocity-curvature power law calculation**

Dagmar S. Fraser, Massimiliano Di Luca, and Jennifer L. Cook

*Correspondence: d.s.fraser@bham.ac.uk (Dagmar S. Fraser).

**Abstract**

The 'one-third power law', relating velocity to curvature is among the most established kinematic invariances in bodily movements. Despite being heralded amongst the “kinematic laws of nature” (Flash 2021, p. 4), there is no consensus on its origin, common reporting practice, or vetted analytical protocol. Many legacy elements of analytical protocols in the literature are suboptimal, such as noise amplification from repeated differentiation, biases arising from filtering, log transformation distortion, and injudicious linear regression, all of which undermine power law calculations. Recent findings of power law divergences in clinical populations have highlighted the need for improved protocols. This article reviews prior power law calculation protocols, identifies suboptimal practices, before proposing solutions grounded in the kinematics literature. We evaluate these protocols via two simple criteria: firstly, they must avoid spurious confirmation of the law, secondly, they must confirm the law when it is present. Ultimately, we synthesise these vetted solutions into a modular protocol which we make freely available to the scientific community. The protocol’s modularity accommodates future analytical advances and permits re-use in broader kinematic science applications. We propose that adoption of this protocol will eliminate artificial confirmation of the law and facilitate more sensitive quantification of recently noted power law divergences, which are associated with neurochemical disturbances arising from dopaminergic drugs, and in conditions such as Parkinson's and autism. 

Keywords Two-thirds power law · One-third power law · Kinematics · Noise · Filtering · Regression

**Description**

This repository contains protocols to calculate the Velocity Gain Factor and Beta exponent of the velocity-curvature one-third power law.  Legacy calculations are presented side by side with vetted analysis choices extracted from the wider literature.  This repository supports the preprint available here - https://osf.io/preprints/psyarxiv/vfq3d .

**How to use**

Clone the repository.  Or launch in MATLAB Online [![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=dagmarfraser/velocity-curvature-power-law-protocol&project=https://github.com/dagmarfraser/velocity-curvature-power-law-protocol/blob/main/velocity-curvature-power-law-protocol.prj)

For synthetic data - Run _PowerLawSynthetic.m_ after navigating to the _src_ subfolder of the main _PowerLawToolChainEBRsubmit_ folder.
Edit the variable _paramChoice_ to choose between **1** Maoz et al. 2005 or **2** Schaal and Sternad 2005 re-implementaions.
Revised graphs for v2 of the paper may be obtained by running _RevisedSynthetic.m_
Figures will be saved in the _figures_ subfolder.  

For empirical data - Run _PowerLawEmpirical.m_ after navigating to the _src_ subfolder of the main _PowerLawToolChainEBRsubmit_ folder.
Revised graphs for v2 of the paper may be obtained by running _RevisedLawEmpirical.m_
This repositiory contains data from the available N=14 (of 40) participants of Zarandi et al. 2023 
Available here https://github.com/lucaoneto/IJCNN2022_Ellipses/tree/main/data
Figure will be saved in the _figures_ subfolder.  

Additional graphs from V2 of the paper _BetaEBR.m, shapesHuhEBR.m & shapesHuhAndBetaEBR.m_.

**Acknowledgments**

This project has received funding from the European Union’s Horizon 2020 Research and Innovation Programme under ERC-2017-STG Grant Agreement No 757583
Massimiliano Di Luca is partially supported by the Engineering and Physical Sciences Research Council 

**Bibliography** [1] Flash T (2021) Brain Representations of Motion Generation and Perception: Space-Time Geometries and the Arts. In: Flash T, Berthoz A (eds) Space-Time Geometries for Motion and Perception in the Brain and the Arts. Springer International Publishing, Cham, pp 3–34. https://doi.org/10.1007/978-3-030-57227-3_1


