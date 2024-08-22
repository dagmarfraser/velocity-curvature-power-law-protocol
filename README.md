# velocity-curvature-power-law-protocol

The velocity-curvature power law in review

Dagmar S. Fraser, Massimiliano Di Luca, and Jennifer L. Cook

*Correspondence: d.s.fraser@bham.ac.uk (Dagmar S. Fraser).

**Abstract**

Bodily movements exhibit kinematic invariances, with the “one-third power law” relating velocity to curvature amongst the most established. Despite being heralded amongst the “kinematic laws of nature” (Flash 2021, p. 4), there is no consensus on its origin, common reporting practice, or vetted analytical protocol. Many legacy elements of analytical protocols in the literature are suboptimal, such as noise amplification from repeated differentiation, biases arising from filtering, log transformation distortion, and injudicious linear regression, all of which undermine power law calculations. This article reviews prior power law calculation protocols, identifies suboptimal practices, before proposing solutions grounded in the kinematics literature and related fields of enquiry. Ultimately, we synthesise these solutions into a vetted, modular protocol which we make freely available to the scientific community. The protocol’s modularity accommodates future analytical advances and permits re-use of modules useful in broader kinematic science applications. We propose that adoption of this protocol will eliminate spurious confirmation of the law and enable more sensitive quantification of recently noted power law divergences. These divergences have been linked to neurochemical disturbances arising from ingestion of dopaminergic drugs, and in neurological conditions such as Parkinson’s and autism. 


Keywords: two-thirds power law; one-third power law, kinematics, noise, filters, regression.

**Description**

This repository contains protocols to calculate the Velocity Gain Factor and Beta exponent of the velocity-curvature one-third power law.  Legacy calculations are presented side by side with vetted analysis choices extracted from the wider literature.

**How to use**

Clone the repository.  

For synthetic data - Run _PowerLawSynthetic.m_ after navigating to the _src_ subfolder of the main _PowerLawToolChainEBRsubmit_ folder.
Edit the variable _paramChoice_ to choose between **1** Maoz et al. 2005 or **2** Schaal and Sternad 2005 replications.

Figures will be saved in the _figures_ subfolder.  

For empirical data - Run _PowerLawEmpirical.m_ after navigating to the _src_ subfolder of the main _PowerLawToolChainEBRsubmit_ folder.
This repositiory contains data from the available N=14 (of 40) participants of Zarandi et al. 2023 
Available here https://github.com/lucaoneto/IJCNN2022_Ellipses/tree/main/data

Figure will be saved in the _figures_ subfolder.  

Flash T (2021) Brain Representations of Motion Generation and Perception: Space-Time Geometries and the Arts. In: Flash T, Berthoz A (eds) Space-Time Geometries for Motion and Perception in the Brain and the Arts. Springer International Publishing, Cham, pp 3–34. https://doi.org/10.1007/978-3-030-57227-3_1




