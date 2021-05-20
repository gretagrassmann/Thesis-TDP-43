# Aggregation of the C-terminal fragment of the TAR DNA-binding Protein 43 (TDP-43) in relation to the Amyotrophic lateral sclerosis (ALS)
## Master degree's thesis in Applied Physics, University of Bologna Physics and Astronomy Department
### Advisors: Prof. Armando Bazzani (University of Bologna Physics and Astronomy Department) and Prof. Giancarlo Ruocco (Center for Life Nano-& Neuro-Science, Fondazione Istituto Italiano di Tecnologia (IIT), Rome).
### Co-advisors: Dr. Edoardo Milanetti and Dr. Claudia Testi (Center for Life Nano-& Neuro-Science, Fondazione Istituto Italiano di Tecnologia (IIT), Rome)

### Table of contents:
  * [Introduction](#introduction)
  * [Molecular dynamics simulations](#molecular-dynamics-simulations)
    + [Simulation analysis.py](#simulation-analysis.py)
    + [Cluster.py](#cluster.py)
  * [2D Zernike expansion](#2D-zernike-expansion)  
    + [Complete.py](#complete.py)
    + [My functions.py](#my-functions.py)
    + [Pca comparison.py](#pca-comparison.py)
    + [Screening plot.py](#screening-plot.py)
- [DISCLAIMER](#disclaimer)  
- [CONTACTS](#contacts)




### Introduction
The Amyotrophic Lateral Sclerosis (ALS) is a neurodegenerative disease specifically affecting cortical and spinal motor neurons. ALS is related to the formation in the cytoplasm of these neurons of protein aggregates.
The human TDP-43, a RNA/DNA binding protein involved
in RNA-related metabolism, is a major component of these pathological inclusions.<br />
While the deposition of the phosphorylated full-length TDP-43 in spinal-cord cells has been widely studied, it has been shown that the brain cortex presents
accumulation of phosphorylated C-terminal fragments (CTFs).
The CTF corresponds to a portion of the full protein including only the last 195-206 residues, and the full understanding of the mechanisms behind its aggregation is a still open and exciting challenge.<br />
We would like to begin our project from what is known in literature to date: CTFs are composed by the disordered C-terminal domain (CTD) and a fragment of RRM2, a folded domain of known structure. The latter could be fundamental importance for the protein's aggregation, since after the TDP-43 proteolysis it partially unfolds and exposes the aggregation prone beta-strands. These beta-strands should be at the core of the aggregation, because they are able to form steric zippers between different CTFs that then, following a typical model for amyloid fibril structure formation, give rise to amyloid structures.<br />
The aim of this work is to unveil these mechanisms: even if CTFs are not a primary cause of ALS, they are a hallmark of TDP-43 related neurodegeneration in the brain.<br />
With the hope of achieving such a result, we will follow three main steps:
1. **Molecular Dynamic (MD) simulations** to observe the conformation at equilibrium of the two kinds of CTFs that have been observed (corresponding to a cleavage at two different sites). In particular, we are hoping to see the formation of well-organized structural arrangements in the sites where the aggregation between different RRM2s should happen; in this case we would be able to study the complementarity between them with a promising approach based on the Zernike polynomial formalism.
2. **2D Zernike polynomial expansion** is a new method (developed at the IIT in 2020) for assessing whether and where two proteins can interact with each other to form a complex. In our case we are going to apply it to the 3D structures, obtained with the MD simulations, of two RRMs.
    As a next step, we will apply Zernike to design an aptamer able to interfere with the aggregation by binding to the site where -would the site be free- the two RRMs would connect. In this case, despite the increased time and computational cost, we could use the Zernike 3D formalism. Thanks to 3D Zernike we would be able to better describe a complementarity that involves extended regions of molecular surface (as in the case of the "flat" aptamers).
3. Finally, to test our theoretical models and computational methods' results, we will perform experimental measurements on cells in vitro, employing **Brillouin scattering**. Brillouin microscopy can probe the viscoelastic properties of biological samples: since the aggregates are characterised by a more solid consistency compared to the surrounding cytoplasm, they are clearly visible in this kind of microscopy. The ideal conclusion of this thesis would be to see a decrease of the aggregates' number and dimensions following the insertion of the aptamer.

A more in-depth explanation can be found !!! METTERE TESI QUANDO SARA' FINITA!!!


### Molecular dynamics simulations
The MD simulations were perfomed with the **GROMACS** software. In particular, we studied:
* The whole RRM2 domain of TDP-43.
* Fragment A, corresponding to the fragment of RRM2 that can be found in the CTFs after the cleavage ar residue 208.
* Fragemnt B, corresponding to the fragment of RRM2 that can be found in the CTFs after the cleavage ar residue 219.<br />

Each simulation was analyzed, using the **GROMACS** commands and the following codes (readable in this repository):
* ***simulation_analysis.py***
* ***cluster.py***

#### Simulation analysis.py
This code is implemented to study the equilibration phase of the simulation and the resulting trajectory. It can be used to visualize
* The *potential energy* minimization during equilibration.
* The *temperature* equilibration during thermalization.
* The *pressure* equilibration during pressurization.
* The evolution of *density* during equilibration. <br />

And
* The *Root Mean Square Deviation* (RMSD)  as  a  function  of  time  of  the evolving structure  respect  to  its  equilibrated  system and respect to its crystal structure,
as well as its distribution.
* The evolution of the *radius of gyration*. <br />
After the trajectory has been studied, the analysis of the possible molecule's conformation can be performed.<br />
With **GROMACS** we implement a Principal Component Analysis (PCA) of the trajectory. We can then visualize the results with this code, in particular:
* The *two-dimensional projection* of the sampled conformations in the subspace spanned by the first two eigenvector.
* The *Explained Variance Ratio* (EVR) for all the eigenvectors resulting from the PCA. <br />

As a next step, we can look for the most representative molecule's conformations, using the code described in the next Section.

#### Cluster.py
This code is implemented to find the most representative conformations, by doing a *k*-cluster analysis of the PCA of the trajectory of a molecule.<br />
The best value *K* for *k* is given by the minimization of the silhouette coefficient.<br />
The code allows for the visualization, for the selected *K* clusters, of both the silhouette plots for each cluster and the *K*-clustering of the scatter plot
of the two-dimensional projection of the sampled conformations (their centroids -labeled by the numbered white circle- are depicted as well).<br />
The points in the dataset (each point corresponds to the configuration of the molecule at a particular time-step of the trajectory evolution) closer to these centroids
are taken as the most representative configurations of the molecule. <br />
To study these configurations we can employ the Zernike method.

### 2D Zernike expansion
All the parameters that control this part of the analysis are defined in ***configuration.txt***, which includes:
* *Rs_select*: the radius of the sphere used to determine the patch whose roughness we are going to study.
* *R_zernike*: the radius of the sphere used to determine the patch whose shape we will reconstruct with the Zernike method.
* *alpha*: a list of the possible values of the parameter determining the ratio between *Rs_select* and *R_c*.<br />
*R_c* is used to determine on how many points we are going to build the patches to implement the Zernike method: given a patch of radius *Rs_select*, the points closer then *R_c* to its center will not be considered as centers for the patches of radius *R_zernike* that we are going to reconstruct with the Zernike method*. Since it is defined as *R_c=Rs_select\*mean_cosine\*alpha*, where *mean_cosine* is the average value of the cosines between the vectors normal to the patch, *R_c* is smaller the more rough the patch is. 
* *fragment*: the name of the directory corresponding to the studied fragment.
* *cluster*: the number of the studied centroid (i.e., one of the representative configurations).
* *step*: every how many points we build a patch whose roughness we want to study.
* *verso*: orientatio of the cone used to build the Zernike patches.<br />
1 compared to 1 -> shape similarity <br />
1 compared to -1 -> shape complementarity
* *Npoint*: every how many points we build a patch whose Zernike reconstruction we want to study.

#### Complete.py
This code is divided in four steps:
1. Calculation of the Zernike descriptors of all the possible points in the surface, defined with *R_zernike*.
2. Calculation of the roughness (i.e., the mean cosine value) for all these possible patches.
3. Evaluation of the points to consider for the Zernike method for all the values selected in *alpha*.
4. Calculation of the Zernike coefficients for each screening obtained in the previous step.

#### My functions.py
Here are defined some of the functions used in ***complete.py***.

#### Pca comparison.py
* This code starts by doing a PCA of the Zernike coefficients obtained in step 1. of ***complete.py***, and calulating the centroid *c_tot* and the inertia *in_tot*
of the founded cluster.<br />
* Then, for each value *alpha'* contained in *alpha*:
  * It projects on the resulting first two principal components, the Zernike coefficients obtained in step 4. of ***complete.py***.
  * It calculates the centroid *c_alpha'* and the inertia *in_alpha'* of this cluster.
  * It calcuates the value of the loss function *Loss(alpha')=(c_tot-c_alpha')\*(in_tot-in_alpha')\*n(alpha')*, where *n(alpha')* is the number of points found in the step 3. of ***complete.py***.
* Finally, it plots *Loss(alpha')* as a function of the values in *alpha*.<br />

The value *ALPHA* that results in the lowest value of *Loss(alpha')* is the best one.


#### Screening plot.py
This code takes the Zernike representation obtained with *ALPHA* and produces a 3D representation of the total points of the surface (in blu) and the one that are considered
as centers to build the patches for the Zernike representation.


## DISCLAIMER
The development of the 2D Zernike expansion and the corresponding code is not part of my work, but was passed to me by the original authors.


## CONTACTS
Please address any question to:
  * gretagrassmann0@gmail.com
