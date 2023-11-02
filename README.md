# Fibromyalgia_FACS_B_Cell_Gating

Fibromyalgia syndrome (FMS) is a condition that results in chronic pain and a multitude of other symptoms such as fatigue, depression and cognitive disorders, responsible for impacting a patient’s quality of life.  Recently, the passive transfer of IgG isolated from FMS patients into mice resulted in the transfer of multiple FMS symptoms suggesting the involvement of B cells in the pathology of FMS. [1]

The B cell gating code was used to gate and analyse B cell subpopulations in FMS patients. The FACS data for this code was obtained from a study conducted by the Alan Edwards Centre for Research on Pain in McGill University [2]. 

The code includes compensation, clean up and transformation of the raw FACS data. The strategy uses the package flowCore [3]. For B cell subpopulation analysis, a CD27/IgD quadGate was applied. CD38 gates were applied to the CD27+/IgD+, CD27+/IgD-, and CD27-/IgD- subpopulations. CD10 gating was performed on the CD27-/IgD+ subpopulation, with further CD21 gating performed on the CD10- subpopulation.

This study was conducted under the guidance of Dr. Rachael Bashford-Rogers at the Department of Biochemistry, University of Oxford.






References:

[1] Goebel, Andreas et al. “Passive transfer of fibromyalgia symptoms from patients to mice.” J Clin. Invest , vol. 131,13 (2021): e144201. doi:10.1172/JCI144201 

[2] Vivek et al. “Unbiased immune profiling reveals a natural killer cell-peripheral nerve axis in fibromyalgia.” Pain vol. 163,7 (2022): e821-e836. doi:10.1097/j.pain.0000000000002498 

[3] Ellis B, Haaland P, Hahne F, Le Meur N, Gopalakrishnan N, Spidlen J, Jiang M, Finak G (2023). flowCore: flowCore: Basic structures for flow cytometry data. doi:10.18129/B9.bioc.flowCore, R package version 2.14.0, https://bioconductor.org/packages/flowCore.
