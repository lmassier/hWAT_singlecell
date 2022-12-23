# A Meta Perspective On The Cellular Landscape of Human White Adipose Tissue
**Lucas Massier, Jutta Jalkanen, Merve Elmastas, Jiawei Zhong, Tongtong Wang, Pamela A. Nono Nankam, Scott Frendo-Cumbo, Jesper Bäckdahl, Narmadha Subramanian, Takuya Sekine, Alastair Kerr, Tzu Pin Tseng, Jurga Laurencikiene, Marcus Buggert, Magda Lourda, Karolina Kublickiene, Nayanika Bhalla, Alma Andersson, Armand Valsesia, Arne Astrup, Ellen E. Blaak, Patrik L. Ståhl, Nathalie Viguerie, Dominique Langin, Christian Wolfrum, Matthias Blüher, Mikael Rydén, Niklas Mejhert**



## Abstract
Adipose single-cell studies have uncovered important biological findings. However, most reports are based on small cohort sizes and there is no cellular consensus nomenclature of human white adipose tissue (WAT). Herein, we performed a comprehensive meta-analysis of publicly available and newly generated results from human WAT single-cell, single-nucleus, and spatial transcriptomic data. Our high-resolution map is based on 401,320  cells from three depots of 14  cohorts including 103 samples from 83 subjects and allowed us to robustly identify >60 subpopulations of adipocytes, fibroblast and adipogenic progenitors, vascular, and immune cells. Using these results, we deconvolved bulk transcriptomic data from eight independent cohorts (864 subjects) and identified associations between cell populations and body weight, insulin resistance, dyslipidemia, fat cell volume, and lipolysis. Altogether, our meta-map of human subcutaneous, omental and perivascular WAT provides a rich resource that defines cell states and how they are affected by health and metabolic disease. 

![Cohort Overview](/images/Cohort_Overview.PNG)

## Description
This repo contains scripts to replicate results of our current study, as well as the source data and CellTypist models trained on this data set.  
Link to publication:  
Link to snSeq data:  



## Content

* `CellTypist`: CellTypist models trained on integrated data
* `FACS`: Raw data for included FACS analysis  
* `Source Data`: Source data for figure panels in the manuscript
* `Spatial_Deconvolution`: Code to replicate spatial deconvolution
* `scripts`: contains scripts used to analyze data
  * `bulk deconvolution`: depots specific and clinical RNA seq deconvolution
  * `Cell Types`: Sublcustering
  * `Integration`: Tested integration methods
  * `Network`: Marker gene comparison before integration



## Contact
For questions regarding data analysis, please write to lucas.massier@ki.se or jiawei.zhong@ki.se  
Corresponding authors: Niklas Mejhert (niklas.mejhert@ki.se) and Mikael Rydén (mikael.ryden@ki.se)
