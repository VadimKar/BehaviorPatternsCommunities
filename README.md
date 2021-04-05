# BehaviorPatternsCommunities
Supplemental material for model simulation and fitting in Karatayev et al "Behavior can regulate larger-scale community patterns"



Behavior_2-12-2021.R contains all required functions.
CICAalldat12all.rds and NENZdatFull.rds contain data used for fitting in California and New Zealand, respectively.

For other uses of this data, please refer to use policies given in primary sources referenced in Appendix A (Kushner et al. 2013 and the PISCO monitoring program in California).



Each data file contains the following colums:

SiteN: Site number (in California, numbers <100 denote NPS KFM sites, numbers >99 denote PISCO sites)

Yr: Sampling year (null in New Zealand)

QuadN: Quadrat Number - the sample ID within each site. For California NPS sites, subsequent quadrat numbers denote spatially adjacent samples

z: Sample depth

Uq: Urchin density in sample

JObs: number of juvenile kelp in sample (as defined in Appendix A and primary sources); not used due to limited coverage

Obs: number of adult kelp in sample

AdArea: area of sample in square meters

N: (relative) modeled nitrate availability (CA only)

L: relative light availability, adjusted for depth

E: relative wave stress (in NZ quantified as wind fetch)

U: site-level urchin density

R: rugosity of sample substrate (available only for a few CA samples)

P: Predator density at site
