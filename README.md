# Wonder drug project (Malaria/ivermectin model)
The following details the calculations/scripts used to obtain results on the impact of ivermectin on reducing malaria incidence in 41 sub-Saharan African countries. 

MAIN script IVMMalScriptRANDsamplingJan30_2025

Estimate transmission rates: findAlphas, calcAB2025

Estimate Skill Scores: Valdiate2019, calcSkillScore2025

Estimate the effect of ivermectin on mosquitoes/parasites: calibrateIvermectinEffectsSAMP2025

Compute baseline malaria predictions: baselineSAMP2025

Compute malaria prections under ivermectin intervention scenarios: IVMinterventionSAMP2025

Setup for the calculation of Disability adjusted life-years: setupForDALYsSAMP2025


Zip-file of Matlab code [MalariaIvermectinProject.zip](https://github.com/mathguypi314/GSISmodel/blob/main/MalariaIvermectinProject.zip)

# GSISmodel
code for parameter fitting / AIC calculation of gSIS and SIS models

The code fits the differential equation

di/dt = beta\*(1-i)\*i - (2\*m\'+1)/m\*i\*m - b\*i

to incidence data on gonorrhea in the US, where beta & b are constant and m is the mean residual waiting-time (& its derivative m').

The code also fits the special case when m = constant = 1/gamma (i.e. the classical SIS model).

The quality of fits is compared by calculating AIC for each model. 

gSIS version from medRxiv preprint
[gSIS project.zip](https://github.com/mathguypi314/GSISmodel/files/8800303/gSIS.project.zip)

gSIS version from publication in Infectious Disease Modeling Journal
[gSISexactsolnsOct102023.m](https://github.com/mathguypi314/GSISmodel/blob/main/gSISexactsolnsOct102023.m)
