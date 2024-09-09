# GSISmodel
code for parameter fitting / AIC calculation of gSIS and SIS models

The code fits the differential equation

di/dt = beta\*(1-i)\*i - (2\*m\'+1)/m\*i\*m - b\*i

to incidence data on gonorrhea in the US, where beta & b are constant and m is the mean residual waiting-time (& its derivative m').

The code also fits the special case when m = constant = 1/gamma (i.e. the classical SIS model).

The quality of fits is compared by calculating AIC for each model. 

[gSIS project.zip](https://github.com/mathguypi314/GSISmodel/files/8800303/gSIS.project.zip)

[gSISexactsolnsOct102023](https://github.com/mathguypi314/GSISmodel/files/8800303/gSISexactsolnsOct102023.m)
