import numpy as np
from model import simulate
from analysis import set_parameter

#48h time, with 12 mn interval between each discretized observation
times = np.linspace(0, 48, 241)

params = {
    "ka": 1.0,          #Central compartment drug absorbtion rate
    "CL": 2.0,          #Drug clearance
    "V": 20.0,          #Drug gistribution volume
    "Imax": 0.8,        #Maximal Biomarker production inhibition
    "IC50": 0.2,        #Half maximal inhibitory concentration (Drug concentration inducing half inhibition effect)
    "kprod": 10.0,      #Biomarker production rate
    "kdeg": 0.1,        #Biomarker degradation rate
    "dose" : 50.0       #Drug dose
}

# First test : dose = 0mg
# Expected behavior With B(0) = B0 : 
# Constant B = B0, Ag = Ac = 0 mg
params_1 = set_parameter(params, "dose", 0.0)
result_1 = simulate(params=params_1, times=times)

Ag_1, Ac_1, B_1 = result_1.y
B0 = params_1["kprod"]/params_1["kdeg"]

#Test hypotheses and print success
print("Testing B = B0 with D = 0 mg")
np.testing.assert_allclose(B_1, B0, rtol=0, atol=1e-3)
print("B is constantly close to B0\n")

print("Testing Ag + Ac = D with D = 0 mg")
np.testing.assert_allclose(Ag_1 + Ac_1, 0.0, rtol=0, atol=1e-3)
print("Ag + Ac is constantly close to 0\n")


# Second test : Imax = 0
# Expected behavior : 
# Constant B = B0
    
params_2 = set_parameter(params, "Imax", 0.0)
result_2 = simulate(params=params_2, times=times)

Ag_2, Ac_2, B_2 = result_2.y
B0 = params_2["kprod"]/params_2["kdeg"]

#Test hypothesis and print success
print("Testing B = B0 with Imax = 0%")
np.testing.assert_allclose(B_2, B0, rtol=0, atol=1e-3)
print("B is constantly close to B0\n")



# Third test : CL = 0
# Expected behavior : 
# Constant Ac + Ag = D

params_3 = set_parameter(params, "CL", 0.0)
result_3 = simulate(params=params_3, times=times)

Ag_3, Ac_3, B_3 = result_3.y

#Test hypothesis and print success
print("Testing Ag + Ac = D with CL = 0 L/h")
np.testing.assert_allclose(Ag_3 + Ac_3, params_3["dose"], rtol=0, atol=1e-3)
print("Ag + Ac is constantly close to D\n")




