import numpy as np
import matplotlib.pyplot as plt
from deepdiff import DeepDiff
from pprint import pprint
from rafiki.rafiki import Rafiki, RafikiInputParameters, RafikiOutputParameters
from yeti.yeti import Yeti, YetiInputParameters
import yeti.yeti
from injector import TripletInjector

# Using IPA and LOX densities respecitivey
FuelDensity__kg_per_m3 = 786
OxidiserDensity__kg_per_m3 = 1141

InputParameters = RafikiInputParameters(	fuel="Isopropanol",
											oxidiser="LOX",
											fuel_concentration=99,
											oxidiser_concentration=100,
											oxidiser_fuel_mass_ratio=1.2,
											peak_thrust__N=2000,
											chamber_pressure__bar=10,
											ambient_pressure__bar=1	)

FineTunedRafikiOutputParameters = RafikiOutputParameters(	epsilon=2.193,
															gamma=1.214,
															Isp__s=214.838,
															Tc__K=2668.905,
															M__kg_per_kmol=18.940,
															Pr_exit=0.412,
															mu_exit__Pa_s=8.792e-05,
															k_exit__W_per_m_K=0.419,
															R__J_per_kg_K=438.958,
															m_dot__kg_per_s=0.949,
															mf_dot__kg_per_s=0.431,
															mo_dot__kg_per_s=0.518,
															engineGeometry__mm=((-212.3, 50), (-62.3, 50), (0, 22.5), (44.9, 34)))

CoolingParameters = YetiInputParameters(	chamberPressure__bar=InputParameters.chamber_pressure__bar.value,
											rafikiOutput=FineTunedRafikiOutputParameters,
											coolant=yeti.yeti.Water,
											wallMaterial=yeti.yeti.Al7049,
											wallThickness__mm=1,
											channelHeight__mm=3,
											channelWidth__mm=2,
											channelCount=40,
											blockageRatio=None,
											coolantInletTemperature__degC=25,
											coolantInletPressure__bar=25,
											coolantMassFlowRate__kg_per_s=5	)

# TODO: Review element count
Injector = TripletInjector(dischargeCoefficient=0.65, elementCount=8, pressureDrop__bar=6)

def calculateInjectorGeometry(mf_dot__kg_per_s, mo_dot__kg_per_s, fuelDensity__kg_per_m3, oxdiserDensity__kg_per_m3, injector):
	# Source: http://athena.leedsrocketry.co.uk/doku.php?id=louie_injector_desgin
	fuelArea__mm2 = (mf_dot__kg_per_s / (injector.dischargeCoefficient * np.sqrt(2 * fuelDensity__kg_per_m3 * injector.pressureDrop__Pa))) * 10**6
	oxidiserArea__mm2 = (mo_dot__kg_per_s / (injector.dischargeCoefficient * np.sqrt(2 * oxdiserDensity__kg_per_m3 * injector.pressureDrop__Pa))) * 10**6

	injector.fuelOrificeDiameter__mm = np.sqrt((2 * fuelArea__mm2) / (injector.elementCount * np.pi))
	injector.oxidiserOrificeDiameter__mm = 2 * np.sqrt(oxidiserArea__mm2 / (injector.elementCount * np.pi))

def displayFineTunedParameters(rafikiOutput):
	print("\nRAFIKI OUTPUT VS FINETUNED VALUES: ")
	pprint(DeepDiff(rafikiOutput[0].outputParameters, FineTunedRafikiOutputParameters, ignore_order=True, ignore_numeric_type_changes=True, significant_digits=3)["values_changed"])

def calculateBurnTime(outputParameters):
	# Work out the maximum burn time according to the tank volumes and denisties of the fuel and oxidiser according to AEL's ICD
	# (we've assumed the LOX tank volume will be no less than that of the NO2 tank)
	maximumFuelBurnTime__s = ((3.5 * 0.786) / outputParameters.mf_dot__kg_per_s.value)
	maximumOxidiserBurnTime__s = ((10 * 1.141) / outputParameters.mo_dot__kg_per_s.value)
	maximumBurnTime__s = min(maximumFuelBurnTime__s, maximumOxidiserBurnTime__s)

	return maximumBurnTime__s

def calculateCoolantPressureDrop(rafikiOutput, inputParameters, cooling_data):
	"""
	Pressure drop assumptions:
		- Inlet and outlet distribution rings are pipes curved along a half-circular arc of the
			same radius as the inner chamber diameter + the wall thickness.
		- The coolant takes a 90 degree turn into the cooling channels, resulting in
			a K-value of 1.3 (see config for reference).
		- The same number of delivery and return pipes as number of inlet and outlet
			pipe.
	"""

	T = inputParameters.coolantInletTemperature__K
	p = inputParameters.coolantInletPressure__Pa

	formattedEngineGeometry__m = np.array(rafikiOutput.engineGeometry__mm.value).transpose() / 1000

	xs__m = formattedEngineGeometry__m[0]
	rs__m = formattedEngineGeometry__m[1]
	
	cool_delivery_len__m = 1 # length of coolant delivery pipe [m]
	cool_delivery_d__m = 0.025 # diameter of coolant delivery pipe [m]
	cool_CD = 0.7 # coolant coefficient of discharge
	cool_inlet_pipe_d__m = 0.015 # coolant inlet pipe diameter [m]
	cool_inlet_pipe_len__m = 0.033 # coolant inlet pipe length [m]
	cool_inlet_ring_r__m = 0.049 # outer radius of the coolant inlet ring [m]
	cool_inlet_ring_w__m = 0.01 # inlet ring width (modeled as a rectangular pipe) [m]
	cool_inlet_ring_h__m = 0.008 # inlet ring height (modeled as a rectangular pipe) [m]
	cool_outlet_pipe_d__m = 0.015 # coolant outlet pipefice diameter [m]
	cool_outlet_pipe_len__m = 0.033 # coolant outlet pipe length [m]
	cool_outlet_pipe_n = 4 # number of coolant outlet pipe
	cool_outlet_ring_r__m = 0.038 # outer radius of the coolant outlet ring [m]
	cool_outlet_ring_w__m = 0.025 # outlet ring width (modeled as a rectangular pipe) [m]
	cool_outlet_ring_h__m = 0.01 # outlet ring height (modeled as a rectangular pipe) [m]
	cool_return_len__m = 1 # length of return pipe [m]
	cool_return_d__m = 0.019 # diameter of return pipe [m]
	cool_halo_inlet_len__m = 0.02 # length of halo inlet pipe [m]
	cool_halo_outlet_len__m = 0.02 # length of halo outlet pipe [m]
	cool_halo_pipe_d__m = 0.011 # diameter of halo inlet pipe [m]
	cool_halo_ring_d__m = 0.181 # diameter of halo ring [m]
	cool_halo_tube_d__m = 0.013 # diameter of halo tube [m]
	channel_width__m = (2/1000)

	# build the pressure drop circuit
	pressure_circuit = PressureDropCircuit(inputParameters.coolant.laminarRe, inputParameters.coolant.turbulentRe)
	# delivery pipe
	pressure_circuit.add_straight_pipes_circ(r__m=cool_delivery_d__m/2,
													n_pipes=1,
													len__m=cool_delivery_len__m,
													mdot__kg_s=inputParameters.coolantMassFlowRate__kg_per_s,
													rho__kg_m3=inputParameters.coolant.density__kg_per_m3(T, p),
													mu__Pas=inputParameters.coolant.viscosity__Pa_s(T, p),
													label="Delivery pipe")
	# halo inlet pipe 
	pressure_circuit.add_straight_pipes_circ(r__m=cool_halo_pipe_d__m/2,
													n_pipes=1,
													len__m=cool_halo_inlet_len__m,
													mdot__kg_s=inputParameters.coolantMassFlowRate__kg_per_s,
													rho__kg_m3=inputParameters.coolant.density__kg_per_m3(T, p),
													mu__Pas=inputParameters.coolant.viscosity__Pa_s(T, p),
													label="Halo inlet pipe")
										
	# halo distribution ring
	pressure_circuit.add_circular_pipe_circular(r__m=cool_halo_ring_d__m/2,
														rp__m=cool_halo_tube_d__m/2,
														percentage_travelled=0.5,
														mdot__kg_s=inputParameters.coolantMassFlowRate__kg_per_s,
														rho__kg_m3=inputParameters.coolant.density__kg_per_m3(T, p),
														mu__Pas=inputParameters.coolant.viscosity__Pa_s(T, p),
														n=1,
														label="Halo distribution inlet")

	# halo turn into cooling jacket inlet
	pressure_circuit.add_90deg_turn_circular(r__m=cool_inlet_pipe_d__m/2,
													mdot__kg_s=inputParameters.coolantMassFlowRate__kg_per_s/4,
													rho__kg_m3=inputParameters.coolant.density__kg_per_m3(T, p),
													n=3,
													label="Halo turn into jacket inlet")

	# inlet pipe
	pressure_circuit.add_straight_pipes_circ(r__m=cool_inlet_pipe_d__m/2,
													n_pipes=4,
													len__m=cool_inlet_pipe_len__m,
													mdot__kg_s=inputParameters.coolantMassFlowRate__kg_per_s/4,
													rho__kg_m3=inputParameters.coolant.density__kg_per_m3(T, p),
													mu__Pas=inputParameters.coolant.viscosity__Pa_s(T, p),
													label="Inlet pipe")

	# inlet distribution ring
	pressure_circuit.add_circular_pipe_rectangular_box(r__m=rs__m[0]+inputParameters.wallThickness__m,
															w__m=cool_inlet_ring_w__m,
															h__m=cool_inlet_ring_h__m,
															percentage_travelled=0.125,
															mdot__kg_s=inputParameters.coolantMassFlowRate__kg_per_s/4,
															rho__kg_m3=inputParameters.coolant.density__kg_per_m3(T, p),
															mu__Pas=inputParameters.coolant.viscosity__Pa_s(T, p),
															n=4,
															label="Inlet distribution")
	# turn into the cooling channels
	pressure_circuit.add_90deg_turn_annulus(r_inner__m=rs__m[0]+inputParameters.wallThickness__m,
													r_outer__m=cool_inlet_ring_r__m,
													mdot__kg_s=inputParameters.coolantMassFlowRate__kg_per_s,
													rho__kg_m3=inputParameters.coolant.density__kg_per_m3(T, p),
													label="90 deg turn into channels")
	# area change into cooling channels
	A_in__m2 = np.pi*(pow(cool_inlet_ring_r__m,2)-pow(rs__m[0]+inputParameters.wallThickness__m,2))
	A_out__m2 = inputParameters.channelHeight__m*channel_width__m*inputParameters.channelCount
	Dh_in__m = (cool_inlet_ring_r__m*2)-((rs__m[0]+inputParameters.wallThickness__m)*2)
	pressure_circuit.add_square_area_reduction(A_in__m2=A_in__m2,
													A_out__m2=A_out__m2,
													Dh_in__m=Dh_in__m,
													mdot__kg_s=inputParameters.coolantMassFlowRate__kg_per_s,
													rho__kg_m3=inputParameters.coolant.density__kg_per_m3(T, p),
													mu__Pas=inputParameters.coolant.viscosity__Pa_s(T, p),
													label="Area change into channels")

	pressure_circuit.add_manual_pressure_drop(abs(cooling_data["p_coolant"][-1]-cooling_data["p_coolant"][0]), label="Channels friction")


	# area change out of cooling channels
	A_in__m2 = inputParameters.channelHeight__m*channel_width__m*inputParameters.channelCount
	A_out__m2 = np.pi*(pow(cool_outlet_ring_r__m,2)-pow(rs__m[-1]+inputParameters.wallThickness__m,2))
	Dh_in__m = ((2*inputParameters.channelHeight__m*channel_width__m)/(inputParameters.channelHeight__m+channel_width__m))*inputParameters.channelCount
	pressure_circuit.add_square_area_expansion(A_in__m2=A_in__m2,
													A_out__m2=A_out__m2,
													Dh_in__m=Dh_in__m,
													mdot__kg_s=inputParameters.coolantMassFlowRate__kg_per_s,
													rho__kg_m3=inputParameters.coolant.density__kg_per_m3(T, p),
													mu__Pas=inputParameters.coolant.viscosity__Pa_s(T, p),
													label="Area change out of channels")
	# turn out of the cooling channels
	pressure_circuit.add_90deg_turn_annulus(r_inner__m=rs__m[-1]+inputParameters.wallThickness__m,
													r_outer__m=cool_outlet_ring_r__m,
													mdot__kg_s=inputParameters.coolantMassFlowRate__kg_per_s,
													rho__kg_m3=inputParameters.coolant.density__kg_per_m3(T, p),
													label="90 deg turn out of channels")
	# outlet distribution ring
	pressure_circuit.add_circular_pipe_rectangular_box(r__m=rs__m[-1]+inputParameters.wallThickness__m,
															w__m=cool_outlet_ring_w__m,
															h__m=cool_outlet_ring_h__m,
															percentage_travelled=0.125,
															mdot__kg_s=inputParameters.coolantMassFlowRate__kg_per_s/4,
															rho__kg_m3=inputParameters.coolant.density__kg_per_m3(T, p),
															mu__Pas=inputParameters.coolant.viscosity__Pa_s(T, p),
															n=4,
															label="Oulet distribution")

	# outlet pipe
	pressure_circuit.add_straight_pipes_circ(r__m=cool_outlet_pipe_d__m/2,
													n_pipes=cool_outlet_pipe_n,
													len__m=cool_outlet_pipe_len__m,
													mdot__kg_s=inputParameters.coolantMassFlowRate__kg_per_s/4,
													rho__kg_m3=inputParameters.coolant.density__kg_per_m3(T, p),
													mu__Pas=inputParameters.coolant.viscosity__Pa_s(T, p),
													label="Outlet pipe")

	# cooling jacket outlet turn into halo
	pressure_circuit.add_90deg_turn_circular(r__m=cool_outlet_pipe_d__m/2,
													mdot__kg_s=inputParameters.coolantMassFlowRate__kg_per_s/4,
													rho__kg_m3=inputParameters.coolant.density__kg_per_m3(T, p),
													n=3,
													label="Jacket outlet into halo")
													
	# halo distribution ring
	pressure_circuit.add_circular_pipe_circular(r__m=cool_halo_ring_d__m/2,
														rp__m=cool_halo_tube_d__m/2,
														percentage_travelled=0.5,
														mdot__kg_s=inputParameters.coolantMassFlowRate__kg_per_s,
														rho__kg_m3=inputParameters.coolant.density__kg_per_m3(T, p),
														mu__Pas=inputParameters.coolant.viscosity__Pa_s(T, p),
														n=1,
														label="Halo distribution outlet")
														
	# halo outlet pipe 
	pressure_circuit.add_straight_pipes_circ(r__m=cool_halo_pipe_d__m/2,
													n_pipes=1,
													len__m=cool_halo_outlet_len__m,
													mdot__kg_s=inputParameters.coolantMassFlowRate__kg_per_s,
													rho__kg_m3=inputParameters.coolant.density__kg_per_m3(T, p),
													mu__Pas=inputParameters.coolant.viscosity__Pa_s(T, p),
													label="Halo outlet pipe")

	# return pipe
	pressure_circuit.add_straight_pipes_circ(r__m=cool_return_d__m/2,
													n_pipes=1,
													len__m=cool_return_len__m,
													mdot__kg_s=inputParameters.coolantMassFlowRate__kg_per_s,
													rho__kg_m3=inputParameters.coolant.density__kg_per_m3(T, p),
													mu__Pas=inputParameters.coolant.viscosity__Pa_s(T, p),
													label="Return pipe")
	# into bucket
	pressure_circuit.add_orifice(r__m=cool_return_d__m/2,
										n_ori=1,
										CD=cool_CD,
										mdot__kg_s=inputParameters.coolantMassFlowRate__kg_per_s,
										rho__kg_m3=inputParameters.coolant.density__kg_per_m3(T, p),
										label="Into bucket")

	# complete pressure drop calcs
	pressure_circuit.calc_total_pressure_drop()

def main():
	# Generate the finetuned engine parameters
	rafiki = Rafiki(InputParameters, displayBranding=False)
	rafikiOutput = rafiki.performConicalNozzleSizing()
	
	displayFineTunedParameters(rafikiOutput)

	# Calculate the burn time
	maximumBurnTime__s = calculateBurnTime(FineTunedRafikiOutputParameters)

	# Perform the cooling calculations
	yeti = Yeti(CoolingParameters, displayBranding=False)
	yetiOutput = yeti.performCoolingCalculations()

	# Perform the injector orifice calculations
	calculateInjectorGeometry(	FineTunedRafikiOutputParameters.mf_dot__kg_per_s.value,
								FineTunedRafikiOutputParameters.mo_dot__kg_per_s.value,
								FuelDensity__kg_per_m3,
								OxidiserDensity__kg_per_m3,
								Injector	)

	# Perform the pressure drop calculations
	calculateCoolantPressureDrop(FineTunedRafikiOutputParameters, CoolingParameters, yetiOutput)

	# TODO: Write this function...
	#calculateInjectorPressureDrop()

	# Display output
	print("\n*** RESULTS ***")
	print("PERFORMANCE:")
	print(f"\tSPECIFIC IMPULSE: %.2f s" % FineTunedRafikiOutputParameters.Isp__s.value)
	print(f"\tBURN TIME: %.0f s" % maximumBurnTime__s)

	print("\nDIMENSIONS:")
	print(f"\tCHAMBER COORDINATES: %s MM" % str(FineTunedRafikiOutputParameters.engineGeometry__mm.value))
	print(f"\tFUEL ORIFICE DIAMETER: %.1f MM" % Injector.fuelOrificeDiameter__mm)
	print(f"\tOXIDISER ORIFICE DIAMETER: %.1f MM" % Injector.oxidiserOrificeDiameter__mm)

	print("\nREQUIREMENTS:")
	print(f"\tFUEL MASS FLOW RATE: %.2f kg/s" % FineTunedRafikiOutputParameters.mf_dot__kg_per_s.value)
	print(f"\tOXIDISER MASS FLOW RATE: %.2f kg/s" % FineTunedRafikiOutputParameters.mf_dot__kg_per_s.value)
	print("***************")

if __name__ == "__main__":
	main()