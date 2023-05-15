import os
import sys
import inspect
import numpy as np
from deepdiff import DeepDiff
from pprint import pprint
from rafiki.rafiki import Rafiki, RafikiInputParameters, RafikiOutputParameters
from yeti.yeti import Yeti, YetiInputParameters
from yeti.pressure_drop_circuit import PressureDropCircuit
import yeti.yeti
from injector import TripletInjector

FuelDensity__kg_per_m3 = 786
FuelViscosity__Pa_s = 0.00237
OxidiserDensity__kg_per_m3 = 1141
OxidiserViscosity__Pa_s=0.000197

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

Injector = TripletInjector(dischargeCoefficient=0.65, elementCount=14, pressureDrop__bar=6)

def calculateInjectorGeometry(fuelMassFlowRate__kg_per_s, oxidiserMassFlowRate__kg_per_s, fuelDensity__kg_per_m3, oxdiserDensity__kg_per_m3, injector):
	# Source: http://athena.leedsrocketry.co.uk/doku.php?id=louie_injector_desgin
	fuelArea__mm2 = (fuelMassFlowRate__kg_per_s / (injector.dischargeCoefficient * np.sqrt(2 * fuelDensity__kg_per_m3 * injector.pressureDrop__Pa))) * 10**6
	oxidiserArea__mm2 = (oxidiserMassFlowRate__kg_per_s / (injector.dischargeCoefficient * np.sqrt(2 * oxdiserDensity__kg_per_m3 * injector.pressureDrop__Pa))) * 10**6

	injector.fuelOrificeDiameter__mm = np.sqrt((2 * fuelArea__mm2) / (injector.elementCount * np.pi))
	injector.oxidiserOrificeDiameter__mm = 2 * np.sqrt(oxidiserArea__mm2 / (injector.elementCount * np.pi))

def displayFineTunedParameters(rafikiOutput):
	print("\nRAFIKI OUTPUT VS FINETUNED VALUES: ")
	pprint(DeepDiff(rafikiOutput[0].outputParameters, FineTunedRafikiOutputParameters, ignore_order=True, ignore_numeric_type_changes=True, significant_digits=3)["values_changed"])

def calculateBurnTime(outputParameters):
	# Work out the maximum burn time according to the tank volumes and denisties of the fuel and oxidiser according to AEL"s ICD
	# (we've assumed the LOX tank volume will be no less than that of the NO2 tank)
	maximumFuelBurnTime__s = ((3.5 * 0.786) / outputParameters.mf_dot__kg_per_s.value)
	maximumOxidiserBurnTime__s = ((10 * 1.141) / outputParameters.mo_dot__kg_per_s.value)
	maximumBurnTime__s = min(maximumFuelBurnTime__s, maximumOxidiserBurnTime__s)

	return maximumBurnTime__s

def calculateCoolantPressureDrop(coolingParameters, cooling_data):
	"""
	Pressure drop assumptions:
		- Inlet and outlet distribution rings are pipes curved along a half-circular arc of the
			same radius as the inner chamber diameter + the wall thickness.
		- The coolant takes a 90 degree turn into the cooling channels, resulting in
			a K-value of 1.3 (see config for reference).
		- The same number of delivery and return pipes as number of inlet and outlet
			pipe.
	"""

	T = coolingParameters.coolantInletTemperature__K
	p = coolingParameters.coolantInletPressure__Pa
	
	cool_delivery_len__m = 1 # length of coolant delivery pipe [m]
	cool_delivery_d__m = 0.025 # diameter of coolant delivery pipe [m]
	cool_CD = 0.7 # coolant coefficient of discharge
	cool_inlet_pipe_d__m = 0.012 # coolant inlet pipe diameter [m]
	cool_inlet_pipe_len__m = 0.039 # coolant inlet pipe length [m]
	cool_inlet_ring_r__m = 0.0435 # outer radius of the coolant inlet ring [m]
	cool_inlet_ring_w__m = 0.0212 # inlet ring width (modeled as a rectangular pipe) [m]
	cool_inlet_ring_h__m = 0.021 # inlet ring height (modeled as a rectangular pipe) [m]
	cool_outlet_pipe_d__m = 0.012 # coolant outlet pipefice diameter [m]
	cool_outlet_pipe_len__m = 0.0325 # coolant outlet pipe length [m]
	cool_outlet_pipe_n = 4 # number of coolant outlet pipe
	cool_outlet_ring_r__m = 0.0545 # outer radius of the coolant outlet ring [m]
	cool_outlet_ring_w__m = 0.016 # outlet ring width (modeled as a rectangular pipe) [m]
	cool_outlet_ring_h__m = 0.0065 # outlet ring height (modeled as a rectangular pipe) [m]
	cool_return_len__m = 1 # length of return pipe [m]
	cool_return_d__m = 0.019 # diameter of return pipe [m]
	cool_halo_inlet_len__m = 0.02 # length of halo inlet pipe [m]
	cool_halo_outlet_len__m = 0.02 # length of halo outlet pipe [m]
	cool_halo_pipe_d__m = 0.0127 # diameter of halo inlet pipe [m]
	cool_halo_ring_d__m = 0.205 # diameter of halo ring [m]
	cool_halo_tube_d__m = 0.019 # diameter of halo tube [m]
	channel_width__m = (2/1000)

	print("deseeeseesntity:  ", coolingParameters.coolant.density__kg_per_m3(T, p))
	print("viscosityyyyyyy:  ", coolingParameters.coolant.viscosity__Pa_s(T, p))

	# build the pressure drop circuit
	pressure_circuit = PressureDropCircuit( coolingParameters.coolantMassFlowRate__kg_per_s, coolingParameters.coolant.density__kg_per_m3(T, p), coolingParameters.coolant.viscosity__Pa_s(T, p))

	# delivery pipe
	pressure_circuit.add_straight_pipes_circ(	r__m=cool_delivery_d__m/2,
												n_pipes=1,
												len__m=cool_delivery_len__m,
												label="Delivery pipe"	)
	# halo inlet pipe 
	pressure_circuit.add_straight_pipes_circ(	r__m=cool_halo_pipe_d__m/2,
												n_pipes=1,
												len__m=cool_halo_inlet_len__m,
												label="Halo inlet pipe"	)
										
	# halo distribution ring
	pressure_circuit.add_circular_pipe_circular(	r__m=cool_halo_ring_d__m/2,
													rp__m=cool_halo_tube_d__m/2,
													percentage_travelled=0.5,
													n_pipes=1,
													label="Halo distribution inlet"	)

	# halo turn into cooling jacket inlet
	pressure_circuit.add_90deg_turn_circular(	r__m=cool_inlet_pipe_d__m/2,
												n_pipes=4,
												label="Halo turn into jacket inlet"	)

	# inlet pipe
	pressure_circuit.add_straight_pipes_circ(	r__m=cool_inlet_pipe_d__m/2,
												n_pipes=4,
												len__m=cool_inlet_pipe_len__m,
												label="Inlet pipe"	)

	# inlet distribution ring
	pressure_circuit.add_circular_pipe_rectangular_box(	r__m=coolingParameters.rs__m[0]+coolingParameters.wallThickness__m,
														w__m=cool_inlet_ring_w__m,
														h__m=cool_inlet_ring_h__m,
														percentage_travelled=0.125,
														n_pipes=4,
														label="Inlet distribution"	)
	# turn into the cooling channels
	pressure_circuit.add_90deg_turn_annulus(	r_inner__m=coolingParameters.rs__m[0]+coolingParameters.wallThickness__m,
												r_outer__m=cool_inlet_ring_r__m,
												label="90 deg turn into channels"	)
	# area change into cooling channels
	A_in__m2 = np.pi*(pow(cool_inlet_ring_r__m,2)-pow(coolingParameters.rs__m[0]+coolingParameters.wallThickness__m,2))
	A_out__m2 = coolingParameters.channelHeight__m*channel_width__m*coolingParameters.channelCount
	Dh_in__m = (cool_inlet_ring_r__m*2)-((coolingParameters.rs__m[0]+coolingParameters.wallThickness__m)*2)

	pressure_circuit.add_square_area_reduction(	A_in__m2=A_in__m2,
												A_out__m2=A_out__m2,
												Dh_in__m=Dh_in__m,
												label="Area change into channels"	)

	pressure_circuit.add_manual_pressure_drop(abs(cooling_data["p_coolant"][-1]-cooling_data["p_coolant"][0]), label="Channels friction")

	# area change out of cooling channels
	A_in__m2 = coolingParameters.channelHeight__m*channel_width__m*coolingParameters.channelCount
	A_out__m2 = np.pi*(pow(cool_outlet_ring_r__m,2)-pow(coolingParameters.rs__m[-1]+coolingParameters.wallThickness__m,2))
	Dh_in__m = ((2*coolingParameters.channelHeight__m*channel_width__m)/(coolingParameters.channelHeight__m+channel_width__m))*coolingParameters.channelCount

	pressure_circuit.add_square_area_expansion(	A_in__m2=A_in__m2,
												A_out__m2=A_out__m2,
												Dh_in__m=Dh_in__m,
												label="Area change out of channels"	)
	# turn out of the cooling channels
	pressure_circuit.add_90deg_turn_annulus(	r_inner__m=coolingParameters.rs__m[-1]+coolingParameters.wallThickness__m,
												r_outer__m=cool_outlet_ring_r__m,
												label="90 deg turn out of channels"	)
	# outlet distribution ring
	pressure_circuit.add_circular_pipe_rectangular_box(	r__m=coolingParameters.rs__m[-1]+coolingParameters.wallThickness__m,
														w__m=cool_outlet_ring_w__m,
														h__m=cool_outlet_ring_h__m,
														percentage_travelled=0.125,
														n_pipes=4,
														label="Oulet distribution"	)

	# outlet pipe
	pressure_circuit.add_straight_pipes_circ(	r__m=cool_outlet_pipe_d__m/2,
												n_pipes=cool_outlet_pipe_n,
												len__m=cool_outlet_pipe_len__m,
												label="Outlet pipe"	)

	# cooling jacket outlet turn into halo
	pressure_circuit.add_90deg_turn_circular(	r__m=cool_outlet_pipe_d__m/2,
												n_pipes=4,
												label="Jacket outlet into halo"	)
													
	# halo distribution ring
	pressure_circuit.add_circular_pipe_circular(	r__m=cool_halo_ring_d__m/2,
													rp__m=cool_halo_tube_d__m/2,
													percentage_travelled=0.5,
													n_pipes=1,
													label="Halo distribution outlet"	)
														
	# halo outlet pipe 
	pressure_circuit.add_straight_pipes_circ(	r__m=cool_halo_pipe_d__m/2,
												n_pipes=1,
												len__m=cool_halo_outlet_len__m,
												label="Halo outlet pipe"	)

	# return pipe
	pressure_circuit.add_straight_pipes_circ(	r__m=cool_return_d__m/2,
												n_pipes=1,
												len__m=cool_return_len__m,
												label="Return pipe"	)
	# into bucket
	pressure_circuit.add_orifice(	r__m=cool_return_d__m/2,
									n_ori=1,
									CD=cool_CD,
									label="Into bucket"	)
	
	# complete pressure drop calcs
	totalPressureDrop__bar = pressure_circuit.calc_total_pressure_drop() / 100000

	return pressure_circuit, totalPressureDrop__bar

def calculateInjectorManifoldPressureDrop(fuelMassFlowRate__kg_per_s, oxidiserMassFlowRate__kg_per_s, fuelDensity__kg_per_m3, fuelViscosity__Pa_s, oxidiserDensity__kg_per_m3, oxidiserViscosity__Pa_s):
	fuelCircuit = PressureDropCircuit(fuelMassFlowRate__kg_per_s, fuelDensity__kg_per_m3, fuelViscosity__Pa_s)
	oxidiserCircuit = PressureDropCircuit(oxidiserMassFlowRate__kg_per_s, oxidiserDensity__kg_per_m3, oxidiserViscosity__Pa_s)
	
	# TODO: Write all the stuff here...
	
	fuelPressureDrop__bar = fuelCircuit.calc_total_pressure_drop() / 100000
	oxidiserPressureDrop__bar = oxidiserCircuit.calc_total_pressure_drop() / 100000

	return 0, 0 # fuelPressureDrop__bar, oxidiserPressureDrop__bar

def main():
	# Generate the finetuned engine parameters
	rafiki = Rafiki(InputParameters, displayBranding=False)
	rafikiOutput = rafiki.performConicalNozzleSizing()
	
	displayFineTunedParameters(rafikiOutput)

	# Calculate the burn time
	maximumBurnTime__s = calculateBurnTime(FineTunedRafikiOutputParameters)

	# Save Rafiki data
	rootdir = os.path.dirname(os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe()))))
	data_path = rootdir + "/10-bar-sizing/data/"
	if not os.path.exists(data_path):
		os.makedirs(data_path)
	rafiki.generateCSV(rafikiOutput, rootdir+"10bar_raw.csv")
	rafikiOutput[0].outputParameters = FineTunedRafikiOutputParameters
	rafiki.generateCSV(rafikiOutput, rootdir+"10bar_fine_tuned.csv")

	# Perform the cooling calculations
	yeti = Yeti(CoolingParameters, displayBranding=False)
	yetiCoolingOutput = yeti.performCoolingCalculations()

	# Perform the injector orifice calculations
	calculateInjectorGeometry(	FineTunedRafikiOutputParameters.mf_dot__kg_per_s.value,
								FineTunedRafikiOutputParameters.mo_dot__kg_per_s.value,
								FuelDensity__kg_per_m3,
								OxidiserDensity__kg_per_m3,
								Injector	)

	# Perform the pressure drop calculations
	coolantPressureDropCircuit, coolantPressureDrop__bar = calculateCoolantPressureDrop(CoolingParameters, yetiCoolingOutput)
	fuelPressureDrop__bar, oxidiserPressureDrop__bar = calculateInjectorManifoldPressureDrop(	FineTunedRafikiOutputParameters.mf_dot__kg_per_s.value,
																								FineTunedRafikiOutputParameters.mo_dot__kg_per_s.value,
																								FuelDensity__kg_per_m3,
																								FuelViscosity__Pa_s,
																								OxidiserDensity__kg_per_m3,
																								OxidiserViscosity__Pa_s	)

	fuelPressureDrop__bar += (Injector.pressureDrop__Pa / 100000)
	oxidiserPressureDrop__bar += (Injector.pressureDrop__Pa / 100000)

	# Display output
	print("\n*** RESULTS ***")
	print("PERFORMANCE:")
	print(f"\tSPECIFIC IMPULSE: %.2f s" % FineTunedRafikiOutputParameters.Isp__s.value)
	print(f"\tBURN TIME: %.1f s" % maximumBurnTime__s)

	print("\nDIMENSIONS:")
	print(f"\tCHAMBER COORDINATES: %s MM" % str(FineTunedRafikiOutputParameters.engineGeometry__mm.value))
	print(f"\tFUEL ORIFICE DIAMETER: %.1f MM" % Injector.fuelOrificeDiameter__mm)
	print(f"\tOXIDISER ORIFICE DIAMETER: %.1f MM" % Injector.oxidiserOrificeDiameter__mm)

	print("\nREQUIREMENTS:")
	print(f"\tFUEL MASS FLOW RATE: %.2f kg/s" % FineTunedRafikiOutputParameters.mf_dot__kg_per_s.value)
	print(f"\tOXIDISER MASS FLOW RATE: %.2f kg/s" % FineTunedRafikiOutputParameters.mo_dot__kg_per_s.value)
	print(f"\tCOOLANT PRESSURE DROP: %.2f bar" % coolantPressureDrop__bar)
	print(f"\tFUEL PRESSURE DROP: %.2f bar" % fuelPressureDrop__bar)
	print(f"\tOXIDISER PRESSURE DROP: %.2f bar" % oxidiserPressureDrop__bar)
	print(f"\tLONGDITUDINAL THERMAL EXPANSION: %.2f MM" % yetiCoolingOutput["longditudinal_thermal_expansion__mm"])
	print("***************")
	
	# All plots
	yeti.plot_main_results(yetiCoolingOutput)
	yeti.plot_aux_results(yetiCoolingOutput)
	coolantPressureDropCircuit.plot_pressure_drop()
	

if __name__ == "__main__":
	main()
