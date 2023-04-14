import numpy as np
import cusfbamboo as bam
from CoolProp.CoolProp import PropsSI
import matplotlib.pyplot as plt
from rafiki.rafiki import Rafiki, RafikiInputParameters
from yeti.yeti import Yeti, YetiInputParameters
import yeti.yeti

RequiredBurnTime__s = 5
OFRatios = np.linspace(1, 2, 15)

StudyInputParameters = RafikiInputParameters(	fuel="Isopropanol",
												oxidiser="LOX",
												fuel_concentration=99,
												oxidiser_concentration=100,
												oxidiser_fuel_mass_ratio=OFRatios,
												peak_thrust__N=2000,
												chamber_pressure__bar=10,
												ambient_pressure__bar=1	)

StudyCoolingParameters = YetiInputParameters(	chamberPressure__bar=StudyInputParameters.chamber_pressure__bar.value,
												rafikiOutput=None,
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
		
def main():
	rafiki = Rafiki(StudyInputParameters, displayBranding=False)
	rafikiOutput = rafiki.performConicalNozzleSizing()

	runCount = len(rafikiOutput)
	maximumAverageWallTemperatures__K = np.full(runCount, np.nan)
	specificImpulses__s = np.full(runCount, np.nan)
	maximumTangentialStresses__MPa = np.full(runCount, np.nan)

	for i in range(0, runCount):
		print(f"\nPROCESSING RUN %i..." % i)
		run = rafikiOutput[i]

		# If the required fuel mass flow rate for this run would be greater than AEL's test rig can provide, move on
		if run.outputParameters.mf_dot__kg_per_s.value > 0.5:
			print(f"IGNORING RUN DUE TO EXCESSIVE FUEL MASS FLOW RATE")
			continue

		# If the required fuel mass flow rate for this run would be greater than AEL's test rig can provide, move on
		# (we've assumed the LOX max mass flow rate will be no less than that of the NO2 supply)
		if run.outputParameters.mo_dot__kg_per_s.value > 2.3:
			print(f"IGNORING RUN DUE TO EXCESSIVE OXIDISER MASS FLOW RATE")
			continue

		# Work out the maximum burn time according to the tank volumes and denisties of the fuel and oxidiser according to AEL's ICD
		# (we've assumed the LOX tank volume will be no less than that of the NO2 tank)
		maximumFuelBurnTime__s = ((3.5 * 0.786) / run.outputParameters.mf_dot__kg_per_s.value)
		maximumOxidiserBurnTime__s = ((10 * 1.141) / run.outputParameters.mo_dot__kg_per_s.value)
		maximumBurnTime__s = min(maximumFuelBurnTime__s, maximumOxidiserBurnTime__s)

		# If the implied burn time for this run would be less than the requirement, move on
		if maximumBurnTime__s < RequiredBurnTime__s:
			print(f"IGNORING RUN DUE TO INSUFFICIENT BURN TIME")
			continue

		# Rafiki output parameters need updating for each run	
		StudyCoolingParameters.rafikiOutput = run.outputParameters

		# Perform the cooling calculations
		yeti = Yeti(StudyCoolingParameters, displayBranding=False)
		yetiOutput = yeti.performCoolingCalculations()

		calculatedTemperatures__K = np.array(yetiOutput["T"])
		tangentialStresses__MPa = np.array(yetiOutput["sigma_t_max"]) / 10**6

		wallCoolantContactTemperature_K = max(calculatedTemperatures__K[:, 1])
		wallExhaustContactTemperature_K = max(calculatedTemperatures__K[:, 2])

		maximumAverageWallTemperature__K = (wallCoolantContactTemperature_K + wallExhaustContactTemperature_K) / 2
		
		maximumTangentialStress__MPa = max(tangentialStresses__MPa[:, 0])
		
		maximumAverageWallTemperatures__K[i] = maximumAverageWallTemperature__K
		specificImpulses__s[i] = run.outputParameters.Isp__s.value
		maximumTangentialStresses__MPa[i] = maximumTangentialStress__MPa

		print(f"\tO/F MASS RATIO: %f" % run.inputParameters.oxidiser_fuel_mass_ratio.value)
		print(f"\tFUEL MASS FLOW RATE: %f kg/s" % run.outputParameters.mf_dot__kg_per_s.value)
		print(f"\tOXIDISER MASS FLOW RATE: %f kg/s" % run.outputParameters.mf_dot__kg_per_s.value)
		print(f"\tBURN TIME: %f s" % maximumBurnTime__s)
		print(f"\tMAXIMUM AVERAGE WALL TEMPERATURE: %f K" % maximumAverageWallTemperatures__K[i])
		print(f"\tSPECIFIC IMPULSE: %f s" % specificImpulses__s[i])
		print(f"\tMAXIMUM TANGENTIAL STRESS: %f MPa" % maximumTangentialStresses__MPa[i])
	
	print("\nPLOTTING...")

	fig, ax1 = plt.subplots()
	ax2 = ax1.twinx()
	ax3 = ax1.twinx()

	fig.subplots_adjust(right=0.75)
	ax3.spines.right.set_position(("axes", 1.2))

	p1, = ax1.plot(OFRatios, maximumAverageWallTemperatures__K, "r", label="Maximum Average Wall Temperature")
	p2, = ax2.plot(OFRatios, specificImpulses__s, "g", label="Specific Impulse")
	p3, = ax3.plot(OFRatios, maximumTangentialStresses__MPa, "b-", label="Maximum Tangential Stress")

	ax1.yaxis.label.set_color(p1.get_color())
	ax2.yaxis.label.set_color(p2.get_color())
	ax3.yaxis.label.set_color(p3.get_color())

	tkw = dict(size=4, width=1.5)
	ax1.tick_params(axis='x', **tkw)
	ax1.tick_params(axis='y', colors=p1.get_color(), **tkw)
	ax2.tick_params(axis='y', colors=p2.get_color(), **tkw)
	ax3.tick_params(axis='y', colors=p3.get_color(), **tkw)

	ax1.set_xlabel("O/F Mass Ratio")

	ax1.set_ylabel("Temperature [K]")
	ax2.set_ylabel("Specific Impulse [s]")
	ax3.set_ylabel("Stress [MPa]")

	ax1.legend(handles=[p1, p2, p3])

	plt.show()

if __name__ == "__main__":
	main()