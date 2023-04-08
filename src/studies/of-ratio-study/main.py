from rafiki.rafiki import Rafiki, RafikiInputParameters, RafikiOutputParameters
from yeti.yeti import Yeti, YetiInputParameters, Coolant, PressureDropCircuit
import numpy as np
import cusfbamboo as bam
from CoolProp.CoolProp import PropsSI
import matplotlib.pyplot as plt

# TODO: Move the following back to Yeti:
# ***
Water = Coolant(    name="Water",
                    density__kg_per_m3=lambda T, p: PropsSI("DMASS", "T", T, "P", p, "WATER"),
                    viscosity__Pa_s=lambda T, p: PropsSI("VISCOSITY", "T", T, "P", p, "WATER"),
                    conductivity__W_per_m_K=lambda T, p: PropsSI("CONDUCTIVITY", "T", T, "P", p, "WATER"),
                    specificHeatCapacity__J_per_kg_per_K=lambda T, p: PropsSI("CPMASS", "T", T, "P", p, "WATER"),
                    prandtl=lambda T, p: PropsSI("PRANDTL", "T", T, "P", p, "WATER"),
                    laminarRe=2300,
                    turbulentRe=3500    )

C106Copper = bam.Material(E = 109e9, poisson = 0.34, alpha = 17.9e-6, k = 380.0)
Al7049 = bam.Material(E = 70e9, poisson = 0.33, alpha = 23.4e-6, k = 160.0)
Al2219 = bam.Material(E = 73e9, poisson = 0.33, alpha = 22.7e-6, k = 125.0)
StainlessSteel = bam.Material(E = 193e9, poisson = 0.29, alpha = 16e-6, k = 14.0)
#***

RequiredBurnTime__s = 3
OFRatios = np.linspace(0.8, 2, 25)

WallThickness__mm = 1
ChannelHeight__mm = 3.2
ChannelWidth__mm = 2
ChannelCount = 30

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
												coolant=Water,
												wallMaterial=Al7049,
												wallThickness__mm=WallThickness__mm,
												channelHeight__mm=ChannelHeight__mm,
												channelCount=ChannelCount,
												blockageRatio=None,
												coolantInletTemperature__degC=25,
												coolantInletPressure__bar=23,
												coolantMassFlowRate__kg_per_s=5	)

class BlockageRatioGenerator():
	def __init__(self, engineGeometry__mm, wallThickness__mm, channelHeight__mm, channelWidth__mm, channelCount):
		self.wallThickness__mm = wallThickness__mm
		self.channelHeight__mm = channelHeight__mm
		self.channelWidth__mm = channelWidth__mm
		self.channelCount = channelCount

		self.channelRadius__mm = ChannelWidth__mm / 2
		self.xs__mm, self.rs__mm = self.formatEngineGeometry(engineGeometry__mm)
	
	def formatEngineGeometry(self, engineGeometry__mm):
		formattedEngineGeometry__mm = np.array(engineGeometry__mm.value).transpose()

		return formattedEngineGeometry__mm[0], formattedEngineGeometry__mm[1]
	
	# Good luck understand this one...
	# Essentially, we're calculating the flow area in a channel formed by cutting a vertical slot in a cyclinder with a ballnose end mill
	def calculateChannelCSA(self, finExternalRadius__mm):
		theta = 2 * np.arcsin((self.channelWidth__mm / 2) / finExternalRadius__mm)
		sectorArea__mm2 = (theta * np.pi * finExternalRadius__mm**2) / (2 * np.pi)
		triangleArea__mm2 = (self.channelWidth__mm / 2) * np.sqrt((finExternalRadius__mm**2) - ((self.channelWidth__mm**2) / 4))
		segmentArea__mm2 = sectorArea__mm2 - triangleArea__mm2

		rectangleHeight__mm = self.channelHeight__mm - self.channelRadius__mm
		rectangleArea__mm2 = rectangleHeight__mm * self.channelWidth__mm

		semicircleArea__mm2 = (np.pi * self.channelRadius__mm**2) / 2

		channelArea__mm2 = segmentArea__mm2 + rectangleArea__mm2 + semicircleArea__mm2

		return channelArea__mm2
	
	# TODO: Consider moving the following function to Yeti
	def generateBlockageRatio(self):
		def f(x__m):
			x__mm = x__m * 1000
			r__mm = np.interp(x__mm, self.xs__mm, self.rs__mm)

			r_inner__mm = r__mm + (self.wallThickness__mm)
			r_outer__mm = r_inner__mm + (self.channelHeight__mm)
			
			total_area__mm2 = (pow(r_outer__mm, 2) * np.pi) - (pow(r_inner__mm, 2) * np.pi)
			channel_area__mm2 = self.calculateChannelCSA(r__mm)
			total_channel_area__mm2 = channel_area__mm2 * self.channelCount
			blockage_ratio = (total_area__mm2 - total_channel_area__mm2) / total_area__mm2

			#print(channel_area__mm2, blockage_ratio)

			return blockage_ratio
		
		return f
		
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

		# Blockage ratio function needs recalculating for each run as it depends on the run geometry
		blockageRatioGenerator = BlockageRatioGenerator(	engineGeometry__mm=run.outputParameters.engineGeometry__mm,
															wallThickness__mm=WallThickness__mm,
															channelHeight__mm=ChannelHeight__mm,
															channelWidth__mm=ChannelWidth__mm, 
															channelCount=ChannelCount	)
		StudyCoolingParameters.blockageRatio = blockageRatioGenerator.generateBlockageRatio()
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