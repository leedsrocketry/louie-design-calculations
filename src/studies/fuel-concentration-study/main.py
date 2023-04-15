import numpy as np
import matplotlib.pyplot as plt
from rafiki.rafiki import Rafiki, RafikiInputParameters

FuelConcentrations = np.linspace(20, 100, 30)

StudyInputParameters = RafikiInputParameters(	fuel="Isopropanol",
												oxidiser="LOX",
												fuel_concentration=FuelConcentrations,
												oxidiser_concentration=100,
												oxidiser_fuel_mass_ratio=1.2,
												peak_thrust__N=2000,
												chamber_pressure__bar=10,
												ambient_pressure__bar=1	)

RowCount = 3
ColumnCount = 4

class Plot():
	def __init__(self, fieldName, yLabel, xLabel, colour):
		self.fieldName = fieldName
		self.xLabel = xLabel
		self.yLabel = yLabel
		self.colour = colour

XLabel = "Fuel Concentration (%)"

Plots = [	Plot("epsilon", "Epsilon", XLabel, "r"),
			Plot("gamma", "Gamma", XLabel, "r"),
			Plot("Isp__s", "Specific Impulse (s)", XLabel, "r"),
			Plot("Tc__K", "Combustion Chamber Temperature (K)", XLabel, "r"),
			Plot("m_dot__kg_per_s", "Total Mass Flow Rate (kg/s)", XLabel, "b"),
			Plot("mf_dot__kg_per_s", "Fuel Mass Flow Rate (kg/s)", XLabel, "b"),
			Plot("mo_dot__kg_per_s", "Oxidiser Mass Flow Rate (kg/s)", XLabel, "b"),
			Plot("De__mm", "Exit Diameter (mm)", XLabel, "g"),
			Plot("Dt__mm", "Throat Diameter (mm)", XLabel, "g"),
			Plot("Dc__mm", "Chamber Diameter (mm)", XLabel, "g"),
			Plot("Lc__mm", "Chamber Length (mm)", XLabel, "g")	]

def plotOutputs(outputs):
	fig = plt.figure()
	fig.suptitle("IPA-LOX Fuel Concentration Parametric Study")
	fig.tight_layout()
	fig.subplots_adjust(wspace=0.4)
	fig.subplots_adjust(hspace=0.5)
	fig.subplots_adjust(top=0.85)

	for i in range(0, len(outputs)):
		plot = Plots[i]
		values = outputs[i]

		# TODO: See of I can suimp;y this
		ax = fig.add_subplot(RowCount, ColumnCount, (i + 1))
		ax.plot(FuelConcentrations, values, plot.colour)

		plt.grid()
		plt.xlabel(plot.xLabel)
		plt.ylabel(plot.yLabel)

	plt.show()

def main():
	rafiki = Rafiki(StudyInputParameters, displayBranding=False)
	rafikiOutput = rafiki.performConicalNozzleSizing()
	
	values = np.empty((len(Plots), len(rafikiOutput)))

	for x in range(0, len(rafikiOutput)):
		outputParameters = rafikiOutput[x].outputParameters

		for y in range(0, len(Plots)):
			plot = Plots[y]
			values[y][x] = getattr(outputParameters, plot.fieldName).value
	
	plotOutputs(values)

if __name__ == "__main__":
	main()