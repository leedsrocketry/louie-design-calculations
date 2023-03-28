from rafiki.rafiki import Rafiki
from yeti.yeti import Yeti
from deepdiff import DeepDiff
import numpy as np
from pprint import pprint

RafikiDesignInputParameters = ("Isopropanol", "LOX", 99.0, 100.0, 1.5, 2000, 10, 1)

FineTunedRafikiOutput = {
	'Fuel': np.array(['Isopropanol'], dtype='<U11'),
	'Oxidiser': np.array(['LOX'], dtype='<U3'),
	'Fuel Concentration [%]': np.array([99.]),
	'Oxidiser Concentration [%]': np.array([100.]),
	'Chamber Pressure [bar]': np.array([10]),
	'Epsilon': np.array([2.31527158]),
	'Gamma': np.array([1.15192012]),
	'Combustion Chamber Temperature [K]': np.array([3092.32065359]),
	'Molecular Weight [kg/kmol]': np.array([21.03667293]),
	'Prandtl Number of Exhaust': np.array([0.37877791]),
	'Viscosity of Exhaust [Pa s]': np.array([8.77097014e-05]),
	'Thermal Conductivity of Exhaust [W/m/K]': np.array([0.58941253]),
	'Total Mass Flow Rate [kg/s]': np.array([0.90483487]),
	'Fuel Mass Flow Rate [kg/s]': np.array([0.36193395]),
	'Oxidiser Mass Flow Rate [kg/s]': np.array([0.54290092]),
	'Throat Diameter [mm]': np.array([45]),
	'Exit Diameter [mm]': np.array([67.92992062]), 
	'Chamber Diameter [mm]': np.array([100]),
	'Chamber Length [mm]': np.array([150]),
}


def main():
	print("PERFORMING SIZING...")
	sizer = Rafiki(*RafikiDesignInputParameters, displayBranding=False)
	rafikiOutput = sizer.performConicalNozzleSizing(displayOutput=False)

	print("\nRAFIKI OUTPUT VS FINE-TUNED RAFIKI OUTPUT CHANGES: ")
	pprint(DeepDiff(rafikiOutput, FineTunedRafikiOutput, ignore_order=True, significant_digits=3)["type_changes"])

	# TODO: Make Yeti like Rafiki...
	'''
	print("\nRUNNING YETI...")
	cooler = Yeti(".", yetiInputFilePath)
	cooler.run_cooling_calcs()
	cooler.plot_data()
	'''

if __name__ == "__main__":
	main()
