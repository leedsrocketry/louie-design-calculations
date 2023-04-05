from rafiki.rafiki import Rafiki, InputParameters, OutputParameters
from yeti.yeti import Yeti
from deepdiff import DeepDiff
import numpy as np
from pprint import pprint

RafikiInputParameters = InputParameters(	fuel="Isopropanol",
											oxidiser="LOX",
											fuel_concentration=99,
											oxidiser_concentration=100,
											oxidiser_fuel_mass_ratio=1.5,
											peak_thrust__N=2000,
											chamber_pressure__bar=10,
											ambient_pressure__bar=1	)

# TODO: Put teh correct values in here
FineTunedRafikiOutputParameters = OutputParameters(	epsilon=2.31527158,
													gamma=1.15192012,
													Isp__s=0,
													Tc__K=3092.32065359,
													M__kg_per_kmol=0,
													Pr_exit=0,
													mu_exit__Pa_s=0,
													k_exit__W_per_m_K=0,
													R__J_per_kg_K=0,
													m_dot__kg_per_s=0,
													mf_dot__kg_per_s=0,
													mo_dot__kg_per_s=0,
													At__mm2=0,
													Ae__mm2=0,
													Dc__mm=0,
													Lc__mm=0,
													Dt__mm=0,
													De__mm=0  )


def main():
	print("PERFORMING SIZING...")
	sizer = Rafiki(RafikiInputParameters, displayBranding=False)
	rafikiOutput = sizer.performConicalNozzleSizing()
	sizer.displayOutput(rafikiOutput)

	print("\nRAFIKI OUTPUT VS FINE-TUNED RAFIKI OUTPUT CHANGES: ")
	pprint(DeepDiff(rafikiOutput[0].outputParameters, FineTunedRafikiOutputParameters, ignore_order=True, ignore_numeric_type_changes=True, significant_digits=3)["values_changed"])

	# TODO: Make Yeti like Rafiki...
	'''
	print("\nRUNNING YETI...")
	cooler = Yeti(".", yetiInputFilePath)
	cooler.run_cooling_calcs()
	cooler.plot_data()
	'''

if __name__ == "__main__":
	main()
