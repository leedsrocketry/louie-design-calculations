from rafiki.rafiki import Rafiki
from yeti.yeti import Yeti
from deepdiff import DeepDiff
import numpy as np

RafikiDesignInputParameters = ("Isopropanol", "LOX", 99.0, 100.0, 1.5, 2000, np.linspace(10,20,10), 1)

def fineTuneRafikiOutput(rafikiOutput):
	return {}

def main():
	print("PERFORMING SIZING...")

	sizer = Rafiki(*RafikiDesignInputParameters, displayBranding=False)
	rafikiOutput = sizer.performConicalNozzleSizing()
	fineTunedRafikiOutput = fineTuneRafikiOutput(rafikiOutput)
	
	print("\nRAFIKI OUTPUT VS FINE-TUNED RAFIKI OUTPUT DIFFERENCES: ")
	print(DeepDiff(rafikiOutput, fineTunedRafikiOutput))

	# TODO: Make Yeti like Rafiki...
	'''
	print("\nRUNNING YETI...")
	cooler = Yeti(".", yetiInputFilePath)
	cooler.run_cooling_calcs()
	cooler.plot_data()
	'''

if __name__ == "__main__":
	main()
