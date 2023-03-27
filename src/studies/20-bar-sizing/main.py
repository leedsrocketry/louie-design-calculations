from rafiki.rafiki import Rafiki
from yeti.yeti import Yeti
from deepdiff import DeepDiff
from csv import DictReader
import pathlib

RafikiInputFile = "rafiki-input.yaml"
RafikiOutputFile = "rafiki-output.csv"
FineTunedRafikiOutputFile = "fine-tuned-rafiki-output.csv"

def main():
	workingDirectory = str(pathlib.Path(__file__).parent.resolve())
	rafikiInputFilePath = f"%s/%s" % (workingDirectory, RafikiInputFile)
	rafikiOutputFilePath = f"%s/%s" % (workingDirectory, RafikiOutputFile)
	fineTunedRafikiOutputFilePath = f"%s/%s" % (workingDirectory, FineTunedRafikiOutputFile)

	print("RUNNING RAFIKI:")
	sizer = Rafiki(".", rafikiInputFilePath, rafikiOutputFilePath, param_study_mode=False)
	sizer.run_sizing_conical_noz(save_data=True)
	
	with open(rafikiOutputFilePath, 'r') as f:
		dictReader = DictReader(f)
		rafikiOutput = list(dictReader)
	
	with open(fineTunedRafikiOutputFilePath, 'r') as f:
		dictReader = DictReader(f)
		fineTunedRafikiOutput = list(dictReader)
	
	print("\nRAFIKI OUTPUT VS FINE-TUNED RAFIKI OUTPUT DIFFERENCES: ")
	print(DeepDiff(rafikiOutput, fineTunedRafikiOutput))

	# TODO: Run yeti...

if __name__ == "__main__":
	main()
