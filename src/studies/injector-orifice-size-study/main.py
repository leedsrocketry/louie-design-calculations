import numpy as np
from injector import TripletInjector
from prettytable import *

ElementCounts = np.arange(8, 25, 1)

Injector = TripletInjector(dischargeCoefficient=0.65, elementCount=ElementCounts, pressureDrop__bar=6)

# Using IPA and LOX densities respecitivey
FuelDensity__kg_per_m3 = 786
OxidiserDensity__kg_per_m3 = 1141

FuelMassFlowRate__kg_per_s=0.431,
OxidiserMassFlowRate__kg_per_s=0.518,

def calculateInjectorGeometry(mf_dot__kg_per_s, mo_dot__kg_per_s, fuelDensity__kg_per_m3, oxdiserDensity__kg_per_m3, injector):
	# Source: http://athena.leedsrocketry.co.uk/doku.php?id=louie_injector_desgin
	fuelArea__mm2 = (mf_dot__kg_per_s / (injector.dischargeCoefficient * np.sqrt(2 * fuelDensity__kg_per_m3 * injector.pressureDrop__Pa))) * 10**6
	oxidiserArea__mm2 = (mo_dot__kg_per_s / (injector.dischargeCoefficient * np.sqrt(2 * oxdiserDensity__kg_per_m3 * injector.pressureDrop__Pa))) * 10**6

	injector.fuelOrificeDiameter__mm = np.sqrt((2 * fuelArea__mm2) / (injector.elementCount * np.pi))
	injector.oxidiserOrificeDiameter__mm = 2 * np.sqrt(oxidiserArea__mm2 / (injector.elementCount * np.pi))

def PrintStyledTable(table):
		table.hrules = NONE
		table.align = "l"
		table.float_format = ".1" # Watch out! We are rounding the output to 1 dp!
		
		print(table)

def main():
	calculateInjectorGeometry(	FuelMassFlowRate__kg_per_s,
								OxidiserMassFlowRate__kg_per_s,
								FuelDensity__kg_per_m3,
								OxidiserDensity__kg_per_m3,
								Injector	)
	
	table = PrettyTable()
	table.field_names = ["Element Count", "Fuel Orifice Diameter (mm)", "Oxidiser Orifice Diameter (mm)"]

	for i in range(0, len(Injector.fuelOrificeDiameter__mm)):
		table.add_row([ElementCounts[i], Injector.fuelOrificeDiameter__mm[i], Injector.oxidiserOrificeDiameter__mm[i]])

	PrintStyledTable(table)

if __name__ == "__main__":
	main()