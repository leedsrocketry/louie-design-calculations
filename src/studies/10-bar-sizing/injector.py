class TripletInjector():
	def __init__(self, dischargeCoefficient, elementCount, pressureDrop__bar, fuelOrificeDiameter__mm=None, oxidiserOrificeDiameter__mm=None):
		self.dischargeCoefficient = dischargeCoefficient
		self.elementCount = elementCount
		self.pressureDrop__Pa = pressureDrop__bar * 100000
		self.fuelOrificeDiameter__mm = fuelOrificeDiameter__mm
		self.oxidiserOrificeDiameter__mm = oxidiserOrificeDiameter__mm