import ephem

class LCO(ephem.Observer):
    def __init__(self):
        super(ephem.Observer, self).__init__()
        self.name = 'Las Campanas Observatory'
        self.lon  = '-70:42:03.06'
        self.lat  = '-29:00:38.65'
        self.elevation = 2285

class Gemini(ephem.Observer):
    def __init__(self):
        super(ephem.Observer, self).__init__()
        self.name = 'Cerro Pachon (Gemini)'
        self.lon  = '-70.73659'
        self.lat  = '-30.24073'
        self.elevation = 2722

locations = dict(LCO         = LCO,
                 Gemini      = Gemini)
                 #basic_ep      = BasicKernelEP,
                 #periodic      = PeriodicKernel,
                 #quasiperiodic = QuasiPeriodicKernel,
                 #quasiperiodic_ep = QuasiPeriodicKernelEP)
