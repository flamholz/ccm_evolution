import numpy as np

DEFAULT_HCARB_RATIO = 100.0

# Lengths in micron units
CELL_SA = 6  # um^2, BNID 101792
CELL_V = 1.5 # um^3, BNID 114924
SA_V_RATIO = CELL_SA / CELL_V 
C_PER_CELL = 1e10 # BNID 103010

# Empirical membrane permeability values [um/s], see Mangan & Flamholz, PNAS 2016 for refs.
C_PERM = 0.3*1e4
# Permeability of HCO3- is mostly due to H2CO3, see Mangan & Flamholz, PNAS 2016.
H_PERM = 3e-3*1e4*np.power(10, 3.2-7.1)

# Default values for C and H permeabilities.
DEFAULT_ALPHA = C_PERM * SA_V_RATIO
DEFAULT_BETA = H_PERM * SA_V_RATIO

# Concentrations in uM units, assume near neutral pH.
# External CO2 roughly equal to Henry's law equilibrium at 25 C. 
# See relevant chapter of Cell Biology by the Numbers for CO2 solubility.
DEFAULT_COUT = 15    # uM
DEFAULT_PH = 7.1     # pH chosen so Keq = 10, see below.

# Effective pKa between CO2 and HCO3-. 
PKA_EFF = 6.1 

# Calculate the Keq as a function of the pH.
def Keq_pH(p, pKa=PKA_EFF):
    return np.power(10, p-pKa)

# By assuming a single Keq in & out of the cell we are assuming pH 
# equilibrium across the cell membrane. This greatly simplifies the math.
# See SI of Mangan & Flamholz, PNAS 2016 for details on this point. 
DEFAULT_KEQ = Keq_pH(DEFAULT_PH)

# Default value for rubisco carboxylation rate constant assumes 1 uM of 
# a relatively fast rubisco. See accompanying Mathematica file for notes.
DEFAULT_GAMMA = 1

# Default value for CO2 hydration is roughly the spontaneous rate constant at 25 C.
DEFAULT_DELTA = 0.01 # /s

# By default we assume that capacity for H-carboxylation scales with rubisco carboxylation capacity.
DEFAULT_GAMMA_OMEGA_RATIO = 100

# By default there is no drive H uptake.
DEFAULT_CHI = 0

class AutotrophModel(object):
    """Model of dual-limitation by CO2 and HCO3- in autotrophy."""
    
    def __init__(self, 
                 a=DEFAULT_ALPHA, b=DEFAULT_BETA,
                 g=DEFAULT_GAMMA,
                 d=DEFAULT_DELTA, 
                 p=np.NaN, 
                 o=np.NaN, 
                 x=DEFAULT_CHI,
                 c_out=DEFAULT_COUT,
                 k_eq=DEFAULT_KEQ, 
                 hcarb_ratio=DEFAULT_HCARB_RATIO,
                 cell_surface_area=CELL_SA,
                 cell_volume=CELL_V):
        """Initialize the model with all parameters.
        
        See accompanying Mathematica file for model description and steady-state solution.
    
        Args:
            a: alpha, CO2 permeability of membrane accounting for SA/V ratio [/s].
            b: beta, effective HCO3- permeability of cell membrane accounting for SA/V ratio [/s].
            g: gamma, first order rate constant for rubisco [/s]
            d: delta, first order rate constant for CO2 hydration (catalyzed or spont) [/s].
            p: phi, first order rate constant for HCO3- dehydration (catalyzed or spont.) [/s].
                If not set, calculated from delta and k_eq via Haldane relation. 
            o: omega, first order rate constant for H-dependent carboxylation [/s].
                If not set, as a constant fraction of rubisco rate constant. 
            x: chi, first order rate constant for HCO3- uptake [/s].
            c_out: extracellular CO2 concentration [uM].
            k_eq: H/C ratio at equilibrium [unitless].
            hcarb_ratio: ratio of rubisco-derived and H-carboxylation derived C in biomass.  
                used to calculate biomass production from rubisco and H-carb fluxes.
            cell_surface_area: surface area in um^2.
            cell_volume: cell volume in um^3 units.
        """
        self.a = a
        self.b = b
        self.g = g
        self.d = d
        self.p = p
        self.o = o
        self.x = x

        self.c_out = c_out
        self.k_eq = k_eq
        self.h_out = k_eq*c_out
        
        self.hcarb_ratio = hcarb_ratio
        
        self.cell_vol = cell_volume
        self.cell_sa = cell_surface_area
        self.flux_conversion_factor = cell_volume*1e-15*1e-6*6.02e23
        
        if np.isnan(self.p):
            self.p = self.d / self.k_eq
        if np.isnan(self.o):
            self.o = self.g / DEFAULT_GAMMA_OMEGA_RATIO
        
    def C_in(self):
        # short names, local vars
        a, b, g, d = self.a, self.b, self.g, self.d
        p, o, x = self.p, self.o, self.x
        c_out, k_eq = self.c_out, self.k_eq
        
        num = c_out*(k_eq*p*(b+x) + a*(b+p+o))
        den = b*(g+d) + g*p + (g+d)*o + a*(b+p+o)
        return num/den
    
    def H_in(self):
        # short names, local vars
        a, b, g, d = self.a, self.b, self.g, self.d
        p, o, x = self.p, self.o, self.x
        c_out, k_eq = self.c_out, self.k_eq
        
        num = c_out*(a*d + k_eq*(a+g+d)*(b+x))
        den = b*(g+d) + g*p + (g+d)*o + a*(b+p+o)
        return num/den
        
    def rubisco_flux(self):
        # rubisco rate is just gamma * C_in
        return self.g*self.C_in()
    
    def rubisco_flux_C_per_s(self):
        return self.rubisco_flux()*self.flux_conversion_factor 
    
    def hcarb_flux(self):
        # h-carboxylation rate is just omega * H_in
        return self.o * self.H_in()
    
    def hcarb_flux_C_per_s(self):
        # h-carboxylation rate is just omega * H_in
        return self.hcarb_flux()*self.flux_conversion_factor 
    
    def biomass_flux(self):
        return np.minimum(self.rubisco_flux(), self.hcarb_ratio*self.hcarb_flux())
    
    def biomass_flux_C_per_s(self):
        return self.flux_conversion_factor*self.biomass_flux()
    
    def doubling_time_hr(self):
        return C_PER_CELL/(3600*self.biomass_flux_C_per_s())
    
    def growth_rate_hr(self):
        return np.log(2)/self.doubling_time_hr()
    
    def C_leakage(self):
        return -self.a*(self.c_out-self.C_in())
    
    def C_leakage_C_per_s(self):
        return -self.a*(self.c_out-self.C_in())*self.flux_conversion_factor 
    