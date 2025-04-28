parms.total <- c(             # Kinesis parms
  T0 = 14,                    # ℃, Optimal temperature for movement
  P0 = 0.8,                   # -, Optimal functional response of food for movement
  sigma_T = 5,                # Spread of the function for temperature relative to T0
  sigma_P = 0.3,              # Spread of the function for food relative to P0
  tau_Temp = 2,               # Weighting factors of f(T)
  tau_Food = 10,              # Weighting factors of f(P)
  H1 = 0.75,                  # Height of function for temperature relative to T0/P0
  H2 = 0.9,                   # Height of function for importance of random component
  SSmax = 3,                  # Maximum swimming speed
  
  # Functional response of food - multi-species prey types
  # below is used to 0-100m for meso-zooplankton of juvenile and adult
  kf_early_larvae_micro_zoo = 0.6, # mmolC/m^3, Half-saturation constant (fish early larvae to microzooplankton)
  kf_late_larvae_micro_zoo = 0.1, # mmolC/m^3, Half-saturation constant (fish late larvae to microzooplankton)
  kf_juvenile_micro_zoo = 0.04,    # mmolC/m^3, Half-saturation constant (fish juvenile to microzooplankton)
  kf_adult_micro_zoo = 0.4,      # mmolC/m^3, Half-saturation constant (fish adult to microzooplankton)
  
  kf_early_larvae_meso_zoo = 0.6,  # mmolC/m^3, Half-saturation constant (fish early larvae to mesozooplankton)
  kf_late_larvae_meso_zoo = 0.06,  # mmolC/m^3, Half-saturation constant (fish late larvae to mesozooplankton)
  kf_juvenile_meso_zoo = 0.02,     # mmolC/m^3, Half-saturation constant (fish juvenile to mesozooplankton)
  kf_adult_meso_zoo = 0.3,       # mmolC/m^3, Half-saturation constant (fish adult to mesozooplankton)
  
#  kf_early_larvae_micro_zoo = 0.6, # mmolC/m^3, Half-saturation constant (fish early larvae to microzooplankton)
#  kf_late_larvae_micro_zoo = 0.1, # mmolC/m^3, Half-saturation constant (fish late larvae to microzooplankton)
#  kf_juvenile_micro_zoo = 0.07,    # mmolC/m^3, Half-saturation constant (fish juvenile to microzooplankton)
#  kf_adult_micro_zoo = 0.35,      # mmolC/m^3, Half-saturation constant (fish adult to microzooplankton)
  
#  kf_early_larvae_meso_zoo = 0.6,  # mmolC/m^3, Half-saturation constant (fish early larvae to mesozooplankton)
#  kf_late_larvae_meso_zoo = 0.1,  # mmolC/m^3, Half-saturation constant (fish late larvae to mesozooplankton)
#  kf_juvenile_meso_zoo = 0.07,     # mmolC/m^3, Half-saturation constant (fish juvenile to mesozooplankton)
#  kf_adult_meso_zoo = 0.35,       # mmolC/m^3, Half-saturation constant (fish adult to mesozooplankton)
  
#  kf_early_larvae_micro_zoo = 0.1, # mmolC/m^3, Half-saturation constant (fish early larvae to microzooplankton)
#  kf_late_larvae_micro_zoo = 0.06, # mmolC/m^3, Half-saturation constant (fish late larvae to microzooplankton)
#  kf_juvenile_micro_zoo = 0.168,    # mmolC/m^3, Half-saturation constant (fish juvenile to microzooplankton)
#  kf_adult_micro_zoo = 0.168,      # mmolC/m^3, Half-saturation constant (fish adult to microzooplankton)
  
#  kf_early_larvae_meso_zoo = 0.1,  # mmolC/m^3, Half-saturation constant (fish early larvae to mesozooplankton)
#  kf_late_larvae_meso_zoo = 0.06,  # mmolC/m^3, Half-saturation constant (fish late larvae to mesozooplankton)
#  kf_juvenile_meso_zoo = 0.168,     # mmolC/m^3, Half-saturation constant (fish juvenile to mesozooplankton)
#  kf_adult_meso_zoo = 0.168,       # mmolC/m^3, Half-saturation constant (fish adult to mesozooplankton)
  
#  kf_earlylarvae = 0.05,      # mmolC/m^3, This coefficient should be adjusted by larvae weight-length relationship
#  kf_latelarvae = 0.8,        # mmolC/m^3, This coefficient should be adjusted by larvae weight-length relationship
#  kf_juvenile = 0.052,        # mmolC/m^3, This coefficient should be adjusted by juvenile weight-at-age data
#  kf_adult = 0.025,           # mmolC/m^3, This coefficient should be adjusted by adult weight-at-age and condition index data
#  kf = 0.035,
#  kf=0.168,
  
  # DEB parms
  # initial values
  E0 = 0.5,
  V0 = 0.000315,
  
  V_early_larvae = 0.00017,   # cm^-3, structural volume at first (exogenous) feeding       # L_early_Larvae = 0.36
  V_late_larvae = 0.0037,     # cm^-3, structural volume at late larvae (generated from TL) # L_late_larvae = 1
  V_juvenile = 0.0755,        # cm^-3, structural volume at metamorphosis                   # L_juvenile = L_metamophosis = 2.5
  V_adult = 1.6585,             # cm^-3, structural volume at puberty                       # L_adult = 7.6 cm
  del_larvae = 0.154,         # shape coefficient for larvae
  del_adult = 0.169,          # shape coefficient for adult
  T_A = 9800,                 # K, Arrhenius temperature for anchovy
  T_ref = 289,                # K, (i.e., 16 degree) reference temperature from Pethybridge et al., (2013)
  T_L = 278,                  # K, lower boundary of tolerance range [[K]] parameters from Joey Volwater (master thesis, Royal Netherlands Institute for Sea Research, 2017)
  T_H = 305,                  # K, upper boundary of tolerance range [[K]] parameters from Joey Volwater (master thesis, Royal Netherlands Institute for Sea Research, 2017),
  T_AL = 50000,               # K, rate of decrease of LOWER boundary [[K]] parameters from Joey Volwater (master thesis, Royal Netherlands Institute for Sea Research, 2017),
  T_AH = 100000,              # K, rate of decrease of UPPER boundary [[K]] parameters from Joey Volwater (master thesis, Royal Netherlands Institute for Sea Research, 2017),
  p_xm = 325,                 # J cm^-2 d^-1, Max. feeding rate
  ae = 0.71,                  # -, Assimilation efficiency
  E_G = 4000,                 # J cm^-3, Volume specific cost for structure
  E_m = 2700,                 # J cm^-3, Max. storage density
  p_M = 48,                   # J cm^-3 d^-1, Volume specific maintenance costs
  kap = 0.71,                 # -, Fraction of energy allocated to growth
  TR = 16,                    # degree, Temperature trigger for gamete production
  uv = 19.9 * 1000,           # J/g, Energy density of 1 g structure
  uE = 35.2 * 1000,           # J/g, Energy density of 1 g reserve
  uG = 33.4 * 1000,           # J/g, Energy density of 1 g gametes
  dv = 0.23,                  # g/cm3, Density of structure
  Eegg = 0.15,                # J/egg, energy content of 1 egg
  RFbatch = 625.3,              # eggs/cm3, relative batch fecundity
  gamma_conversion = 4.1,      # -, Wdw to Www conversion
  
  # IBM mortality
  mort_egg_yolksac = 0.5,
  mort_early_larvae = 0.286,
  mort_late_larvae = 0.057,
  mort_juvenile = 0.004,
  mort_adult = 0.003,
  
  max_day_egg_yolksac = 6,
  max_day_early_larvae = 40,
  max_day_late_larvae = 80,
  max_day_juvenile = 280,
  max_day_adult = 1825
)

# colnames for DEB model
name.col <- c(
  'E',
  'V',
  'ER',
  'EGAM',
  'PA',
  'PG',
  'PC',
  'PM',
  'PM2',
  'PR2',
  'PM3',
  'PJ',
  'PR',
  'TL',
  'Wdw',
  'Www',
  'WGON',
  'Ed',
  'Kful',
  'GSI',
  'gam_dev_period',
  'spawn_time',
  'spawn_month',
  'spawn_batch_energy',
  'eggs',
  'life_stage',
  'growth_weight',
  'growth_length'
)