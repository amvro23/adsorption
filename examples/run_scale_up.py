from adsorption import Adsorbent_ScaleUp


scale_up = Adsorbent_ScaleUp()

# =================
# Experiment
# =================
ads_cap_mg_g = 125
exp_ads_time_min = 500

# =================
# Pilot unit
# =================
total_FR_mL_min = 2500 
pilot_ads_time_min = 2000

print('')
scale_up.exp_unit(ads_capacity=ads_cap_mg_g, exp_ads_time=exp_ads_time_min)
scale_up.pilot_unit(total_flow_rate=total_FR_mL_min, pilot_ads_time=pilot_ads_time_min)


