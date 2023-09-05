import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from sklearn.metrics import mean_squared_error,r2_score
import datetime
#plotyear = 2016
#icemeltstart = datetime.datetime(2016,6,26)
#icemeltend = datetime.datetime(2016,8,29)


plotyear = 2009
icemeltstart = datetime.datetime(2009,7,5)
icemeltend = datetime.datetime(2009,8,16)

#load_reindexed_dat
ice_albedo = 0.25
aws_gaps = pd.read_csv('/home/shl/OneDrive/projects/aws_processing_v1.0/data_v1.0/gem_database/2022/zac_l_2008_2022_QC_final.csv', index_col = 0, parse_dates = True)

aws_gaps['SW_net'] = aws_gaps.dsr_corr - aws_gaps.usr_corr
aws_gaps['LW_net'] = aws_gaps.dlr - aws_gaps.ulr
aws_gaps['ice_ablation'] = aws_gaps.ice_ablation.values
aws_interpolated = aws_gaps.interpolate()
#aws_interpolated = aws_interpolated
albedo_day = aws_interpolated['albedo'].loc[str(plotyear)].resample('D').median()
#albedo_day.plot()
# The coefficients to calibrate the melt model
roughness_div = 100

aws = aws_interpolated

roughness_w = 0.0001 
zw1= str(roughness_w)
# Calculate turbulent heat fluxes
rho_0 = 1.29 # kg/m3 air density at P_0
cp = 1005 #J/(Kg K), specific heat of dry air at constant pressure
kappa = 0.41 # von karmans constant
P_0 = 101325 # Pa standard atm pressure 
z = 2.7 #m instrument height
Ls = 2.849*10**6 # J/Kg latent heat of sublimation
Lv = 2.514*10**6 # J/kg latent heat of vaporization
e0 = 611 # Pa Vapour pressure at the melting surface

transfer_coeff_heat = 1.5/1000 # Only for the bulk method
transfer_coeff_vapor = 1.5/1000 # Only for the bulk method

T2 = aws['t_u'] # temperature at 2m
T0 = aws['t_surf']
# calculating vapour pressure at 2m e2
RH = aws['rh_u_corr']
e2 = 611*np.exp(17.27*T2/(243.04+T2))*RH/100 

P = aws['p_u']*100 # atm pressure from hPa to Pa
u2 = aws['wspd']# wind speed at 2m

z0w = roughness_w # roughness parameters for logarithmic profiles of wind 
z0t = z0w/roughness_div # roughness parameters for logarithmic profiles of temp
z0e = z0t # roughness parameters for logarithmic profiles of vapour pressure

L = T2*0+Lv
L[(e2-e0 < 0)]= Ls
L[(e2-e0 >= 0) & (T0<0) ]= Ls


sensible = P/P_0*cp*(kappa**2)*rho_0*u2*T2/(np.log(z/z0w)*np.log(z/z0t))
#latent = L*(0.623/P_0)*(kappa**2/(np.log(z/z0w)*np.log(z/z0e)))*u2*(e2-e0)

latent = 0.623*L*kappa**2*(rho_0/P_0)*(u2*(e2-e0)/(np.log(z/z0w)*np.log(z/z0e)))


#L = T2*0+Lv
#L[(e2-e0 < 0)]= Ls
#L[(e2-e0 >= 0) & (T0<0) ]= Ls
#Ch=transfer_coeff_heat
#Ce=transfer_coeff_vapor
#
#sensible=0.0129*Ch*P*u2*(T2) #Eh
#latent=0.622*rho_0*L*Ce*u2*(e2-e0)/P_0 #Ee
aws = aws.assign(sensible=sensible)
aws = aws.assign(latent=latent)


# Calculate the energy available for melt 
Emelt = aws['SW_net']+aws['LW_net']+aws['sensible']+aws['latent']
Emelt[Emelt<0] = 0

# Calculate the amount of ice that can be melted with the energy available for melt
Lf = 334000 # J/kg
melt = Emelt*3600*(1)/(Lf*1000) # Qmelt/Lf * number of seconds per timestep / density of water = m WE/ timestep
aws['melt, z0w='+str(roughness_w)] = melt
aws['Emelt, z0w='+str(roughness_w)] = Emelt

roughness_w = 0.001
zw2= str(roughness_w)
# Calculate turbulent heat fluxes
rho_0 = 1.29 # kg/m3 air density at P_0
cp = 1005 #J/(Kg K), specific heat of dry air at constant pressure
kappa = 0.41 # von karmans constant
P_0 = 101325 # Pa standard atm pressure 
z = 2.7 #m instrument height
Ls = 2.849*10**6 # J/Kg latent heat of sublimation
Lv = 2.514*10**6 # J/kg latent heat of vaporization
e0 = 611 # Pa Vapour pressure at the melting surface

transfer_coeff_heat = 1.5/1000 # Only for the bulk method
transfer_coeff_vapor = 1.5/1000 # Only for the bulk method

T2 = aws['t_u'] # temperature at 2m
T0 = aws['t_surf']
# calculating vapour pressure at 2m e2
RH = aws['rh_u_corr']
e2 = 611*np.exp(17.27*T2/(243.04+T2))*RH/100 

P = aws['p_u']*100 # atm pressure from hPa to Pa
u2 = aws['wspd']# wind speed at 2m

z0w = roughness_w # roughness parameters for logarithmic profiles of wind 
z0t = z0w/roughness_div # roughness parameters for logarithmic profiles of temp
z0e = z0t # roughness parameters for logarithmic profiles of vapour pressure

L = T2*0+Lv
L[(e2-e0 < 0)]= Ls
L[(e2-e0 >= 0) & (T0<0) ]= Ls


sensible = P/P_0*cp*(kappa**2)*rho_0*u2*T2/(np.log(z/z0w)*np.log(z/z0t))
#latent = L*(0.623/P_0)*(kappa**2/(np.log(z/z0w)*np.log(z/z0e)))*u2*(e2-e0)

latent = 0.623*L*kappa**2*(rho_0/P_0)*(u2*(e2-e0)/(np.log(z/z0w)*np.log(z/z0e)))


#L = T2*0+Lv
#L[(e2-e0 < 0)]= Ls
#L[(e2-e0 >= 0) & (T0<0) ]= Ls
#Ch=transfer_coeff_heat
#Ce=transfer_coeff_vapor
#
#sensible=0.0129*Ch*P*u2*(T2) #Eh
#latent=0.622*rho_0*L*Ce*u2*(e2-e0)/P_0 #Ee
aws = aws.assign(sensible=sensible)
aws = aws.assign(latent=latent)


# Calculate the energy available for melt 
Emelt = aws['SW_net']+aws['LW_net']+aws['sensible']+aws['latent']
Emelt[Emelt<0] = 0

# Calculate the amount of ice that can be melted with the energy available for melt
Lf = 334000 # J/kg
melt = Emelt*3600*(1)/(Lf*1000) # Qmelt/Lf * number of seconds per timestep / density of water = m WE/ timestep
aws['melt, z0w='+str(roughness_w)] = melt
aws['Emelt, z0w='+str(roughness_w)] = Emelt

roughness_w = 0.01
zw3= str(roughness_w)
# Calculate turbulent heat fluxes
rho_0 = 1.29 # kg/m3 air density at P_0
cp = 1005 #J/(Kg K), specific heat of dry air at constant pressure
kappa = 0.41 # von karmans constant
P_0 = 101325 # Pa standard atm pressure 
z = 2.7 #m instrument height
Ls = 2.849*10**6 # J/Kg latent heat of sublimation
Lv = 2.514*10**6 # J/kg latent heat of vaporization
e0 = 611 # Pa Vapour pressure at the melting surface

transfer_coeff_heat = 1.5/1000 # Only for the bulk method
transfer_coeff_vapor = 1.5/1000 # Only for the bulk method

T2 = aws['t_u'] # temperature at 2m
T0 = aws['t_surf']
# calculating vapour pressure at 2m e2
RH = aws['rh_u_corr']
e2 = 611*np.exp(17.27*T2/(243.04+T2))*RH/100 

P = aws['p_u']*100 # atm pressure from hPa to Pa
u2 = aws['wspd']# wind speed at 2m

z0w = roughness_w # roughness parameters for logarithmic profiles of wind 
z0t = z0w/roughness_div # roughness parameters for logarithmic profiles of temp
z0e = z0t # roughness parameters for logarithmic profiles of vapour pressure

L = T2*0+Lv
L[(e2-e0 < 0)]= Ls
L[(e2-e0 >= 0) & (T0<0) ]= Ls


sensible = P/P_0*cp*(kappa**2)*rho_0*u2*T2/(np.log(z/z0w)*np.log(z/z0t))
#latent = L*(0.623/P_0)*(kappa**2/(np.log(z/z0w)*np.log(z/z0e)))*u2*(e2-e0)

latent = 0.623*L*kappa**2*(rho_0/P_0)*(u2*(e2-e0)/(np.log(z/z0w)*np.log(z/z0e)))


#L = T2*0+Lv
#L[(e2-e0 < 0)]= Ls
#L[(e2-e0 >= 0) & (T0<0) ]= Ls
#Ch=transfer_coeff_heat
#Ce=transfer_coeff_vapor
#
#sensible=0.0129*Ch*P*u2*(T2) #Eh
#latent=0.622*rho_0*L*Ce*u2*(e2-e0)/P_0 #Ee
aws = aws.assign(sensible=sensible)
aws = aws.assign(latent=latent)


# Calculate the energy available for melt 
Emelt = aws['SW_net']+aws['LW_net']+aws['sensible']+aws['latent']
Emelt[Emelt<0] = 0

# Calculate the amount of ice that can be melted with the energy available for melt
Lf = 334000 # J/kg
melt = Emelt*3600*(1)/(Lf*1000) # Qmelt/Lf * number of seconds per timestep / density of water = m WE/ timestep
aws['melt, z0w='+str(roughness_w)] = melt
aws['Emelt, z0w='+str(roughness_w)] = Emelt

datapath = '/home/shl/OneDrive/projects/aws_processing_v1.0/data_v1.0/gem_database/2022/'
sr50 = pd.read_csv(datapath+'preQC/zac_l_hour_SR50_stake_height_ice.csv', parse_dates = True, index_col=0)
for year_int in [plotyear]:
    
    year = str(year_int)
    #print(year)
    first = icemeltstart #albedo_day[year][albedo_day<=ice_albedo].index[0]
    second = icemeltend #str(year)+'-September'
    aws_year = aws[first:second] 
    
    observed_ice_ablation = aws_year['ice_ablation']*0.917
    
    observed_ice_ablation = observed_ice_ablation-observed_ice_ablation[0]
    
    fig,ax = plt.subplots(1,1,figsize = (3.5,5))
    ax1 = ax.twinx()
    observed_ice_ablation.rolling(24,center=True).mean().plot(ax = ax, label = 'Observed', linewidth = 3)
    aws_year['melt, z0w='+zw1].cumsum().plot(ax = ax, label = 'z0w = '+zw1, linestyle = '--')
    aws_year['melt, z0w='+zw2].cumsum().plot(ax = ax, label = 'z0w = '+zw2, linestyle = '--')
    aws_year['melt, z0w='+zw3].cumsum().plot(ax = ax, label = 'z0w = '+zw3, linestyle = '--')
    albedo_day.plot(ax=ax1, alpha = 0.5, color = 'gray')
    #aws_year['SR_in'].resample('D').sum().plot(ax = ax1)
    #aws_gaps_year['SR_in'].resample('D').sum().plot(ax=ax1)
    
    ax.legend(loc = 4)
    ax1.set_ylim(0,1)
    ax1.set_ylabel('Albedo', color = 'gray')
    ax1.tick_params(axis='y', labelcolor='gray')
    ax1.spines['right'].set_color('gray')
    ax.set_xlabel('')
    ax.set_ylabel('Ice ablation, m w.e.q.')
    ax1.set_xlim(first,second)
    fig.tight_layout()
    
    observed_ice_melt = observed_ice_ablation.resample('D').mean().diff()
    observed_ice_melt_nans = observed_ice_melt.where(observed_ice_melt>=0.02, np.nan )
    observed_ice_melt_clean = observed_ice_melt_nans.dropna()
    modelled_ice_melt = aws_year['melt, z0w='+zw2].resample('D').sum()
    modelled_ice_melt_clean = modelled_ice_melt.where(observed_ice_melt>=0.02, np.nan )
    modelled_ice_melt_clean = modelled_ice_melt_clean.dropna()
    
    
    
    r2 = 1-(((observed_ice_melt_clean-modelled_ice_melt_clean)**2).sum())/(((observed_ice_melt_clean-observed_ice_melt_clean.mean())**2).sum())
    r2 = np.round(r2*100)/100
    #r2 = r2_score(observed_ice_melt_clean,modelled_ice_melt_clean) 
    rmse = np.sqrt(mean_squared_error(observed_ice_melt_clean, modelled_ice_melt_clean))
    rmse = np.round(rmse*100)/100
    
    fig2, ax2 = plt.subplots(1,1,figsize = (3.5,3))
    ax2.scatter(observed_ice_melt_nans,modelled_ice_melt)
    ax2.plot([0, 0.09], [0, 0.09], color = 'gray', alpha = 0.5)
    ax2.set_ylim([0.0,0.09])
    ax2.set_xlim([0.0,0.09])
    #ax2.text(0.001,0.08,'r² = '+str(r2))
    #ax2.text(0.001,0.075,'rmse = '+str(rmse))
    ax2.set_ylabel('Modelled, m w.eq./day')
    ax2.set_xlabel('Observed, m w.eq./day')
    fig2.tight_layout()
    
    #v1 = sr50.loc[first:second]-sr50.loc[first].mean()
    #v1.plot(ax = ax)
    fig.savefig('figures/point_melt_model/zac_l_model_vs_obs_'+year+'.png', dpi = 300)

    #first = albedo_day[year][albedo_day<=ice_albedo].index[0]
    #second = str(year)+'-September'
    #aws_year = aws[first:second] 
    #aws_gaps_year = aws_gaps[first:second]
    ##print(aws_year.keys())
    #observed_ice_ablation = aws_year['ice_ablation']*0.917
    #observed_ice_ablation = observed_ice_ablation-observed_ice_ablation[0]
    #
    #fig,ax = plt.subplots(1,1)
    #ax1 = ax.twinx()
    #observed_ice_ablation.rolling(24,center=True).mean().plot(ax = ax1, label = 'Observed ice melt', linewidth = 3)
    #aws_year['SW_net'].resample('D').sum().plot(ax = ax, label = 'SR_net', linestyle = '--')
    #aws_year['LW_net'].resample('D').sum().plot(ax = ax, label = 'LR_net', linestyle = '--')
    #aws_year['sensible'].resample('D').sum().plot(ax = ax, label = 'Sensible', linestyle = '--')
    #aws_year['latent'].resample('D').sum().plot(ax = ax, label = 'Latent', linestyle = '--')
    #
    #
    #ax.legend()
    #ax.set_ylabel('Energy for melt, W/m2')
    #
    #fig.savefig('figures/point_melt_model/energy_balance_'+year+'.png', dpi = 300)

    fig2.savefig('figures/point_melt_model/zac_l_scatter_'+year+'.png', dpi = 300)

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
#load_reindexed_dat
ice_albedo = 0.25
aws_gaps = pd.read_csv('/home/shl/OneDrive/projects/aws_processing_v1.0/data_v1.0/gem_database/2022/zac_u_2008_2022_QC_final.csv', index_col = 0, parse_dates = True)
#aws_l = pd.read_csv('/home/shl/OneDrive/projects/aws_processing_v1.0/data_v1.0/gem_database/2022/zac_l_2008_2022_QC_final.csv', index_col = 0, parse_dates = True)

aws_gaps['SW_net'] = aws_gaps.SR_in - aws_gaps.SR_out
aws_gaps['LW_net'] = aws_gaps.LR_in - aws_gaps.LR_out
aws_gaps['ice_ablation'] = aws_gaps.Ablation.values
aws_interpolated = aws_gaps.interpolate()

albedo_day = aws_gaps['albedo'].resample('D').median()

# The coefficients to calibrate the melt model
roughness_div = 100

aws = aws_interpolated

roughness_w = 0.0001 
# Calculate turbulent heat fluxes
rho_0 = 1.29 # kg/m3 air density at P_0
cp = 1005 #J/(Kg K), specific heat of dry air at constant pressure
kappa = 0.41 # von karmans constant
P_0 = 101325 # Pa standard atm pressure 
z = 2.7 #m instrument height
Ls = 2.849*10**6 # J/Kg latent heat of sublimation
Lv = 2.514*10**6 # J/kg latent heat of vaporization
e0 = 611 # Pa Vapour pressure at the melting surface

transfer_coeff_heat = 1.5/1000 # Only for the bulk method
transfer_coeff_vapor = 1.5/1000 # Only for the bulk method

T2 = aws['t_u'] # temperature at 2m
T0 = aws['t_surf']
# calculating vapour pressure at 2m e2
RH = aws['rh_u_corr']
e2 = 611*np.exp(17.27*T2/(243.04+T2))*RH/100 

P = aws['p_u']*100 # atm pressure from hPa to Pa
u2 = aws['wspd']# wind speed at 2m

z0w = roughness_w # roughness parameters for logarithmic profiles of wind 
z0t = z0w/roughness_div # roughness parameters for logarithmic profiles of temp
z0e = z0t # roughness parameters for logarithmic profiles of vapour pressure

L = T2*0+Lv
L[(e2-e0 < 0)]= Ls
L[(e2-e0 >= 0) & (T0<0) ]= Ls


sensible = P/P_0*cp*(kappa**2)*rho_0*u2*T2/(np.log(z/z0w)*np.log(z/z0t))
#latent = L*(0.623/P_0)*(kappa**2/(np.log(z/z0w)*np.log(z/z0e)))*u2*(e2-e0)

latent = 0.623*L*kappa**2*(rho_0/P_0)*(u2*(e2-e0)/(np.log(z/z0w)*np.log(z/z0e)))


#L = T2*0+Lv
#L[(e2-e0 < 0)]= Ls
#L[(e2-e0 >= 0) & (T0<0) ]= Ls
#Ch=transfer_coeff_heat
#Ce=transfer_coeff_vapor
#
#sensible=0.0129*Ch*P*u2*(T2) #Eh
#latent=0.622*rho_0*L*Ce*u2*(e2-e0)/P_0 #Ee
aws = aws.assign(sensible=sensible)
aws = aws.assign(latent=latent)


# Calculate the energy available for melt 
Emelt = aws['SW_net']+aws['LW_net']+aws['sensible']+aws['latent']
Emelt[Emelt<0] = 0

# Calculate the amount of ice that can be melted with the energy available for melt
Lf = 334000 # J/kg
melt = Emelt*3600*(1)/(Lf*1000) # Qmelt/Lf * number of seconds per timestep / density of water = m WE/ timestep
aws['melt, z0w='+str(roughness_w)] = melt
aws['Emelt, z0w='+str(roughness_w)] = Emelt

roughness_w = 0.001
# Calculate turbulent heat fluxes
rho_0 = 1.29 # kg/m3 air density at P_0
cp = 1005 #J/(Kg K), specific heat of dry air at constant pressure
kappa = 0.41 # von karmans constant
P_0 = 101325 # Pa standard atm pressure 
z = 2.7 #m instrument height
Ls = 2.849*10**6 # J/Kg latent heat of sublimation
Lv = 2.514*10**6 # J/kg latent heat of vaporization
e0 = 611 # Pa Vapour pressure at the melting surface

transfer_coeff_heat = 1.5/1000 # Only for the bulk method
transfer_coeff_vapor = 1.5/1000 # Only for the bulk method

T2 = aws['t_u'] # temperature at 2m
T0 = aws['t_surf']
# calculating vapour pressure at 2m e2
RH = aws['rh_u_corr']
e2 = 611*np.exp(17.27*T2/(243.04+T2))*RH/100 

P = aws['p_u']*100 # atm pressure from hPa to Pa
u2 = aws['wspd']# wind speed at 2m

z0w = roughness_w # roughness parameters for logarithmic profiles of wind 
z0t = z0w/roughness_div # roughness parameters for logarithmic profiles of temp
z0e = z0t # roughness parameters for logarithmic profiles of vapour pressure

L = T2*0+Lv
L[(e2-e0 < 0)]= Ls
L[(e2-e0 >= 0) & (T0<0) ]= Ls


sensible = P/P_0*cp*(kappa**2)*rho_0*u2*T2/(np.log(z/z0w)*np.log(z/z0t))
#latent = L*(0.623/P_0)*(kappa**2/(np.log(z/z0w)*np.log(z/z0e)))*u2*(e2-e0)

latent = 0.623*L*kappa**2*(rho_0/P_0)*(u2*(e2-e0)/(np.log(z/z0w)*np.log(z/z0e)))


#L = T2*0+Lv
#L[(e2-e0 < 0)]= Ls
#L[(e2-e0 >= 0) & (T0<0) ]= Ls
#Ch=transfer_coeff_heat
#Ce=transfer_coeff_vapor
#
#sensible=0.0129*Ch*P*u2*(T2) #Eh
#latent=0.622*rho_0*L*Ce*u2*(e2-e0)/P_0 #Ee
aws = aws.assign(sensible=sensible)
aws = aws.assign(latent=latent)


# Calculate the energy available for melt 
Emelt = aws['SW_net']+aws['LW_net']+aws['sensible']+aws['latent']
Emelt[Emelt<0] = 0

# Calculate the amount of ice that can be melted with the energy available for melt
Lf = 334000 # J/kg
melt = Emelt*3600*(1)/(Lf*1000) # Qmelt/Lf * number of seconds per timestep / density of water = m WE/ timestep
aws['melt, z0w='+str(roughness_w)] = melt
aws['Emelt, z0w='+str(roughness_w)] = Emelt

roughness_w = 0.01
# Calculate turbulent heat fluxes
rho_0 = 1.29 # kg/m3 air density at P_0
cp = 1005 #J/(Kg K), specific heat of dry air at constant pressure
kappa = 0.41 # von karmans constant
P_0 = 101325 # Pa standard atm pressure 
z = 2.7 #m instrument height
Ls = 2.849*10**6 # J/Kg latent heat of sublimation
Lv = 2.514*10**6 # J/kg latent heat of vaporization
e0 = 611 # Pa Vapour pressure at the melting surface

transfer_coeff_heat = 1.5/1000 # Only for the bulk method
transfer_coeff_vapor = 1.5/1000 # Only for the bulk method

T2 = aws['t_u'] # temperature at 2m
T0 = aws['t_surf']
# calculating vapour pressure at 2m e2
RH = aws['rh_u_corr']
e2 = 611*np.exp(17.27*T2/(243.04+T2))*RH/100 

P = aws['p_u']*100 # atm pressure from hPa to Pa
u2 = aws['wspd']# wind speed at 2m

z0w = roughness_w # roughness parameters for logarithmic profiles of wind 
z0t = z0w/roughness_div # roughness parameters for logarithmic profiles of temp
z0e = z0t # roughness parameters for logarithmic profiles of vapour pressure

L = T2*0+Lv
L[(e2-e0 < 0)]= Ls
L[(e2-e0 >= 0) & (T0<0) ]= Ls


sensible = P/P_0*cp*(kappa**2)*rho_0*u2*T2/(np.log(z/z0w)*np.log(z/z0t))
#latent = L*(0.623/P_0)*(kappa**2/(np.log(z/z0w)*np.log(z/z0e)))*u2*(e2-e0)

latent = 0.623*L*kappa**2*(rho_0/P_0)*(u2*(e2-e0)/(np.log(z/z0w)*np.log(z/z0e)))


#L = T2*0+Lv
#L[(e2-e0 < 0)]= Ls
#L[(e2-e0 >= 0) & (T0<0) ]= Ls
#Ch=transfer_coeff_heat
#Ce=transfer_coeff_vapor
#
#sensible=0.0129*Ch*P*u2*(T2) #Eh
#latent=0.622*rho_0*L*Ce*u2*(e2-e0)/P_0 #Ee
aws = aws.assign(sensible=sensible)
aws = aws.assign(latent=latent)


# Calculate the energy available for melt 
Emelt = aws['SW_net']+aws['LW_net']+aws['sensible']+aws['latent']
Emelt[Emelt<0] = 0

# Calculate the amount of ice that can be melted with the energy available for melt
Lf = 334000 # J/kg
melt = Emelt*3600*(1)/(Lf*1000) # Qmelt/Lf * number of seconds per timestep / density of water = m WE/ timestep
aws['melt, z0w='+str(roughness_w)] = melt
aws['Emelt, z0w='+str(roughness_w)] = Emelt

for year_int in [2016,2019]:
    print(year)
    year = str(year_int)
    #print(year)
    first = icemeltstart #albedo_day[year][albedo_day<=ice_albedo].index[0]
    second = icemeltend #str(year)+'-September'
    aws_year = aws[first:second] 
    
    observed_ice_ablation = aws_year['ice_ablation']*0.917
    
    observed_ice_ablation = observed_ice_ablation-observed_ice_ablation[0]
    
    fig,ax = plt.subplots(1,1,figsize = (3.5,5))
    ax1 = ax.twinx()
    observed_ice_ablation.rolling(24,center=True).mean().plot(ax = ax, label = 'Observed', linewidth = 3)
    aws_year['melt, z0w='+zw1].cumsum().plot(ax = ax, label = 'z0w = '+zw1, linestyle = '--')
    aws_year['melt, z0w='+zw2].cumsum().plot(ax = ax, label = 'z0w = '+zw2, linestyle = '--')
    aws_year['melt, z0w='+zw3].cumsum().plot(ax = ax, label = 'z0w = '+zw3, linestyle = '--')
    albedo_day.plot(ax=ax1, alpha = 0.5, color = 'gray')
    #aws_year['SR_in'].resample('D').sum().plot(ax = ax1)
    #aws_gaps_year['SR_in'].resample('D').sum().plot(ax=ax1)
    
    ax.legend(loc = 4)
    ax1.set_ylim(0,1)
    ax1.set_ylabel('Albedo', color = 'gray')
    ax1.tick_params(axis='y', labelcolor='gray')
    ax1.spines['right'].set_color('gray')
    ax.set_xlabel('')
    ax.set_ylabel('Ice ablation, m w.e.q.')
    ax1.set_xlim(first,second)
    fig.tight_layout()
    
    observed_ice_melt = observed_ice_ablation.resample('D').mean().diff()
    observed_ice_melt_nans = observed_ice_melt.where(observed_ice_melt>=0.02, np.nan )
    observed_ice_melt_clean = observed_ice_melt_nans.dropna()
    modelled_ice_melt = aws_year['melt, z0w='+zw2].resample('D').sum()
    modelled_ice_melt_clean = modelled_ice_melt.where(observed_ice_melt>=0.02, np.nan )
    modelled_ice_melt_clean = modelled_ice_melt_clean.dropna()
    
    
    
    r2 = 1-(((observed_ice_melt_clean-modelled_ice_melt_clean)**2).sum())/(((observed_ice_melt_clean-observed_ice_melt_clean.mean())**2).sum())
    r2 = np.round(r2*100)/100
    #r2 = r2_score(observed_ice_melt_clean,modelled_ice_melt_clean) 
    rmse = np.sqrt(mean_squared_error(observed_ice_melt_clean, modelled_ice_melt_clean))
    rmse = np.round(rmse*100)/100
    
    fig2, ax2 = plt.subplots(1,1,figsize = (3.5,3))
    ax2.scatter(observed_ice_melt_nans,modelled_ice_melt)
    ax2.plot([0, 0.09], [0, 0.09], color = 'gray', alpha = 0.5)
    ax2.set_ylim([0.0,0.09])
    ax2.set_xlim([0.0,0.09])
    #ax2.text(0.001,0.08,'r² = '+str(r2))
    #ax2.text(0.001,0.075,'rmse = '+str(rmse))
    ax2.set_ylabel('Modelled, m w.eq./day')
    ax2.set_xlabel('Observed, m w.eq./day')
    fig2.tight_layout()
    
    fig.savefig('figures/point_melt_model/zac_u_model_vs_obs_'+year+'.png', dpi = 300)

    #first = albedo_day[year][albedo_day<=ice_albedo].index[0]
    #second = str(year)+'-September'
    #aws_year = aws[first:second] 
    #aws_gaps_year = aws_gaps[first:second]
    ##print(aws_year.keys())
    #observed_ice_ablation = aws_year['ice_ablation']*0.917
    #observed_ice_ablation = observed_ice_ablation-observed_ice_ablation[0]
    #
    #fig,ax = plt.subplots(1,1)
    #ax1 = ax.twinx()
    #observed_ice_ablation.rolling(24,center=True).mean().plot(ax = ax1, label = 'Observed ice melt', linewidth = 3)
    #aws_year['SW_net'].resample('D').sum().plot(ax = ax, label = 'SR_net', linestyle = '--')
    #aws_year['LW_net'].resample('D').sum().plot(ax = ax, label = 'LR_net', linestyle = '--')
    #aws_year['sensible'].resample('D').sum().plot(ax = ax, label = 'Sensible', linestyle = '--')
    #aws_year['latent'].resample('D').sum().plot(ax = ax, label = 'Latent', linestyle = '--')
    #
    #
    #ax.legend()
    #ax.set_ylabel('Energy for melt, W/m2')
    #
    #fig.savefig('figures/point_melt_model/energy_balance_'+year+'.png', dpi = 300)

# Calculate turbulent heat fluxes
rho_0 = 1.29 # kg/m3 air density at P_0
cp = 1005 #J/(Kg K), specific heat of dry air at constant pressure
kappa = 0.41 # von karmans constant
P_0 = 101325 # Pa standard atm pressure 
z = 2.7 #m instrument height
Ls = 2.849*10**6 # J/Kg latent heat of sublimation
Lv = 2.514*10**6 # J/kg latent heat of vaporization
e0 = 611 # Pa Vapour pressure at the melting surface

transfer_coeff_heat = 1.5/1000 # Only for the bulk method
transfer_coeff_vapor = 1.5/1000 # Only for the bulk method

T2 = aws['t_u'] # temperature at 2m
T0 = aws['t_surf']
# calculating vapour pressure at 2m e2
RH = aws['rh_u_corr']
e2 = 611*np.exp(17.27*T2/(243.04+T2))*RH/100 

P = aws['p_u']*100 # atm pressure from hPa to Pa
u2 = aws['wspd']# wind speed at 2m

z0w = roughness_w # roughness parameters for logarithmic profiles of wind 
z0t = z0w/roughness_div # roughness parameters for logarithmic profiles of temp
z0e = z0t # roughness parameters for logarithmic profiles of vapour pressure

L = T2*0+Lv
L[(e2-e0 < 0)]= Ls
L[(e2-e0 >= 0) & (T0<0) ]= Ls


sensible = P/P_0*cp*(kappa**2)*rho_0*u2*T2/(np.log(z/z0w)*np.log(z/z0t))
#latent = L*(0.623/P_0)*(kappa**2/(np.log(z/z0w)*np.log(z/z0e)))*u2*(e2-e0)

latent = 0.623*L*kappa**2*(rho_0/P_0)*(u2*(e2-e0)/(np.log(z/z0w)*np.log(z/z0e)))


#L = T2*0+Lv
#L[(e2-e0 < 0)]= Ls
#L[(e2-e0 >= 0) & (T0<0) ]= Ls
#Ch=transfer_coeff_heat
#Ce=transfer_coeff_vapor
#
#sensible=0.0129*Ch*P*u2*(T2) #Eh
#latent=0.622*rho_0*L*Ce*u2*(e2-e0)/P_0 #Ee
aws = aws.assign(sensible=sensible)
aws = aws.assign(latent=latent)


# Calculate the energy available for melt 
Emelt = aws['SW_net']+aws['LW_net']+aws['sensible']+aws['latent']
Emelt[Emelt<0] = 0

# Calculate the amount of ice that can be melted with the energy available for melt
Lf = 334000 # J/kg
melt = Emelt*3600*(1)/(Lf*1000) # Qmelt/Lf * number of seconds per timestep / density of water = m WE/ timestep
aws['melt, z0w='+str(roughness_w)] = melt
aws['Emelt, z0w='+str(roughness_w)] = Emelt

first = albedo_day[year][albedo_day<=ice_albedo].index[0]
second = str(year)+'-September'
aws_year = aws[first:second] 
aws_gaps_year = aws_gaps[first:second]
#print(aws_year.keys())
observed_ice_ablation = aws_year['ice_ablation']*0.917
observed_ice_ablation = observed_ice_ablation-observed_ice_ablation[0]

fig,ax = plt.subplots(1,1)
ax1 = ax.twinx()
observed_ice_ablation.rolling(24,center=True).mean().plot(ax = ax1, label = 'Observed ice melt', linewidth = 3)
aws_year['SW_net'].resample('D').sum().plot(ax = ax, label = 'SR_net', linestyle = '--')
aws_year['LW_net'].resample('D').sum().plot(ax = ax, label = 'LR_net', linestyle = '--')
aws_year['sensible'].resample('D').sum().plot(ax = ax, label = 'Sensible', linestyle = '--')
aws_year['latent'].resample('D').sum().plot(ax = ax, label = 'Latent', linestyle = '--')


ax.legend()
ax.set_ylabel('Energy for melt, W/m2')

first = icemeltstart #albedo_day[year][albedo_day<=ice_albedo].index[0]
second = icemeltend #str(year)+'-September'
aws_year = aws[first:second] 

observed_ice_ablation = aws_year['ice_ablation']*0.917

observed_ice_ablation = observed_ice_ablation-observed_ice_ablation[0]

fig,ax = plt.subplots(1,1,figsize = (3.5,5))
ax1 = ax.twinx()
observed_ice_ablation.rolling(24,center=True).mean().plot(ax = ax, label = 'Observed', linewidth = 3)
aws_year['melt, z0w='+zw1].cumsum().plot(ax = ax, label = 'z0w = '+zw1, linestyle = '--')
aws_year['melt, z0w='+zw2].cumsum().plot(ax = ax, label = 'z0w = '+zw2, linestyle = '--')
aws_year['melt, z0w='+zw3].cumsum().plot(ax = ax, label = 'z0w = '+zw3, linestyle = '--')
albedo_day.plot(ax=ax1, alpha = 0.5, color = 'gray')
#aws_year['SR_in'].resample('D').sum().plot(ax = ax1)
#aws_gaps_year['SR_in'].resample('D').sum().plot(ax=ax1)

ax.legend(loc = 4)
ax1.set_ylim(0,1)
ax1.set_ylabel('Albedo', color = 'gray')
ax1.tick_params(axis='y', labelcolor='gray')
ax1.spines['right'].set_color('gray')
ax.set_xlabel('')
ax.set_ylabel('Ice ablation, m w.e.q.')
ax1.set_xlim(first,second)
fig.tight_layout()

observed_ice_melt = observed_ice_ablation.resample('D').mean().diff()
observed_ice_melt_nans = observed_ice_melt.where(observed_ice_melt>=0.02, np.nan )
observed_ice_melt_clean = observed_ice_melt_nans.dropna()
modelled_ice_melt = aws_year['melt, z0w='+zw2].resample('D').sum()
modelled_ice_melt_clean = modelled_ice_melt.where(observed_ice_melt>=0.02, np.nan )
modelled_ice_melt_clean = modelled_ice_melt_clean.dropna()



r2 = 1-(((observed_ice_melt_clean-modelled_ice_melt_clean)**2).sum())/(((observed_ice_melt_clean-observed_ice_melt_clean.mean())**2).sum())
r2 = np.round(r2*100)/100
#r2 = r2_score(observed_ice_melt_clean,modelled_ice_melt_clean) 
rmse = np.sqrt(mean_squared_error(observed_ice_melt_clean, modelled_ice_melt_clean))
rmse = np.round(rmse*100)/100

fig2, ax2 = plt.subplots(1,1,figsize = (3.5,3))
ax2.scatter(observed_ice_melt_nans,modelled_ice_melt)
ax2.plot([0, 0.09], [0, 0.09], color = 'gray', alpha = 0.5)
ax2.set_ylim([0.0,0.09])
ax2.set_xlim([0.0,0.09])
#ax2.text(0.001,0.08,'r² = '+str(r2))
#ax2.text(0.001,0.075,'rmse = '+str(rmse))
ax2.set_ylabel('Modelled, m w.eq./day')
ax2.set_xlabel('Observed, m w.eq./day')
fig2.tight_layout()

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
#load_reindexed_dat
ice_albedo = 0.4
aws_gaps = pd.read_csv('/home/shl/OneDrive/projects/aws_processing_v1.0/data_v1.0/gem_database/2022/zac_l_2008_2022_QC_final.csv', index_col = 0, parse_dates = True)

aws_gaps['SW_net'] = aws_gaps.SR_in - aws_gaps.SR_out
aws_gaps['LW_net'] = aws_gaps.LR_in - aws_gaps.LR_out
aws_gaps['ice_ablation'] = aws_gaps.Z_pt.values
aws_interpolated = aws_gaps.interpolate()

albedo_day = aws_gaps['albedo'].resample('D').median()

# The coefficients to calibrate the melt model
roughness_div = 100


roughness_w = 0.001
aws = aws_interpolated
# Calculate turbulent heat fluxes
rho_0 = 1.29 # kg/m3 air density at P_0
cp = 1005 #J/(Kg K), specific heat of dry air at constant pressure
kappa = 0.41 # von karmans constant
P_0 = 101325 # Pa standard atm pressure 
z = 2.7 #m instrument height
Ls = 2.849*10**6 # J/Kg latent heat of sublimation
Lv = 2.514*10**6 # J/kg latent heat of vaporization
e0 = 611 # Pa Vapour pressure at the melting surface

transfer_coeff_heat = 1.5/1000 # Only for the bulk method
transfer_coeff_vapor = 1.5/1000 # Only for the bulk method

T2 = aws['t_u'] # temperature at 2m
T0 = aws['t_surf']
# calculating vapour pressure at 2m e2
RH = aws['rh_u_corr']
e2 = 611*np.exp(17.27*T2/(243.04+T2))*RH/100 

P = aws['p_u']*100 # atm pressure from hPa to Pa
u2 = aws['wspd']# wind speed at 2m

z0w = roughness_w # roughness parameters for logarithmic profiles of wind 
z0t = z0w/roughness_div # roughness parameters for logarithmic profiles of temp
z0e = z0t # roughness parameters for logarithmic profiles of vapour pressure

L = T2*0+Lv
L[(e2-e0 < 0)]= Ls
L[(e2-e0 >= 0) & (T0<0) ]= Ls


sensible = P/P_0*cp*(kappa**2)*rho_0*u2*T2/(np.log(z/z0w)*np.log(z/z0t))
#latent = L*(0.623/P_0)*(kappa**2/(np.log(z/z0w)*np.log(z/z0e)))*u2*(e2-e0)

latent = 0.623*L*kappa**2*(rho_0/P_0)*(u2*(e2-e0)/(np.log(z/z0w)*np.log(z/z0e)))


#L = T2*0+Lv
#L[(e2-e0 < 0)]= Ls
#L[(e2-e0 >= 0) & (T0<0) ]= Ls
#Ch=transfer_coeff_heat
#Ce=transfer_coeff_vapor
#
#sensible=0.0129*Ch*P*u2*(T2) #Eh
#latent=0.622*rho_0*L*Ce*u2*(e2-e0)/P_0 #Ee
aws = aws.assign(sensible=sensible)
aws = aws.assign(latent=latent)


# Calculate the energy available for melt 
Emelt = aws['SW_net']+aws['LW_net']+aws['sensible']+aws['latent']
Emelt[Emelt<0] = 0

# Calculate the amount of ice that can be melted with the energy available for melt
Lf = 334000 # J/kg
melt = Emelt*3600*(1)/(Lf*1000) # Qmelt/Lf * number of seconds per timestep / density of water = m WE/ timestep
aws['melt, z0w='+str(roughness_w)] = melt
aws['Emelt, z0w='+str(roughness_w)] = Emelt
melt_interpolated = aws['melt, z0w=0.001']

roughness_w = 0.001
aws = aws_gaps
# Calculate turbulent heat fluxes
rho_0 = 1.29 # kg/m3 air density at P_0
cp = 1005 #J/(Kg K), specific heat of dry air at constant pressure
kappa = 0.41 # von karmans constant
P_0 = 101325 # Pa standard atm pressure 
z = 2.7 #m instrument height
Ls = 2.849*10**6 # J/Kg latent heat of sublimation
Lv = 2.514*10**6 # J/kg latent heat of vaporization
e0 = 611 # Pa Vapour pressure at the melting surface

transfer_coeff_heat = 1.5/1000 # Only for the bulk method
transfer_coeff_vapor = 1.5/1000 # Only for the bulk method

T2 = aws['t_u'] # temperature at 2m
T0 = aws['t_surf']
# calculating vapour pressure at 2m e2
RH = aws['rh_u_corr']
e2 = 611*np.exp(17.27*T2/(243.04+T2))*RH/100 

P = aws['p_u']*100 # atm pressure from hPa to Pa
u2 = aws['wspd']# wind speed at 2m

z0w = roughness_w # roughness parameters for logarithmic profiles of wind 
z0t = z0w/roughness_div # roughness parameters for logarithmic profiles of temp
z0e = z0t # roughness parameters for logarithmic profiles of vapour pressure

L = T2*0+Lv
L[(e2-e0 < 0)]= Ls
L[(e2-e0 >= 0) & (T0<0) ]= Ls


sensible = P/P_0*cp*(kappa**2)*rho_0*u2*T2/(np.log(z/z0w)*np.log(z/z0t))
#latent = L*(0.623/P_0)*(kappa**2/(np.log(z/z0w)*np.log(z/z0e)))*u2*(e2-e0)

latent = 0.623*L*kappa**2*(rho_0/P_0)*(u2*(e2-e0)/(np.log(z/z0w)*np.log(z/z0e)))


#L = T2*0+Lv
#L[(e2-e0 < 0)]= Ls
#L[(e2-e0 >= 0) & (T0<0) ]= Ls
#Ch=transfer_coeff_heat
#Ce=transfer_coeff_vapor
#
#sensible=0.0129*Ch*P*u2*(T2) #Eh
#latent=0.622*rho_0*L*Ce*u2*(e2-e0)/P_0 #Ee
aws = aws.assign(sensible=sensible)
aws = aws.assign(latent=latent)


# Calculate the energy available for melt 
Emelt = aws['SW_net']+aws['LW_net']+aws['sensible']+aws['latent']
Emelt[Emelt<0] = 0

# Calculate the amount of ice that can be melted with the energy available for melt
Lf = 334000 # J/kg
melt = Emelt*3600*(1)/(Lf*1000) # Qmelt/Lf * number of seconds per timestep / density of water = m WE/ timestep
aws['melt, z0w='+str(roughness_w)] = melt
aws['Emelt, z0w='+str(roughness_w)] = Emelt
melt_gaps = aws['melt, z0w=0.001']


fig,ax = plt.subplots(1,1)
#melt.plot(ax=ax)
melt_gaps.cumsum().plot(ax=ax, label = 'Raw data') 
melt_interpolated.cumsum().plot(ax=ax, label = 'Interpolated data')
ax.legend()
ax.set_ylabel('Cummulated modelled melt')
