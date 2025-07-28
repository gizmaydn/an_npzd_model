
library(deSolve)


npzd <- function(time, variables, parameters) {
  with(as.list(c(variables, parameters)), {
    
    light <- 1 + 0.5 * abs(sin( pi * time / 365))  # no unit, it is a relative value 
    light_lim = light/ (ksLight + light)
    
    
    temp  <- 2 + 18 * abs(sin( pi * (time) / 365)) #celcius
    
    
    temp_lim <- exp(0.063 * temp) #eppley function
    
    p_resp = Ipn * pp  # phyto respiration
    z_resp = Izn * zoo  # zoo respiration 
    n_upt = rmax * light_lim * nut / (kn + nut) * temp_lim * pp # nutrient uptake
    mineralization = Idn * detr  # remineralization
    graz = gmax * (1 - exp(-Iv*pp)) *  zoo  * temp_lim # grazing
    z_mort = Izd* zoo # zoo mortality
    p_mort = Ipd *  pp # phyto mortality
    
    dpp = n_upt - p_mort - graz - p_resp
    dzoo = graz - z_mort - z_resp
    ddetr = p_mort + z_mort - mineralization
    dnut =  p_resp + z_resp + mineralization - n_upt
    
    
    return(list(c(dpp, dzoo, ddetr, dnut)))
    
  })
}

#parameters

parameters <- c(rmax = 1.2, # nutrient uptake
                kn = 0.1, # half saturation constant
                Ipn = 0.01, # phyto respiration
                Ipd = 0.02, # phyto mortality
                Idn = 0.001, # remineralization
                gmax = 0.2, # maximum grazing
                Iv = 1.8, # Ivlev constant for grazing of zoo on pp
                Izn = 0.01, # zoo respiration
                Izd = 0.02 ,# zoo mortality
                ksLight = 20
)

#initial conditionss

variables <- c(pp =0.01, zoo=0.01, detr=0.99, nut= 0.99 )

#similation period

simtime <- seq(0, 730, by = 1)


#simulate
output <- as.data.frame(ode(func = npzd, y = variables, parms = parameters, times = simtime))

#mass conservation check
output$total_mass <- output$pp + output$zoo + output$detr + output$nut

# plot the simulation

matplot(output[, -1], type = "l", xlab = "Time (Days)", ylab = "Concentration",
        col = c("green", "blue", "brown", "red","orange"), lty = c(1, 1, 1, 1))

legend("topright", c("Phytoplankton", "Zooplankton", "Detritus", "Nutrient","mass"),
       lty = c(1, 1, 1, 1), col = c("green", "blue", "brown", "red","orange"), box.lwd = 0)


#light <- 1 + 0.5 * abs(sin( pi * simtime / 365))  # Seasonal light variation (1-year period)
#light_lim = light/ (ksLight + light)
#ksLight = 2
#plot(light)

#temp  <- 2 + 18 * abs(sin( pi * (simtime) / 360)) 
#temp_lim <- exp(0.063 * temp)
#plot(temp)
