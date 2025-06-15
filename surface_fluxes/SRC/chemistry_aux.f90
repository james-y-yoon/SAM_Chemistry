module chemistry_aux
! auxiliary functions for chemistry 
  
   use chemistry_params, only: OH_night, OH_day_peak

   implicit none

   real :: pi=3.1415927

   real :: tsunrise = 10.     ! hour UTC hardcoded for March 15
   real :: tzenith  = 15.     ! hour UTC hardcoded for March 15 

   real :: NO2_photolysis_night = 0  ! /s
   real :: NO2_photolysis_peak =  0.01 !  /s
   
 CONTAINS

   ! still not used
   real function get_OH(time)  ! time in days 
     real :: t_solar_peak = 0.167 ! days
     real :: time
     get_OH = MAX(OH_night, OH_day_peak*cos((time - t_solar_peak)*2.*pi))
     ! get OH concentration in molecules/cm3 as a function of hour
     !get_OH = truncated_cos(t, tsunrise, tzenith, OH_night, OH_day_peak)

   end function get_OH
    
   !still not used
   real function get_NO2_photolysis_rate(t)
     real :: t  ! time in hours
     !  get NO2 photolysis rate constant as a function of hour  
       get_NO2_photolysis_rate = truncated_cos(t, tsunrise, tzenith, &
            NO2_photolysis_night, NO2_photolysis_peak)

   end function get_NO2_photolysis_rate    
   
   real function truncated_cos(t, tsunrise, tzenith, Vnight, Vpeak)

     ! compute truncated cosine curve based on time of day
     !  tsunrise is time of sunrise in hours
     !  tzenith is time of local noon in hours
     !  V night is the desired value whenever sun is down     
     !  V peak  is the desired value at local noon     
     
     real :: t, tsunrise, tzenith, Vnight, Vpeak, A
         
     tsunrise = tsunrise * pi/12.
     tzenith = tzenith * pi/12

     A = (Vnight - Vpeak)/(cos(tsunrise-tzenith) - 1)
     truncated_cos = MAX(A * cos(t*pi/12. - tzenith) + Vpeak - A, Vnight)     

     
   end function truncated_cos

   ! Following function not in use
   real function deposition_loss_rate(conc, rate_constant)
     ! Compute deposition loss rate based on concentration
     real::  conc, rate_constant
     deposition_loss_rate = conc * rate_constant

   end function deposition_loss_rate  


   

   
end module chemistry_aux
