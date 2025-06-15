module iepox_aqueous

implicit none  
  
contains 

subroutine iepox_aqueous_tendencies(nzm, T, p, Rp, Wl, H, NO3, SO4, HSO4, &
  IEPOXg, IEPOXaq, PRODUCTaq, TETROLg, TETROLaq, &
  IEPOX_NO3g, IEPOX_NO3aq, IEPOX_SO4g, IEPOX_SO4aq, &
  IEPOXg_tend, IEPOXaq_tend, PRODUCTaq_tend, TETROLg_tend, TETROLaq_tend, &
  IEPOX_NO3g_tend, IEPOX_NO3aq_tend, IEPOX_SO4g_tend, IEPOX_SO4aq_tend)

  implicit none
  
  integer, intent(in) :: nzm   ! number of vertical grid levels  
  real, dimension(1:nzm), intent(in) :: T   !  K 
  real, dimension(1:nzm), intent(in) :: p   ! atm
  real, dimension(1:nzm), intent(in) :: Rp  ! droplet radius (m)
  real, dimension(1:nzm), intent(in) :: Wl  ! liquid water content vol H2O/vol air
  real, dimension(1:nzm), intent(in) :: H   !  M of H+
  real, dimension(1:nzm), intent(in) :: NO3   !  M 
  real, dimension(1:nzm), intent(in) :: SO4   !  M 
  real, dimension(1:nzm), intent(in) :: HSO4   !  M 

  real, dimension(1:nzm), intent(in) :: IEPOXg, IEPOXaq
  real, dimension(1:nzm), intent(in) :: PRODUCTaq
  real, dimension(1:nzm), intent(in) :: TETROLg, TETROLaq
  real, dimension(1:nzm), intent(in) :: IEPOX_NO3g, IEPOX_NO3aq
  real, dimension(1:nzm), intent(in) :: IEPOX_SO4g, IEPOX_SO4aq

  real, dimension(1:nzm), intent(out) :: IEPOXg_tend, IEPOXaq_tend
  real, dimension(1:nzm), intent(out) :: PRODUCTaq_tend
  real, dimension(1:nzm), intent(out) :: TETROLg_tend, TETROLaq_tend
  real, dimension(1:nzm), intent(out) :: IEPOX_NO3g_tend, IEPOX_NO3aq_tend
  real, dimension(1:nzm), intent(out) :: IEPOX_SO4g_tend, IEPOX_SO4aq_tend

  integer :: k

  ! constants
  real :: R = 0.082057 !L*atm/k*mol
  real :: Dg = 1e-5  !  m^2/s - diffusitivity of IEPOX in air (0.1 cm^2/s)
  real :: ah = 0.25 !value is from eddingsaas - H/55.5 #no units - mol/L of H3O+/55.5 mol/L water - mol fraction of H+ in water - activity 
  real :: a = 1 !no unit - molec enering liquid/molec collisions with surface 

  real :: Daq = 1e-9 !m^2/s - diffusion of A in water 
  real :: Ha = 1.7e8  ! M/atm - henry's law constant - value from gaston 

  real :: kh = 1.2e-3 !M-1s-1 - rate constant of H
  real :: kno3 = 2.e-4 !M-1s-1 rate constant of nucleohpile NO3-  
  real :: kso4 = 1.e-4 !M-1s-1 - PLACEHOLDER - rate constant of nucleophile SO4-2 
  real :: khso4 = 7.3e-4 !M-1s-1 rate constant of general acids  

  real :: Ma = 0.1191390 !kg/mol - mass of IEPOX
  real :: Mt = 0.13615     !kg/mol - mass of tetrol
  real :: Ms = 0.216123   ! kg/mol mass of sulfate product
  real :: Mn = 0.18215    ! kg/mol estimate of nitrate product

  real :: Htet = 1e8 !M/atm tetrol estimate 
  real :: Hs = 1e16 !M/atm sulfate product estimate 
  real :: Hn = 1e7 !M/atm nitrate product estimate 
  
  real, dimension(1:nzm) :: kaq  ! s-1 full rate constant
  real, dimension(1:nzm) :: faqH, faqNO3, faqSO4  ! rate fractions
  real, dimension(1:nzm) :: qq   ! m/(1/s/m^2/s)**0.5 = no units
  real, dimension(1:nzm) :: Q   ! no units
  real, dimension(1:nzm) :: w  ! m/s   mean speed
  real, dimension(1:nzm) :: Kmt  ! s-1 mass transfer coefficient
  real :: blah

  kaq = (kh*H) + (kno3*NO3*ah) + (kso4*SO4*ah) + (khso4*HSO4) ! full rate constant
     ! add epps to avoid divide by zero????
  faqH = (kh*H)/kaq
  faqNO3 = (kno3*NO3*ah)/kaq
  faqSO4 = ((kso4*SO4*ah) + (khso4*HSO4))/kaq
  qq = Rp*((kaq/Daq)**0.5)
  Q = 3.0*((1.0/tanh(qq))/qq)-(1.0/(qq**2))
  w = ((8*R*T)/(3.1415*Ma))**0.5 ! m/s book equ for mean speed
  Kmt = ((Rp**2/(3*Dg)) + (4*Rp/(3*w*a))) **(-1) ! s-1 mass transfer coefficient

  !write (*,*) ah, kaq, faqH, faqNO3, faqSO4, qq, Q, w
  !write (*,*) Kmt(1), IEPOXg(1), IEPOXaq(1)
 
  !   ODE terms

  IEPOXg_tend = -(Kmt*Wl*IEPOXg) + ((1/Ha)*Kmt*IEPOXaq*Wl)
  IEPOXaq_tend = ((Kmt*IEPOXg)/(R*T))-((Kmt*IEPOXaq)/(Ha*R*T))-(Q*kaq*IEPOXaq)
  PRODUCTaq_tend = (Q*kaq*IEPOXaq)
  TETROLg_tend =  ((Kmt*TETROLaq*Wl)/Htet)-(Kmt*TETROLg*Wl)
  TETROLaq_tend = (faqH*Q*kaq*IEPOXaq)-(Kmt*TETROLaq/(R*T*Htet))+(Kmt*TETROLg/(R*T))
  IEPOX_NO3g_tend = ((Kmt*IEPOX_NO3aq*Wl)/Hn)-(Kmt*IEPOX_NO3g*Wl)
  IEPOX_NO3aq_tend = (faqNO3*Q*kaq*IEPOXaq)-(Kmt*IEPOX_NO3aq/(R*T*Hn))+(Kmt*IEPOX_NO3g/(R*T))
  IEPOX_SO4g_tend = ((Kmt*IEPOX_SO4aq*Wl)/Hs)-(Kmt*IEPOX_SO4g*Wl)
  IEPOX_SO4aq_tend = (faqSO4*Q*kaq*IEPOXaq)-(Kmt*IEPOX_SO4aq/(R*T*Hs))+(Kmt*IEPOX_SO4g/(R*T))

end subroutine iepox_aqueous_tendencies

end module iepox_aqueous
