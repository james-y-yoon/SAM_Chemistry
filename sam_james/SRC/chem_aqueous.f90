module chem_aqueous

  implicit none  

  ! Aqueous chemistry variables, one aqueous and one gas for each species
  integer, parameter :: naqchem_fields = 8
  integer :: iTETROL = 1
  integer :: iIEPOX_NO3 = 2
  integer :: iIEPOX_SO4 = 3
  integer :: iIEPOX = 4  
  integer :: iIP1NIT = 5
  integer :: iTRIOL = 6
  integer :: iHNO3 = 7
  integer :: iIPDINIT = 8
    
  real,  parameter, dimension(naqchem_fields) :: molwt = (/ &
       136.15, 182.15, 216.123, 119.14, &
        79.8, 118.13, 63.01, 133.3/) ! g/mol corresponding to above species  
  
  logical, dimension(naqchem_fields), public :: flag_aqchemvar_out3D  ! which aqueous chem array to  output

  logical, dimension(naqchem_fields), public :: flag_aqchemgasvar_out3D  ! which aqueous chem product array to output

  character(LEN=7), parameter, dimension(naqchem_fields) :: aq_species_names = (/ &
       'TETROLa', 'IPXNO3a', 'IPXSO4a', 'IEPOXa ', &
        'IP1NITa', 'TRIOLa ', 'HNO3a  ', 'IPDiNta'/)
  character(LEN=7), parameter, dimension(naqchem_fields) :: aq_gasprod_species_names = (/ &
       'TETROLg', 'IPXNO3g', 'IPXSO4g', 'IEPOXg ', &
        'IP1NITg', 'TRIOLg ', 'HNO3g  ', 'IPDiNtg'/)

  contains
  
subroutine iepox_aqueous_tendencies(nzm, T, p, Rp, Wl, H, NO3, SO4, HSO4, &
  gas_input, aq_input, gas_tend, aq_tend, print_debug_output)

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

  real, dimension(1:nzm, 1:naqchem_fields), intent(in) :: gas_input, aq_input
  real, dimension(1:nzm, 1:naqchem_fields), intent(out) :: gas_tend, aq_tend  

  logical, intent(in) :: print_debug_output
  logical :: add_text 
  
  real, dimension(1:nzm) :: IEPOXg, IEPOXaq
  real, dimension(1:nzm) :: TETROLg, TETROLaq
  real, dimension(1:nzm) :: IEPOX_NO3g, IEPOX_NO3aq
  real, dimension(1:nzm) :: IEPOX_SO4g, IEPOX_SO4aq

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


  kaq = (kh*H) + (kno3*NO3*ah) + (kso4*SO4*ah) + (khso4*HSO4) ! full rate constant
     ! add epps to avoid divide by zero????
  faqH = (kh*H)/kaq
  faqNO3 = (kno3*NO3*ah)/kaq
  faqSO4 = ((kso4*SO4*ah) + (khso4*HSO4))/kaq
  qq = Rp*((kaq/Daq)**0.5)
  Q = 3.0*(((1.0/tanh(qq))/qq)-(1.0/(qq**2)))
  w = ((8*R*T)/(3.1415*Ma))**0.5 ! m/s book equ for mean speed
  Kmt = ((Rp**2/(3*Dg)) + (4*Rp/(3*w*a))) **(-1) ! s-1 mass transfer coefficient

  !write (*,*) ah, kaq, faqH, faqNO3, faqSO4, qq, Q, w
  !write (*,*) Kmt(1), IEPOXg(1), IEPOXaq(1)


  IEPOXg = gas_input(:, iIEPOX)
  IEPOXaq = aq_input(:, iIEPOX)
  TETROLg = gas_input(:, iTETROL)
  TETROLaq = aq_input(:, iTETROL)
  IEPOX_NO3g = gas_input(:, iIEPOX_NO3)
  IEPOX_NO3aq = aq_input(:, iIEPOX_NO3)
  IEPOX_SO4g = gas_input(:, iIEPOX_SO4)
  IEPOX_SO4aq = aq_input(:, iIEPOX_SO4)
  
  
  !   ODE terms

  gas_tend(:, iIEPOX) = -(Kmt*Wl*IEPOXg) + ((1/Ha)*Kmt*IEPOXaq*Wl)
  aq_tend(:, iIEPOX) = ((Kmt*IEPOXg)/(R*T))-((Kmt*IEPOXaq)/(Ha*R*T))-(Q*kaq*IEPOXaq)
  gas_tend(:, iTETROL) = ((Kmt*TETROLaq*Wl)/Htet)-(Kmt*TETROLg*Wl)
  aq_tend(:, iTETROL) = (faqH*Q*kaq*IEPOXaq)-(Kmt*TETROLaq/(R*T*Htet))+(Kmt*TETROLg/(R*T))
  gas_tend(:, iIEPOX_NO3) = ((Kmt*IEPOX_NO3aq*Wl)/Hn)-(Kmt*IEPOX_NO3g*Wl)
  aq_tend(:, iIEPOX_NO3) = (faqNO3*Q*kaq*IEPOXaq)-(Kmt*IEPOX_NO3aq/(R*T*Hn))+(Kmt*IEPOX_NO3g/(R*T))
  gas_tend(:, iIEPOX_SO4) =  ((Kmt*IEPOX_SO4aq*Wl)/Hs)-(Kmt*IEPOX_SO4g*Wl)
  aq_tend(:, iIEPOX_SO4) = (faqSO4*Q*kaq*IEPOXaq)-(Kmt*IEPOX_SO4aq/(R*T*Hs))+(Kmt*IEPOX_SO4g/(R*T))

  if (print_debug_output) then
     add_text=.true.
     do k = 1,nzm
        
        if (Wl(k).gt.1.e-10) then
           if (add_text) then
              write (*,*) 'k, Wl, IEPOXg, IEPOXg+, IEPOXa, IEPOXa+, TETROLg, TETROLg+, TETROLa, TETROLa+'
              add_text = .false.
           end if   
              
              write(*,*) k, Wl(k), gas_input(k,iIEPOX), gas_tend(k,iIEPOX), aq_input(k,iIEPOX), aq_tend(k,iIEPOX) !
           !  gas_input(k,iTETROL), aq_input(k,iTETROL), gas_tend(k,iTETROL), aq_tend(k,iTETROL)
        end if
     end do   
  end if              

  
end subroutine iepox_aqueous_tendencies

subroutine isop1nit_aqueous_tendencies(nzm, T, p, Rp, Wl, gas_input, aq_input, gas_tend, aq_tend, print_debug_output)

  implicit none
  
  integer, intent(in) :: nzm   ! number of vertical grid levels  
  real, dimension(1:nzm), intent(in) :: T   !  K 
  real, dimension(1:nzm), intent(in) :: p   ! atm
  real, dimension(1:nzm), intent(in) :: Rp  ! droplet radius (m)
  real, dimension(1:nzm), intent(in) :: Wl  ! liquid water content vol H2O/vol air

  real, dimension(1:nzm, 1:naqchem_fields), intent(in) :: gas_input, aq_input
  real, dimension(1:nzm, 1:naqchem_fields), intent(out) :: gas_tend, aq_tend  

  logical, intent(in) :: print_debug_output
  logical :: add_text 
  
  real, dimension(1:nzm) :: IP1NITg, IP1NITaq
  real, dimension(1:nzm) :: TRIOLg, TRIOLaq
  real, dimension(1:nzm) :: HNO3g, HNO3aq
  real, dimension(1:nzm) :: IPDINITg, IPDINITaq
  integer :: k

  ! constants
  real :: R = 0.082057 !L*atm/k*mol
  real :: Dg = 1e-5  !  m^2/s - diffusitivity of IEPOX in air (0.1 cm^2/s)
  real :: ah = 0.25 !value is from eddingsaas - H/55.5 #no units - mol/L of H3O+/55.5 mol/L water - mol fraction of H+ in water - activity 
  real :: a = 1 !no unit - molec enering liquid/molec collisions with surface 

  real :: Daq = 1e-9 !m^2/s - diffusion of A in water 
  real :: Ha = 1.7e8  ! M/atm - henry's law constant - value from gaston 

  real :: H_ip1 = 1e8 !M/atm tetrol estimate 
  real :: H_triol = 1e8 !M/atm sulfate product estimate 
  real :: H_hno3 = 1e8 !M/atm nitrate product estimate
  real :: H_dinit = 1e8

  real :: M_ip1  ! molecular weight in kg/mol
  real :: M_dinit  ! molecular weight in kg/mol
  
  real, dimension(1:nzm) :: kaq_isop1  ! s-1 full rate constant
    real, dimension(1:nzm) :: kaq_dinit  ! s-1 full rate constant
  real, dimension(1:nzm) :: qq1, qq2   ! m/(1/s/m^2/s)**0.5 = no units
    
  real, dimension(1:nzm) :: Q_isop1, Q_dinit   ! no units
  real, dimension(1:nzm) :: w1, w2  ! m/s   mean speed
  real, dimension(1:nzm) :: Kmt_isop1, Kmt_dinit  ! s-1 mass transfer coefficient

  M_ip1 = molwt(iIP1NIT)/1000. ! converting g/mol to kg/mol
  M_dinit = molwt(iIPDINIT)/1000.

  kaq_isop1 = 1.
  kaq_dinit = 1.
  qq1 = Rp*((kaq_isop1/Daq)**0.5)
  qq2 = Rp*((kaq_dinit/Daq)**0.5)
  Q_isop1 = 3.0*(((1.0/tanh(qq1))/qq1)-(1.0/(qq1**2)))
  Q_dinit = 3.0*(((1.0/tanh(qq2))/qq2)-(1.0/(qq2**2)))
  w1 = ((8*R*T)/(3.1415*M_ip1))**0.5 ! m/s book equ for mean speed
  w2 = ((8*R*T)/(3.1415*M_dinit))**0.5 ! m/s book equ for mean speed
  Kmt_isop1 = ((Rp**2/(3*Dg)) + (4*Rp/(3*w1*a))) **(-1) ! s-1 mass transfer coefficient
  Kmt_dinit = ((Rp**2/(3*Dg)) + (4*Rp/(3*w2*a))) **(-1) ! s-1 mass transfer coefficient


  IP1NITg = gas_input(:, iIP1NIT)
  IP1NITaq = aq_input(:, iIP1NIT)
  TRIOLg = gas_input(:, iTRIOL)
  TRIOLaq = aq_input(:, iTRIOL )
  HNO3g = gas_input(:, iHNO3)
  HNO3aq = aq_input(:, iHNO3)
  
  !   ODE terms

  gas_tend(:, iIP1NIT) = -(Kmt_isop1*Wl*IP1NITg) + ((1/H_ip1)*Kmt_isop1*IP1NITaq*Wl)
  aq_tend(:, iIP1NIT) = ((Kmt_isop1*IP1NITg)/(R*T))-((Kmt_isop1*IP1NITaq)/(H_ip1*R*T))-(Q_isop1*kaq_isop1*IP1NITaq)
  gas_tend(:, iTRIOL) = ((Kmt_isop1*TRIOLaq*Wl)/H_triol)-(Kmt_isop1*TRIOLg*Wl)
  aq_tend(:, iTRIOL) = (Q_isop1*kaq_isop1*IP1NITaq)-(Kmt_isop1*TRIOLaq/(R*T*H_triol))+(Kmt_isop1*TRIOLg/(R*T))
  gas_tend(:, iHNO3) = ((Kmt_isop1*HNO3aq*Wl)/H_hno3)-(Kmt_isop1*HNO3g*Wl)
  aq_tend(:, iHNO3) = (Q_isop1*kaq_isop1*IP1NITaq)-(Kmt_isop1*HNO3aq/(R*T*H_hno3))+(Kmt_isop1*HNO3g/(R*T))

  gas_tend(:, iIPDINIT) = -(Kmt_dinit*Wl*IPDINITg) + ((1/H_dinit)*Kmt_dinit*IPDINITaq*Wl)
  aq_tend(:, iIPDINIT) = ((Kmt_dinit*IPDINITg)/(R*T))-((Kmt_dinit*IPDINITaq)/(H_dinit*R*T))-(Q_dinit*kaq_dinit*IPDINITaq)
  aq_tend(:, iHNO3) = aq_tend(:, iHNO3) + 2*(Q_dinit*kaq_dinit*IPDINITaq)
  aq_tend(:, iTETROL) =  aq_tend(:, iTETROL) + (Q_dinit*kaq_dinit*IPDINITaq)
  
!  if (print_debug_output) then
!     add_text=.true.
!     do k = 1,nzm
        
!        if (Wl(k).gt.1.e-10) then
!           if (add_text) then
!              write (*,*) 'k, Wl, IEPOXg, IEPOXg+, IEPOXa, IEPOXa+, TETROLg, TETROLg+, TETROLa, TETROLa+'
!              add_text = .false.
!           end if   
              
!              write(*,*) k, Wl(k), gas_input(k,iIEPOX), gas_tend(k,iIEPOX), aq_input(k,iIEPOX), aq_tend(k,iIEPOX) !
!           !  gas_input(k,iTETROL), aq_input(k,iTETROL), gas_tend(k,iTETROL), aq_tend(k,iTETROL)
!        end if
!     end do   
!  end if              

  
end subroutine isop1nit_aqueous_tendencies

end module chem_aqueous
