program test_iepox_parcel

  use iepox_aqueous

  implicit none
  integer :: nzm
  real :: dt, tt

  real, dimension(1) :: T   !  K 
real, dimension(1) :: p   ! atm
real, dimension(1) :: Rp  ! droplet radius (m)
real, dimension(1) :: Wl  ! liquid water content vol H2O/vol air
real, dimension(1) :: H   !  M of H+
real, dimension(1) :: NO3   !  M 
real, dimension(1) :: SO4   !  M 
real, dimension(1) :: HSO4   !  M 

real, dimension(1) :: IEPOXg, IEPOXaq
real, dimension(1) :: PRODUCTaq
real, dimension(1) :: TETROLg, TETROLaq
real, dimension(1) :: IEPOX_NO3g, IEPOX_NO3aq
real, dimension(1) :: IEPOX_SO4g, IEPOX_SO4aq

real, dimension(1) :: IEPOXg_tend, IEPOXaq_tend
real, dimension(1) :: PRODUCTaq_tend
real, dimension(1) :: TETROLg_tend, TETROLaq_tend
real, dimension(1) :: IEPOX_NO3g_tend, IEPOX_NO3aq_tend
real, dimension(1) :: IEPOX_SO4g_tend, IEPOX_SO4aq_tend

nzm = 1

T = 280.0
p = 0.85
Rp = 3.e-5 
Wl = 8.e-7 !m^3/m^3 

H = 1e-5 !M of H+ (pye et al) (placeholder values for all other concentrations)
NO3 = 1e-7 !M of nucleophile NO3-
SO4 = 1e-7 !M of nucleophile SO4-2 
HSO4 = 1e-10 !M of general acid HSO4-

IEPOXg = 0
IEPOXaq = 0
PRODUCTaq = 0
TETROLg=0
TETROLaq=0
IEPOX_NO3g=0
IEPOX_NO3aq=0
IEPOX_SO4g=0
IEPOX_SO4aq=0


! mimic IEPOX model initialization
IEPOXg = 1e-9
IEPOXaq = 0.0

dt = 0.01

! do timesteps

print '(120A)', 'time             IEPOXg         IEPOXaq     PRODUCTaq     TETROLg     TETROLaq     IEPOX_NO3g     IEPOX_NO3aq    IEPOX_SO4g    IEPOX_SO4aq'

do tt = 0, 1000, dt

   call iepox_aqueous_tendencies(nzm, T, p, Rp, Wl, H, NO3, SO4, HSO4, &
     IEPOXg, IEPOXaq, PRODUCTaq, TETROLg, TETROLaq, &
     IEPOX_NO3g, IEPOX_NO3aq, IEPOX_SO4g, IEPOX_SO4aq, &
     IEPOXg_tend, IEPOXaq_tend, PRODUCTaq_tend, TETROLg_tend, TETROLaq_tend, &
     IEPOX_NO3g_tend, IEPOX_NO3aq_tend, IEPOX_SO4g_tend, IEPOX_SO4aq_tend)

  
   IEPOXg = IEPOXg + IEPOXg_tend * dt
   IEPOXaq = IEPOXaq + IEPOXaq_tend * dt
   PRODUCTaq = PRODUCTaq + PRODUCTaq_tend * dt
   TETROLg = TETROLg + TETROLg_tend * dt
   TETROLaq = TETROLaq + TETROLaq_tend * dt
   IEPOX_NO3g = IEPOX_NO3g + IEPOX_NO3g_tend * dt
   IEPOX_NO3aq = IEPOX_NO3aq + IEPOX_NO3aq_tend * dt
   IEPOX_SO4g = IEPOX_SO4g + IEPOX_SO4g_tend * dt
   IEPOX_SO4aq = IEPOX_SO4aq + IEPOX_SO4aq_tend * dt


   print  '(F8.2, 4x, 9(E14.4))', tt, IEPOXg, IEPOXaq, PRODUCTaq, TETROLg, TETROLaq, &
        IEPOX_NO3g, IEPOX_NO3aq, IEPOX_SO4g, IEPOX_SO4aq

end do   

end program test_iepox_parcel  
