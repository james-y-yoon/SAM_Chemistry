MODULE chem_isoprene_rate_functions
  USE chem_isoprene_Parameters, only: NREACT
  IMPLICIT NONE
  
CONTAINS

subroutine compute_rate_constants(nzm, TempK, M, qv0, rate_const)
    integer, intent(in) :: nzm
    real, intent(in) :: TempK(nzm), M(nzm), qv0(nzm)   ! temperature and air molecules/cm3 profiles
    real, intent(inout) :: rate_const(nzm, NREACT)

    rate_const(:,1) = k_T(2.10E-11,465.,TempK)
    rate_const(:,2) = KRO2HO2(TempK)
    rate_const(:,3) = KRO2NO(TempK)
    rate_const(:,7) = k_T(3.3e12,-8660.,TempK)
    rate_const(:,8) = KRO2HO2(TempK)
    rate_const(:,9) = KRO2NO(TempK)
    rate_const(:,12) = KRO2HO2(TempK)
    rate_const(:,13) = KRO2NO(TempK)
    rate_const(:,17) = k_T(8.e-12,380.,TempK)
    rate_const(:,18) = KFPAN(M, TempK)
    rate_const(:,19) = KBPAN(M, TempK)
    rate_const(:,20) = KRO2NO(TempK)
    rate_const(:,22) = k_T(1.05E-11,465.,TempK)
    rate_const(:,23) = k_T(2.6E-12,610.,TempK)
    rate_const(:,24) = k_T(1.05e-11,465.,TempK)
    rate_const(:,25) = k_T(1.05e-11,465.,TempK)
    rate_const(:,26) = k_T(2.07E-12,-1400.,TempK)
    rate_const(:,27) = k_T(1.7E-12,-940.,TempK)
    rate_const(:,28) = KHO2O3(TempK)
    rate_const(:,29) = k_T(3.45E-12,270.,TempK)
    rate_const(:,30) = KOHNO2(M, TempK)
    rate_const(:,31) = KMT11(M, TempK)
    rate_const(:,32) = k_T(1.8E-11,110.,TempK)
    rate_const(:,33) = KMT09(M, TempK)
    rate_const(:,34) = k_T(3.2E-13,690.,TempK)
    rate_const(:,35) = k_T(4.8E-11,250.,TempK)
    rate_const(:,36) = k_T(5.4E-12,135.,TempK)
    rate_const(:,37) = KHO2HO2(M, TempK, qv0)
    rate_const(:,38) = KMT05(M)
    rate_const(:,39) = k_T(1.85E-12,-1690.,TempK)
    rate_const(:,40) = KRO2HO2ROOH(TempK)
    rate_const(:,41) = KRO2NO(TempK)
    ! rate_const(:,4) = constant rate coefficient
    ! rate_const(:,5) = constant rate coefficient
    ! rate_const(:,6) = constant rate coefficient
end subroutine compute_rate_constants

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!! FUNCTIONS !!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

ELEMENTAL REAL FUNCTION KRO2HO2(T)
! T in K
!    use isoprene_k_T
    REAL, intent(IN) :: T
    REAL :: C1, C2
    C1 = 1.4e-13
    C2 = 1330.
    KRO2HO2 = k_T(C1, C2, T)
END FUNCTION KRO2HO2

ELEMENTAL REAL FUNCTION KRO2HO2ROOH(T)
! T in K
!    use isoprene_k_T
    REAL, intent(IN) :: T
    REAL :: C1, C2
    C1 = 2.91e-13
    C2 = 1300.
    KRO2HO2ROOH = k_T(C1, C2, T)
END FUNCTION KRO2HO2ROOH

ELEMENTAL REAL FUNCTION KHO2O3(T)
    ! T in K
    REAL, intent(IN) :: T
    REAL :: C1, C2
    C1 = 2.03e-16
    C2 = 693.
    KHO2O3 = k_T(C1, C2, T) * (T/300.)**4.57
END FUNCTION KHO2O3

ELEMENTAL REAL FUNCTION KOHNO2(M, T)
    ! T in K
    ! M in air molecules/cm3

    REAL, intent(IN) :: M, T
    REAL :: K80, K8I, KR8, FC8, NC8, F8

    K80 = 3.2e-30*M*(T/300)**-4.5
    K8I = 3.0e-11
    KR8 = K80/K8I
    FC8 = 0.41
    NC8 = 0.75-1.27*(log10(FC8))
    F8 = 10**(log10(FC8)/(1+(log10(KR8)/NC8)**2))
    KOHNO2 = (K80*K8I)*F8/(K80+K8I)
END FUNCTION KOHNO2

ELEMENTAL REAL FUNCTION KMT11(M, T)
    ! T in K
    ! M in air molecules/cm3

    REAL, intent(IN) :: M, T
    REAL :: K1, K3, K4, K2

    K1 = 2.40e-14*exp(460./T)
    K3 = 6.50e-34*exp(1335./T)
    K4 = 2.70e-17*exp(2199./T)
    K2 = (K3*M)/(1+(K3*M/K4))
    KMT11 = K1 + K2
END FUNCTION KMT11

ELEMENTAL REAL FUNCTION KHO2HO2(M, T, qv)
    ! T in K
    ! M in air molecules/cm3

    REAL, intent(IN) :: M, T, qv
    REAL :: unit_f, H2O, f0, k1, k2

    H2O = qv / 1000 * (28.97 / 18.02) * M ! this quantity is in [molecules cm-3]
    f0 = 1 + 1.4e-21 * H2O * exp(2200./T)
    k1 = 2.2e-13 * exp(600./T) * f0
    k2 = 1.9e-33 * exp(980./T) * M * f0
    KHO2HO2 = k1 + k2
END FUNCTION KHO2HO2

ELEMENTAL REAL FUNCTION KMT09(M, T)
    ! T in K
    ! M in air molecules/cm3

    REAL, intent(IN) :: M, T
    REAL :: K90, K9I, KR9, FC9, NC9, F9

    K90 = 1.4e-31*M*(T/300.)**(-3.1)
    K9I = 4.0e-12
    KR9 = K90/K9I
    FC9 = 0.4
    NC9 = 0.75-1.27*(log10(FC9))
    F9 = 10**(log10(FC9)/(1+(log10(KR9)/NC9)**2))
    KMT09 = (K90*K9I)*F9/(K90+K9I)
END FUNCTION KMT09

ELEMENTAL REAL FUNCTION KMT05(M)
    ! M in air molecules/cm3
    REAL, intent(IN) :: M

    KMT05 = 1.44e-13 * (1 + M/4.2E19)
END FUNCTION KMT05

ELEMENTAL REAL FUNCTION KRO2NO(T)
    REAL, intent(IN) :: T
    REAL :: C1, C2
    C1 = 2.7e-12
    C2 = 360.
    KRO2NO = k_T(C1, C2, T)
END FUNCTION KRO2NO 

ELEMENTAL REAL FUNCTION k_T(C1, C2, T)
! T in Kelvin
    REAL, intent(IN) :: C1, C2, T
    k_T =  C1 * EXP(C2/T)
END FUNCTION k_T

ELEMENTAL REAL FUNCTION KFPAN(M, T)
! M in air molecules/cm3
! T in Kelvin
    REAL, intent(IN) :: M, T
    REAL :: KC0, KC1, KRC, FCC, NC, FC
    KC0 = 3.28e-28 * M * (T/300)**(-6.87)
    KC1 = 1.25e-11 * (T/300)**(1.105)
    KRC = KC0/KC1
    FCC = 0.30
    NC = 0.75 - 1.27*LOG10(FCC)
    FC = 10**(LOG10(FCC)/(1 + (LOG10(KRC)/NC)**2))
    KFPAN = (KC0 * KC1) * FC/(KC0 + KC1)

END FUNCTION KFPAN

ELEMENTAL REAL FUNCTION KBPAN(M, T)
! M in air molecules/cm3
! T in Kelvin
    REAL, intent(IN) :: M, T
    REAL :: KD0, KD1, KRD, FCD, NCD, FD
    KD0 = 1.10e-5*M*EXP(-10100/T)
    KD1 = 1.90e17*EXP(-14100/T)
    KRD = KD0/KD1
    FCD = 0.30
    NCD = 0.75 - 1.27*LOG10(FCD)
    FD = 10**(LOG10(FCD)/(1 + (LOG10(KRD)/NCD)**2))
    KBPAN = (KD0 * KD1) * FD/(KD0 + KD1)
END FUNCTION KBPAN

END MODULE chem_isoprene_rate_functions
