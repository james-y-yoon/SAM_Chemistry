subroutine nudging()
	
use vars
use params
use microphysics, only: micro_field, index_water_vapor
implicit none

real coef, coef1
integer i,j,k
	
!bloss: nudging that adjusts to inversion height
real :: nudge_ramp(nzm)
real :: pii

real :: itau_transient, itau_transient_aloft
real :: itauz_t(nzm), itauz_q(nzm)
real :: zramp_transient(nzm), itau_transient_z(nzm)
real :: tmp_zinv_obs
real :: dist_zinv

! tendency for CGILS qfloor
!   This is a strictly positive nudging that prevents the above-inversion
!   air from becoming unrealistically dry (mostly due to horizontal advective
!   drying that is applied above the inversion).
real dqdt_qfloor

real, external :: get_inversion_height

call t_startf ('nudging')

tnudge = 0.
qnudge = 0.
unudge = 0.
vnudge = 0.

if(donudging_transient) then
    if(masterproc) write(*,*) 'nudge_d0, day, nudge_df = ', transient_nudging_start, day, transient_nudging_end 
 !bloss: Define a vertically-uniform inverse nudging timescale
  !   that can be switched on/off during the simulation.
  !   The nudging of temperature and moisture will always 
  !   be at least as strong as itau_transient.
  if((day.gt.transient_nudging_start).AND.(day.lt.transient_nudging_end)) then
    donudging_tq = .true.
    itau_transient = 1./tau_transient_nudging
    itau_transient_aloft = 1./600. ! ten minute nudging timescale above inversion
    
    !smooth start/finish to transient nudging
    itau_transient = itau_transient &
         *0.5*(1.-cos(pi*MAX(0., MIN(1., (day-transient_nudging_start)/transient_nudging_ramp ) ) ) ) &
         *0.5*(1.-cos(pi*MAX(0., MIN(1., (transient_nudging_end-day)/transient_nudging_ramp ) ) ) ) 
!!$    if(masterproc) write(*,*) 'tau0, tauf = ', (day-transient_nudging_start)/transient_nudging_ramp, &
!!$         (transient_nudging_end-day)/transient_nudging_ramp 
!!$    if(masterproc) write(*,*) 'fac0, facf = ', &
!!$         0.5*(1.-cos(pi*MAX(0., MIN(1., (day-transient_nudging_start)/transient_nudging_ramp ) ) ) ), &
!!$         0.5*(1.-cos(pi*MAX(0., MIN(1., (transient_nudging_end-day)/transient_nudging_ramp ) ) ) ) 

    ! turn off nudging within 50m of observed inversion, so that model can adopt
    !   its own inversion structure if close to correct inversion height)
    tmp_zinv_obs = get_inversion_height(nzm,z,pres,tg0+gamaz,tg0,qg0) 

    zramp_transient(:) = 1.
    do k = 1,nzm
      dist_zinv = ABS(z(k)-tmp_zinv_obs) 
      IF (dist_zinv.lt.50.) then
        zramp_transient(k) = 0.
      elseif (dist_zinv.lt.100.) then
        zramp_transient(k) = 0.5*(1. - cos(pi*(dist_zinv-50.)/50.))
      end IF
      if(z(k).lt.tmp_zinv_obs) then
        itau_transient_z(k) = itau_transient*zramp_transient(k)
      else
        ! different, shorter nudging time above inversion
        itau_transient_z(k) = itau_transient_aloft*zramp_transient(k)  
      end if
    end do

    if(icycle.eq.1.AND.mod(nstep,10).eq.1) then
      if(masterproc) write(*,*) 'Transient nudging on with itau(day) = ', itau_transient
    end if
  else
    itau_transient = 0.
    itau_transient_z(:) = 0.
  end if
end if


coef = 1./tauls

if(donudging_uv) then
  if(nudge_to_sounding_winds) then
    do k=1,nzm
      if(z(k).ge.nudging_uv_z1.and.z(k).le.nudging_uv_z2) then
        unudge(k)=unudge(k) - (u0(k)-usounding0(k))*coef
        vnudge(k)=vnudge(k) - (v0(k)-vsounding0(k))*coef
        do j=1,ny
          do i=1,nx
             dudt(i,j,k,na)=dudt(i,j,k,na)-(u0(k)-usounding0(k))*coef
             dvdt(i,j,k,na)=dvdt(i,j,k,na)-(v0(k)-vsounding0(k))*coef
          end do
        end do
      end if
    end do
  else
    do k=1,nzm
      if(z(k).ge.nudging_uv_z1.and.z(k).le.nudging_uv_z2) then
        unudge(k)=unudge(k) - (u0(k)-ug0(k))*coef
        vnudge(k)=vnudge(k) - (v0(k)-vg0(k))*coef
        do j=1,ny
          do i=1,nx
             dudt(i,j,k,na)=dudt(i,j,k,na)-(u0(k)-ug0(k))*coef
             dvdt(i,j,k,na)=dvdt(i,j,k,na)-(v0(k)-vg0(k))*coef
          end do
        end do
      end if
    end do
  end if
endif

if(donudging_tq.or.donudging_t.or.donudging_q) then

  ! set up nudging heights automatically by tracking inversion height
  if(dovariable_tauz) then
    nudging_t_z1 = get_inversion_height(nzm,z,pres,t0,tabs0,q0) &
         + variable_tauz_offset_above_inversion
    nudging_t_z1 = MAX(nudging_t_z1, variable_tauz_minimum_height)
    nudging_t_zramp = variable_tauz_thickness_of_onset

    nudging_q_z1 = nudging_t_z1
    nudging_q_zramp = nudging_t_zramp
  end if

  ! vertically-varying nudging amplitude for temperature
  nudge_ramp(:) = 0.

  pii = acos(-1.)
  do k=1,nzm
    if(z(k).ge.nudging_t_z1.and.z(k).le.nudging_t_z2) then
      !nudging will be applied
      nudge_ramp(k) = 1.

      if(z(k).lt.nudging_t_z1+nudging_t_zramp) then
        ! gradual onset of nudging between z1 and z1+zramp
        nudge_ramp(k) = 0.5*(1-cos(pii*(z(k)-nudging_t_z1) &
             / nudging_t_zramp ) )
      elseif(z(k).gt.nudging_t_z2 - nudging_t_zramp) then
        ! gradual falloff of nudging between z2 and z2-zramp
        nudge_ramp(k) = 0.5*(1-cos(pii*(nudging_t_z2 - z(k)) &
             / nudging_t_zramp ) )
      end if
    end if
  end do

  !bloss: Use itauz_t(1:nzm) to keep track of vertically-varying inverse nudging timescale
  coef = 1./tautqls
  itauz_t(:) = coef*nudge_ramp(:)

  if(donudging_transient) then
    !bloss: Apply vertically-uniform nudging where stronger than standard nudging
    do k = 1,nzm
      itauz_t(k) = MAX(itauz_t(k), itau_transient_z(k))
    end do
  end if

  ! vertically-varying nudging amplitude for moisture
  nudge_ramp(:) = 0.

  pii = acos(-1.)
  do k=1,nzm
    if(z(k).ge.nudging_q_z1.and.z(k).le.nudging_q_z2) then
      !nudging will be applied
      nudge_ramp(k) = 1.

      if(z(k).lt.nudging_q_z1+nudging_q_zramp) then
        ! gradual onset of nudging between z1 and z1+zramp
        nudge_ramp(k) = 0.5*(1-cos(pii*(z(k)-nudging_q_z1) &
             / nudging_q_zramp ) )
      elseif(z(k).gt.nudging_q_z2 - nudging_q_zramp) then
        ! gradual falloff of nudging between z2 and z2-zramp
        nudge_ramp(k) = 0.5*(1-cos(pii*(nudging_q_z2 - z(k)) &
             / nudging_q_zramp ) )
      end if
    end if
  end do

  !bloss: Use itauz_q(1:nzm) to keep track of vertically-varying inverse nudging timescale
  coef = 1./tautqls
  itauz_q(:) = coef*nudge_ramp(:)

  if(donudging_transient) then
    !bloss: Apply vertically-uniform nudging where stronger than standard nudging
    do k = 1,nzm
      itauz_q(k) = MAX(itauz_q(k), itau_transient_z(k))
    end do
  end if

end if


if(donudging_tq.or.donudging_t) then
    do k=1,nzm
!bloss: Apply nudging for all levels -- itauz_t(k) will be zero where no nudging is applied.
        tnudge(k)=tnudge(k) -(t0(k)-tg0(k)-gamaz(k))*itauz_t(k)

        do j=1,ny
          do i=1,nx
            t(i,j,k)=t(i,j,k)-(t0(k)-tg0(k)-gamaz(k))*dtn*itauz_t(k)
          end do
        end do
    end do
endif

if(donudging_tq.or.donudging_q) then
    do k=1,nzm
!bloss: Apply nudging for all levels -- itauz_t(k) will be zero where no nudging is applied.
      qnudge(k)=qnudge(k) -(q0(k)-qg0(k))*itauz_q(k)
      do j=1,ny
        do i=1,nx
          micro_field(i,j,k,index_water_vapor)=micro_field(i,j,k,index_water_vapor)-(q0(k)-qg0(k))*dtn*itauz_q(k)
        end do
      end do
    end do
endif

!bloss: option for CGILS cases to avoid over-drying of sounding above BL
if(doenforce_cgils_qfloor) then
  do k = 1,nzm
    if((z(k).lt.ztop_qfloor).AND.(q0(k).lt.qfloor)) then
      dqdt_qfloor = (qfloor-q0(k))/tau_qfloor
      qnudge(k) = qnudge(k) + dqdt_qfloor
      do j=1,ny
        do i=1,nx
          micro_field(i,j,k,index_water_vapor) &
               = micro_field(i,j,k,index_water_vapor) + dtn*dqdt_qfloor
        end do
      end do
    end if
  end do
end if

call t_stopf('nudging')

end subroutine nudging
