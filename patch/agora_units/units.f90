subroutine units(scale_l,scale_t,scale_d,scale_v,scale_nH,scale_T2)
  use amr_commons
  use hydro_commons
  use cooling_module
  real(dp)::scale_nH,scale_T2,scale_t,scale_v,scale_d,scale_l
  !-----------------------------------------------------------------------
  ! Conversion factors from user units into cgs units
  ! For gravity runs, make sure that G=1 in user units.
  !-----------------------------------------------------------------------

  !Units are now in kpc, Gyr, km/s and time is 0.98 Gyr

  ! scale_d converts mass density from user units into g/cc
  !scale_d =6.80714757d-23  !1d9 msun
  scale_d=1.573391528d-26   !2.32e5 msun --> 0.977 Gyr etc
  if(cosmo) scale_d = omega_m * rhoc *(h0/100.)**2 / aexp**3

  ! scale_l converts distance from user units into cm
  scale_l = 3.086568025d21
  if(cosmo) scale_l = aexp * boxlen_ini * 3.08d24 / (h0/100)

  ! scale_t converts time from user units into seconds
  scale_t = 1.0/sqrt((6.67d-8)*scale_d)  
  if(cosmo) scale_t = aexp**2 / (h0*1d5/3.08d24)

  ! scale_v convert velocity in user units into cm/s
  scale_v = scale_l / scale_t

  ! scale_T2 converts (P/rho) in user unit into (T/mu) in Kelvin
  scale_T2 = mH/kB * scale_v**2

  ! scale_nH converts rho in user units into nH in H/cc
  scale_nH = X/mH * scale_d

end subroutine units
