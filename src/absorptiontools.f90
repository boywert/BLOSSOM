module absorptiontools
  use common_vars
  implicit none

contains
  subroutine n_cal(zcoll,n,r,rho)
    implicit none
    integer :: i,n
    real(kind=8) :: rho(0:max_size), r(0:max_size)
    real(kind=8) :: Delta_c, Delta_tophat,rho_crit_z,rho0,volume,mass,den_bar,Mxtil,zeta_x,zcoll
#ifdef USERHO178
    Delta_tophat = 18.*pi**2 
    rho_crit_z = rho_crit_0*(lambda0+Omega_0*(zcoll+1.)**3. &    
         +(1.-lambda0-Omega_0)*(1.+zcoll)**2.)
    Delta_c = (etaSUS/etaTIS)**3.*Delta_tophat
    rho0 = 6.*pi**2.*(b_T/5.)**3. * Mttil**2. * rho_crit_z

    open(3,file=trim(den_profile_file)) !contains dimensionless TIS density
    ! profile in the first 2 columns     
    i = 1
    r(0) = 0.0 
    rho(0) = 1.0d0
    !read the TIS profile in

    do while(r(i-1) .le. zeta_t)
       read(3,*) r(i), rho(i)
       i = i + 1
    enddo
    n = i-1
    close(3)

    mass = 0.
    volume = 0.

    do i= 1,n
       volume = 4./3.*pi*r(i)**3.
       mass = mass + rho0*(rho(i)+rho(i-1))/2. * 4./3.*pi*(r(i)**3.-r(i-1)**3.)
       den_bar = mass/volume/rho_crit_z
       if(den_bar <= Delta_tophat) goto 128
    end do
 128   n=i
    zeta_x = r(n)
    Mxtil = (Delta_tophat/Delta_c)*(Mttil/zeta_t**3.)*zeta_x**3. 
    !set new values 
    Mttil = Mxtil
    zeta_t = zeta_x
    if(rank==0) print*, "MXtil:",Mttil,"zeta_X:",zeta_t
#endif





    open(3,file=trim(den_profile_file)) !contains dimensionless TIS density
    ! profile in the first 2 columns     
    i = 1
    r(0) = 0.0 
    rho(0) = 1.0d0
    !read the TIS profile in

    do while(r(i-1) .le. zeta_t)
       read(3,*) r(i), rho(i)
       i = i + 1
    enddo
    n = i-1
    close(3)
    return
  end subroutine n_cal
  
  function find_T_k(M0,zcoll)
    implicit none
    real(kind=8) :: rho_crit_z, Delta_c,nu_z
    real(kind=8) :: r_t, r0, rho0, sigma_V, find_T_k, M0, zcoll, M_J



    nu_z = nu0/(1.d0+zcoll)
    Delta_c = 18.*pi**2 

    rho_crit_z = rho_crit_0*(lambda0+Omega_0*(zcoll+1.)**3. &    
         +(1.-lambda0-Omega_0)*(1.+zcoll)**2.)

#ifndef USERHO178
    Delta_c = (etaSUS/etaTIS)**3.*Delta_c
#endif
    
    M_J = 5.7d3*(Omega_0*h**2/0.15)**(-0.5)* &     
         (Omega_b*h**2/0.02)**(-0.6)*((1.0d0+zcoll)/10.)**1.5 !M_J in M_sol 

    r_t  = (3.*M0*M_sol/(4.*pi*Delta_c*rho_crit_z))**(1./3.) 

    r0  = r_t/zeta_t
    rho0 = zeta_t**3/3./Mttil*Delta_c*rho_crit_z
    sigma_V = sqrt(4.*pi*grav_const*rho0*r0**2)

    find_T_k = amass/boltzk*sigma_V**2 !assuming fully neutral
    return
  end function find_T_k


  subroutine tau_cal(M0,zcoll,impact_param,n,r,rho,sigma_V,areatau0,tau0)
    use omp_lib
    implicit none
    integer :: size
    real(kind=8) :: rho(0:max_size), r(0:max_size), n_HI(0:max_size)
    real(kind=8) :: T_S(0:max_size), y_c(0:max_size), y_H, y_e, ne
    real(kind=8) :: f(0:max_size)
    real(kind=8) :: tau_tot(0:max_size), check, D_A, V_A, theta, omegabh2
    real(kind=8) :: rho_crit_z
    real(kind=8) :: diff, x, T_CMB_z
    real(kind=8) :: delta_nu, nu_z, kappa_coeff
    real(kind=8) :: r_t, r0, rho0, sigma_V, T_k, M0, zcoll, dnum, M_J
    real(kind=8) :: Delta_c,tau0,areatau0,absorption,impact_param

    real(kind=8) :: shape, kappa1, delTb_bol, bol_abs_nu
    real(kind=8) :: X1(8), neutr_fraction, delOmega, bol_absorp
    !real(kind=8) :: lgtau(0:100),area_tau0(0:100), tau_current
    real(kind=8),allocatable :: tau(:,:)
    real(kind=4) :: start_time, stop_time
    integer :: area_flag
    integer :: i, n, j, q, niter, k, upper_r

    call cpu_time(start_time)
    nu_z = nu0/(1.d0+zcoll)
    Delta_c = 18.*pi**2 

    rho_crit_z = rho_crit_0*(lambda0+Omega_0*(zcoll+1.)**3. &    
         +(1.-lambda0-Omega_0)*(1.+zcoll)**2.)
#ifndef USERHO178
    Delta_c = (etaSUS/etaTIS)**3.*Delta_c
#endif

    M_J = 5.7d3*(Omega_0*h**2/0.15)**(-0.5)* &     
         (Omega_b*h**2/0.02)**(-0.6)*((1.0d0+zcoll)/10.)**1.5 !M_J in M_sol 

    r_t  = (3.*M0*M_sol/(4.*pi*Delta_c*rho_crit_z))**(1./3.) 
    
    r0  = r_t/zeta_t
    rho0 = zeta_t**3/3./Mttil*Delta_c*rho_crit_z
    sigma_V = sqrt(4.*pi*grav_const*rho0*r0**2)

    T_k = amass/boltzk*sigma_V**2 !assuming fully neutral

    neutr_fraction = nH
    delta_nu=sqrt(2.*mu*pi)*nu0*sigma_V/c !thermal width of the 21-cm line  

    T_CMB_z = 2.73d0*(1.0d0+zcoll)    
    kappa1 = 3.0d0*c**2./32./pi*A10*T_star
    kappa_coeff = kappa1/nu0**2./delta_nu
    do i = 0,n
       n_HI(i)=neutr_fraction*Omega_b/Omega_0*rho(i)*rho0/amass
       y_H = 7.345d-5*n_HI(i)*T_k**(8.781-4.7755*log10(T_k) &        
            +1.075*(log10(T_k))**2-0.091*log10(T_k)**3) !combined fit
       y_e =0.0d0
       y_c(i) = y_H+y_e
       T_S(i) = (T_CMB_z+y_c(i)*T_k)/(1.0d0+y_c(i))
    end do
    delTb_bol =0.0d0
    do i = 0,n    
       f(i)=kappa1/delta_nu*n_HI(i)/nu0**2/T_S(i)
    end do
    allocate(tau(0:max_size,-max_size:max_size))
    tau(n,n)   = 0.0d0
    tau_tot(n) = tau(n,n)
    
    !$omp parallel private(j,diff) 
    !$omp do
    do i=n-1,0,-1
       tau(i,-n)=0.0d0
       do j=-n,-i-1
          !Trapezoidal rule
          diff = sqrt(r(-j)*r(-j)-r(i)*r(i)) &
               - sqrt(r(-j-1)*r(-j-1)-r(i)*r(i))
          tau(i,j+1) = tau(i,j) + (f(-j)+f(-j-1))*diff/2.0d0
       enddo
       tau(i,i) = tau(i,-i)
       do j=i+1,n 
          tau(i,j)=2.0d0*tau(i,i) - tau(i,-j)
       end do
       tau_tot(i) = 2.0d0*tau(i,i)*r0 !total optical depth at r_i
    enddo
    !$omp end do
    !$omp end parallel
    deallocate(tau)
    do i=1,n
       if(impact_param/r0 < r(i)) then
          goto 127
       endif
    end do

127 tau0 = (tau_tot(i)-tau_tot(i-1))/(r(i) - r(i-1))*(impact_param/r0-r(i-1)) + tau_tot(i-1) 


    areatau0=0.0d0
    do i=1,n
       areatau0 = areatau0 + (tau_tot(i-1)*r(i-1)+tau_tot(i)*r(i)) * (r(i)-r(i-1))
    end do
    areatau0 =  areatau0/zeta_t**2


    !print*,T_k,exp(-1*tau0)
    !absorption = 1.0d0 - exp(-1.d0*tau0)

    ! do j=0,100
    !    tau_current=10**lgtau(j)
    !    do i=0,n-1
    !       !print*,'t_c',tau_current,'t_t-',tau_tot(i-1),'t_t',tau_tot(i)
    !       if(tau_current.lt.tau_tot(i-1).and. tau_current.gt.tau_tot(i)) then
    !          area_tau0(j)=(r(i)/zeta_t)**2*3.1415*(r_t/3.086e24)**2 
    !          go to 111
    !       end if
    !    end do
    !111    continue
    ! end do
    return
  end subroutine tau_cal


  subroutine lorentz_term_cal(n,f0,delta_f,gamma,lorentz_term)
    implicit none
    integer :: n,i
    real(kind=8) :: delta_f,gamma,f0
    real(kind=8) :: f(-n:n)
    real(kind=8) :: lorentz_term(-n:n)

    do i= -n,n
       f(i) = real(i,8)*delta_f+f0
       lorentz_term(i) = gamma/(4.*pi*f(i)**2. + (gamma/2.)**2.)
       !print*, f(i), lorentz_term(i)
    end do
    return
  end subroutine lorentz_term_cal


  subroutine doppler_term_cal(n,f0,delta_f,sigma,doppler_term)
    implicit none
    integer :: n,i
    real(kind=8) :: delta_f,sigma,sigma_f,f0
    real(kind=8) :: f(-n:n), doppler_term(-n:n)
    do i=-n,n
       f(i) = real(i,8)*delta_f+f0
       sigma_f = sigma*f0/c*sqrt(2./3.)
       doppler_term(i) = exp(-1.*(f(i)/sigma_f)**2.)/sqrt(pi)/sigma
       !print*, f(i), doppler_term(i) 
    end do
    return
  end subroutine doppler_term_cal

  subroutine voigt_term_cal(n,f0,delta_f,gamma,sigma,voigt_term)
    implicit none
    integer :: n,i,j
    real(kind=8) :: delta_f,sigma,gamma,f0
    real(kind=8) :: f(-n:n), doppler_term(-n:n),lorentz_term(-n:n)
    real(kind=8) :: voigt_term(-n:n),sum

    call doppler_term_cal(n,f0,delta_f,sigma,doppler_term)
    !print*, doppler_term(-n:n)
    call lorentz_term_cal(n,f0,delta_f,gamma,lorentz_term)
    !print*, lorentz_term

    do i=-n,n
       f(i) = real(i,8)*delta_f + f0
       voigt_term(i) = 0.
       do j=-n,n
          if(i-j >= -n .and. i-j <= n) then
             voigt_term(i) = voigt_term(i) + lorentz_term(i-j)*doppler_term(j)
          endif
       end do
    end do
    sum = 0.
    do i =-n,n
       sum = sum + voigt_term(i)*delta_f
    end do
    voigt_term(-n:n) = voigt_term(-n:n)/sum
    do i =-n,n
       print*, f(i), voigt_term(i)
    end do
    return
  end subroutine voigt_term_cal

  function voigt_FWHM(gamma,sigma,nu0)
    implicit none
    real(kind=8) :: gamma, sigma, voigt_FWHM,nu0
    real(kind=8) :: f_L, f_D, nu_th
    nu_th = nu0/c*sqrt(2./3.)*sigma
    f_D = sqrt(alog(2.)*2.)*2.*nu_th
    f_L = gamma/2./pi
    voigt_FWHM = 0.5346*f_L + sqrt(0.2166*f_L**2. + f_D**2.)
    !print*,f_D,f_L,voigt_FWHM
  end function voigt_FWHM


  function absorption_line(nu_min,nu_max,nu_binsize,totalbin,nu_center,fwhm,tau)
    implicit none
    integer :: totalbin,left_block,right_block
    real(kind=8) :: nu_min, nu_max, nu_binsize, nu_center, fwhm, tau
    real(kind=8) :: absorption_line(totalbin)
    real(kind=8) :: left_edge,right_edge
    
    if((nu_max-nu_min)/nu_binsize .ne. totalbin) call abort("Bins don't match")
    absorption_line(1:totalbin) = 0.d0
    
    left_edge = nu_center-fwhm/2.0d0
    if(left_edge .le. nu_min) then
       left_block = 1
       absorption_line(left_block) =  1.d0 - exp(-tau)
    else if (left_edge > nu_min .and. left_edge < nu_max) then
       left_block = ceiling((left_edge - nu_min)/nu_binsize)
       absorption_line(left_block) =  &
            (1.d0 - exp(-tau))*(real(left_block,8)*nu_binsize + nu_min - left_edge)/nu_binsize
    else
       call abort("Halo is not in the box")
    end if

    right_edge = nu_center+fwhm/2.0d0
    if(right_edge .ge. nu_max) then
       right_block = totalbin
       absorption_line(right_block) =  1.d0 - exp(-tau)
    else if (right_edge > nu_min .and. right_edge < nu_max) then
       right_block = ceiling((right_edge - nu_min)/nu_binsize)
       absorption_line(right_block) = &
            (1.d0 - exp(-tau))*(nu_binsize-(real(right_block,8)*nu_binsize + nu_min - right_edge))/nu_binsize
    else
       call abort("Halo is not in the box")
    end if   

    absorption_line(left_block+1:right_block-1) = 1.d0 - exp(-tau)
    return
  end function absorption_line

  subroutine save_absorptionline(absorpline,n,nu_min,nu_binsize,z_s,filename)
    implicit none
    integer :: n,i
    real(kind=8) :: nu_min, nu_binsize
    real(kind=8) :: absorpline(n)
    character(len=100) :: z_s
    character(len=256) :: filename
    filename = adjustl(filename)
    print*, "Print absorpline to ",(result_path)//z_s(1:len_trim(z_s))//'/'//trim(filename)
    z_s = adjustl(z_s)
    call system("rm -f "//trim(result_path)//z_s(1:len_trim(z_s))//'/'//trim(filename))
    open (unit=53,file=trim(result_path)//z_s(1:len_trim(z_s))//'/'//trim(filename),STATUS = 'NEW')
    do i=0, n
       write(53,*) (real(i,8)-0.5)*real(nu_binsize,8)+nu_min, 1.0-absorpline(i)
    end do
    close(53)
  end subroutine save_absorptionline

end module absorptiontools

