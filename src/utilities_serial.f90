!## ALL ROUTINES ARE SERIAL
module utilities_serial
  implicit none
  
contains
  
  subroutine correlation_function_unitcube(rho,n_point, pos_in, xi_out, count)
    use common_vars
    use omp_lib
    use mpitools
    implicit none
    integer(kind=4) :: n_point,n_block,blockx,blocky,blockz,block,ix,iy,iz
    integer(kind=4), allocatable :: np_block(:)
    real(kind=4), allocatable :: block_pos(:,:)
    real(kind=4) :: pos_in(1:3,1:n_point)
    real(kind=4) :: xi_out(1:n_bin_xi)
    real(kind=4) :: rad(-2:n_bin_xi+2) 
    integer(kind=8) :: count(1:n_bin_xi),group(1:n_bin_xi)

    integer(kind=8) :: i, j, index, indx, m1, m2, k, m
    real(kind=4) :: dx,dy,dz,dist, r2, avgcount,didx,r3min,rho
    real(kind=4) :: x1,y1,z1,dr 

    real(kind=4) :: dmass,x,y,z,vel,lum,absorp
    integer(kind=8),allocatable :: subcount(:,:),subgroup(:,:)
    real(kind=4) :: smalldr,start_time,stop_time
    real(kind=8) :: sumtest
    call cpu_time(start_time)
    !if(n_point .gt. maxpoint_serial_xi) call abort
    smalldr = 1.0/n_smallcell
    !print*, "all"
    allocate(np_block(0:n_smallcell**3-1))
    allocate(block_pos(1:3,0:n_smallcell**3-1))
    np_block(0:n_smallcell**3-1) = 0
    !print*, "put points in blocks"
    do blockx=0,n_smallcell-1
       do blocky=0,n_smallcell-1
          do blockz=0,n_smallcell-1
             block = blockz + blocky*n_smallcell + blockx*n_smallcell**2
             block_pos(1:3,block) = (/ blockx, blocky, blockz  /) * smalldr
          end do
       end do    
    end do
    do m1=1,n_point
       blockx = pos_in(1,m1)/smalldr
       if(blockx == n_smallcell) blockx = n_smallcell-1
       blocky = pos_in(2,m1)/smalldr
       if(blocky == n_smallcell) blocky = n_smallcell-1
       blockz = pos_in(3,m1)/smalldr
       if(blockz == n_smallcell) blockz = n_smallcell-1
       block = blockz + blocky*n_smallcell + blockx*n_smallcell**2
       np_block(block) = np_block(block) + 1
    end do

    allocate(subcount(1:n_bin_xi,0:omp_workthreads-1))
    allocate(subgroup(1:n_bin_xi,0:omp_workthreads-1))

    dr=(max_radius_xi)/n_bin_xi


    do i=1,n_bin_xi
       rad(i)=(real(i)-.5)*dr
    end do
    subcount(1:n_bin_xi,0:omp_workthreads-1) = 0
    subgroup(1:n_bin_xi,0:omp_workthreads-1) = 0
    !print*, "start counting"
    !$omp parallel private(x1,y1,z1,m2,dx,dy,dz,r2,dist,indx,index) 
    !$omp do
    do m1=0,n_smallcell**3-1
       !print*,pos_in(1:3,m1)
       !if(mod(m1,10000) == 0) print*,m1,omp_get_thread_num(),block_pos(1:3,m1)
       x1=block_pos(1,m1)
       y1=block_pos(2,m1)
       z1=block_pos(3,m1)

       do m2=0,m1-1
          dx=abs(x1-block_pos(1,m2))
          dy=abs(y1-block_pos(2,m2))
          dz=abs(z1-block_pos(3,m2))
          r2=dx**2+dy**2+dz**2
          dist=(r2)**(0.5)
          index=ceiling(dist/dr)
          !if(index<1 .or. index > n_bin_xi)  print*,index
          !print*, m1,m2,index

          subcount(index,omp_get_thread_num())=2+subcount(index,omp_get_thread_num())
          subgroup(index,omp_get_thread_num())=2*np_block(m1)*np_block(m2)+subgroup(index,omp_get_thread_num())
          !          count(index)=count(index)+mass(m2)
          !    
       enddo
    end do
    !$omp end do
    !$omp end parallel

    !print*,"combine OMP"

    count(1:n_bin_xi) = 0
    group(1:n_bin_xi) = 0
    do i=0,omp_workthreads-1
       !print*, subcount(1:n_bin_xi,i)
       count(1:n_bin_xi) = count(1:n_bin_xi)+subcount(1:n_bin_xi,i)
       group(1:n_bin_xi) = group(1:n_bin_xi)+subgroup(1:n_bin_xi,i)
    end do
    deallocate(subcount)
    deallocate(subgroup)

    sumtest = 0.
    !print*,"calculate xi"
    do i=1,n_bin_xi
       sumtest = sumtest + real(count(i))*(smalldr**6)
       xi_out(i) = real(group(i)) / (real(count(i))*(rho**2*smalldr**6)) -1.
       !xi_out(i) = xi_out(i)*count(i)*(rho**2*dr**6)
    end do
#ifdef printoutputxi
    print*, "rho",rho
    print*,"n1*n2 in bins"
    print*, group
    print*, "count cells in bins"
    print*, count
    print*, "[n1*n2] = n1*n2/count"
    print*, group/count
    print*,"xi"
    print*, xi_out(1:n_bin_xi)/real(count(1:n_bin_xi)*(rho**2*smalldr**6))
    print*,"rho*rho* xi* dV1dV2"
    print*, xi_out(1:n_bin_xi)
    print*,"sum dV1dV2 in cell"
    print*,real(count(1:n_bin_xi)*(smalldr**6))
    print*, "all dV1dV2",sumtest
#endif !printoutputxi
    call cpu_time(stop_time)
    !print*,'sumtest',sumtest
    !print*, rank,"use ",stop_time-start_time, "to calculate"
    return
  end subroutine correlation_function_unitcube


  subroutine correlation_function_unity(rho,n_point, pos_in, xi_out)
    use common_vars
    use omp_lib
    implicit none
    integer(kind=4) :: n_point
    real(kind=4) :: pos_in(1:3,1:n_point)
    real(kind=4) :: xi_out(0:n_bin_xi)
    real(kind=4) :: rad(-2:n_bin_xi+2) 
    integer(kind=8) :: count(-2:n_bin_xi+3)

    integer(kind=8) :: i, j, index, indx, m1, m2, k, m
    real(kind=4) :: dx,dy,dz,distlog, r2, avgcount,didx,r3min,rho
    real(kind=4) :: x1,y1,z1,r1log,r2log,r3log,drlog
    real(kind=4) :: rmaxlog,rminlog   

    real(kind=4) :: dmass,x,y,z,vel,lum,absorp
    integer(kind=8) ::count1(-2:n_bin_xi+3)
    integer(kind=8),allocatable :: subcount(:,:)
    real(kind=4) :: intg(-2:n_bin_xi+2)

    if(n_point .gt. maxpoint_serial_xi) call abort
    allocate(subcount(-2:n_bin_xi+3,0:omp_workthreads-1))

    rminlog=log10(min_radius_xi)
    rmaxlog=log10(max_radius_xi)

    drlog=(rmaxlog-rminlog)/n_bin_xi

    r1log=rminlog-drlog
    r2log=rminlog-2*drlog
    r3log=rminlog-3*drlog

    do i=-2,n_bin_xi+2
       rad(i)=10.**(rminlog+i*drlog)
    end do
    subcount(-2:n_bin_xi+3,0:omp_workthreads-1) = 0

    !$omp parallel private(x1,y1,z1,m2,dx,dy,dz,r2,distlog,indx,index) 
    !$omp do
    do m1=1,n_point
       !print*,pos_in(1:3,m1)
       x1=pos_in(1,m1)
       y1=pos_in(2,m1)
       z1=pos_in(3,m1)
       !if(x1 <= 0. .or. x1 > 1.) print*, 'invalid',pos_in(1:3,m1)
       !if(y1 <= 0. .or. y1 > 1.) print*, 'invalid',pos_in(1:3,m1)
       !if(z1 <= 0. .or. z1 > 1.) print*, 'invalid',pos_in(1:3,m1)
       do m2=1,m1-1
          dx=abs(x1-pos_in(1,m2))
          if(dx.gt.0.5) dx=1.-dx
          dy=abs(y1-pos_in(2,m2))
          if(dy.gt.0.5) dy=1.-dy
          dz=abs(z1-pos_in(3,m2))
          if(dz.gt.0.5) dz=1.-dz
          r2=dx**2+dy**2+dz**2
          distlog=0.5*log10(r2)
          indx=int((distlog-r3log)/drlog)
          index=min(max(indx,0)-2,n_bin_xi+3)
          !if(index>nb+3 .or. index < -2)  print*,index

          subcount(index,omp_get_thread_num())=2+subcount(index,omp_get_thread_num())
          !          count(index)=count(index)+mass(m2)
          !    
       enddo
    end do
    !$omp end do
    !$omp end parallel

    ! Compute cumulative counts.
    count(-2:n_bin_xi+3) = 0
    do i=0,omp_workthreads-1
       count(-2:n_bin_xi+3) = count(-2:n_bin_xi+3)+subcount(-2:n_bin_xi+3,i)
    end do
    deallocate(subcount)
    !print*,count(0:n_bin_xi)

    do i=0,n_bin_xi+2
       count(i)=count(i)+count(i-1)
    end do


    !rho=float(n_point)
    r3min=rad(-2)**3

    intg(-2)=0.
    do i=-1,n_bin_xi+2
       avgcount=float(count(i))/n_point
       !intg(i)=avgcount/(4.*pi*rho)-(rad(i)**3-r3min)/3.
       intg(i)=avgcount/(4.*pi*rho)
    end do

    do i=0,n_bin_xi
       !didx=(intg(i-2)-intg(i+2)+8.*(intg(i+1)-intg(i-1)))/(12.*drlog)
       !didx=(intg(i-2)-intg(i+2)+2.*(intg(i+1)-intg(i-1)))/(8.*drlog)
       didx = (intg(i+1)-intg(i-1))/(2.*drlog)
       !xi_out(i)=didx/(rad(i)**3*log(1.0d01)) 
       xi_out(i)=didx/(rad(i)**3*log(1.0d01)) - 1.
       !       print*,'check corr.',didx,rad(i),corr(i),drlog,intg
       !       pause
    end do

    !print*, xi_out(0:n_bin_xi)
    return
  end subroutine correlation_function_unity

end module utilities_serial
