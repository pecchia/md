module simulations

  use constants 
  use parameters
  use list
  use boxes
  use forces
  use dynamics
  use clock
  use lyap_n
  use vector
  implicit none
  private

  public :: transform_units
  public :: init_positions_fcc
  public :: init_velocities
  public :: init_seed
  public :: nve_sim
  public :: write_coords
  public :: compute_g
  

  contains

  ! ---------------------------------------------------------
  ! Initialize positions of particles as in fcc lattice
  subroutine init_positions_fcc()
    integer :: i,err
    real(dp) :: rndr, ax, ay, az, aa(3)
    real(dp) :: bb(3,4)
    integer :: ii,jj,kk
    integer :: l,k,m

    if (4*Nx*Ny*Nz .ne. Natoms) then
       STOP 'Error: Nx*Ny*Nz != Natoms'
    endif

    ax = Lx/Nx; ay = Ly/Ny; az = Lz/Nz
    bb(:,1) = (/        ax/4.0_dp,        ay/4.0_dp, az/4.0_dp /)
    bb(:,2) = (/ 3.0_dp*ax/4.0_dp, 3.0_dp*ay/4.0_dp, az/4.0_dp /)
    bb(:,3) = (/ ax/4.0_dp, 3.0_dp*ay/4.0_dp, 3.0_dp*az/4.0_dp /)
    bb(:,4) = (/ 3.0_dp*ax/4.0_dp, ay/4.0_dp, 3.0_dp*az/4.0_dp /)

    i = 1

    do m = 1, Nz
      do k = 1, Ny
         do l= 1, Nx

            aa(1) = (l-1)*Lx/Nx
            aa(2) = (k-1)*Ly/Ny
            aa(3) = (m-1)*Lz/Nz
            ! Initialize positions in the box
            x(:,i) = aa(:) + bb(:,1)
            call boxind(x(:,i),ii,jj,kk)
            !write(*,'(3f16.6,3i5)') x(:,i),ii,jj,kk
            call add(boxlists(ii,jj,kk),i)
            i = i + 1

            x(:,i) = aa(:) + bb(:,2)
            call boxind(x(:,i),ii,jj,kk)
            !write(*,'(3f16.6,3i5)') x(:,i),ii,jj,kk
            call add(boxlists(ii,jj,kk),i)
            i = i + 1

            x(:,i) = aa(:) + bb(:,3)
            call boxind(x(:,i),ii,jj,kk)
            !write(*,'(3f16.6,3i5)') x(:,i),ii,jj,kk
            call add(boxlists(ii,jj,kk),i)
            i = i + 1

            x(:,i) = aa(:) + bb(:,4)
            call boxind(x(:,i),ii,jj,kk)
            !write(*,'(3f16.6,3i5)') x(:,i),ii,jj,kk
            call add(boxlists(ii,jj,kk),i)
            i = i + 1

         enddo
      enddo
    enddo

  end subroutine init_positions_fcc

  ! ---------------------------------------------------------
  ! according to the MB distribution
  subroutine init_velocities()

    integer :: i

    do i = 1, Natoms

      !print*,'set velocities'
      ! Initialize random MB velocities.
      v(1,i) = maxwell_boltzmann()
      v(2,i) = maxwell_boltzmann()
      v(3,i) = maxwell_boltzmann()
    end do

  end subroutine init_velocities
    

  ! ---------------------------------------------------------
  subroutine init_lyapunov()    

   integer :: i, j, k

   do i=1,6*Natoms
     do j=1,Natoms
        do k=1,3
           if (i.le.3*Natoms) then
              if (i+j-k==4*j-3) then
                 dx(k,j,i)=0.0_dp
              else
                 dx(k,j,i)=0.0_dp
              end if
           else
              dx(k,j,i)=0.0_dp
           end if
           
        end do

     end do

  end do
          
   do i=1,6*Natoms
     do j=1,Natoms
        do k=1,3
           if (i.gt.3*Natoms) then
              if (i-3*Natoms+j-k==4*j-3) then
                 dv(k,j,i)=0.0_dp
              else
                 dv(k,j,i)=0.0_dp
              end if
          

           else
              dv(k,j,i)=0.0_dp
           end if
           
        end do

     end do

  end do

  end subroutine init_lyapunov
  
  ! ---------------------------------------------------------
  ! In ogni direzione P = (m/2PikT)^1/2 exp(-m v^2 / 2kT)
  ! P0=(m/2PikT)^1/2
  ! P/P0 = exp[-m v^2/ 2kT]
  ! => v^2 = +/- sqrt( - m/2kT * ln(P/P0) )
  function maxwell_boltzmann() result(vv)
    real(dp) :: vv
    real(dp) :: rndr, b, eps

    ! Define machine precision
    ! This is used to change the interval [0,1) into (0,1]
    eps = mach()

    call random_number(rndr)
    b= cos(2.0_dp*Pi*rndr)
    do
       call random_number(rndr)

       vv = b*sqrt( - 2.0_dp * kb * Temp/Mass * log(rndr+eps) )
       exit
    enddo

  end function maxwell_boltzmann

  ! ---------------------------------------------------------
  subroutine init_seed(seed_in)
    integer, intent(in), optional :: seed_in

    integer, dimension(:),allocatable :: seed
    integer :: j

    call random_seed(size=j)
    allocate(seed(j))

    if (present(seed_in)) then
      seed=seed_in
    else
      call system_clock(j)
      seed=j
    end if
    call random_seed(put=seed)
    deallocate(seed)

  end subroutine init_seed


  ! ---------------------------------------------------------
  ! Perform an NVE simulation
  subroutine nve_sim()
    real(dp), dimension(:,:), allocatable :: xf,vf
    real(dp), dimension(:,:,:), allocatable :: etaf
    real(dp), dimension(:,:), allocatable :: Vxlf,Vylf,Vzlf
    real(dp), dimension(:,:), allocatable :: Xlf,Ylf,Zlf


    real(dp), dimension(:,:,:), allocatable :: Vxlf2,Vylf2,Vzlf2
    real(dp), dimension(:,:,:), allocatable :: Xlf2,Ylf2,Zlf2

    real(dp), dimension(:,:,:),allocatable :: dxf,dvf
     real(dp), dimension(6*Natoms,6*Natoms) ::d
    real(dp), dimension(:,:),allocatable :: uv,vu
    real(dp), dimension(:,:,:),allocatable :: inf_v

    
    real(dp), dimension(:), allocatable :: Lyapunov,L_plus,L_minus


    !real(dp),dimension(:,:),allocatable :: Nor1,Nor2,Nor3


    real(dp), dimension(:),allocatable :: kine
    
    
    integer :: nstep1,nstep2, err
    integer :: n,ng, iter , pq,i,j
    character(3) :: ind

  
    real(dp) :: U      ! Potential energy
    real(dp) :: K      ! Kinetic energy
   
    real(dp) :: K0   !kin
    real(dp) :: virial ! virial = - Sum_i Sum_j (r_ij * F_ij)
    real(dp) :: P      ! Pressure (from virial)
    real(dp) :: lambda ! scaling parameter 
    real(dp) :: Kav,Uav,Pav ! Averaged quantities
    real(dp) :: R(3),R0(3)  ! For diffusivity 
    real(dp), allocatable :: dR(:)  ! For diffusivity 
    real(dp) :: R2, R02, Diff, Rcm(3)
    ! write(*,*)
    
    
   
    pq=0
    Pav = 0.0_dp
    Uav = 0.0_dp
    Kav = 0.0_dp
   
    Diff= 0.0_dp 
    Rcm = 0.0_dp
    K0=(3.d0*Natoms*kb*Temp/2.d0)

    nstep1=nint(tinit/dt)
    nstep2=nint(tsim/dt)

    ng=nstep2/ngram

    allocate(dR(3*Natoms),stat=err)
    allocate(etaf(3,Natoms,5),stat=err)
    
    allocate(xf(3,Natoms),stat=err)
    allocate(vf(3,Natoms),stat=err)

    allocate(dxf(3,Natoms,6*Natoms),stat=err)
    allocate(dvf(3,Natoms,6*Natoms),stat=err)

    allocate(uv(6*Natoms,6*Natoms),stat=err)
    allocate(vu(6*Natoms,6*Natoms),stat=err)
    allocate(inf_v(6*Natoms,6*Natoms,ng),stat=err)

    allocate(Lyapunov(6*Natoms),stat=err)
    allocate(L_plus(6*Natoms),stat=err)
    allocate(L_minus(6*Natoms),stat=err)

    allocate(kine(nstep1),stat=err)
    
    
    if (err /= 0) STOP 'ALLOCATION ERROR'
    

    call init_lj(Natoms,eps,sigma,Rc)
    
    !!$    call init_ly(Natoms,eps,sigma,Rc)
    
    write(*,*) 'Target T=',Temp,'K=',(3.d0*Natoms*kb*Temp/2.d0)

    call set_clock()
    !write(ind,'(i3.3)') n
    !open(101,file='coord'//ind//'.xyz')

    if (print_xyz) then
      open(101,file='data/coords.xyz')
      write(101,'(i0)') Natoms
      write(101,*) 'Frame',0
      call write_xyz(101)
    end if
    open(102,file='data/R2.dat')
    open(103,file='data/kin.dat')
    open(113,file='data/e_tot.dat')
    open(123,file='data/U.dat')
    open(133,file='data/P.dat')
    
    open(104,file='data/lyap_pl_200_2.dat')
    open(106,file='data/lyap_mi_200_2.dat')
    open(105,file='data/xv_x300.dat')
    open(115,file='data/xv_y300.dat')
    open(125,file='data/xv_z300.dat')

    ! init time
    write(*,*) 'Warm up phase:'
    do n=1,nstep1

       if (print_xyz .and. mod(n,xyz_interval) == 0) then
          write(101,'(i0)') Natoms
          write(101,*) 'Frame',n
          call write_xyz(101)
       endif

       !call verlet(x,v,xf,vf,U,virial,dt,lj)
       call verlet_nh2_5(x,v,U,virial,dt,lj,K,K0,eta,etaf,xf,vf,dx,dv,dxf,dvf)

       call update_boxes(x,xf)

       v = vf
       x = xf

       !K = kinetic()
       P = (Natoms*kb*Temp + virial/3.d0)/Vol

       if (mod(n,xyz_interval) == 0) then
         write(*,'(i6,a,i6,3x,4(a3,ES14.6,2x))') n,'/',nstep1,'Ek=',K,'U=',U,'E=',K+U,'P=',P
       end if

       !if (scaling) then
       !   lambda = sqrt((3.d0*Natoms*kb*Temp/2.d0)/K)
       !   v = vf * lambda
       !endif

    end do

    ! simulation time
    R0 = x(:,Natoms/2)
    R02 = dot_product(R0,R0)
    dR = 0.0_dp   
    
    write(*,*) 'Simulation phase:'
    ! init time
    d=newshape2(dx,dv)
    
    do n=1,nstep2
     
       if (print_xyz .and. mod(n,xyz_interval) == 0) then
          write(101,'(i0)') Natoms
          write(101,*) 'Frame',n
          call write_xyz(101)
       endif
       
       
       !call verlet(x,v,xf,vf,U,virial,dt,lj,K)
    
       !call verlet_nh1(x,v,U,virial,dt,lj,K,K0,eta,etaf,xf,vf,dx,dv,dxf,dvf)
       !call verlet_nh2(x,v,U,virial,dt,lj,K,K0,eta,etaf,xf,vf,dx,dv,dxf,dvf)
       !call verlet_nh1_5(x,v,U,virial,dt,lj,K,K0,eta,etaf,xf,vf,dx,dv,dxf,dvf)
       call verlet_nh2_5(x,v,U,virial,dt,lj,K,K0,eta,etaf,xf,vf,dx,dv,dxf,dvf)

       d=newshape2(dxf,dvf)

       ! GRAMS-SHMIDT ORTHOGONALIZATION 
       if (mod(n,ngram)==0) then
          pq=pq+1
          
          call grams( dxf,dvf,vu,uv)

          do i=1,6*Natoms
             do j=1,6*Natoms
                if(j.ne.i) then
                   if (dot_product(uv(:,i),uv(:,j)).ge. 0.8) then
                      write(*,*) 'ortornormal error'
                      write(*,*) dot_product(uv(:,i),uv(:,j))
                      write(*,*) sqrt(dot_product(uv(:,j),uv(:,j)))
                      write(*,*) sqrt( dot_product(uv(:,i),uv(:,i)))
                      stop
                   end if

                end if
             end do
          end do
          inf_v(:,:,pq)=vu(:,:)
       end if
       
       call update_boxes(x,xf)
 
       x = xf     
       v = vf
 
       !K = kinetic()  
       P = (Natoms*kb*Temp + virial/3.d0)/Vol

       write(*,'(i6,3x,3(a3,ES14.6,2x))') n,'Ek=',K,'U=',U,'P=',P
       write(103,*) n*dt, K
       write(113,*) n*dt, U+K
       write(123,*) n*dt, U
       write(133,*) n*dt, P
       kine(n)=K
      
       !if (scaling) then 
       !   lambda = sqrt((3.d0*Natoms*kb*Temp/2.d0)/K)
       !   v = vf * lambda
       !endif

       eta=etaf

       dx=dxf
       dv=dvf
   
       do i=1,Natoms
          write(105,*) x(1,i),v(1,i)
          write(115,*) x(2,i),v(2,i)
          write(125,*) x(3,i),v(3,i)
       end do

    end do

    call lyap_numbers(inf_v,Lyapunov,ng,nstep1)

    do i=1,6*Natoms
      if (Lyapunov(i) .gt.0) then
        L_plus(i)=Lyapunov(i)
      else if(Lyapunov(i) .lt.0) then
        L_minus(i-96)=Lyapunov(i)
      end if   
    end do

    do i=1,(6*Natoms)/2
       write(104,*) i,L_plus (i)
       write(106,*) (6*Natoms)/2- i, L_minus (i)
       write(*,*) Lyapunov(i)
    end do


    close(101)
    close(102)
    close(103)
    close(113)
    close(123)
    close(133)
    close(104)
    close(105)
    close(115)
    close(125)
    close(106)


    deallocate(Lyapunov)

    write(*,*)
    write(*,*) 'K0=', K0
    write(*,*) 'scarto valore K0', abs(K0-sum(kine)/nstep1)
    write(*,*) 'Average K' , sum(kine)/nstep1
    write(*,*) 'scarto=',sqrt( (sum((kine-(sum(kine)/nstep1))**2))/nstep1)
    write(*,*) 'pot=', U
    write(*,'(a16,3x)',advance='NO') 'Simulation time:'
    call write_clock()

!!$    write(*,*)
!!$    write(*,*) 'Averaged quantities:'
!!$    write(*,'(9x,3(a3,ES14.6,2x))') 'Ek=',Kav,'U=',Uav,'P=',Pav 
!!$    write(*,*) 'Rcm=',Rcm
!!$    write(*,*) 'Diffusivity=',Diff/dt,'nm^2/fs'
!!$    !write(*,*) ' K0=' K0
!!$   
    deallocate(xf,vf,etaf)
    deallocate(dxf,dvf)
 

    open(102,file='av.dat')
    write(102,'(9x,3(a3,ES14.6,2x))') 'Ek=',Kav,'U=',Uav,'P=',Pav 
    write(102,*) 'Rcm=',Rcm
    write(102,*) 'Diffusivity=',Diff/dt,'nm^2/fs'
    
    close(102)

  end subroutine nve_sim
 
  function kinetic() result(Ek)
     real(dp) :: Ek
     integer :: n
     
     Ek =0.0_dp
     do n = 1, Natoms
        Ek = Ek + dot_product(v(:,n),v(:,n))   
     end do
     Ek = Ek*Mass*0.5_dp 

  end function kinetic

  subroutine write_xyz(id)
    integer :: id

    integer :: n, ii, jj, kk
 
    do n = 1, Natoms      
       !call boxind(x(:,n),ii,jj,kk)
       write(id,'(a,3(f12.6),3(ES20.8))') 'He  ', x(:,n)*10.0_dp  !,v(:,n)
    enddo

  end subroutine write_xyz

  subroutine write_coords(id)
    integer :: id

    integer :: n
 
    write(id,'(1X,L2,I7,3E23.15)') .true.,Natoms,Lx/sigma,Ly/sigma,Lz/sigma
    do n = 1, Natoms      
       write(id,'(1X,3E23.15)') x(:,n)/sigma
    enddo
    do n = 1, Natoms      
       write(id,'(1X,3E23.15)') v(:,n)*sqrt(Mass/eps)
    enddo
    do n = 1, Natoms      
       write(id,'(1X,3E23.15)') v(:,n)/dt*(sigma*Mass/eps)
    enddo

  end subroutine write_coords


  
  subroutine compute_g()
    real(dp) :: rij(3), g(3), r, a, b 
    integer :: ii, jj, kk, ci, cj, ck, u,v,w
    integer :: m, l
    type(TNode), pointer :: it  
     
    real(dp),dimension(:), allocatable :: gg
    integer :: Nk, basket, err
    
    Nk = aint(Rc/dr)
    allocate(gg(Nk),stat=err)
    if(err.ne.0) STOP 'ALLOCATION ERROR gg'

    gg = 0.0_dp

    do m = 1, Natoms

       ! cerca la scatola ci,cj,ck di i
       call boxind(x(:,m),ci,cj,ck)          
       
       do w=-1,1
         do v=-1,1 
           do u=-1,1
              ii = ci + u
              jj = cj + v     
              kk = ck + w  
                    
              ! g e' vettore supercella
              call folding(ii,jj,kk,g)

              ! Iterates over atoms in box (ii,jj,kk)
              it => boxlists(ii,jj,kk)%start
         
              do while (associated(it))     
             
                  l = it%val

                  if (l .eq. m) then 
                      it => it%next
                      cycle
                  endif

                  ! segno corretto rij = rj - ri
                  rij(:) = x(:,l)-x(:,m)+g(:)
           
                  r = sqrt(dot_product(rij,rij))
                
                  if (r<Rc) then
                    basket = aint(r/dr) + 1
                    gg(basket) = gg(basket) + 1
                  endif

                  it => it%next

              end do

            enddo
          enddo
        enddo

     end do
      

     open(101,file='g.dat')
     do m = 1, Nk
        r = m*dr
        write(101,*) r, gg(m)*Vol/(4.d0*Natoms*Natoms*Pi*r*r*dr)  
     enddo
 
     close(101)

   end subroutine compute_g




end module simulations

