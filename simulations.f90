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
  public :: init_positions_rnd
  public :: init_velocities
  public :: init_velocities_rnd
  public :: init_positions_file
  public :: init_velocities_couette
  public :: init_seed
  public :: nve_sim
  public :: write_coords
  public :: compute_g
  public :: init_lyapunov
  

  contains

  subroutine init_positions_rnd()
    integer :: i, ii, jj
    call random_number(x(1,:))
    x(1,:) = x(1,:)*Lx
    call random_number(x(2,:))
    x(2,:) = x(2,:)*Ly

    do i = 1, Natoms
       call boxind(x(:,i),ii,jj)
       call add(boxlists(ii,jj),i)
    end do 
  end subroutine init_positions_rnd

  subroutine init_velocities_rnd()
 
    call random_number(v(1,:))
    call random_number(v(2,:)) 

  end subroutine init_velocities_rnd

  ! ---------------------------------------------------------
  ! Initialize positions of particles as in fcc lattice
  subroutine init_positions_fcc()
    integer :: i,err
    real(dp) :: rndr, ax, ay, aa(2)
    real(dp) :: bb(2,2)
    integer :: ii,jj
    integer :: l,k

    
    ax = Lx/Nx; ay = Ly/Ny

    bb(:,1) = (/        ax/4.0_dp,        ay/4.0_dp       /)
    bb(:,2) = (/ 3.0_dp*ax/4.0_dp,  3.0_dp*ay/4.0_dp        /)

    i = 1

      do k = 1, Ny
         do l= 1, Nx

            aa(1) = (l-1)*Lx/Nx
            aa(2) = (k-1)*Ly/Ny
            ! Initialize positions in the box
            x(:,i) = aa(:) + bb(:,1)
            call boxind(x(:,i),ii,jj)
            !write(*,'(3f16.6,3i5)') x(:,i),ii,jj
            call add(boxlists(ii,jj),i)
            i = i + 1

            x(:,i) = aa(:) + bb(:,2)
            call boxind(x(:,i),ii,jj)
            !write(*,'(3f16.6,3i5)') x(:,i),ii,jj
            call add(boxlists(ii,jj),i)
            i = i + 1

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
    end do

  end subroutine init_velocities
    

  ! ---------------------------------------------------------
  ! Init velocities for the couette BC
subroutine init_velocities_couette()

    integer :: i
    real(dp) :: upperlayers, lowerlayers

    upperlayers = maxval(x(2,:)) - Ly/Ny !first two layers of the fluid
    lowerlayers = minval(x(2,:)) + Ly/Ny

    do i = 1, Natoms

       if(x(2,i) .ge. upperlayers) then
          v(1,i) = v_drift
          v(2,i) = 0.d0

        elseif(x(2,i) .le. lowerlayers) then
          v(1,i) = 0.d0
          v(2,i) = 0.d0          

       else
          v(1,i) = maxwell_boltzmann()
          v(2,i) = maxwell_boltzmann()
       end if
    end do
    print*,'init velocities couette:'
    print*,'vmax=',maxval(v(1,:)),maxval(v(2,:))
    print*,'vmin=',minval(v(1,:)),minval(v(2,:))
    print*,'done'

  end subroutine init_velocities_couette

   subroutine init_positions_file
    real(dp) :: dummy
    character(5) :: string
    integer :: i, ii, jj
    
    print*,"Reading file "//trim(file_name)

    open(100, file=trim(file_name))
    read(100, *)
    read(100, *)
    do i=1,Natoms

      read(100, * ) string, x(1,i), x(2,i), dummy, v(1,i), v(2,i), dummy
      call boxind(x(:,i),ii,jj)
      call add(boxlists(ii,jj),i)

    end do
    
    
   end subroutine init_positions_file

  
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

       vv = b*sqrt( - kb * Temp/Mass * log(rndr+eps) )
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

    real(dp), dimension(:,:,:),allocatable :: dxf,dvf
    real(dp), dimension(:,:),allocatable :: uu, vv
    real(dp), dimension(:,:,:),allocatable :: inf_v
    
    real(dp), dimension(:), allocatable :: Lyapunov,L_plus,L_minus
    real(dp), dimension(:),allocatable :: kine
    
    integer :: nstep1, nstep2, nstepgram, ng, err
    integer :: n, iter, pq, i, j, kminus, kplus
    character(3) :: ind

    real(dp) :: U           ! Potential energy
    real(dp) :: K           ! Kinetic energy
    real(dp) :: K0          ! kin ?
    real(dp) :: virial      ! virial = - Sum_i Sum_j (r_ij * F_ij)
    real(dp) :: P           ! Pressure (from virial)
    real(dp) :: lambda      ! scaling parameter
    real(dp) :: Kav,Uav,Pav ! Averaged quantities
    real(dp) :: R(2),R0(2)  ! For diffusivity 
    real(dp), allocatable :: dR(:)  ! For diffusivity 
    real(dp) :: R2, R02, Diff, Rcm(2)
  
    Pav = 0.0_dp
    Uav = 0.0_dp
    Kav = 0.0_dp
   
    Diff= 0.0_dp 
    Rcm = 0.0_dp
    !Lennard-Jones time units:
    !LJTU = sqrt(Mass/eps)*sigma
    !print*,'LJ TIME UNITS:',LJTU

    nstep1=nint(tinit/dt)
    nstep2=nint(tsim/dt)
    nstepgram = nint(tgram/dt)
    if (do_lyapunov) then
      ng = tsim/tgram
      pq = 0
      allocate(dxf(2,Natoms,4*Natoms),stat=err)
      allocate(dvf(2,Natoms,4*Natoms),stat=err)
 
      allocate(vv(4*Natoms,4*Natoms),stat=err)
      allocate(uu(4*Natoms,4*Natoms),stat=err)
 
      allocate(Lyapunov(4*Natoms),stat=err)
      Lyapunov = 0.0_dp
      allocate(L_plus(4*Natoms),stat=err)
      allocate(L_minus(4*Natoms),stat=err)
      allocate(inf_v(4*Natoms,4*Natoms,ng),stat=err)
    end if

    allocate(dR(2*Natoms),stat=err)
    allocate(xf(2,Natoms),stat=err)
    allocate(vf(2,Natoms),stat=err)
    allocate(kine(nstep2),stat=err)
    allocate(etaf, mold=eta)
    
    if (err /= 0) STOP 'ALLOCATION ERROR'
    
    call init_lj(Natoms,eps,sigma,Rc)
    
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
  
    
    !open(104,file='data/lyap_pl_200_2.dat')
    !open(106,file='data/lyap_mi_200_2.dat')
    !open(105,file='data/xv_x300.dat')
    !open(115,file='data/xv_y300.dat')
    !open(125,file='data/xv_z300.dat')

    ! init time
    write(*,*) 'Warm up phase:',nstep1,'steps'
    ! Target mean Kinetic energy:
    K0=(Natoms*kb*Temp)
    write(*,*) 'Target T=',Temp,'K=',K0


    do n=1,nstep1

       if (print_xyz .and. mod(n,print_interval) == 0) then
          write(101,'(i0)') Natoms
          write(101,*) 'Frame',n
          call write_xyz(101)
       endif
      !
      call verlet(x,v,xf,vf,U,virial,dt,lj,K)
      call sym4(x,v,xf,vf,U,virial,dt,lj,K)
      call update_boxes(x,xf)


       !P = (Natoms*kb*Temp + virial/2.d0)/Area

       if (mod(n,print_interval) == 0) then
         write(*,'(i8,a,i8,3x,4(a3,ES14.6,2x))') n,'/',nstep1,'Ek=',K,'U=',U,'E=',K+U
       end if

       if (scaling) then
          lambda = sqrt((Natoms*kb*Temp)/K)
          v = vf * lambda
       else 
          v = vf
       endif
       x = xf
       eta = etaf
       
    end do

    ! simulation time
    R0 = x(:,Natoms/2)
    R02 = dot_product(R0,R0)
    dR = 0.0_dp   
   
    write(*,*) '**********************************************************************' 
    write(*,*) 'Simulation phase:',nstep2,'steps'
    ! init time

    do n=1, nstep2
     
       if (print_xyz .and. mod(n,print_interval) == 0) then
          write(101,'(i0)') Natoms
          write(101,*) 'Frame',n
          call write_xyz(101)
       endif
       
       if (do_lyapunov) then     
          call verlet_ly(x,v,xf,vf,U,virial,dt,lj_ly,K,dx,dv,dxf,dvf)
          !call sym4_ly(x,v,xf,vf,U,virial,dt,lj_ly,K,dx,dv,dxf,dvf) 
       else   
          call verlet(x,v,xf,vf,U,virial,dt,lj,K)
          !call sym4(x,v,xf,vf,U,virial,dt,lj,K) 
       end if   
       
       if (do_lyapunov) then
          ! Ly = sum_i log(|dx|)/tfin
          ! => |dx| = exp(Ly*tfin)
          ! RE-ORTHOGONALIZATION procedure
          if (mod(n,nstepgram)==0) then
             call message_clock("QR orthogonalization") 
             call qr(dxf, dvf, vv, lyapunov)
             call write_clock()
             !print*,'lyap_max=',maxval(abs(lyapunov))   

             !call grams(dxf, dvf, vv, uu)
 
             if (algorithm /= LyapunovAlgorithm%QR) then 
                pq=pq+1
                inf_v(:,:,pq) = vv(:,:)
             end if  
          end if
       end if

       call update_boxes(x,xf)
 
       !P = (Natoms*kb*Temp + virial/2.d0)/Area
      

       if (mod(n,print_interval) == 0) then
         write(*,'(i8,a,i8,3x,4(a3,ES14.6,2x))') n,'/',nstep2,'Ek=',K,'U=',U,'E=',K+U
       end if 
       !write(103,*) n*dt, K
       !write(113,*) n*dt, U+K
       !write(123,*) n*dt, U
       !write(133,*) n*dt, P
       kine(n)=K
       Uav = Uav + U/nstep2

       if (scaling) then 
          lambda = sqrt((Natoms*kb*Temp)/K)
          v = vf * lambda
       else 
          v = vf
       endif

       x = xf     
       eta=etaf
       if (do_lyapunov) then
         dx=dxf
         dv=dvf
       end if

      dR = dR + reshape(v ,(/ 2*Natoms /))
      dR = dR + reshape(vf,(/ 2*Natoms /))
      R2 = dot_product(dR,dR)*dt*dt/Natoms/4.0_dp
      write(102,*) n*dt, R2

      Pav = Pav + P/nstep2
      Kav = Kav + K/nstep2
      Uav = Uav + U/nstep2
      R = x(:,Natoms/2) - R0
      Rcm = Rcm + R/nstep2
      R2 = dot_product(R,R)
      Diff = Diff + R2/nstep2
   
    end do

    if (do_lyapunov) then
      print*,'COMPUTATION OF LYAPUNOV SPECTRUM'
      select case(algorithm) 
      case(LyapunovAlgorithm%volumes)     
        print*,'Volumes algorithm is used'     
        call lyap_numbers(inf_v, Lyapunov, ng)
      case(LyapunovAlgorithm%lengths)
        print*,'Lengths algorithm is used'     
        call lyap_numbers2(inf_v, Lyapunov, ng)
      case(LyapunovAlgorithm%QR)
        print*,'QR algorithm is used'     
      end select      

      kplus=0; kminus=0;
      do i=1,4*Natoms
        if (Lyapunov(i).gt.0) then
          kplus=kplus+1    
          L_plus(kplus)=Lyapunov(i)
        else if(Lyapunov(i).lt.0) then
          kminus=kminus+1    
          L_minus(kminus)=Lyapunov(i)
        end if   
      end do
  
      do i=1,2*Natoms
         !write(104,*) i, L_plus(i)*LJTU
         write(*,'(a,i4,a,f10.5)') 'l(',i,')=',L_plus(i)*LJTU
         write(*,'(a,i4,a,f10.5)') 'l(',i,')=',L_minus(i)*LJTU
      end do

      ! Kaplan-Yorke dimension
      !lj1=0.d0
      !jj=0.d0
      !do i=1,6*Natoms
      !  if (lj1 .ge. 0.d0) then
      !     lj1=lj1+lyapunov(i)
      !     jj=0.d0+i
      !  end if
      !end do
     
      !if (nint(jj) .eq. (6*Natoms)) then
      !   Dfrac=(jj-6)
      !else
      !   Dfrac=(jj-6)+lj1/abs(lyapunov(int(jj+1)))
      !end if
    end if

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
   


    write(*,*)
    write(*,*) 'K0 =', K0
    write(*,*) '<K> =' , sum(kine)/nstep2
    write(*,*) '|<K> - K0| =', abs(K0-sum(kine)/nstep2)
    write(*,*) 'sqrt(<K^2>-<K>^2) =',sqrt( (sum((kine-(sum(kine)/nstep2))**2))/nstep2)
    write(*,*) '<U> =', Uav
    write(*,'(a16,3x)',advance='NO') 'Simulation time:'
    call write_clock()

!!$    write(*,*)
!!$    write(*,*) 'Averaged quantities:'
!!$    write(*,'(9x,3(a3,ES14.6,2x))') 'Ek=',Kav,'U=',Uav,'P=',Pav 
!!$    write(*,*) 'Rcm=',Rcm
!!$    write(*,*) 'Diffusivity=',Diff/dt,'nm^2/fs'
!!$    !write(*,*) ' K0=' K0
    if (do_lyapunov) then 
      deallocate(Lyapunov,L_plus,L_minus)
      deallocate(dxf,dvf,uu,vv)
    end if  
    deallocate(xf,vf,etaf)

    !open(102,file='av.dat')
    !write(102,'(9x,3(a3,ES14.6,2x))') 'Ek=',Kav,'U=',Uav,'P=',Pav 
    !write(102,*) 'Rcm=',Rcm
    !write(102,*) 'Diffusivity=',Diff/dt,'nm^2/fs'
    !close(102)

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

    integer :: n, ii, jj
 
    do n = 1, Natoms      
       !call boxind(x(:,n),ii,jj)
      write(id,'(a,4(f12.6),3(ES20.8))') 'He  ', x(:,n), 0.0_dp ,v(:,n), 0.0_dp
    enddo

  end subroutine write_xyz

  subroutine write_coords(id)
    integer :: id

    integer :: n
 
    write(id,'(1X,L2,I7,3E23.15)') .true.,Natoms,Lx/sigma,Ly/sigma
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

  ! ---------------------------------------------------------
  subroutine init_lyapunov()    

    integer :: i, j, k

    call create_xvly()

    ! dx(2*Natoms,4*Natoms) -> dx(l,m)  = delta(l,m)
    ! (1 0 0 0 . . .  )
    ! (0 1 0 0 . . .  )
    ! (0 0 1 0 .....  )
    !
    ! Mapping to dx(k,j,i) into dx(l,m):
    !   l=2(j-1)+k, m=i
    !   
    !  2j-2 + k == i ==>  3*j-2 == i + j - k

    dx = 0.0_dp
    dv = 0.0_dp

    do i = 1, 4*Natoms
      do j = 1, Natoms
         do k = 1, 2
            if (i.le.2*Natoms) then
               ! 4*j-3    
               if (i+j-k == 3*j-2) then
                  dx(k,j,i)=D0
               end if
            else
               if (i+j-k-2*Natoms == 3*j-2) then
                  dv(k,j,i)=D0
               end if
            end if
         end do
       end do
    end do
 
   end subroutine init_lyapunov
   
   ! ---------------------------------------------------------
  
  subroutine compute_g()
    real(dp) :: rij(2), g(2), r, a, b 
    integer :: ii, jj, ci, cj, u,v
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
       call boxind(x(:,m),ci,cj)          
    
         do v=-1,1 
           do u=-1,1
              ii = ci + u
              jj = cj + v
                    
              ! g e' vettore supercella
              call folding(ii,jj,g)

              ! Iterates over atoms in box (ii,jj)
              it => boxlists(ii,jj)%start
         
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


     end do
      

     open(101,file='data/g.dat')
     do m = 1, Nk
        r = m*dr
        write(101,*) r, gg(m)*Area/(2.d0*Natoms*Natoms*Pi*r*r*dr) 
     enddo
 
     close(101)

   end subroutine compute_g




end module simulations

