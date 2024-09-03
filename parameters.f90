module parameters
 use constants 
 implicit none

 integer :: Natoms       ! Number of atoms 
 integer :: Nx, Ny       ! Number of atoms in each direction
 real(dp) :: Lx, Ly      ! Box sizes  [nm]
 real(dp) :: Area
 real(dp) :: Rc          ! Cutoff radius [nm]
 real(dp) :: Temp        ! Temperature [K]
 real(dp) :: tinit       ! initialization time [fs]
 real(dp) :: tsim        ! simulation time [fs]
 real(dp) :: dt          ! time step [fs] 
 integer :: Nsteps       ! number of time steps

 real(dp) :: eps         ! LJ Energy [eV]
 real(dp) :: sigma       ! LJ sigma [nm]

 real(dp) :: Mass        ! in AMU (1822.886 me)
 real(dp) :: dr          ! step in sampling g(r)
 real(dp) :: v_drift     ! velocity of boundary particles
 
 real(dp) :: tgram
 integer :: Ngram        ! Gram-S steps
 real(dp):: D0=1.0_dp    ! initial lenght tangent vectors
 real(dp):: Fe=0.0_dp    ! External field

 real(dp):: Q=5.0_dp     ! Nose-Hoover mass
 
 logical :: scaling=.false.     ! velocity rescaling
 logical :: nose_hoover=.false. ! nose-hoover thermostat
 logical :: print_xyz=.false.   ! if xyz should be printed
 integer :: print_interval=1    ! xyz print interval


 character(50) :: file_name


! Stuff for Lyapunov
 real(dp) :: LJTU = 1.0_dp
 logical :: Dtemp
 logical :: shnn
 logical :: do_lyapunov

 type TLyapunovEnum
   integer :: volumes = 1
   integer :: lengths = 2
   integer :: QR = 3
 end type TLyapunovEnum  

  type(TLyapunovEnum), parameter, public :: LyapunovAlgorithm = TLyapunovEnum()
 
 integer :: algorithm = LyapunovAlgorithm%QR

 real(dp) :: Tempf
 real(dp) :: tsh
 integer(dp) :: partsint



 real(dp), dimension(:,:), allocatable :: x
 real(dp), dimension(:,:), allocatable :: v
 real(dp), dimension(:,:,:), allocatable :: eta ! Nose Hoover
 real(dp), dimension(:,:,:), allocatable :: dx
 real(dp), dimension(:,:,:), allocatable :: dv
 
 contains

 subroutine create_xv() !agg
   integer :: err
   allocate(x(2,Natoms),stat=err)
   allocate(v(2,Natoms),stat=err)
   if (err /= 0) STOP 'ALLOCATION ERROR x or v or eta'
 end subroutine create_xv

 subroutine create_eta() !agg
   integer :: err
   allocate(eta(2,Natoms,5), stat=err)
   if (err /= 0) STOP 'ALLOCATION ERROR x or v or eta'
 end subroutine create_eta

 subroutine destroy_xv()
   if (allocated(x)) deallocate(x)
   if (allocated(v)) deallocate(v)
 end subroutine destroy_xv

 subroutine destroy_eta()
   if (allocated(eta)) deallocate(eta)
 end subroutine destroy_eta

 subroutine create_xvly() !agg
   integer :: err
   allocate(dx(2,Natoms,4*Natoms),stat=err)
   allocate(dv(2,Natoms,4*Natoms),stat=err)
   if (err /= 0) STOP 'ALLOCATION ERROR xl or vl or etal'
 end subroutine create_xvly

 subroutine destroy_xvly()
   if (allocated(dx)) deallocate(dx)
   if (allocated(dv)) deallocate(dv)
 end subroutine destroy_xvly

 subroutine transform_units()
   
     ! transform particle mass from AMU to 
     ! something such that
     ! a = F/M  ([F]=eV/nm [a]=nm/fs^2)
     ! [M2F] = eV * (fs/nm)^2
   
     Mass = Mass * M2F
     Area = Lx*Ly
 end subroutine transform_units


end module parameters
