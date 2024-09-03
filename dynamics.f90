module dynamics
  use constants, only : dp
  use parameters, only : Natoms, Mass, Q
  use boxes, only : update_boxes
  implicit none
  private

  public :: verlet          ! + ext forces
  public :: verlet_ly       ! + ext forces
  !public :: verlet_nh15     ! + ext forces
  !public :: verlet_nh15_ly  ! + ext forces

  !public :: verlet_nh1_ly   
  !public :: verlet_nh2_ly
  !public :: verlet_nh25_ly
  
  public :: sym4
  public :: sym4_ly

  interface  
    subroutine Tforces(x,F,UU,virial)
      use precision, only : dp
      real(dp), dimension(:,:), intent(in) :: x
      real(dp), dimension(:,:), intent(out) :: F
      real(dp), intent(out) :: UU
      real(dp), intent(out) :: virial
    end subroutine TForces

    subroutine Tforces_ly(x,dx,F,lF,UU,virial)
      use precision, only : dp
      real(dp), dimension(:,:), intent(in) :: x
      real(dp), dimension(:,:,:), intent(in) :: dx
      real(dp), dimension(:,:), intent(out) :: F
      real(dp), dimension(:,:,:), intent(out) :: lF
      real(dp), intent(out) :: UU
      real(dp), intent(out) :: virial
    end subroutine Tforces_ly
  end interface

contains

  ! ---------------------------------------------------------------------------------
  ! BASIC VERLET ALGORITHM 
  !
  !Aggiungere Subroutine con un algoritmo Simplettico migliore
  !
  subroutine verlet(x,v,x1,v1,U,virial,dt, forces, K)
    real(dp), dimension(:,:), intent(in) :: x, v
    real(dp), dimension(:,:), intent(out) :: x1, v1
    real(dp), intent(out) :: U
    real(dp), intent(out) :: virial
    real(dp), intent(in) :: dt
    procedure(Tforces) :: forces
    real(dp), intent(out) :: K

    real(dp), dimension(:,:), allocatable :: F, F1
    integer :: err

    allocate(F(2,Natoms), stat=err)
    allocate(F1(2,Natoms), stat=err)
    if (err.ne.0) STOP 'ALLOCATION ERROR'


    call forces(x,F,U,virial)
    x1 = x + v * dt + 0.5_dp*dt*dt*F/Mass

    call forces(x1,F1,U,virial)
    v1 = v + 0.5_dp * (F+F1)/Mass * dt

    K=kinetic(v1)

    deallocate(F,F1)  

  end subroutine verlet
  
  subroutine verlet_ly(x,v,x1,v1,U,virial,dt,forces,K,dx,dv,dxf,dvf)
    real(dp), dimension(:,:), intent(in) :: x, v
    real(dp), dimension(:,:), intent(out) :: x1, v1
    real(dp), dimension(:,:,:), intent(in) :: dx,dv
    real(dp), dimension(:,:,:), intent(out) :: dxf,dvf  
    real(dp), intent(out) :: U
    real(dp), intent(out) :: virial
    real(dp), intent(in) :: dt
    procedure(Tforces_ly) :: forces
    real(dp), intent(out) :: K

    real(dp), dimension(:,:), allocatable :: F, F1
    real(dp), dimension(:,:,:), allocatable :: lF, lF1
    integer :: err

    allocate(F(3,Natoms), stat=err)
    allocate(F1(3,Natoms), stat=err)
    allocate(lF(3,Natoms,6*Natoms), stat=err)
    allocate(lF1(3,Natoms,6*Natoms), stat=err)
    
    if (err.ne.0) STOP 'ALLOCATION ERROR'


    call forces(x,dx,F,lF,U,virial)
    x1 = x + v * dt + 0.5_dp*dt*dt*F/Mass
    dxf = dx + dv*dt + 0.5_dp*dt*dt*lF/Mass

    call forces(x1,dxf,F1,lF1,U,virial) 
    v1 = v + 0.5_dp * (F+F1)/Mass * dt
    dvf= dv+ 0.5_dp*dt*(lF+lF1)/Mass
    
    K=kinetic(v1)

    deallocate(F,F1, lF, lF1)  

  end subroutine verlet_ly
  ! ---------------------------------------------------------------------------------
  ! Symplectic integration of 4th order
  subroutine sym4(x,v,x1,v1,U,virial,dt,forces,K)
    real(dp), dimension(:,:), intent(in) :: x, v
    real(dp), dimension(:,:), intent(out) :: x1, v1
    real(dp), intent(out) :: U
    real(dp), intent(out) :: virial
    real(dp), intent(in) :: dt
    procedure(TForces) :: forces
    real(dp), intent(out) :: K

    real(dp), dimension(:,:), allocatable :: F
    integer :: err
    real(dp) :: dt1
    real(dp), parameter :: alpha = 1.0_dp/(2.0_dp - 2.0_dp**(1.0/3.0))
    real(dp), parameter :: beta = -2.0_dp**(1.0/3.0)*alpha

    real(dp), parameter :: a1=0.205177661542290_dp 
    real(dp), parameter :: a2=0.403021281604210_dp 
    real(dp), parameter :: a3=-0.12092087633891_dp 
    real(dp), parameter :: a4=0.512721933192410_dp 
    real(dp), parameter :: a5=0.0_dp 

    real(dp), parameter :: b1=0.061758858135626_dp
    real(dp), parameter :: b2=0.33897802655364_dp
    real(dp), parameter :: b3=0.61479130717558_dp
    real(dp), parameter :: b4=-0.14054801465937_dp 
    real(dp), parameter :: b5=0.12501982279453_dp 

    allocate(F(2,Natoms), stat=err)
    if (err.ne.0) STOP 'ALLOCATION ERROR'

    ! See S.K.Gray et al., J.Chem.Phys. 101, 4062 (1994)
    call forces(x,F,U,virial)
    v1 = v + F/Mass * b1 * dt 
    x1 = x + v1 * a1 * dt

    call forces(x1,F,U,virial)
    v1 = v1 + F/Mass * b2 * dt 
    x1 = x1 + v1 * a2 * dt
    
    call forces(x1,F,U,virial)
    v1 = v1 + F/Mass * b3 * dt 
    x1 = x1 + v1 * a3 * dt

    call forces(x1,F,U,virial)
    v1 = v1 + F/Mass * b4 * dt 
    x1 = x1 + v1 * a4 * dt

    call forces(x1,F,U,virial)
    v1 = v1 + F/Mass * b5 * dt 
    
    K=kinetic(v1)
    deallocate(F)

  end subroutine sym4
  
  ! ---------------------------------------------------------------------------------
  ! Symplectic integration of 4th order + prop in tangent space
  subroutine sym4_ly(x,v,xf,vf,U,virial,dt,forces,K,dx,dv,dxf,dvf)
    real(dp), dimension(:,:), intent(in) :: x, v
    real(dp), dimension(:,:,:), intent(in) :: dx, dv
    real(dp), dimension(:,:), intent(out) :: xf, vf
    real(dp), dimension(:,:,:), intent(out) :: dxf, dvf
    real(dp), intent(out) :: U
    real(dp), intent(out) :: virial
    real(dp), intent(in) :: dt
    procedure(TForces_ly) :: forces
    real(dp), intent(out) :: K

    real(dp), dimension(:,:), allocatable :: F
    real(dp), dimension(:,:,:), allocatable :: lF

    integer :: err
    real(dp) :: dt1
    real(dp), parameter :: alpha = 1.0_dp/(2.0_dp - 2.0_dp**(1.0/3.0))
    real(dp), parameter :: beta = -2.0_dp**(1.0/3.0)*alpha

    real(dp), parameter :: a1=0.205177661542290_dp 
    real(dp), parameter :: a2=0.403021281604210_dp 
    real(dp), parameter :: a3=-0.12092087633891_dp 
    real(dp), parameter :: a4=0.512721933192410_dp 
    real(dp), parameter :: a5=0.0_dp 

    real(dp), parameter :: b1=0.061758858135626_dp
    real(dp), parameter :: b2=0.33897802655364_dp
    real(dp), parameter :: b3=0.61479130717558_dp
    real(dp), parameter :: b4=-0.14054801465937_dp 
    real(dp), parameter :: b5=0.12501982279453_dp 

    allocate(F(2,Natoms), stat=err)
    allocate(lF(2,Natoms,4*Natoms), stat=err)
    if (err.ne.0) STOP 'ALLOCATION ERROR'

    ! See S.K.Gray et al., J.Chem.Phys. 101, 4062 (1994)
    call forces(x,dx,F,lF,U,virial)
    vf = v + F/Mass * b1 * dt 
    xf = x + vf * a1 * dt

    dxf = dx + dv * dt
    dvf = dv + dt * lF/Mass

    call forces(xf,dxf,F,lF,U,virial)
    vf = vf + F/Mass * b2 * dt 
    xf = xf + vf * a2 * dt

    dxf = dxf + dvf * dt
    dvf = dvf + dt * lF/Mass

    call forces(xf,dxf,F,lF,U,virial)
    vf = vf + F/Mass * b3 * dt 
    xf = xf + vf * a3 * dt

    dxf = dxf + dvf * dt
    dvf = dvf + dt * lF/Mass

    call forces(xf,dxf,F,lF,U,virial)
    vf = vf + F/Mass * b4 * dt 
    xf = xf + vf * a4 * dt

    dxf = dxf + dvf * dt
    dvf = dvf + dt * lF/Mass

    call forces(xf,dxf,F,lF,U,virial)
    vf = vf + F/Mass * b5 * dt 
    
    dxf = dxf + dvf * dt
    dvf = dvf + dt * lF/Mass
    
    K=kinetic(vf)
    deallocate(F)
    deallocate(lF)

  end subroutine sym4_ly

  function kinetic(v) result(Ek)
     real(dp), dimension(:,:), intent(in):: v
     real(dp) :: Ek
     integer :: n
     
     Ek =0.d0
     do n = 1, Natoms
        Ek = Ek + dot_product(v(:,n),v(:,n))   
     end do
     Ek = Ek*Mass*0.5_dp

  end function kinetic

end module dynamics
