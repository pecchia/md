Module forces
  use constants
  use list
  use boxes
  implicit none
  private

  public :: init_lj
  public :: lj
  public :: lj_ly

  type Tpar
     integer :: Natoms
     real(dp) :: eps
     real(dp) :: sigma
     real(dp) :: Rc
     real(dp) :: Fe
  end type

  type(Tpar) :: par
  real(dp) :: sg2, ra2, rc2, Fc, Fa, Uc, Ua

  contains

  subroutine init_lj(Natoms,eps,sigma,Rc)
     integer :: Natoms
     real(dp) :: eps, sigma, Rc
     real(dp) :: rm2, rm6, rm12

     par%Natoms = Natoms
     par%eps = eps
     par%sigma = sigma
     par%Rc = Rc

     ! Compute LJ forces and energy at cutoff
     sg2 = par%sigma*par%sigma
     rc2 = par%Rc*par%Rc
     
     rm2 = sg2/rc2
     rm6 = rm2*rm2*rm2
     rm12 = rm6*rm6
     Fc = 24.0_dp*rm2*(2.0_dp*rm12-rm6)
     Uc = 4.0_dp*(rm12-rm6)

     ! Compute LJ forces and energy at Rmin= sigma/2
     ra2 = sg2/4.0_dp
     rm2 = sg2/ra2
     rm6 = rm2*rm2*rm2
     rm12 = rm6*rm6
     Fa = 24.0_dp*rm2*(2.0_dp*rm12-rm6)
     Ua = 4.0_dp*(rm12-rm6)

     write(*,*) 'Uc=',Uc*par%eps,'Fc=',Fc*par%eps
     write(*,*) 'Um=',Ua*par%eps,'Fm=',Fa*par%eps

  end subroutine init_lj

  ! Lennard jones forces
  subroutine lj(x, F, UU, virial)
    real(dp), dimension(:,:), intent(in) :: x
    real(dp), dimension(:,:), intent(out) :: F

    real(dp), intent(out) :: UU
    real(dp), intent(out) :: virial

    real(dp) :: rij(2), g(2),  Fm(2)
    real(dp) :: r2, rm1, rm2, rm6, rm12, tmp
    integer :: ii, jj, ci, cj, u, err
    integer :: m, l, Natoms, k, p
    type(TNode), pointer :: it

    !integer, external :: omp_get_thread_num

    Natoms = par%Natoms
    ! Virial should be corrected due to cutoff potential
    UU = 0.0_dp
    virial = 0.0_dp
    !$OMP PARALLEL DO DEFAULT(PRIVATE), &
    !$OMP SHARED(Natoms,map,boxlists,Fa,Fc,Ua,Uc,ra2,rc2,sg2,x,F) &
    !$OMP&   REDUCTION( + : UU, virial)
    do m = 1, Natoms

       ! cerca la scatola ci,cj di m
       call boxind(x(:,m),ci,cj)

       Fm = 0.0_dp

       do u = 1, 9

         ii = ci + map(1,u)
         jj = cj + map(2,u)

         ! controlla se la scatola e' una copia periodica
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

             rij(:) = x(:,l) + g(:) - x(:,m)

             r2 = dot_product(rij,rij)

             if (r2 .lt. rc2) then
               rm2 = sg2/r2
               rm6 = rm2*rm2*rm2
               rm12 = rm6*rm6
               tmp = 24.0_dp*rm2*(2.0_dp*rm12-rm6)

               Fm(:) = Fm(:) - (tmp-Fc) * rij(:)
               UU = UU + (4.0_dp*(rm12-rm6)-Uc)
               
               !isotropic pressure = trace of stress tensor
               virial = virial + dot_product(rij,Fm(:))         
             endif

             it => it%next

          end do

        end do
        !print*,omp_get_thread_num(),':',Fm
        F(:,m) = Fm(:)
     end do
     !$OMP END PARALLEL DO

     F = F * par%eps/sg2
     UU = UU * 0.5_dp * par%eps
     virial = virial * 0.5_dp * par%eps/sg2

  end subroutine lj

  ! Lennad-Jones forces + tangent space for lyapunov
  subroutine lj_ly(x,dx,F,lF,UU,virial)
    real(dp), dimension(:,:), intent(in) :: x
    real(dp), dimension(:,:,:), intent(in) :: dx
    real(dp), dimension(:,:), intent(out) :: F
    real(dp), dimension(:,:,:), intent(out) :: lF

    real(dp), intent(out) :: UU
    real(dp), intent(out) :: virial

    real(dp) :: rij(2), g(2), r2, rm1, rm2, rm6, rm12, tmp
    real(dp) :: drij(2), Fm(2), ltmp
    integer :: ii, jj, ci, cj, u, h
    integer :: m, l, Natoms
    type(TNode), pointer :: it
    !integer, external :: omp_get_thread_num

    Natoms = par%Natoms
    ! Virial should be corrected due to cutoff potential
    UU = 0.0_dp
    virial = 0.0_dp
    lF=0.0_dp
    !$OMP PARALLEL DO DEFAULT(PRIVATE), &
    !$OMP& SHARED(Natoms,map,boxlists,Fa,Fc,Ua,Uc,ra2,rc2,sg2,x,F,lF,dx) &
    !$OMP& REDUCTION( + : UU, virial)
    do m = 1, Natoms

       ! cerca la scatola ci,cj di m
       call boxind(x(:,m),ci,cj)

       Fm = 0.0_dp

       do u = 1, 9 

         ii = ci + map(1,u)
         jj = cj + map(2,u)
         ! controlla se la scatola una copia periodica
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

             rij(:) = x(:,l) + g(:) - x(:,m)

             r2 = dot_product(rij,rij)

             if (r2 .lt. rc2) then

               rm2 = sg2/r2
               rm6 = rm2*rm2*rm2
               rm12 = rm6*rm6
               tmp = 24.0_dp*rm2*(2.0_dp*rm12-rm6)

               Fm(:) = Fm(:) - (tmp-Fc) * rij(:)
               UU = UU + (4.0_dp*(rm12-rm6)-Uc)
               virial = virial + dot_product(rij,Fm(:))

               ltmp = 24.0_dp*(8.0_dp*rm6-28.0_dp*rm12)*rm2*rm2

               do h=1,4*Natoms
                 drij(:) = dx(:,l,h) - dx(:,m,h)
                 lF(:,m,h)=lF(:,m,h)-(tmp-Fc)*sg2*drij(:)- &
                     & ltmp*rij(:)*dot_product(drij(:),rij(:))
               end do

             endif

             it => it%next

          end do

        end do
        !print*,omp_get_thread_num(),':',Fm
        F(:,m) = Fm(:)
     end do
     !$OMP END PARALLEL DO

     lF = lF * par%eps/(sg2*sg2)
     F = F * par%eps/sg2
     UU = UU * 0.5_dp * par%eps
     virial = virial * 0.5_dp * par%eps/sg2

  end subroutine lj_ly


end module forces
