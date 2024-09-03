module boxes
  use constants
  use list
  implicit none
  private
  !!!!!!!!!!!!!!!!!!
  !----------------!
  !!!!!!!!!!!!!!!!!!
  public :: boxlists
  public :: create_boxes
  public :: destroy_boxes
  public :: update_boxes
  public :: boxind
  public :: boxinfo
  public :: folding
  public :: map
  public :: init_map

  type(Tlist), dimension(:,:), allocatable :: boxlists
  integer :: map(2,9)


  ! PRIVATE STUFF:
  integer :: Nx, Ny
  real(dp) :: lx, ly
  real(dp) :: LLx, LLy

  contains

  !               x
  ! |---|---|---|---|---|---|---|---|---|---|
  ! 0 lx          4                         LLx

  subroutine create_boxes(x,y,cutoff)
    real(dp) :: x, y, cutoff
    integer :: err


    LLx=x
    LLy=y

    Nx = floor(LLx/cutoff) ! CON IL FLOOR MI ASSICURO CHE
    Ny = floor(LLy/cutoff) ! CUTOFF SIA MAGGIORE DI lx ...

    lx = LLx/Nx
    ly = LLy/Ny

    ! Algebra modulare, funziona bene partendo da 0
    allocate(boxlists(0:Nx-1,0:Ny-1), stat=err)

    if (err.ne.0) STOP 'ERROR ALLOCATION Boxlist'

  end subroutine create_boxes

  ! ----------------------------------------------------
  subroutine destroy_boxes()
     integer i,j

     do i=0,Nx-1
      do j=0,Ny-1
         call delete_list(boxlists(i,j))
      enddo
     enddo

     deallocate(boxlists)
  end subroutine destroy_boxes
  ! ----------------------------------------------------

  subroutine update_boxes(x,x1)
    real(dp), dimension(:,:) :: x,x1

    integer :: i, natoms, ii,jj,ii1,jj1
    integer :: ii2,jj2
    type(TNode), pointer :: node
    real(dp) :: g(2)

    natoms = size(x,2)

    do i = 1, natoms
       call boxind(x(:,i),ii,jj)
       call boxind(x1(:,i),ii1,jj1)

       if (ii.ne.ii1 .or. jj.ne.jj1) then ! ESEGUO SOLO SE
                                                         ! ESCO AD UNA BOX

          call folding(ii1,jj1,g)

          x1(:,i) = x1(:,i) - g(:) ! SOTTRAGGO LUNGHEZZA DI CUI DEVO RIENTRARE

          ! ---------------------------------------------------
          ! RIMUOVE PARTICELLA i DAL BOX
          call remove(boxlists(ii,jj),node,i)
          if (.not.associated(node)) then
              print*, i,'not found in',ii,jj
              STOP 'error'
          endif
          ! AGGIUNGE PARTICELLE I AL BOX NUOVA
          call add(boxlists(ii1,jj1),node)
       end if
    end do

  end subroutine update_boxes


  !//////////////////////////////////////////////////////////////
  subroutine boxinfo()

    write(*,*) 'Simulation box properties:'
    write(*,*) 'Lx=', LLx,  'Nx=', Nx
    write(*,*) 'Ly=', LLy,  'Ny=', Ny

  end subroutine boxinfo

  !//////////////////////////////////////////////////////////////
  subroutine boxind(r,i,j)
     real(dp) :: r(2)

     integer :: i,j

     i = floor(r(1)/lx)
     j = floor(r(2)/ly)

  end subroutine boxind

  !//////////////////////////////////////////////////////////////
  !
  !   ---||---|---|---|---|---||-o-|---|---|---|---|---|
  !      0                   LLx
  !                             x(t+dt)
  !    x(t+dt)=x(t+dt) + g
  subroutine folding(ii,jj,g)

     integer, intent(inout) :: ii, jj
     real(dp), intent(out) :: g(2)

     integer :: m

     g=0.d0

     ! Versione che puo' foldare di piu' celle
     if (ii.lt.0) then
        m=(ii-Nx+1)/Nx
        g(1)=LLx*m
        ii=mod(ii+16*Nx,Nx)
     else if (ii.gt.Nx-1) then
        m=ii/Nx
        g(1)=LLx*m
        ii=mod(ii,Nx)
     endif

     if (jj.lt.0) then
        m=(jj-Ny+1)/Ny
        g(2)=LLy*m
        jj=mod(jj+16*Ny,Ny)
     else if (jj.gt.Ny-1) then
        m=jj/Ny
        g(2)=LLy*m
        jj=mod(jj,Ny)
     endif


     if (ii<0 .or. jj<0 .or. ii>Nx-1 .or. jj>Ny-1) then
       print*, "FOLDING ERROR", ii, jj
       stop
     end if

  end subroutine folding

  subroutine init_map()
    integer :: u, v, i
    i = 1
      do v=-1, +1
        do u=-1, +1
          map(:,i) = (/u,v/)
          i = i + 1
        end do
      end do
  end subroutine init_map

end module boxes
