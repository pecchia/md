Module vector
  use constants
  use parameters
  implicit none

  public :: grams
  public :: qr

contains



  subroutine grams(q,p,v,u)
    real(dp),dimension(2,Natoms,4*Natoms),intent(inout) :: q, p
    real(dp),dimension(4*Natoms,4*Natoms),intent(out):: V
    real(dp),dimension(4*Natoms,4*Natoms),intent(out):: U

    integer:: i,j
    
    !Combina q-p da 2 array (3, Natoms, 6*Natoms) -> (6Natoms, 6Natoms)
    V=newshape2(q,p)
    
    U(:,1) = V(:,1)/sqrt(dot_product(V(:,1),V(:,1)))

    do i = 2, 4*Natoms
      U(:,i) = V(:,i)
      do j = 1,i-1
        U(:,i) = U(:,i) - U(:,j)*dot_product( U(:,j),U(:,i) )/dot_product(U(:,j),U(:,j))
      end do
      U(:,i) = U(:,i)/sqrt(dot_product(U(:,i),U(:,i)))
    end do

    ! Orthonormalization check 
    !do i=1,6*Natoms
    !   do j=1,6*Natoms
    !      if (j.ne.i) then
    !         if (abs(dot_product(u(:,i),u(:,j))).ge.0.0001) then
    !            write(*,*) 'ortornormal error',i,j
    !            stop
    !         end if
    !      end if
    !   end do
    !end do
   
    !riporta u agli array q p 
    call oldshape(U,q,p)

  end subroutine grams

   
  subroutine qr(q,p,V,lyap_ex)
    real(dp),dimension(2,Natoms,4*Natoms),intent(inout) :: q,p
    real(dp),dimension(4*Natoms,4*Natoms),intent(out):: V
    real(dp),dimension(4*Natoms),intent(inout) :: lyap_ex
    
    integer:: i,j,info,lwork,n,nb
    real(dp),dimension(:,:), allocatable :: RR
    real(dp),allocatable:: tau(:),work(:)
    integer, external :: ilaenv

    n=4*Natoms
    lwork=n
    allocate(RR(n,n))
    allocate(tau(n),work(lwork))

    V=newshape2(q,p)
    RR = V
   
    call dgeqr2p(n,n,RR,n,tau,work,info)

    do i=1,4*Natoms
      if (RR(i,i)> 0.0_dp) then
        lyap_ex(i) = lyap_ex(i) + log(RR(i,i))/tsim
      end if    
    end do

    ! Compute optimal block size 
    nb = ilaenv(1, 'DORGQR', ' ', n, n, n, -1)
    lwork = nb * n
    deallocate(work)
    allocate(work(lwork))
 
    call dorgqr(n,n,n,RR,n,tau,work,lwork,info)
       
    call oldshape(RR,q,p)
 
    deallocate(RR)
    deallocate(work,tau)

  end subroutine qr

  

  function Projection(a,b) result(c)
    real(dp),dimension(4*Natoms),intent(in) :: a,b
      real(dp),dimension(4*Natoms) :: c
      
      c=dot_product(a,b)*a/dot_product(a,a)
      write(*,*)'oooooooooooooooooooooooo'
      write(*,*) dot_product(a,a),dot_product(a,b)
      write(*,*)'oooooooooooooooooooooooo'

  end function Projection



  

  function newshape(a,b) result(d)
    integer:: i,k,j
    real(dp),dimension(2,Natoms,4*Natoms),intent(in) :: a,b
    real(dp),dimension(4,Natoms,4*Natoms) :: c
    real(dp),dimension(4*Natoms,4*Natoms) :: d
     
    do i=1,4*Natoms
       c(2:3,:,i)=a(1:2,:,i)
       c(2:3,:,i)=b(1:2,:,i)
    end do

    do i=1,4*Natoms
      d(:,i)=reshape(c(:,:,i),(/4*Natoms/))
    end do

  end function newshape

   
  function newshape2(a,b) result(d)
    real(dp),dimension(2,Natoms,4*Natoms),intent(in) :: a,b
    real(dp),dimension(4*Natoms,4*Natoms) :: d
     
    integer:: i,k,j
    real(dp),dimension(2*Natoms,4*Natoms) :: e,f
        
    do i=1,4*Natoms
      e(:,i)=reshape(a(:,:,i),(/2*Natoms/))
      f(:,i)=reshape(b(:,:,i),(/2*Natoms/))
    end do  

    do i=1,4*Natoms
      d(1:2*Natoms,i)=e(1:2*Natoms,i)
      d(2*Natoms+1:4*Natoms,i)=f(1:2*Natoms,i)
    end do

   end function newshape2


   subroutine oldshape(a,b,c)
     real(dp),dimension(4*Natoms,4*Natoms),intent(in)::a
     real(dp),dimension(2,Natoms,4*Natoms),intent(out)::b,c

     real(dp),dimension(2*Natoms,4*Natoms):: e, f
     integer:: i,k

     do i=1,4*Natoms
        e(1:2*Natoms,i) = a(1:2*Natoms,i) 
        f(1:2*Natoms,i) = a(2*Natoms+1:4*Natoms,i) 
        b(:,:,i)=reshape(e(:,i),(/2,Natoms/))
        c(:,:,i)=reshape(f(:,i),(/2,Natoms/))
     end do

   end subroutine oldshape
     

 end module vector

 

   


        



     
  
