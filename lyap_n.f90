module lyap_n
  use constants
  use parameters
  implicit none
  private

  public:: lyap_numbers
  public:: lyap_numbers2
  public:: write_Lyap

contains

  subroutine lyap_numbers(a,c,n)
    integer, intent(in) :: n 
    real(dp),dimension(4*Natoms,4*Natoms,n),intent(in):: a
    real(dp),dimension(4*Natoms),intent(out):: c
  
    integer:: i,j,pq,h,l
    integer:: err
    real(dp),dimension(:,:,:),allocatable :: G
    real(dp),dimension(:,:),allocatable :: Vol,Vol1,lh
    real(dp),dimension(:),allocatable:: lyapsommati
    
    allocate(G(4*Natoms,4*Natoms,n),stat=err)
    allocate(Vol(4*Natoms,n),stat=err)
    allocate(Vol1(4*Natoms,n),stat=err)
    allocate(lh(4*Natoms,n),stat=err)
    allocate(lyapsommati(4*Natoms),stat=err)
    if (err /= 0) STOP 'ALLOCATION ERROR'
    
    G=0.d0
    Vol=0.d0
    Vol1=0.d0
    lyapsommati=0.d0
    c=0.d0

    ! Norma vettore
    open(702,file='LD1.dat')
    open(703,file='LD96.dat')
    open(704,file='LD192.dat')
    open(705,file='vol10.dat')
    open(706,file='vol100.dat')
    open(707,file='vol500.dat')


    do pq=1,n
      write(*,*) pq,n
       
      do i=1,4*Natoms
        if (i==1) then
          Vol(i,pq)=sqrt(dot_product(a(:,i,pq),a(:,i,pq)))
        else
          Vol(i,pq)=sqrt(GramMatrix(a(:,1:i,pq),i,pq))
        end if

        !if (vol(i,pq) .eq. 0) then
        !   write(*,*) "zero volume"
        !   write(*,*) Vol(i-1,pq)
        !   write(*,*) Vol(i,pq)
        !   do h=1,i
        !      write(*,*)'..............................................',i
        !      write(*,*) a(:,h,pq)
        !   end do
        !
        !   do h=1,i
        !      write(*,*)'..............................................',i
        !      do l=1,i
        !         if (l.eq.h) then
        !            write(*,*)'--------norm-------'
        !            write(*,*) dot_product(a(:,h,pq),a(:,l,pq))
        !            write(*,*)'--------norm-------'
        !         else
        !            write(*,*) dot_product(a(:,h,pq),a(:,l,pq))
        !         end if
        !      end do 
        !   end do
        !   stop
        !end if

      end do
    end do
    
    do i=1,4*Natoms
       write(705,*) vol(i,10)
    end do
    do i=1,4*Natoms
       write(706,*) vol(i,100)
    end do
    do i=1,4*Natoms
       write(707,*) vol(i,500)
    end do

    do i=1,4*Natoms
      Vol1(i,:)=log(Vol(i,:))
      if (i==1) then
         lh(i,:) = vol1(i,:) 
      else
         lh(i,:) = vol1(i,:)-vol1(i-1,:)
      end if
    end do

    do i=1,n
       write(702,*) i,lh(1,i)*LJTU*n/(tsim)
       write(703,*) i,lh(96,i)*LJTU*n/(tsim)
       write(704,*) i,lh(192,i)*LJTU*n/(tsim)
    end do

    close(702)
    close(703)
    close(704)
    close(705)
    close(706)
    close(707)

    do i=1,4*Natoms
       lyapsommati(i)=sum(Vol1(i,:))/(tsim)
    end do

    do i=1,6*Natoms
       if (i==1) then
         c(i) = lyapsommati(i) 
       else
         c(i) = lyapsommati(i)-lyapsommati(i-1)
       end if   
    end do

  end subroutine lyap_numbers



  subroutine lyap_numbers2(a,c,n)
    integer, intent(in) :: n
    real(dp),dimension(4*Natoms,4*Natoms,n),intent(in):: a
    real(dp),dimension(4*Natoms),intent(out):: c
  
    integer :: i,pq 
    integer:: err
    real(dp),dimension(:,:),allocatable :: leng
    
   
    allocate(leng(4*Natoms,n),stat=err)

    if (err /= 0) STOP 'ALLOCATION ERROR'

    leng=0.d0
    c=0.d0

    do pq = 1, n
      do i = 1, 4*Natoms
         leng(i,pq) =log(sqrt(dot_product(a(:,i,pq),a(:,i,pq))))
      end do
    end do

    do i=1,4*Natoms
      c(i) = sum(leng(i,:))/tsim
    end do
  
  end subroutine lyap_numbers2

  ! ---------------------------------------------------------
  subroutine write_lyap(Lyapunov,Lp,Lm)
    real(dp),dimension(4*Natoms)::Lyapunov
    real(dp),dimension(Lp):: L_plus
    real(dp),dimension(Lm):: L_minus
    integer :: i,Lp,Lm
    
    open(104,file='lyap_pl.dat')
    open(106,file='lyap_mi.dat')
  
    Lp=0
    Lm=0
    do i=1,4*Natoms
       if (Lyapunov(i) .ge.0) then
          Lp=Lp+1
          L_plus(Lp)=Lyapunov(i)
       else if(Lyapunov(i) .lt.0) then
         Lm=Lm+1
        L_minus(Lm)=Lyapunov(i)
        end if
    end do
  
    do i=1,Lp
       write(104,*) i,L_plus(i)*(2170)
    end do
    
    do i=1,Lm
       write(106,*)  i, L_minus(i)*(2170)
    end do
     
    write(*,*) sum(L_plus)*(2170),sum(L_minus)*(2170),(sum(L_plus)+sum(L_minus))*(2170)
    write(*,*) Lp,Lm
   
    close(104)
    close(106)
       
   end subroutine write_lyap

   ! ---------------------------------------------------------
   function GRAMmatrix(Vectors,n,pq) result(f)
  
     real(dp),dimension(4*natoms,n) :: Vectors
     real(dp),dimension(n,n) :: G_matrix
     real(dp)::f
     integer,intent(in) :: n
     integer::p,q,pq
  
     do p=1,n
       do q=1,n
          G_matrix(p,q)=dot_product(Vectors(:,p),Vectors(:,q))
       end do
     end do
  
     f= finddet(G_matrix,n)
      
    end function GRAMmatrix
  
  
    ! ---------------------------------------------------------
    ! Finds the determinant of a Upper tri-diagonal matrix
    real(dp) function FindDet(matrix, n) 
  
      REAL(dp), DIMENSION(n,n) :: matrix
      INTEGER, INTENT(IN) :: n
      REAL(dp) :: m, temp
      INTEGER :: i, j, k, l
      LOGICAL :: DetExists = .TRUE.
      l = 1
      !Convert to upper triangular form
      DO k = 1, n-1
          IF (matrix(k,k) == 0) THEN
              DetExists = .FALSE.
              DO i = k+1, n
                  IF (matrix(i,k) /= 0) THEN
                      DO j = 1, n
                          temp = matrix(i,j)
                          matrix(i,j)= matrix(k,j)
                          matrix(k,j) = temp
                      END DO
                      DetExists = .TRUE.
                      l=-l
                      EXIT
                  ENDIF
              END DO
              IF (DetExists .EQV. .FALSE.) THEN
                  FindDet = 0
                  return
              END IF
          ENDIF
          DO j = k+1, n
              m = matrix(j,k)/matrix(k,k)
              DO i = k+1, n
                  matrix(j,i) = matrix(j,i) - m*matrix(k,i)
              END DO
          END DO
      END DO
     
      !Calculate determinant by finding product of diagonal elements
      FindDet = l
      DO i = 1, n
          FindDet = FindDet * matrix(i,i)
      END DO
     
   end function FindDet


end module lyap_n

