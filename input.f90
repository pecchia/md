module input
  use parameters
  implicit none
  private

  public :: parse_input
 
  type(TLyapunovEnum), parameter, public :: LyAlgo = TLyapunovEnum()

  contains
  
  subroutine parse_input()

    character(20) :: token, val
    character(256) :: line
    integer :: idx1, idx2, err

    do  
      read(*,fmt='(A)',iostat=err) line
      if (err/=0) then
         exit
      endif   
      idx1 = index(trim(line),"=")
      idx2 = index(trim(line),"!")
      if (idx1==0 .or. idx2==1) cycle

      token = line(1:idx1-1)
      if (idx2>0) then
        val = line(idx1+1:idx2-1)
      else
        val = line(idx1+1:idx1+20)
      end if 
      write(*,'(A10,3x,A)') trim(adjustl(token)), trim(adjustl(val))


      if (trim(val)=="") then
        write(*,*) trim(token)    
        stop "error empty val"
      end if   

      call readvalue(token, val)
    end do   
 
  end subroutine parse_input


  subroutine readvalue(token, val)
    character(20) :: token, val

    if (to_lower(trim(token))=="natoms") then
       read(val,*) Natoms
    else if (to_lower(trim(token))=="lx") then
       read(val,*) Lx
    else if (to_lower(trim(token))=="ly") then
       read(val,*) Ly
    else if (to_lower(trim(token))=="nx") then
       read(val,*) Nx
    else if (to_lower(trim(token))=="ny") then
       read(val,*) Ny
    else if (to_lower(trim(token))=="rc") then
       read(val,*) Rc
   else if (to_lower(trim(token))=="temp") then
       read(val,*) Temp
    else if (to_lower(trim(token))=="mass") then
       read(val,*) Mass
    else if (to_lower(trim(token))=="eps") then
       read(val,*) eps
    else if (to_lower(trim(token))=="sigma") then
       read(val,*) Sigma
    else if (to_lower(trim(token))=="tinit") then
       read(val,*) Tinit
    else if (to_lower(trim(token))=="tsim") then
       read(val,*) Tsim
    else if (to_lower(trim(token))=="dt") then
       read(val,*) dt
    else if (to_lower(trim(token))=="scaling") then
       read(val,*) Scaling
    else if (to_lower(trim(token))=="nose_hoover") then
       read(val,*) nose_hoover
    else if (to_lower(trim(token))=="print_xyz") then
       read(val,*) print_xyz
    else if (to_lower(trim(token))=="print_interval") then
       read(val,*) print_interval 
    else if (to_lower(trim(token))=="file_name") then
       read(val, *) file_name
    else if (to_lower(trim(token))=="nh_mass") then
       read(val,*) Q
    else if (to_lower(trim(token))=="dr") then
       read(val,*) dr
    else if (to_lower(trim(token))=="lyapunov") then
       read(val,*) do_lyapunov
    else if (to_lower(trim(token))=="lyapunovalgo") then
       if (to_lower(trim(val)) == "volumes") then 
         algorithm = LyAlgo%volumes
       else if (to_lower(trim(val)) == "lengths") then    
         algorithm = LyAlgo%lengths
       else if (to_lower(trim(val)) == "qr") then   
         algorithm = LyAlgo%QR
       end if
    else if (to_lower(trim(token))=="tortho") then
       read(val,*) tgram
    else if (to_lower(trim(token))=="d0") then
       read(val,*) D0
    else if (to_lower(trim(token))=="dtemp") then
       read(val,*) Dtemp
    else if (to_lower(trim(token))=="tempf") then
       read(val,*) Tempf
    else if (to_lower(trim(token))=="shnn") then
       read(val,*) shnn
    else if (to_lower(trim(token))=="tsh") then
       read(val,*) tsh
    else if (to_lower(trim(token))=="fext") then
       read(val,*) Fe
    else 
       write(*,*) "Token '",to_lower(trim(token)),"' not recognized"  
    end if

  end subroutine readvalue

  function to_lower(strIn) result(strOut)
    character(*), intent(in) :: strIn
    character(len=len(strIn)) :: strOut

    integer :: i,j 

     do i = 1, len(strIn)
        j = iachar(strIn(i:i))
        if (j>= iachar("A") .and. j<=iachar("Z") ) then
           strOut(i:i) = achar(iachar(strIn(i:i))+32)
        else
           strOut(i:i) = strIn(i:i)
        end if
     end do

  end function to_lower

end module input

