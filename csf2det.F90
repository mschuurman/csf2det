 MODULE globaldata
 
  type kvpair 
   integer          :: ind
   double precision :: val
  end type kvpair

  ! CSF file
  character*144                             :: csf_file
  ! numerical cutoff for reading CSFs
  double precision                          :: csf_min = 0. 
  ! numerical cutoff for printing determinants
  double precision                          :: det_min = 0.
  ! value of Ms for determinants/csfs
  double precision                          :: m_s

  ! norm of the wavefunction in CSF basis
  double precision                          :: csf_norm
  ! number of internal orbitals, including frozen
  integer                                   :: n_intl
  ! maximum number of external orbitals (i.e. mcscf=0, foci=1, soci=2)
  integer                                   :: n_extl
  ! total number of occupied orbitals in a determinant (i.e. nintl+nextl)
  integer                                   :: n_occ
  ! total number of orbitals in basis
  integer                                   :: n_orb

  ! maximum number of determinants (for array allocation)
  integer                                   :: n_det_max
  ! actual number of determinants
  integer                                   :: n_det
  ! record length -- is nintl + 2 * next (i.e. external orbs, plus orb label)
  integer                                   :: rec_len

  ! vector in which to store representation of corresponding determinants
  integer,dimension(:,:),allocatable        :: det_vec
  ! coefficients in determinant basis
  type(kvpair),dimension(:),allocatable     :: det_cf

 end MODULE globaldata

!
! Program: CSF2DET 
!          -- convert a CI expansion in CSFs to one in determinants
!
!
 PROGRAM csf2det 

  call read_cmdline()

  call parse_csf_list()

  call write_det_list()

  call cleanup()

 end PROGRAM csf2det

!******************************************************
! INITIALIZATION ROUTINES
!******************************************************
 subroutine read_cmdline()
  use globaldata, only: csf_file,csf_min,det_min
  implicit none
  integer                   :: n_arg, i_arg
  character*144             :: abuf

  n_arg = iargc()

  if(n_arg.lt.1)stop 'need to specify input file...'

  ! read in input file
  call getarg(1,abuf)
  csf_file = adjustl(abuf)

  i_arg = 2
  do while (i_arg < n_arg)
    call getarg(i_arg,abuf)

    if (trim(adjustl(abuf)) == '-csf_min') then
       i_arg = i_arg + 1
       call getarg(i_arg,abuf)
       abuf = trim(adjustl(abuf))
       read(abuf,*) csf_min

    elseif (trim(adjustl(abuf)) == '-det_min') then
       i_arg = i_arg + 1
       call getarg(i_arg,abuf)
       abuf = trim(adjustl(abuf))
       read(abuf,*) det_min

    else
       stop ' command line argument argument not recognized...'
           
    endif

    i_arg = i_arg + 1
  enddo   

  return
 end subroutine read_cmdline

!
!
!
 subroutine parse_line(line, cf, csf_vec)
  use globaldata
  implicit none
  character*144, intent(inout) :: line
  real, intent(out)            :: cf
  integer, intent(out)         :: csf_vec(rec_len)
  character*20                 :: str
  integer                      :: ni,ne

  line = adjustl(line)
  str  = line(1:index(line,' ')-1)
  read(str,*)cf

  if (abs(cf).lt.csf_min) then
   cf = 0.
   return
  endif

  ! move start of string to first index past coefficient
  line  = line(index(line,' ')+1:)
  ni      = 0
  ne      = 0
  csf_vec = 0
  do while (len(trim(line)).gt.0)
    str = line(1:index(line,' ')-1)
    if(index(str,':').ne.0) then
     ne = ne + 1
     read(str(1:index(str,':')-1),*)csf_vec(n_occ+ne) 
     line = adjustl(line(index(line,' ')+1:))
     str = line(1:index(line,' ')-1)
     read(str,*)csf_vec(ne)
    else
     ni = ni + 1
     read(str,*)csf_vec(n_extl+ni)
    endif
    if (index(line,' ').eq.0) exit
    line = adjustl(line(index(line,' ')+1:))
  enddo

  if (ni.ne.n_intl.or.ne.ne.n_extl)stop 'error parsing CSF file, ni!=n_intl.or.ne!=n_extl'

 end subroutine
 
!
!  Read the ciplcs file:
!      1. determine the size of the orbital spaces
!      2. figure out how many CSFs there are
!
 subroutine parse_csf_list()
  use globaldata
  implicit none
  integer                :: ios
  integer                :: n_csf
  integer,allocatable    :: csf_vec(:)
  integer                :: cfile=99
  integer                :: orb_occ
  real                   :: cf
  character*144          :: line
  character*20           :: str
  real,dimension(0:3)    :: spin = (/0.0, 0.5, -0.5, 0.0/)

  ! open the csf file --
  ! this file is formated such that there is one csf per line. A CSF is specified by:
  ! CF  X X X X X .. ORB: X ORB: X
  ! where CF is the coefficient of the CSF, X is an occupation value (==1,2,3) and
  ! "ORB" is an orbital index for excitations out of the frzn+internal space
  !
  open(unit=cfile,file=trim(csf_file));

  ! get the number of CSFs in the file
  n_csf=0
  do 
   read(cfile, 1000, iostat=ios)line
   if(ios.lt.0)exit

   ! get rid of leading and trailing white space and grab leading coefficient
   ! if coefficient smaller than cutoff, we're done
   line = trim(adjustl(line))
   str = line(1:index(line,' ')-1)
   read(str,*)cf
   if (abs(cf).lt.csf_min) exit
  
   ! else, increment the csf counter
   n_csf = n_csf + 1

   ! if this is the first line, also figure out how many orbital indices
   ! we need to store
   if(n_csf.eq.1) then

    ! advance string to first orbital index
    line = adjustl(line(index(line,' ')+1:))

    n_intl = 0
    n_extl = 0 
    m_s    = 0.
    do while (len(trim(adjustl(line))).gt.0) 
      str = line(1:index(line,' ')-1)
      if(index(str,':').ne.0) then
       n_extl = n_extl + 1
       line   = adjustl(line(index(line,' ')+1:))
       str    = line(1:index(line,' ')-1)
      else
       n_intl = n_intl + 1
      endif
      ! read the orbital occupation
      read(str,*)orb_occ
      m_s = m_s + spin(orb_occ)
      if (index(line,' ').eq.0) exit
      line = adjustl(line(index(line,' ')+1:))
    enddo
   endif
  enddo

  ! initially allocate the determinant array to 5 * CSFs
  n_occ     = n_intl + n_extl
  rec_len   = n_intl + 2*n_extl
  n_det_max = 5 * n_csf
  ! initialize total number of orbitals to number of occupied
  n_orb     = n_occ
  n_det     = 0
  csf_norm  = 0.
  allocate(csf_vec(rec_len))
  allocate(det_vec(rec_len,n_det_max))
  allocate(det_cf(n_det_max))

  ! go back to the beginning of the file
  rewind(cfile)

  do
   read(cfile, 1000, iostat=ios)line
   if(ios.lt.0)exit
   if(ios.gt.0)stop 'error reading csf_file'

   ! read the CSF
   call parse_line(line, cf, csf_vec)
   !print *,'csf_vec=',csf_vec

   ! compute the norm of the csf wfn for diagnostic purposes
   csf_norm = csf_norm + cf**2

   ! if the returned coefficient is zero, we've read all csfs less than cutoff
   if (cf == 0.) exit

   ! check the external orbital index: update n_orb if necessary
   n_orb = max(n_orb, maxval(csf_vec(n_occ:)))

   ! else, convert the csf to a linear combination of determinants
   call unroll_csf(cf, csf_vec)  

  enddo

  csf_norm = sqrt(csf_norm)
  
  close(cfile)

  return
1000 format(144a)
1001 format(a9,i5,a9,i5,a9,i5,a9,i5)
 end subroutine parse_csf_list 

!
!  Project out the contribution of allowed determinants from a CSF
!
!
 subroutine unroll_csf(csf_cf, csf_vec)
  use globaldata
  implicit none
  real, intent(in)    :: csf_cf
  integer, intent(in) :: csf_vec(rec_len)
 
  integer             :: db(4),mz2(2),d1f(2),d2f(2),del(2)
  integer             :: num,denom,sgn
  integer             :: j,k,l,icnt,bt,found
  integer             :: ms2,m2,nalpha,nopen,nloops
  integer             :: bvec(n_occ),aloc(n_occ),oopen(n_occ)
  integer             :: idet(rec_len),det(rec_len),refdet(rec_len)
  double precision    :: cf
  logical             :: zero
  integer             :: ifac
  integer             :: oparity
  double precision    :: eps=1.d-8

  db  = (/ 0,  1, -1, 0 /)
  mz2 = (/ 1, -1 /)
  d1f = (/ 1, -1 /)
  d2f = (/-1,  1 /) 
  del = (/ 1,  0 /)

  ms2   = int(2.*m_s)
  bt    = 0
  nopen = 0
  oopen = 0
  bvec  = 0
  aloc  = 0

  refdet = csf_vec

  do j = 1,n_occ
   bt = bt + db(csf_vec(j)+1)
   bvec(j) = bt
   if(csf_vec(j).eq.1.or.csf_vec(j).eq.2)then
    nopen = nopen + 1
    oopen(nopen) = j
    refdet(j) = 2
   endif
  enddo !do j = 1,nocc
 
  ! set alpha spin counter
  nalpha = (ms2 + nopen)/2
  do j = 1,nalpha
   aloc(j) = j
  enddo
  if(nalpha.gt.0)aloc(nalpha) = aloc(nalpha)-1
  nloops = ifac(nopen)/(ifac(nalpha)*ifac(nopen-nalpha))

  !
  !loop over the allowed determinants
  !
  do j = 1,nloops
   det = refdet

   if(nalpha.gt.0)then
    ! loop over all perumtations of alpha/beta occupations
    icnt = nalpha
    do while(aloc(icnt) .eq. (nopen-nalpha+icnt))
     icnt = icnt - 1
    enddo
    aloc(icnt) = aloc(icnt) + 1
    do k = icnt+1,nalpha
      aloc(k) = aloc(icnt)+k-icnt
    enddo

    ! set selected unpaired electrons to alpha spin
    do k = 1,nalpha
     det(oopen(aloc(k)))=1
    enddo
   endif

   cf = 1.
   sgn = 1
   m2 = 0
   zero = .false.
   do k = 1,n_occ
    select case(csf_vec(k))

      case(1)
        m2 = m2 + mz2(det(k))
        num = bvec(k) + d1f(det(k)) * m2
        if(num == 0)then
         zero = .true.
         exit
        endif
        denom = 2 * bvec(k) 
        cf = cf * num / denom

      case(2)
        m2 = m2 + mz2(det(k))
        num = bvec(k) + 2 + d2f(det(k)) * m2
        if(num == 0)then
         zero = .true.
         exit
        endif
        denom = 2 * (bvec(k) + 2)
        sgn = sgn * (-1)**(bvec(k)+del(det(k)))
        cf = cf * num / denom

      case(3)
        sgn = sgn * (-1)**(bvec(k))

      case default
     end select
    enddo ! do k =1 ,nocc 

    if(.not.zero) then

      ! convert the determinant to multigrid format
      call convert_det(det,idet)
      ! include parity in value of cf -- necessary to compute S^2 consistently
      if(cf.lt.0)stop 'ERROR computing determinant'
      cf = sqrt(cf)*sgn*oparity(idet)

      found = 0
      do k = 1,n_det
       if(all(det_vec(:,k).eq.idet))then
         found = k
         exit
       endif
      enddo

      if(found.ne.0)then
       det_cf(found)%val = det_cf(found)%val + cf*csf_cf
      else
       n_det = n_det + 1
       if(n_det.gt.n_det_max) stop 'ndet .gt. numdet, increase rdet2csf'
       det_cf(n_det)%ind = n_det
       det_cf(n_det)%val = cf*csf_cf
       det_vec(:,n_det)  = idet
      endif

    endif !if(.not.zero) 
   enddo !do j = 1,nloops

  return
1000 format(20(i3))
 end subroutine unroll_csf 


!
! write the determinant list to standard output
!
!
 subroutine write_det_list
  use globaldata
  implicit none
  integer            :: i
  integer            :: ofile=99
  double precision   :: norm,s2

  if(n_det.eq.0)return

  ! sort so determinants printed largest coefficient to smallest
  call sort_det_list()

  ! print determinants
  i = 1 
  do while(abs(det_cf(i)%val).gt.det_min)
   call print_det(det_cf(i)%val,det_vec(:,det_cf(i)%ind))
   i = i + 1
   if(i.gt.n_det)exit
  enddo

  ! compute S^2 and norm
  call compute_norm_s2(norm, s2)
 
  ! not sure how to best handle this.  For the time being, useful to check
  ! wavefunction norm and S^2 to make sure cutoff is not too severe.
  open(unit=ofile,file='csf2det.out',status='replace')
  write(ofile,1000)csf_norm
  write(ofile,1001)norm
  write(ofile,1002)s2
  close(ofile)  

  return
1000 format('norm of wfn in CSF basis: ',f15.10)
1001 format('norm of wfn in det basis: ',f15.10)
1002 format('S^2: ',f8.4)
 end subroutine write_det_list

!
! deallocate dynamic memory arrays
!
 subroutine cleanup()
  use globaldata
  implicit none

  deallocate(det_vec)
  deallocate(det_cf)

  return
 end subroutine cleanup

!***********************************************************************
!*
!*   Utility Functions
!*
!**********************************************************************
!
!
! sorts the coefficient array, and maps a make for the vector list
!
!
 subroutine sort_det_list
  use globaldata
  implicit none
  integer                         :: i
  type(kvpair),allocatable        :: scr(:)

  allocate(scr(int((n_det+1)/2)))
  call mergesort(det_cf(1:n_det),n_det,scr)
  deallocate(scr)

  return
 end subroutine sort_det_list


!
! Convert a determinant in COLUMBUS 'step' notation to MULTIGRID notation
!
 subroutine convert_det(indet,outdet)
  use globaldata, only: rec_len,n_occ
  implicit none
  integer,intent(in)  :: indet(rec_len)
  integer,intent(out) :: outdet(rec_len)
  integer             :: i,ivec(4)

  ivec = (/0, 1, -1, 2 /)
  do i = 1,n_occ
   outdet(i) = ivec(indet(i)+1)
  enddo
  outdet(n_occ+1:rec_len) = indet(n_occ+1:rec_len)

  return
 end subroutine convert_det

!
!  compute effect of S2 operator on a single determinant
!
  subroutine compute_norm_s2(norm,s2)
   use globaldata, only: rec_len,n_occ,n_det,det_cf,det_vec
   implicit none
   double precision,intent(out) :: norm
   double precision, intent(out):: s2
   integer                      :: i,j,k,l
   integer                      :: trial(rec_len)
   integer                      :: par(n_det)
   integer                      :: oparity

   s2 = 0.
   norm = 0.
 
   do i = 1,n_det
    par(i) = oparity(det_vec(:,i))
   enddo

   do i = 1,n_det
    norm = norm + det_cf(i)%val**2
    do j = 1,n_occ
     if(abs(det_vec(j,i)).ne.1)cycle
     s2 = s2 + 0.5 * det_cf(i)%val**2 ! diagonal contribution of the ladder operators

     do k = 1,n_occ
      if(abs(det_vec(k,i)).ne.1) cycle
      s2 = s2 + 0.25 * det_vec(j,i) * det_vec(k,i) * det_cf(i)%val**2  ! contribution from Sz^2

      if(det_vec(j,i).ne.det_vec(k,i)) then
       trial    = det_vec(:,i)
       trial(j) = -trial(j)
       trial(k) = -trial(k)
       
       do l = 1,n_det
        if(any(det_vec(:,l).ne.trial)) cycle
 
        s2 = s2 + 0.5 * det_cf(i)%val * det_cf(l)%val * par(i) * par(l)

       enddo ! end l = 1,ndet
      endif ! if(det(i,j)!=det(i,k)
     enddo ! k = 1,norb
    enddo ! j = 1,norb
   enddo ! i = 1,ndet

   norm = sqrt(norm)

   return
  end subroutine compute_norm_s2

!
! prints a determinant in multigrid format
!
 subroutine print_det(cf,det)
  use globaldata, only: rec_len,n_intl,n_extl,n_occ,n_orb
  implicit none
  double precision,intent(in)       :: cf
  integer,intent(in)                :: det(rec_len)
  integer                           :: n_vrt
  integer                           :: dif(n_extl)
  integer                           :: i,orb_i,ind
  character*4                       :: lstr
  character*3                       :: istr
  character*17                      :: ofmt
  character*(3*n_orb)               :: ovec

  n_vrt = n_orb - n_intl
!  print *,'printing: ',det

  write(lstr,'(i4)')3*n_orb
  ofmt = '(f15.10,1x,a'//trim(adjustl(lstr))//')'

  ovec = ''
  do i = n_extl+1,n_occ
   if(abs(det(i)).eq.1)then
    write(istr,'(sp,i3)')det(i)
   else
    write(istr,'(ss,i3)')det(i)
   endif
   ovec = trim(adjustl(ovec))//adjustr(istr)
  enddo

  do i = 1,n_vrt
   orb_i = n_intl + i
   if(any(det(n_occ+1:)==orb_i)) then
    dif = abs(det(n_occ+1:)-orb_i)
    ind = minloc(dif,dim=1)
    if (abs(det(ind)).eq.1) then
      write(istr,'(sp,i3)')det(ind)
    else
      write(istr,'(ss,i3)')det(ind)
    endif
    ovec = trim(adjustl(ovec))//adjustr(istr)
   else
    ovec = trim(adjustl(ovec))//'  0'
   endif
  enddo

  write(*,ofmt)cf,trim(adjustl(ovec))

  return
 end subroutine print_det


!
! Returns n!
!
 recursive function ifac(n) result(fac)
   integer, intent(in) :: n
   integer             :: fac

   if(n.eq.0)then
    fac = 1
   else
    fac = n * ifac(n-1)
   endif
   return
 end function

!
!  merge two lists in descending order
! 
 subroutine mergearr(a,b,c,na,nb,nc)
   use globaldata
   integer, intent(in)         :: na,nb,nc   ! Normal usage: NA+NB = NC
   type(kvpair), intent(inout) :: a(na)      ! B overlays C(NA+1:NC)
   type(kvpair), intent(in)    :: b(nb)
   type(kvpair), intent(inout) :: c(nc)

   integer :: i,j,k

   i = 1
   j = 1
   k = 1
   do while(i.le.na .and.j.le.nb)
    if(abs(a(i)%val).ge.abs(b(j)%val)) then
       c(k) = a(i)
       i = i + 1
     else
       c(k) = b(j)
       j = j + 1
     endif
     k = k + 1
   enddo
   do while (i.le.na)
      c(k) = a(i)
      i = i + 1
      k = k + 1
   enddo
   do while (j.le.nb)
      c(k) = b(j)
      j = j + 1
      k = k + 1
   enddo
  
   return
 
 end subroutine mergearr
 
!
! sort list by divid and conquer
!
 recursive subroutine mergesort(a,n,t)
   use globaldata 
   integer, intent(in)                :: n
   type(kvpair), intent(inout)        :: a(n)
   type(kvpair), intent(out)          :: t(int((n+1)/2))
   integer                            :: na,nb
   type(kvpair)                        :: dval

   if (n.lt.2) return
   if (n.eq.2) then
    if (abs(a(1)%val).lt.abs(a(2)%val)) then
     dval = a(1)
     a(1) = a(2)
     a(2) = dval
    endif
    return
   endif      
   na=int((n+1)/2)
   nb=n-na
 
   call mergesort(a(1:na),na,t)
   call mergesort(a(na+1:n),nb,t)
 
   if(abs(a(na)%val).lt.abs(a(na+1)%val)) then
      t(1:na)=a(1:na)
      call mergearr(t(1:na),a(na+1:n),a(1:n),na,nb,n)
   endif

   return
 
 end subroutine mergesort

!
!
! get the sign change upon putting a determinant in alpha string notation
!
 integer function oparity(det)
  use globaldata, only: rec_len,n_occ,n_extl
  implicit none
  integer,intent(in)           :: det(rec_len)
  integer                      :: dsort(2*n_occ)
  integer                      :: i,j,na,nb,at,nt
  integer                      :: icnt,nbnd,newbnd
  logical                      :: sorted

  at = count(det(1:n_occ)==1) + count(det(1:n_occ)==2)
  nt = 0
  na = 0
  ! initialize n[beta] to the end of the alpha string
  nb = at 

  ! loop over all internal orbitals 
  do i = n_extl+1,n_occ
   if(any(det(i)==(/1, 2/)))then
    nt = nt + 1;na = na + 1;dsort(nt) = na
   endif
   if(any(det(i)==(/-1, 2/)))then
    nt = nt + 1;nb = nb + 1;dsort(nt) = nb
   endif
  enddo

  ! loop over all external orbitals 
  do i = 1,n_extl
   if(any(det(i)==(/1, 2/)))then
    nt = nt + 1;na = na + 1;dsort(nt) = na
   endif
   if(any(det(i)==(/-1, 2/)))then
    nt = nt + 1;nb = nb + 1;dsort(nt) = nb
   endif
  enddo

  if(na.ne.at.or.nb.ne.nt)stop 'oparity ERROR!'

  oparity = 1
  sorted  = .false.

  do while(.not.sorted)
   sorted = .true.
   do i = 2,nt
     if(dsort(i-1)<dsort(i)) cycle
     j          = dsort(i-1)
     dsort(i-1) = dsort(i)
     dsort(i)   = j
     oparity    = -oparity
     sorted     = .false. 
   enddo
  enddo

  return
 end function oparity
