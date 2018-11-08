 MODULE globaldata
 
  type kvpair 
   integer          :: ind
   double precision :: val
  end type kvpair

  ! CSF file
  character*144                             :: csf_file
  ! numerical cutoff for reading CSFs
  double precision                          :: ccutoff = 1e-5 
  ! numerical cutoff for printing determinants
  double precision                          :: dcutoff = 1e-5
  ! value of Ms for determinants/csfs
  double precision                          :: ms

  ! number of internal orbitals, including frozen
  integer                                   :: n_intl
  ! maximum number of external orbitals (i.e. mcscf=0, foci=1, soci=2)
  integer                                   :: n_extl
  ! total number of occupied orbitals (i.e. nintl+nextl)
  integer                                   :: n_occ

  ! maximum number of determinants 
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
  use globaldata, only: csf_file,ccutoff,dcutoff
  implicit none
  integer                   :: narg
  character*144             :: abuf

  narg = iargc()

  if(narg.lt.1)stop 'need to specify input file...'

  ! read in input file
  call getarg(1,abuf)
  csf_file = adjustl(abuf)

  iarg = 2
  do while iarg < narg
    call getarg(iarg,abuf)

    if trim(adjustl(abuf)) == '-cmin' then
       iarg += 1
       call getarg(iarg,abuf)
       read(trim(adjustl(abuf)),*) ccutoff

    elseif trim(adjustl(abuf)) == '-dmin' then
       iarg += 1
       call getarg(iarg,abuf)
       read(trim(adjustl(abuf)),*) dcutoff

    else
       stop ' command line argument argument not recognized...'
           
    endif

    iarg += 1
  enddo   

  return
 end subroutine read_cmdline

!
!
!
 subroutine parse_line(line, cf, csf_vec)
  implicit none
  character*144, intent(inout) :: line
  real, intent(out)            :: cf
  integer, intent(out)         :: csf_vec(:)
  character*20                 :: str
  integer                      :: ni,ne

  line = trim(adjustl(line))
  str  = line(1:index(line,' ')-1)
  read(str,*)cf

  if (cf.lt.ccutoff) then
   cf = 0.
   return
  endif

  line  = line(index(line,' ')+1:)
  ni = 0
  ne = 0
  do while (len(line).gt.0)
    str = line(1:index(line,' ')-1)
    if(index(str,':').ne.0) then
     ne += 1
     read(str,*)csf_vec(nocc+ne) 
     line = line(index(line,' ')+1:)
     str = line(1:index(line,' ')-1)
     read(str,*)csf_vec(nintl+ne)
    else
     ni += 1
     read(str,*)csf_vec(ni)
    endif
    if (index(line,' ').eq.0) exit
  enddo
 endif

 end subroutine
  
!
!  Read the ciplcs file:
!      1. determine the size of the orbital spaces
!      2. figure out how many CSFs there are
!
 subroutine create_csf_list()
  use globaldata
  implicit none
  integer                :: ios
  integer                :: ncsf
  integer                :: cfile=99
  real                   :: cf
  character*144          :: line
  character*20           :: str
  real                   :: spin = (/0, 0.5, -0.5, 0/)
  ! open the csf file --
  ! this file is formated such that there is one csf per line. A CSF is specified by:
  ! CF  X X X X X .. ORB: X ORB: X
  ! where CF is the coefficient of the CSF, X is an occupation value (==1,2,3) and
  ! "ORB" is an orbital index for excitations out of the frzn+internal space
  !
  open(unit=cfile,file=trim(csf_file));

  ! get the number of CSFs in the file
  ncsf=0
  do 
   read(cfile, 1000, iostat=ios)line
   if(ios.lt.0)exit

   ! get rid of leading and trailing white space and grab leading coefficient
   ! if coefficient smaller than cutoff, we're done
   line = trim(adjustl(line))
   str = line(1:index(line,' ')-1)
   read(str,*)cf
   if (cf.lt.ccutoff) exit

   ! else, increment the csf counter
   ncsf += 1

   ! if this is the first line, also figure out how many orbital indices
   ! we need to store
   if(ncsf.eq.1) then
    nintl = 0
    nextl = 0 
    line = line(index(line,' ')+1:)
    do while (len(line).gt.0) 
      str = line(1:index(line,' ')-1)
      if(index(str,':').ne.0) then
       nextl += 1
       line   = line(index(line,' ')+1:)
      else
       nintl += 1
      endif
      if (index(line,' ').eq.0) exit
      line = line(index(line,' ')+1:)
    enddo
   endif
  enddo

  ! initially allocate the determinant array to 5 * CSFs
  nocc     = nintl + nextl
  ndet_max = 5 * ncsf
  allocate(csf_vec(nintl+2*nextl))
  allocate(det_vec(nintl+2*nextl,ndet_max))

  ! go back to the beginning of the file
  rewind(cfile)

  do
   read(ifile, 1000, iostat=ios)line
   if(ios.lt.0)exit
   if(ios.gt.0)stop 'error reading csf_file'

   call parse_line(line, cf, csf_vec)

   ! if the returned coefficient is zero, we've read all csfs less than cutoff
   if (cf.eq.0) exit

   ! else, convert the csf to a linear combination of determinants
   call unroll_csf(cf, csf_vec)  

  enddo

  close(ifile)

  return
1000 format(144a)
1001 format(a9,i5,a9,i5,a9,i5,a9,i5)
 end subroutine parseInput

!
!  Project out the contribution of allowed determinants from a CSF
!
!
 subroutine unroll_csf(csf_cf, csf_vec)
  use globaldata
  implicit none
  real, intent(in)    :: csf_cf
  integer, intent(in) :: csf_vec
 
  integer             :: db(4),mz2(2),d1f(2),d2f(2),del(2)
  integer             :: num,denom,sgn
  integer             :: i,j,k,l,icnt,bt,found
  integer             :: ms2,m2,nalpha,nopen,nloops
  integer             :: bvec(nocc),aloc(nocc),oopen(nocc)
  integer             :: idet(nocc+nextl),det(nocc+nextl),refdet(nocc+nextl)
  double precision    :: cf
  logical             :: zero
  integer             :: ifac
  integer             :: oparity

! DEBUG --------------
!  integer             :: nchk
!  integer             :: chkdet(rlen,1024)
!  double precision    :: chkcf(1024)
  double precision    :: compute_s2
  double precision    :: eps=1.d-8

  db  = (/ 0,  1, -1, 0 /)
  mz2 = (/ 1, -1 /)
  d1f = (/ 1, -1 /)
  d2f = (/-1,  1 /) 
  del = (/ 1,  0 /)

  ndet ` = 0
  ms2   = int(2.*ms)
  bt    = 0
  nopen = 0
  bvec  = 0
 
  refdet = csf_vec

  do j = 1,nocc
   bt = bt + db(csfvec(j,i)+1)
   bvec(j) = bt
   if(csfvec(j,i).eq.1.or.csfvec(j,i).eq.2)then
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
   do k = 1,nocc
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
      call convertDet(det,idet)
      ! include parity in value of cf -- necessary to compute S^2 consistently
      if(cf.lt.0)stop 'ERROR computing determinant'
      cf = sqrt(cf)*sgn*oparity(idet)

!      nchk = nchk + 1
!      chkdet(:,nchk) = idet
!      chkcf(nchk) = cf

      found = 0
      do k = 1,ndet
       if(all(detvec(:,k).eq.idet))then
         found = k
         exit
       endif
      enddo

      if(found.ne.0)then
       det_cf(found)%val = det_cf(found)%val + cf*csf_cf
      else
       ndet = ndet + 1
       if(ndet.gt.numdet) then
        stat = .false.
        stop 'ndet .gt. numdet, increase rdet2csf'
       endif
       det_cf(ndet)%ind = ndet
       det_cf(ndet)%val = cf*csf_cf
       det_vec(:,ndet)  = idet
      endif
    endif !if(.not.zero) 
   enddo !do j = 1,nloops

   !DEBUG------------ compute S2 for this CSF
!   if(abs(compute_s2(nchk,chkcf,chkdet(:,1:nchk))-ms*(ms+1.)).gt.eps) then
!        write(*,*)'vec=',csfvec(:,i)
!        write(*,*)'nchk=',nchk
!        do j = 1,nchk
!         write(*,*)'coef=',chkcf(j)
!         write(*,*)'chkdet=',chkdet(:,j)
!        enddo
!        write(*,*)'cf=',csfcf(i)
!        write(*,*)'ms(ms+1)=',ms*(ms+1.)
!        write(*,*)'i=',i,' S^2',compute_s2(nchk,chkcf,chkdet(:,1:nchk)) 
!        stop 'S^2 value incorrect for CSF -> DET conversion'
!    endif

  enddo !do i = 1,ncsf

  ! DEBUG ---- comptue S2 for entire wfn
  !num = compute_s2(ndet,detcf(1:ndet)%val,detvec(:,1:ndet))
  !write(*,*)' s2 value = ',num

  return
1000 format(20(i3))
 end subroutine makeDetList


!
! write the determinant list to standard output
!
!
 subroutine write_det_list
  use globaldata
  implicit none
  integer            :: i
  double precision   :: dnorm

  if(ndet.eq.0)return

  call sortDetList()

  dnorm = 0.
  i = 1 
  do while(abs(det_cf(i)%val).gt.d_cutoff)
   call print_det(det_cf(i)%val,det_vec(:,det_cf(i)%ind))
   dnorm = dnorm + det_cf(i)%val**2
   i = i + 1
   if(i.gt.ndet)exit
  enddo

  ! not sure how to best handle this.  For the time being, useful to check
  ! wavefunction norm to make sure cutoff is not too severe.
  open(unit=ifile,file='norm.dat',status='replace')
  write(ifile,'(f16.12)')sqrt(dnorm)
  close(ifile)  

  return
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
 subroutine sortDetList
  use globaldata
  implicit none
  integer                         :: i
  type(kvpair),allocatable        :: scr(:)

  allocate(scr(int((ndet+1)/2)))
  call mergesort(detcf(1:ndet),ndet,scr)
  deallocate(scr)

  return
 end subroutine sortDetList


!
! Convert a determinant in COLUMBUS 'step' notation to MULTIGRID notation
!
 subroutine convertDet(indet,outdet)
  use globaldata, only: rlen,nocc
  implicit none
  integer,intent(in)  :: indet(rlen)
  integer,intent(out) :: outdet(rlen)
  integer             :: i,ivec(4)

  ivec = (/0, 1, -1, 2 /)
  do i = 1,nocc
   outdet(i) = ivec(indet(i)+1)
  enddo
  outdet(nocc+1:rlen) = indet(nocc+1:rlen)

  return
 end subroutine convertDet

!
!  compute effect of S2 operator on a single determinant
!
  double precision function compute_s2(ndet,cf,det)
   use globaldata, only: rlen,nocc
   implicit none
   integer,intent(in)           :: ndet
   double precision, intent(in) :: cf(ndet)
   integer, intent(in)          :: det(rlen,ndet)
   integer                      :: i,j,k,l
   integer                      :: trial(rlen)
   integer                      :: par(ndet)
   integer                      :: oparity
   double precision             :: norm

   compute_s2 = 0.
   norm = 0.
 
   do i = 1,ndet
    par(i) = oparity(det(:,i))
   enddo

   do i = 1,ndet
    norm = norm + cf(i)**2
    do j = 1,nocc
     if(abs(det(j,i)).ne.1)cycle
     compute_s2 = compute_s2 + 0.5 * cf(i)**2 ! diagonal contribution of the ladder operators

     do k = 1,nocc
      if(abs(det(k,i)).ne.1) cycle
      compute_s2 = compute_s2 + 0.25 * det(j,i) * det(k,i) * cf(i)**2  ! contribution from Sz^2

      if(det(j,i).ne.det(k,i)) then
       trial = det(:,i)
       trial(j) = -trial(j)
       trial(k) = -trial(k)
       
       do l = 1,ndet
        if(any(det(:,l).ne.trial)) cycle
 
        compute_s2 = compute_s2 + 0.5 * cf(i) * cf(l) * par(i) * par(l)

       enddo ! end l = 1,ndet
      endif ! if(det(i,j)!=det(i,k)
     enddo ! k = 1,norb
    enddo ! j = 1,norb
   enddo ! i = 1,ndet

   return
  end function compute_s2

!
! prints a determinant in multigrid format
!
 subroutine printDet(cf,det)
  use globaldata, only: rlen,nfrzn,nintl,nocc,next,nvrt,norb,caswfn
  implicit none
  double precision,intent(in)       :: cf
  integer,intent(in)        :: det(rlen)
  integer                           :: i,exind,oind
  character*4                       :: lstr
  character*3                       :: istr
  character*17                      :: ofmt
  character*(3*norb)                :: ovec

  if(caswfn)then
   write(lstr,'(i4)')3*(nfrzn+nintl)
  else
   write(lstr,'(i4)')3*norb
  endif
  ofmt = '(f15.10,1x,a'//trim(adjustl(lstr))//')'

  ovec = ''
  do i = 1,nfrzn
   ovec = trim(adjustl(ovec))//'  2'
  enddo

  do i = 1,nintl
   if(abs(det(next+i)).eq.1)then
    write(istr,'(sp,i3)')det(next+i)
   else
    write(istr,'(ss,i3)')det(next+i)
   endif
   ovec = trim(adjustl(ovec))//adjustr(istr)
  enddo

  if(.not.caswfn) then
   exind=1
   do i = 1,nvrt
    oind = nfrzn + nintl + i
    if(oind.eq.det(nocc+exind)) then
     if(abs(det(exind)).eq.1)then
     write(istr,'(sp,i3)')det(exind)
      else
      write(istr,'(ss,i3)')det(exind)
     endif
     ovec = trim(adjustl(ovec))//adjustr(istr)
     if(exind.lt.next)exind = exind + 1
    else
     ovec = trim(adjustl(ovec))//'  0'
    endif
   enddo
  endif

  write(*,ofmt)cf,trim(adjustl(ovec))

  return
 end subroutine printDet


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
   integer, intent(in)    :: na,nb,nc   ! Normal usage: NA+NB = NC
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
  use globaldata, only: rlen,nfrzn,nintl,next,nocc
  implicit none
  integer,intent(in)   :: det(rlen)
  integer                      :: dsort(2*(nocc+nfrzn))
  integer                      :: i,j,na,nb,at,nt
  integer                      :: icnt,nbnd,newbnd

  at = nfrzn + count(det(1:nocc)==1) + count(det(1:nocc)==2)
  nt = 0
  na = 0
  nb = at 

  ! first do frozen orbitals
  do i = 1,nfrzn
    nt = nt + 1;na = na + 1;dsort(nt) = na
    nt = nt + 1;nb = nb + 1;dsort(nt) = nb
  enddo

  ! then do internal orbitals 
  do i = next+1,nocc
   if(any(det(i)==(/1, 2/)))then
    nt = nt + 1;na = na + 1;dsort(nt) = na
   endif
   if(any(det(i)==(/-1, 2/)))then
    nt = nt + 1;nb = nb + 1;dsort(nt) = nb
   endif
  enddo

  ! next do external orbitals 
  do i = 1,next
   if(any(det(i)==(/1, 2/)))then
    nt = nt + 1;na = na + 1;dsort(nt) = na
   endif
   if(any(det(i)==(/-1, 2/)))then
    nt = nt + 1;nb = nb + 1;dsort(nt) = nb
   endif
  enddo

  if(na.ne.at.or.nb.ne.nt)stop 'oparity ERROR!'

  oparity = 1
  nbnd = nt

  do while(nbnd.gt.1)
   newbnd = 0
   do i = 2,nbnd
    if(dsort(i-1).gt.dsort(i)) then
     j = dsort(i-1)
     dsort(i-1) = dsort(i)
     dsort(i) = j
     oparity = -oparity
     newbnd = i
    endif
   enddo
   nbnd = newbnd
  enddo

  return
 end function oparity



