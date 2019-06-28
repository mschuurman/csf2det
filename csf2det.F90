 MODULE csf2det
 
  type kvpair 
   integer          :: ind
   double precision :: val
  end type kvpair

  type t_region
   character*144    :: region_name
   double precision :: wall_last
   double precision :: cpu_last
   double precision :: wall_total
   double precision :: cpu_total
  end type t_region

  ! CSF file
  character*144                             :: csf_file
  ! numerical cutoff for reading CSFs
  double precision                          :: csf_min = 0. 
  ! numerical cutoff for printing determinants
  double precision                          :: det_min = 0.
  ! value of Ms for determinants/csfs
  double precision                          :: m_s

  ! number of csfs to expand
  integer                                   :: n_csf
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
  ! user defined number of MOs, i.e. pad "n_orb" so that it is nmo_tot
  integer                                   :: nmo_pad = 0

  ! record length -- is nintl + 2 * next (i.e. external orbs, plus orb label)
  integer                                   :: rec_len
  ! vector in which to store the csfs
  integer,dimension(:,:),allocatable        :: csf_vec
  ! coefficients for each csf
  double precision,dimension(:),allocatable :: csf_cf

  ! maximum number of determinants, i.e. size of allocated arrays
  integer                                    :: ndet_max = 0
  ! total number of determinants that we've read/generated
  integer                                    :: ndet_all = 0
  ! number of determinants to be printed to wavefunction file
  integer                                    :: n_det = 0
  ! vector to store the determinants
  integer,dimension(:,:),allocatable         :: det_vec
  ! determinant coefficients: a scalar real value, and a sorting index
  ! to make sorting by magnitude simpler/more efficient
  type(kvpair),dimension(:),allocatable      :: det_cf
  ! the associated phase to bring determinant into alpha-string notation
  ! which is required by superdyson
  integer,dimension(:),allocatable           :: det_phase
  ! norm of determinantal wfn
  double precision                           :: det_norm
  ! <S^2>
  double precision                           :: det_s2

  ! number of timers
  integer                                   :: n_timers = 0
  ! maximum number of timers
  integer,parameter                         :: nt_max = 10
  ! list of timers
  type(t_region),dimension(10)              :: t_list

 contains

!******************************************************
! INITIALIZATION ROUTINES
!******************************************************
 subroutine read_cmdline()
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

    elseif (trim(adjustl(abuf)) == '-n_mos') then
       i_arg = i_arg + 1
       call getarg(i_arg, abuf)
       abuf = trim(adjustl(abuf))
       read(abuf,*) nmo_pad

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
 subroutine parse_line(line, cf, step)
  implicit none
  character*144, intent(inout) :: line
  double precision, intent(out):: cf
  integer,intent(out)          :: step(rec_len)
  character*20                 :: str
  integer                      :: ni,ne

  line = adjustl(line)
  str  = line(1:index(line,' ')-1)
  step = 0
  read(str,*)cf

  if (abs(cf).lt.csf_min) then
   cf = 0.
   return
  endif

  ! move start of string to first index past coefficient
  line  = line(index(line,' ')+1:)
  ni      = 0
  ne      = 0
  do while (len(trim(line)).gt.0)
    str = line(1:index(line,' ')-1)
    if(index(str,':').ne.0) then
     ne = ne + 1
     read(str(1:index(str,':')-1),*)step(n_occ+ne) 
     line = adjustl(line(index(line,' ')+1:))
     str = line(1:index(line,' ')-1)
     read(str,*)step(ne)
    else
     ni = ni + 1
     read(str,*)step(n_extl+ni)
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
  implicit none
  integer                :: ios
  integer                :: cfile=99
  integer                :: orb_occ, icsf
  character*144          :: line
  character*20           :: str
  double precision       :: cf
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
  scan_csf_list: do 
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
    scan_csf: do while (len(trim(adjustl(line))).gt.0) 
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
    enddo scan_csf
   endif
  enddo scan_csf_list

  ! initially allocate the determinant array to 5 * CSFs
  n_occ     = n_intl + n_extl
  rec_len   = n_intl + 2*n_extl
  n_orb       = max(n_occ, nmo_pad)
  csf_norm    = 0.
  allocate(csf_vec(rec_len,n_csf))
  allocate(csf_cf(n_csf))

  ! go back to the beginning of the file
  rewind(cfile)

  icsf = 0
  csf_norm = 0.
  read_csf_list: do
   read(cfile, 1000, iostat=ios)line
   if(ios.lt.0)exit
   if(ios.gt.0)stop 'error reading csf_file'

   ! read the CSF
   icsf = icsf + 1
   call parse_line(line, csf_cf(icsf), csf_vec(:,icsf))
   !print *,'csf_vec=',csf_vec

   ! compute the norm of the csf wfn for diagnostic purposes
   csf_norm = csf_norm + csf_cf(icsf)**2

   ! check the external orbital index: update n_orb if necessary
   n_orb = max(n_orb, maxval(csf_vec(n_occ:,icsf)))

   ! if the returned coefficient is zero, we've read all csfs less than cutoff
   if (icsf == n_csf) exit

  enddo read_csf_list

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
 subroutine expand_csfs()
  implicit none
  integer             :: db(4),mz2(2),d1f(2),d2f(2),del(2)
  integer             :: num, denom, sgn, loc_max, det_cnt
  integer             :: icsf,j,k,l,icnt,bt
  integer             :: ms2,m2,nalpha,nopen,nloops
  integer             :: bvec(n_occ),aloc(n_occ),oopen(n_occ)
  integer             :: ref_det(rec_len), step_det(rec_len), multi_det(rec_len)
  double precision    :: cf
  logical             :: zero,found
  !
  integer,dimension(:,:),allocatable        :: vec_loc
  double precision,dimension(:),allocatable :: cf_loc, phase_loc
  

  !  each csf gets 16 slots, fill in determinants as necessary
  !  combine in serial as a second step, no reduction
  ndet_max = 8
  ndet_all = 0
  allocate(det_vec(rec_len,ndet_max*n_csf))
  allocate(det_cf(ndet_max*n_csf))
  allocate(det_phase(ndet_max*n_csf))

  !$omp parallel default(none) &
  !$omp& shared(n_occ,rec_len,n_csf,m_s,ndet_all,csf_vec,csf_cf,det_vec,det_cf,det_phase)  &
  !$omp& private(db,mz2,d1f,d2f,del) &
  !$omp& private(num,denom,sgn,loc_max,det_cnt,icsf,j,k,l,icnt,bt) &
  !$omp& private(ms2,m2,nalpha,nopen,nloops,bvec,aloc,oopen) &
  !$omp& private(ref_det,step_det,multi_det,vec_loc,cf_loc,phase_loc) &
  !$omp& private(cf,zero,found)

  loc_max = 16
  allocate(vec_loc(rec_len,loc_max*n_csf))
  allocate(cf_loc(loc_max*n_csf))
  allocate(phase_loc(loc_max*n_csf))  

  db      = (/ 0,  1, -1, 0 /)
  mz2     = (/ 1, -1 /)
  d1f     = (/ 1, -1 /)
  d2f     = (/-1,  1 /) 
  del     = (/ 1,  0 /)
  ms2     = int(2.*m_s)
  det_cnt = 0
  
  !$omp do
  loop_csf_array: do icsf = 1,n_csf

   bt      = 0
   nopen   = 0
   oopen   = 0
   bvec    = 0
   aloc    = 0
   ref_det = csf_vec(:,icsf)

   do j = 1,n_occ
    bt = bt + db(csf_vec(j,icsf)+1)
     bvec(j) = bt
    if(csf_vec(j,icsf).eq.1.or.csf_vec(j,icsf).eq.2)then
     nopen = nopen + 1
     oopen(nopen) = j
     ref_det(j) = 2
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
    step_det = ref_det
 
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
      step_det(oopen(aloc(k)))=1
     enddo
    endif
 
    ! determine the projection of the determinant on the CSF
    cf   = 1.
    sgn  = 1
    m2   = 0
    zero = .false.
    do k = 1,n_occ
     select case(csf_vec(k,icsf))
  
       case(1)
         m2 = m2 + mz2(step_det(k))
         num = bvec(k) + d1f(step_det(k)) * m2
         if(num == 0)then
          zero = .true.
          exit
         endif
         denom = 2 * bvec(k) 
         cf = cf * num / denom
 
       case(2)
         m2 = m2 + mz2(step_det(k))
         num = bvec(k) + 2 + d2f(step_det(k)) * m2
         if(num == 0)then
          zero = .true.
          exit
         endif
         denom = 2 * (bvec(k) + 2)
         sgn = sgn * (-1)**(bvec(k)+del(step_det(k)))
         cf = cf * num / denom
 
       case(3)
         sgn = sgn * (-1)**(bvec(k))
 
       case default
      end select
     enddo ! do k =1 ,nocc 
 
     if(.not.zero) then

       ! convert the determinant to multigrid format
        call convert_det(step_det,multi_det)

       ! include parity in value of cf -- necessary to compute S^2 consistently
       if(cf.lt.0)stop 'ERROR computing determinant'
       cf = sqrt(cf)*sgn*csf_cf(icsf)

       ! scan local list to see if vector already stored
       found = .false.
       scan_detlist: do k = 1,det_cnt
        if (any(multi_det.ne.vec_loc(:,k))) cycle scan_detlist
        found = .true.
        exit
       enddo scan_detlist

       if(found) then
        cf_loc(k)          = cf_loc(k) + cf*phase_loc(k)
       else
        det_cnt            = det_cnt + 1
        vec_loc(:,det_cnt) = multi_det
        phase_loc(det_cnt) = oparity(multi_det)
        cf_loc(det_cnt)    = cf*phase_loc(det_cnt)
       endif

     endif !if(.not.zero) 
    enddo !do j = 1,nloops
   enddo loop_csf_array
   !$omp end do

   ! now collapse all determiants in to a unqiue list
   !$omp critical 
   do j = 1,det_cnt
    multi_det = vec_loc(:,j)
 
    ! scan master list to see if vector already stored
    found = .false.
    scan_loclist: do k = 1,ndet_all
     if (any(multi_det.ne.det_vec(:,k))) cycle scan_loclist
     found = .true.
     exit
    enddo scan_loclist

    if(found) then
     det_cf(k)%val = det_cf(k)%val + cf_loc(j)
    else
     ndet_all = ndet_all + 1
     det_vec(:,ndet_all)  = multi_det
     det_cf(ndet_all)%ind = ndet_all
     det_cf(ndet_all)%val = cf_loc(j)
     det_phase(ndet_all)  = phase_loc(j)
    endif
   enddo
   !$omp end critical
   !$omp end parallel

1000 format(20(i3))
1001 format(f18.12)
1002 format(40(i3))
 end subroutine expand_csfs 

!
! write the determinant list to standard output
!
!
 subroutine write_det_list
  implicit none

  if(ndet_all.eq.0)return

  ! sort so determinants printed largest coefficient to smallest
  call sort_det_list()

  ! print determinants
  n_det = 0 
  do while(n_det.lt.ndet_all.and.abs(det_cf(n_det+1)%val).gt.det_min)
   n_det = n_det + 1
   call print_det(det_cf(n_det)%val,det_vec(:,det_cf(n_det)%ind))
  enddo

  return
 end subroutine write_det_list

!
! print diagnostics of the determinatn list
!
 subroutine run_diagnostics()
  implicit none

  ! compute S^2 and norm -- should only include determinants that are printed, hence 'n_print'
  call compute_norm_s2()

 end subroutine run_diagnostics 

!
!
! 
 subroutine print_summary()
  implicit none
  integer                      :: ofile=99

  ! not sure how to best handle this.  For the time being, useful to check
  ! wavefunction norm and S^2 to make sure cutoff is not too severe.
  open(unit=ofile,file='csf2det.out',status='replace')
  write(ofile,"(' Diagnostic Report')")
  write(ofile,"(' -----------------')")
  write(ofile,"(' Norm of wfn in CSF basis: ',f15.10)")csf_norm
  write(ofile,"(' Norm of wfn in det basis: ',f15.10)")det_norm
  write(ofile,"(' <S^2>:                    ',f15.4)")det_s2

  call print_timer(ofile)
  close(ofile)

  return
 end subroutine print_summary


!
! deallocate dynamic memory arrays
!
 subroutine cleanup()
  implicit none

  deallocate(csf_vec)
  deallocate(csf_cf)
  deallocate(det_vec)
  deallocate(det_cf)
  deallocate(det_phase)

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
  implicit none
  integer                         :: i
  type(kvpair),allocatable        :: scr(:)

  allocate(scr(int((ndet_all+1)/2)))
  call mergesort(det_cf(1:ndet_all),ndet_all,scr)
  deallocate(scr)

  return
 end subroutine sort_det_list


!
! Convert a determinant in COLUMBUS 'step' notation to MULTIGRID notation
!
 subroutine convert_det(indet,outdet)
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
  subroutine compute_norm_s2()
   implicit none
   integer                      :: ibra, braloc, iket, ketloc
   integer                      :: iorb, jorb
   integer                      :: bra(rec_len), ket(rec_len)
   integer                      :: trial(rec_len)
   double precision             :: bra_cf, ket_cf
   double precision             :: s2_loc, norm_loc

   !$omp parallel default(none) &
   !$omp& shared(n_det,n_occ,det_vec,det_cf,det_phase)  &
   !$omp& private(iket,ibra,iorb,jorb,braloc,ketloc,trial,bra,ket,bra_cf,ket_cf) &
   !$omp& reduction(+:s2_loc,norm_loc)
   s2_loc   = 0.
   norm_loc = 0.

   !$omp do
   loop_ket: do iket = 1,n_det
    ketloc = det_cf(iket)%ind
    ket    = det_vec(:,ketloc)
    ket_cf = det_cf(iket)%val
    norm_loc = norm_loc + ket_cf**2

    scan_orbi: do iorb = 1,n_occ
     if(abs(ket(iorb)).ne.1)cycle scan_orbi
     s2_loc = s2_loc + 0.5 * ket_cf**2 ! diagonal contribution of the ladder operators

     scan_orbj: do jorb = 1,n_occ
      if(abs(ket(jorb)).ne.1) cycle scan_orbj
      s2_loc = s2_loc + 0.25 * ket(iorb) * ket(jorb) * ket_cf**2 ! contribution from Sz^2

      if(ket(iorb).ne.ket(jorb)) then
       trial       =  ket
       trial(iorb) = -trial(iorb)
       trial(jorb) = -trial(jorb)
       
       scan_bra: do ibra = 1,n_det
        braloc = det_cf(ibra)%ind
        bra    = det_vec(:,braloc)
        bra_cf = det_cf(ibra)%val 

        ! we're looking for trial in the det list..
        if(any(bra.ne.trial)) cycle scan_bra
        s2_loc = s2_loc + 0.5 * ket_cf * bra_cf * det_phase(ketloc) * det_phase(braloc)
        exit ! each determinant in the list is unique, exit loop if we find right one.

       enddo scan_bra ! end l = 1,ndet
      endif ! if(det(i,j)!=det(i,k)
     enddo scan_orbj ! k = 1,norb
    enddo scan_orbi! j = 1,norb
   enddo loop_ket ! i = 1,ndet
  !$omp end do
  !$omp end parallel

   det_norm = sqrt(norm_loc)
   det_s2   = s2_loc

   return
  end subroutine compute_norm_s2

!
! prints a determinant in multigrid format
!
 subroutine print_det(cf,det)
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
   integer, intent(in)         :: na,nb,nc   ! Normal usage: NA+NB = NC
   type(kvpair), intent(inout) :: a(na)      ! B overlays C(NA+1:NC)
   type(kvpair), intent(in)    :: b(nb)
   type(kvpair), intent(inout) :: c(nc)
   integer :: i,j,k

   i = 1
   j = 1
   k = 1
   do while(i.le.na.and.j.le.nb)
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
   integer, intent(in)                :: n
   type(kvpair), intent(inout)        :: a(n)
   type(kvpair), intent(out)          :: t(int((n+1)/2))
   integer                            :: na,nb
   type(kvpair)                       :: dval

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
     if(dsort(i-1) < dsort(i)) cycle
     j          = dsort(i-1)
     dsort(i-1) = dsort(i)
     dsort(i)   = j
     oparity    = -oparity
     sorted     = .false. 
   enddo
  enddo

  return
 end function oparity

!
!
!
 subroutine timer_start(reg_name)
  implicit none
  character(len=*),intent(in) :: reg_name
  character*144               :: label
  double precision            :: tcpu, twall
  integer                     :: cnt, cnt_rate, cnt_max
  integer                     :: it
! integer                     :: timer_pos

  label = trim(adjustl(reg_name))

  it = timer_pos(label)
  call cpu_time(tcpu)
  call system_clock(cnt,cnt_rate,cnt_max)

  twall                = dble(1.*cnt/cnt_rate)
  t_list(it)%cpu_last  = tcpu
  t_list(it)%wall_last = twall

 end subroutine timer_start

!
!
!
 subroutine timer_stop(reg_name)
  implicit none
  character(len=*),intent(in) :: reg_name
  character*144               :: label
  double precision            :: tcpu, twall
  integer                     :: cnt, cnt_rate, cnt_max
  integer                     :: it
!  integer                     :: timer_pos

  label = trim(adjustl(reg_name))

  it = timer_pos(label)

  call cpu_time(tcpu)
  call system_clock(cnt,cnt_rate,cnt_max)

  twall                    = dble(1.*cnt/cnt_rate)
  t_list(it)%cpu_total     = tcpu - t_list(it)%cpu_last
  t_list(it)%wall_total    = twall - t_list(it)%wall_last
  t_list(it)%cpu_last      = tcpu
  t_list(it)%wall_last     = twall
 end subroutine timer_stop

!
!
!
 subroutine print_timer(ofile)
  implicit none
  integer,intent(in)   :: ofile
  integer              :: i
  double precision     :: total_cpu, total_wall
  
  total_cpu = 0.
  total_wall = 0.

  write(ofile,"('')")
  write(ofile,"(' Timing Report')")
  write(ofile,"(' -------------')")
  write(ofile,"(50x,'  CPU time',4x,' Real time (seconds)')")
  do i = 1,n_timers
    write(ofile,"(1x,a35,f24.2,f24.2)")adjustl(t_list(i)%region_name), t_list(i)%cpu_total, t_list(i)%wall_total
    total_cpu  = total_cpu  + t_list(i)%cpu_total
    total_wall = total_wall + t_list(i)%wall_total
  enddo
  write(ofile,"(' -----------------------------------------------------------------------------------')")
  write(ofile,"(1x,a5,30x,f24.2,f24.2)")'Total',total_cpu,total_wall

  return
 end subroutine print_timer
 
!
!
!
 function timer_pos(label)
  implicit none
  character*144       :: label
  integer             :: timer_pos
  logical             :: found
 
  found      = .false.
  timer_pos = 1
  do while (.not.found.and.timer_pos.le.n_timers)
    if(t_list(timer_pos)%region_name == adjustl(label)) then
      found = .true.
      exit
    endif
    timer_pos = timer_pos + 1
  enddo

  if(.not.found) then
   n_timers = n_timers + 1
   t_list(timer_pos)%region_name = adjustl(label)
  endif

 end function timer_pos

end module csf2det

!
! Program: CSF2DET 
!          -- convert a CI expansion in CSFs to one in determinants
!
!
 PROGRAM csf2det_driver
  use csf2det
  implicit none

  call read_cmdline()

  call timer_start('parse_csf_list')
  call parse_csf_list()
  call timer_stop('parse_csf_list')

  call timer_start('expand_csfs')
  call expand_csfs()
  call timer_stop('expand_csfs')

  call timer_start('write_det_list')
  call write_det_list()
  call timer_stop('write_det_list')

  call timer_start('run_diagnostics')
  call run_diagnostics()
  call timer_stop('run_diagnostics')

  call print_summary()
  call cleanup()

 end PROGRAM csf2det_driver

