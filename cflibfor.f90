  !=======================================================================
  !           CORRELATION FUNCION LIBRARY FOR PYTHON WRAPPING
  !=======================================================================
  ! These are Fortran that implement fast counting of pairs in correlation
  ! function calculations and intended to be called from python wrappers
  ! but can calso be called from standard Fortran. Functions are availabe to
  !
  !   * compute (auto)counts in 3D/projected/angular space
  !   * compute (cross)counts in 3D/projected/angular space
  !   * build and count bootstrap samples to derive bootstrap errors
  !   * build sk-ll tables to speed up counting
  !
  ! PYTHON INTERFACE
  !  ---------------------------------------------------------------------
  !  To compile this library and generate cflibfor.so, use
  !     > f2py -c cflibfor.pyf cosmolib.f90 cflibfor.f90
  !  
  !  Note 1 : cflibfor.pyf holds the in/out declaration of variables necessary 
  !           for using this functions from python and only those functions with
  !           proper interfase declarations are exposed in python.
  !  Note 2 : Remember to update these declarations if you intend to add/remove 
  !           parameters to any function exposed or if you want to 
  !           expose another one. Just copy and paste the examples, but if 
  !           needed you can put subroutine() body in a .f90 file and run
  !              > f2py -h subroutine.pyf -m cflibfor subroutine.f90
  !           This will generate a basic .pyf interface that you can fine tune
  !           and paste back in cflibfor.pyf 
  !  Note 3 : There are some optimization options you can try
  !     > f2py -c --opt="-O3 -funroll-loops -march=corei7-avx -ftree-vectorize" [...]
  !     > f2py -c --verbose --opt='-O3 -march=corei7-avx -ftree-vectorize -ftree-vectorizer-verbose=3 -fopt-info-vec-all -fdump-tree-vect' cflibfor.pyf cosmolib.f90 cflibfor.f90
  ! 
  ! Contact Author
  !===============
  ! Emilio Donoso (edonoso@conicet.gov.ar)

  
module mod
use omp_lib
real(kind=8) :: deg2rad = acos(-1.d0)/180.d0
real(kind=8) :: rad2deg = 180.d0/acos(-1.d0)

contains


!==========================================================================================
!   AUXILIARY ROUTINES
!==========================================================================================

subroutine init_random_seed()
  !=======================================================================
  ! NAME
  !  init_random_seed()
  !
  ! PURPOSE
  !  Seed the pseudo-RN generators
  !
  ! DESCRIPTION
  !  Use if you need truly diffrent RN sequences in multiple processors
  
  use iso_fortran_env, only: int64
  implicit none
  integer, allocatable :: seed(:)
  integer :: i, n, un, istat, dt(8), pid
  integer(int64) :: t

  call random_seed(size = n)
  allocate(seed(n))
  ! First try if the OS provides a random number generator
  open(newunit=un, file="/dev/urandom", access="stream", &
       form="unformatted", action="read", status="old", iostat=istat)
  if (istat == 0) then
     read(un) seed
     close(un)
  else
     ! Fallback to XOR:ing the current time and pid. The PID is
     ! useful in case one launches multiple instances of the same
     ! program in parallel.
     call system_clock(t)
     if (t == 0) then
        call date_and_time(values=dt)
        t = (dt(1) - 1970) * 365_int64 * 24 * 60 * 60 * 1000 &
             + dt(2) * 31_int64 * 24 * 60 * 60 * 1000 &
             + dt(3) * 24_int64 * 60 * 60 * 1000 &
             + dt(5) * 60 * 60 * 1000 &
             + dt(6) * 60 * 1000 + dt(7) * 1000 &
             + dt(8)
     end if
     pid = getpid()
     t = ieor(t, int(pid, kind(t)))
     do i = 1, n
        seed(i) = lcg(t)
     end do
  end if
  call random_seed(put=seed)
contains
  ! This simple PRNG might not be good enough for real work, but is
  ! sufficient for seeding a better PRNG.
  function lcg(s)
    integer :: lcg
    integer(int64) :: s
    if (s == 0) then
       s = 104729
    else
       s = modulo(s, 4294967296_int64)
    end if
    s = modulo(s * 279470273_int64, 4294967291_int64)
    lcg = int(modulo(s, int(huge(0), int64)), kind(0))
  end function lcg
end subroutine init_random_seed

subroutine bootstrap(npt,nboot,fixseed,w)
  !=======================================================================
  ! NAME
  !  bootstrap()
  !
  ! PURPOSE
  !  Create bootstrap resampling weights
  !
  ! DESCRIPTION
  !  Given npt particles and nboot desired bootstrap samples
  !  within the area specified by (ral,rau,decl,decu,dcoml,dcomu).
  !
  ! INPUTS
  !  Variable----Type------------Description------------------------------
  !  npt         int             Number of particles
  !  nboot       int             Number of bootstrap samples to return
  !  fixseed     int             If 0, then the PRNG is not initialized
  !                              so every call will give different randoms
  !                              If >0, then the PRNG is seeded with the single
  !                              value of fixseed, so every call will give same 
  !                              randoms.
  !                              ALWAYS use fixseed>0 if you want the SAME
  !                              bootstrap samples on each paralell process.
  !                              Otherwise, every bootstrap() call running on
  !                              multiple processes get different PRN sequences, 
  !                              so w will be different in each spawned process
  !
  ! OUTPUTS
  !  Variable----Type----------------Description--------------------------
  !  w           real*4(npt,nboot)   Bootstrap resampling weights

implicit none
integer      :: npt,nboot,fixseed,n,i,nb,idx
real(kind=4) :: sd
real(kind=4) :: w(nboot,npt)
integer,allocatable :: seed(:)

if (fixseed>0) then
   call random_seed(SIZE=n)
   !print *,n
   allocate(seed(n))
   seed = fixseed
   call random_seed(PUT=seed)
   !  call init_random_seed()
endif

w = 0.0

do nb=1,nboot
   do i=1,npt
      call random_number(sd)
      idx = floor(sd*(npt-1)) + 1
      w(nb,idx) = w(nb,idx) + 1.
   enddo
enddo

end subroutine bootstrap

real(kind=8) function sthmax2(rpmax,rg,rkl)
  ! Auxiliary routine used by counting routines
  implicit none
  real(kind=8)  :: rpmax,rg
  real(kind=8)  :: rkl
  if(rkl==0)then
     sthmax2 = 1.d0
  else
     sthmax2 = rpmax/sqrt(2.0d0*rg*rkl)
     if(sthmax2>1.0d0) sthmax2 = 1.0d0
  endif
end function sthmax2

real(kind=8) function dalp(sthmax2,decg,decr,i,hc1,decl,isele)
  ! Auxiliary routine used by counting routines
  implicit none
  real(kind=8) :: decg,hc1,decl,decr,sthmax2
  real(kind=8) :: cdeccl,cdeccu,cdeccmin,atmp,btmp,sdalp2,sdalp22
  integer      :: i,isele
  
  cdeccl   = cos((decl+(i-1)*hc1)*deg2rad)
  cdeccu   = cos((decl+i*hc1)*deg2rad)
  cdeccmin = min(cdeccu,cdeccl)
  
  if(cdeccmin==0) then
     dalp = 180.0d0
     return
  endif

  if(isele==1) then
     sdalp2 = sthmax2/sqrt(cos(decr)*cdeccmin)
  else
     atmp = (sin((decg-(decl+(i-1)*hc1))*0.5*deg2rad))**2
     btmp = (sin((decg-(decl+i*hc1))*0.5*deg2rad))**2
     atmp = min(atmp,btmp)
     sdalp22 = (sthmax2**2-atmp)/(cos(decr)*cdeccmin)
     if(sdalp22<0) sdalp22=0.
     sdalp2 = sqrt(sdalp22)
  endif

  if(sdalp2>1)then
     dalp = 180.0d0
  else
     !dalp = asind(sdalp2)*2.
     dalp = (asin(sdalp2)*2.0d0)*rad2deg
  endif
end function dalp

real(kind=4) function wfiber(shth2)
  ! Apply genereic SDSS fiber correction to a pair separated by shth2
  ! Note this is based on quick and old fits to the angular correlation of a
  ! spectroscopic sample and that of its parent sample.
  !
  ! It should provide decent results down to scales of ~0.05 Mpc/h, but it is
  ! provided just as an example of the technique. Check Guo et al. 2012 for 
  ! more sophisticated approaches.
  implicit none
  real(kind=8) :: shth2,theta

  theta = sqrt(shth2)*2.0d0*3600.0d0*rad2deg

  wfiber = 1.0
  if(theta<55.0d0.and.theta>0.001d0) wfiber=(1.0d0+30.973409d0*theta**(-0.79937208d0))/ &
                      (4.7734205d0*theta**(-0.47210907d0))
end function wfiber

!==========================================================================================
!   COUNTING ROUTINES - 3D SPHERICAL (PROJECTED)
!==========================================================================================

subroutine skll3d(mxh1,mxh2,mxh3,npt,ra,dec,dc,sbound,sepv,nsepv,sk,ll)
!===============================================================================
! NAME
!  skll3d()
!
! PURPOSE
!  Construct a skip table (SK) and linked list (LL) for a set of particles with 
!  (ra,dec,dcom) coordinates. 
!  The SK table of size (mxh3,mxh2,mxh1) is constructed within the area 
!  specified by sbound=(ramin,ramax,decmin,decmax,dcmin,dcmax).
!
! INPUTS
!  Variable----Type------------Description--------------------------------------
!  mxh1        int             Nr. of DEC cells in SK table
!  mxh2        int             Nr. of RA cells in SK table
!  mxh3        int             Nr. of DCOM cells in SK table
!  npt         int             Number of particles
!  ra          real8(npt)      RA of particles
!  dec         real8(npt)      DEC of particles
!  dc          real8(npt)      Comoving distance of particles
!  sbound      real8(6)        Survey boundaries in RA,DEC,DCOM. Form is 
!                              (ramin,ramax,decmin,decmax,dcmin,dcmax)
!  nseps       int             Number of 3D separation bins
!  seps        real4(nseps+1)  Bins in 3D separation
!  nsepv       int             Number of radial separation bins
!  sepv        real4(nsepv+1)  Bins in radial separation
!
! OUTPUTS
!  Variable---Type-------------------Description--------------------------------
!  sk          int(mxh3,mxh2,mxh1)   Skip table
!  ll          int(np)               Linked list
!
! NOTES  -----------------------------------------------------------------------
!  1. Angles are in degrees
!  2. The boundaries of RA should be 0 - 360

implicit none
integer      :: mxh1,mxh2,mxh3,npt,nsepv
real(kind=8) :: ra(npt),dec(npt),sbound(6),ral,rau,decl,decu,dcoml,dcomu,hc1,hc2,hc3
real(kind=8) :: dc(npt),sepv(nsepv+1),rvmax  !kind=4
integer      :: sk(mxh3,mxh2,mxh1),ll(npt),q(3)
integer      :: nc3,i,i0

sk = 0
!-----------------------------------
! Unpack boundaries
ral = sbound(1) ; rau = sbound(2) ; decl = sbound(3) ; decu = sbound(4)
dcoml = sbound(5) ; dcomu = sbound(6)

!------------------------------------
! Get nc3
rvmax = sepv(nsepv+1)
nc3   = int((dcomu-dcoml)/rvmax)
if(nc3>mxh3) nc3=mxh3

!-------------------------------------
! Get cell dimensions in (dec,ra,dcom)
hc1 = (decu-decl)/float(mxh1)
hc2 = (rau-ral)/float(mxh2)
hc3 = (dcomu-dcoml)/float(nc3)

!------------------------------
! Build SK and LL
do i=1, npt
   ! Find cell coordinates of the i-th particle
   q(1) = int((dec(i)-decl)/hc1)+1
   q(2) = int((ra(i)-ral)/hc2)+1
   q(3) = int((dc(i)-dcoml)/hc3)+1
 
   ! Check coordinates are within range
   if(q(1)>mxh1.or.q(1)<1) cycle
   if(q(2)>mxh2) then
      q(2) = q(2)-mxh2
   else if(q(2)<1) then
      q(2) = q(2)+mxh2
   end if
   if(q(3)>nc3.or.q(3)<1) cycle

   ! Assign LL and SK entry
   ll(i) = sk(q(3),q(2),q(1))
   sk(q(3),q(2),q(1)) = i

   ! Sort in increasing dc(i) so that the sk element has the smallest dc
   if(ll(i)==0) cycle
   if(dc(i)<=dc(ll(i))) cycle
   sk(q(3),q(2),q(1)) = ll(i)
   i0 = ll(i)
   11 if(ll(i0)/=0)then
      if(dc(i)<=dc(ll(i0)))then
         ll(i)  = ll(i0)
         ll(i0) = i
         cycle
      else
         i0 = ll(i0)
         goto 11
      endif
   else
      ll(i0) = i
      ll(i)  = 0
   endif
end do

end subroutine skll3d


subroutine rppi_A(nt,npt,dec,dc,x,y,z,nsepp,sepp,nsepv,sepv,sbound, &
           mxh1,mxh2,mxh3,cntid,logf,sk,ll,aapv)
!===============================================================================
! NAME
!  rppi_A()
!
! DESCRIPTION
!  Count data pairs in projected space for ONE sample of particles
!
! INPUTS
!  Variable---Type-----------Description----------------------------------------
!  nt         int            Number of thread to use (OpenMP)
!  npt        int            Number of particles
!  dec        r8(npt)        DEC of particles [deg]
!  dc         r8(npt)        DCOM of particles [Mpc/h]
!  x,y,z      r8(npt)        X,Y,Z coordinates of particles (see radec2xyz())
!  nsepp      int            Number of projected separation bins
!  sepp       r8(nsepp+1)    Bins in projected separation [Mpc/h]
!  nsepv      int            Number of radial separation bins
!  sepv       r8(nsepv+1)    Bins in radial separation [Mpc/h]
!  sbound     r8(6)          Survey boundaries in RA,DEC,DCOM
!                            Form is (ramin,ramax,decmin,decmax,dcmin,dcmax)
!  mxh1       int            Nr of DEC cells of skip table
!  mxh2       int            Nr of RA cells of skip table
!  mxh3       int            Nr of DCOM cells of skip table
!  cntid      char*2         Two-character string to identify the samples,
!                            e.g. DD, RR, DR,... (used only for ouput messages)
!  logf       char*80        String for the file to log Fortran status messages
!  sk         int(mxh3,mxh2,mxh1)  Skip table (SK) contructed with skll2d()
!  ll         int(npt)       Linked list (LL) contructed with skll2d()
!
! OUTPUTS
!  Variable---Type-------------Description--------------------------------------
!  aapv       r8(nsepv,nsepp)  Counts in radial and projected separation bins
!
! NOTES  -----------------------------------------------------------------------
!  1. RA limits should be set to ramin=0 and ramax=360.
!  2. Remember to update the declarations in cflibfor.pyf if you add/remove 
!     in/out parameters to this exposed subroutine.

implicit none
real(kind=8)  :: dec(npt),sbound(6),ral,rau,decl,decu,dcoml,dcomu,hc1,hc2,hc3,rp2
real(kind=8)  :: x(npt),y(npt),z(npt),dc(npt),rcl,stm2,dltdec,dltra,dci,xi,yi,zi
real(kind=8)  :: sepp(nsepp+1),sepv(nsepv+1),sepp2(nsepp+1),rpmax,rpmax2,rvmax,rv,idsepv
real(kind=8)  :: aapv(nsepv,nsepp)
integer       :: sk(mxh3,mxh2,mxh1),ll(npt),mxh1,mxh2,mxh3,npt,nsepp,nsepv
integer       :: nt,nthr,nc1,nc2,nc3,jq1m,jq2m
integer       :: i,ii,j,jj,iq1,iq2,iq3,jq1,jq2,jq3,jq1min,jq2min,jq2max,jq2t
character     :: cntid*2,logf*80

!------Open log file---------------------------------
open(11,file=logf,action='write',access='append')

! Some useful print statements. Place where desired
!print *, mxh1,mxh2,mxh3
!print *, nc1,nc2,nc3
!print *, hc1, hc2, hc3
!print *, omp_get_num_threads(), omp_get_max_threads(), omp_get_num_procs()

!----------------------------------------------------
! Reset the counts, set (rpmax,rvmax) and square rp bins
aapv   = 0.0d0
rpmax  = sepp(nsepp+1)
rvmax  = sepv(nsepv+1)
sepp2  = sepp*sepp
rpmax2 = sepp2(nsepp+1)
idsepv = 1./(sepv(2)-sepv(1))  !inverse radial bin size

!----------------------------------------------------
! Unpack boundaries and get cell dimensions in (dec,ra,dcom)
ral = sbound(1) ; rau = sbound(2) ; decl = sbound(3) ; decu = sbound(4)
dcoml = sbound(5) ; dcomu = sbound(6)
nc1 = mxh1
nc2 = mxh2
nc3 = int((dcomu-dcoml)/rvmax)
if(nc3>mxh3) nc3 = mxh3
hc1 = (decu-decl)/float(nc1)     !effective nr of DEC cells
hc2 = (rau-ral)/float(nc2)       !effective nr of RA cells
hc3 = (dcomu-dcoml)/float(nc3)   !effective nr of DCOM cells

!----------------------------------------------------
! Set the number of threads
if(nt<=0) then
   nthr = omp_get_num_procs()
else
   nthr = nt
endif
call omp_set_num_threads(nthr)

!----------------------------------------------------
! Count pairs in SK grid
write(*,*) ' '
write(*,fmt='(a,i3,a)') '====  Counting '//cntid//' pairs in ', mxh1, ' DEC strips  ===='
!$omp parallel do reduction(+:aapv) default(shared) &
!$omp& private(iq1,iq2,iq3,jq1,jq2,jq3,rcl,stm2,dltdec,jq1m,jq1min) &
!$omp& private(dltra,jq2m,jq2min,jq2max,jq2t,j,i,xi,yi,zi,dci,rv,rp2,ii,jj) &
!$omp& schedule(guided) if(nthr>1)
do iq1=1,nc1
   !$omp critical
   write(*,fmt="(i4)",advance='no') iq1                             ! for screen
   write(11,*) cntid//' counting in DEC strip > ',iq1,'/',mxh1      ! for disk
   flush(11)
   !$omp end critical
   do iq2=1,nc2
      do iq3=1,nc3
         i = sk(iq3,iq2,iq1) ! index of ith particle 
         do while(i/=0) !---- Loop over ith particles ----
            lp_jq3: do jq3=iq3,iq3+1
               if(jq3>nc3) cycle lp_jq3
               rcl    = (jq3-1)*hc3+dcoml
               stm2   = sthmax2(rpmax,dc(i),rcl)
               dltdec = 2.0d0*asin(stm2)*rad2deg
               jq1m   = int(dltdec/hc1)+1 !dltdec/hc1+1.
               if(jq3==iq3) then
                  jq1min = iq1
               else
                  jq1min = iq1-jq1m
               end if
               lp_jq1: do jq1=jq1min,iq1+jq1m
                  if(jq1>nc1.or.jq1<1) cycle lp_jq1
                  if(jq1==iq1) then
                     dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,1)
                  else
                     dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,0)
                  end if
                  jq2m   = int(dltra/hc2)+1  !dltra/hc2+1.
                  jq2max = iq2+jq2m
                  jq2min = iq2-jq2m
                  if(jq2max-jq2min+1>nc2) jq2max = jq2min-1+nc2
                  lp_jq2: do jq2=jq2min,jq2max
                     if(jq2>nc2) then
                        jq2t = jq2-nc2
                     else if(jq2<1) then
                        jq2t = jq2+nc2
                     else
                        jq2t = jq2
                     end if
                     if(jq3==iq3.and.jq2t<iq2.and.jq1==iq1) cycle lp_jq2
                     if(jq3==iq3.and.jq2t==iq2.and.jq1==iq1) then
                        j = ll(i)
                     else
                        j = sk(jq3,jq2t,jq1)
                     endif
                     xi  = x(i)
                     yi  = y(i)
                     zi  = z(i)
                     dci = dc(i)

                     do while(j/=0)
                        rv = abs(dci-dc(j))
                        if(rv<=rvmax) then
                           rp2 = 4.*dci*dc(j)*((xi-x(j))**2 + (yi-y(j))**2 + (zi-z(j))**2)
                           if(rp2<=rpmax2) then
                               jj = 1 + int(rv*idsepv)       !find rv bin number
                               if(rp2>sepp2(nsepp)) then 
                                  aapv(jj,nsepp) = aapv(jj,nsepp) + 1.0d0
                                  goto 70
                               endif
                               do ii=nsepp-1,1,-1
                                  if(rp2>sepp2(ii)) then
                                     aapv(jj,ii) = aapv(jj,ii) + 1.0d0
                                     goto 70
                                  endif
                               enddo
                           endif !rp2<rpmax2
                        else
                           if(jq3>iq3) cycle lp_jq2  !if(rv>rvmax.and.jq3>iq3)
                        endif
                        70 j = ll(j)
                     end do
                     
                  end do lp_jq2
               end do lp_jq1
            end do lp_jq3
            i = ll(i)
         end do !---- End loop over ith particles ----
      end do !loop iq3
   end do !loop iq2
end do
!$omp end parallel do
close(11)  ! close log
write(*,*) ' '
end subroutine rppi_A


subroutine rppi_A_wg(nt,npt,dec,dc,wei,x,y,z,nsepp,sepp,nsepv,sepv,sbound, &
           mxh1,mxh2,mxh3,wfib,cntid,logf,sk,ll,aapv)
!===============================================================================
! NAME
!  rppi_A_wg()
!
! DESCRIPTION
!  Count weighted data pairs in projected space for ONE sample of particles and
!  optionally applies fiber-collision corrections
!
! INPUTS
!  Variable---Type-----------Description----------------------------------------
!  nt         int            Number of thread to use (OpenMP)
!  npt        int            Number of particles
!  dec        r8(npt)        DEC of particles [deg]
!  dc         r8(npt)        DCOM of particles [Mpc/h]
!  wei        r4(npt)        WEIGHT of particles
!  x,y,z      r8(npt)        X,Y,Z coordinates of particles (see radec2xyz())
!  nsepp      int            Number of projected separation bins
!  sepp       r8(nsepp+1)    Bins in projected separation [Mpc/h]
!  nsepv      int            Number of radial separation bins
!  sepv       r8(nsepv+1)    Bins in radial separation [Mpc/h]
!  sbound     r8(6)          Survey boundaries in RA,DEC,DCOM
!                            Form is (ramin,ramax,decmin,decmax,dcmin,dcmax)
!  mxh1       int            Nr of DEC cells of skip table
!  mxh2       int            Nr of RA cells of skip table
!  mxh3       int            Nr of DCOM cells of skip table
!  wfib       int            If 1, also apply fiber-correction weights
!  cntid      char*2         Two-character string to identify the samples,
!                            e.g. DD, RR, DR,... (used only for ouput messages)
!  logf       char*80        String for the file to log Fortran status messages
!  sk         int(mxh3,mxh2,mxh1)  Skip table (SK) contructed with skll2d()
!  ll         int(npt)       Linked list (LL) contructed with skll2d()
!
! OUTPUTS
!  Variable---Type-------------Description--------------------------------------
!  aapv       r8(nsepv,nsepp)  Counts in radial and projected separation bins
!
! NOTES  -----------------------------------------------------------------------
!  1. RA limits should be set to ramin=0 and ramax=360.
!  2. Remember to update the declarations in cflibfor.pyf if you add/remove 
!     in/out parameters to this exposed subroutine.

implicit none
real(kind=8)  :: dec(npt),sbound(6),ral,rau,decl,decu,dcoml,dcomu,hc1,hc2,hc3,rp2
real(kind=4)  :: wei(npt),wpp,wi
real(kind=8)  :: x(npt),y(npt),z(npt),dc(npt),rcl,stm2,dltdec,dltra,dci,xi,yi,zi
real(kind=8)  :: sepp(nsepp+1),sepv(nsepv+1),sepp2(nsepp+1),rpmax,rpmax2,rvmax,rv,shth2,idsepv
real(kind=8)  :: aapv(nsepv,nsepp)
integer       :: sk(mxh3,mxh2,mxh1),ll(npt),mxh1,mxh2,mxh3,npt,nsepp,nsepv
integer       :: nt,nthr,nc1,nc2,nc3,jq1m,jq2m
integer       :: i,ii,j,jj,iq1,iq2,iq3,jq1,jq2,jq3,jq1min,jq2min,jq2max,jq2t,wfib
character     :: cntid*2,logf*80

!------Open log file---------------------------------
open(11,file=logf,action='write',access='append')

! Some useful print statements. Place where desired
!print *, mxh1,mxh2,mxh3
!print *, nc1,nc2,nc3
!print *, hc1, hc2, hc3
!print *, omp_get_num_threads(), omp_get_max_threads(), omp_get_num_procs()

!----------------------------------------------------
! Reset the counts, set (rpmax,rvmax) and square rp bins
aapv   = 0.0d0
rpmax  = sepp(nsepp+1)
rvmax  = sepv(nsepv+1)
sepp2  = sepp*sepp
rpmax2 = sepp2(nsepp+1)
idsepv = 1./(sepv(2)-sepv(1))  !inverse radial bin size

!----------------------------------------------------
! Unpack boundaries and get cell dimensions in (dec,ra,dcom)
ral = sbound(1) ; rau = sbound(2) ; decl = sbound(3) ; decu = sbound(4)
dcoml = sbound(5) ; dcomu = sbound(6)
nc1 = mxh1
nc2 = mxh2
nc3 = int((dcomu-dcoml)/rvmax)
if(nc3>mxh3) nc3 = mxh3
hc1 = (decu-decl)/float(nc1)     !effective nr of DEC cells
hc2 = (rau-ral)/float(nc2)       !effective nr of RA cells
hc3 = (dcomu-dcoml)/float(nc3)   !effective nr of DCOM cells

!----------------------------------------------------
! Set the number of threads
if(nt<=0) then
   nthr = omp_get_num_procs()
else
   nthr = nt
endif
call omp_set_num_threads(nthr)

!----------------------------------------------------
! Count pairs in SK grid 
write(*,*) ' '
write(*,fmt='(a,i3,a)') '====  Counting '//cntid//' pairs in ', mxh1, ' DEC strips  ===='
!$omp parallel do reduction(+:aapv) default(shared) &
!$omp& private(iq1,iq2,iq3,jq1,jq2,jq3,rcl,stm2,dltdec,jq1m,jq1min) &
!$omp& private(dltra,jq2m,jq2min,jq2max,jq2t,j,i,xi,yi,zi,dci,rv,rp2,ii,jj,wi,wpp,shth2) &
!$omp& schedule(guided) if(nthr>1)
do iq1=1,nc1
   !$omp critical
   write(*,fmt="(i4)",advance='no') iq1                             ! for screen
   write(11,*) cntid//' counting in DEC strip > ',iq1,'/',mxh1      ! for disk
   flush(11)
   !$omp end critical
   do iq2=1,nc2
      do iq3=1,nc3
         i = sk(iq3,iq2,iq1) ! index of ith particle 
         do while(i/=0) !---- Loop over ith particles ----
            lp_jq3: do jq3=iq3,iq3+1
               if(jq3>nc3) cycle lp_jq3
               rcl    = (jq3-1)*hc3+dcoml
               stm2   = sthmax2(rpmax,dc(i),rcl)
               dltdec = 2.0d0*asin(stm2)*rad2deg
               jq1m   = int(dltdec/hc1)+1 !dltdec/hc1+1.
               if(jq3==iq3) then
                  jq1min = iq1
               else
                  jq1min = iq1-jq1m
               end if
               lp_jq1: do jq1=jq1min,iq1+jq1m
                  if(jq1>nc1.or.jq1<1) cycle lp_jq1
                  if(jq1==iq1) then
                     dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,1)
                  else
                     dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,0)
                  end if
                  jq2m   = int(dltra/hc2)+1  !dltra/hc2+1.
                  jq2max = iq2+jq2m
                  jq2min = iq2-jq2m
                  if(jq2max-jq2min+1>nc2) jq2max = jq2min-1+nc2
                  lp_jq2: do jq2=jq2min,jq2max
                     if(jq2>nc2) then
                        jq2t = jq2-nc2
                     else if(jq2<1) then
                        jq2t = jq2+nc2
                     else
                        jq2t = jq2
                     end if
                     if(jq3==iq3.and.jq2t<iq2.and.jq1==iq1) cycle lp_jq2
                     if(jq3==iq3.and.jq2t==iq2.and.jq1==iq1) then
                        j = ll(i)
                     else
                        j = sk(jq3,jq2t,jq1)
                     endif
                     xi  = x(i)
                     yi  = y(i)
                     zi  = z(i)
                     dci = dc(i)
                     wi  = wei(i)
                     
                     do while(j/=0)  !---- Loop over jth particles ----
                        rv = abs(dci-dc(j))
                        if(rv<=rvmax) then
                           rp2 = 4.*dci*dc(j)*((xi-x(j))**2 + (yi-y(j))**2 + (zi-z(j))**2)
                           if(rp2<=rpmax2) then
                               jj = 1 + int(rv*idsepv)       !find rv bin number
                               wpp = wi*wei(j)               !weight by input
                               if(wfib==1) then              !weight by fiber
                                   shth2 = rp2/(4.*dci*dc(j))
                                   wpp   = wpp*wfiber(shth2)
                               endif
                               if(rp2>sepp2(nsepp)) then 
                                  aapv(jj,nsepp) = aapv(jj,nsepp) + wpp
                                  goto 70
                               endif
                               do ii=nsepp-1,1,-1
                                  if(rp2>sepp2(ii)) then
                                     aapv(jj,ii) = aapv(jj,ii) + wpp
                                     goto 70
                                  endif
                               enddo
                           endif !rp2<rpmax2
                        else
                            if(jq3>iq3) cycle lp_jq2  !if(rv>rvmax.and.jq3>iq3)
                        endif
                        70 j = ll(j)
                     end do
                     
                  end do lp_jq2
               end do lp_jq1
            end do lp_jq3
            i = ll(i)
         end do !---- End loop over ith particles ----
      end do  !loop iq3
   end do !loop iq2
end do
!$omp end parallel do
close(11)  ! close log
write(*,*) ' '
end subroutine rppi_A_wg


subroutine rppi_Ab(nt,npt,dec,dc,x,y,z,nsepp,sepp,nsepv,sepv,sbound, &
           mxh1,mxh2,mxh3,nbts,bseed,cntid,logf,sk,ll,aapv,baapv)
!===============================================================================
! NAME
!  rppi_Ab()
!
! DESCRIPTION
!  Count data pairs in projected space for ONE sample of particles. Also count
!  boostraped data pairs
!
! INPUTS
!  Variable---Type-----------Description----------------------------------------
!  nt         int            Number of thread to use (OpenMP)
!  npt        int            Number of particles
!  dec        r8(npt)        DEC of particles [deg]
!  dc         r8(npt)        DCOM of particles [Mpc/h]
!  x,y,z      r8(npt)        X,Y,Z coordinates of particles (see radec2xyz())
!  nsepp      int            Number of projected separation bins
!  sepp       r8(nsepp+1)    Bins in projected separation [Mpc/h]
!  nsepv      int            Number of radial separation bins
!  sepv       r8(nsepv+1)    Bins in radial separation [Mpc/h]
!  sbound     r8(6)          Survey boundaries in RA,DEC,DCOM
!                            Form is (ramin,ramax,decmin,decmax,dcmin,dcmax)
!  mxh1       int            Nr of DEC cells of skip table
!  mxh2       int            Nr of RA cells of skip table
!  mxh3       int            Nr of DCOM cells of skip table
!  nbts       int            Number of boostrap samples
!  bseed      int            Seed for RNG during boostrap resampling
!  cntid      char*2         Two-character string to identify the samples,
!                            e.g. DD, RR, DR,... (used only for ouput messages)
!  logf       char*80        String for the file to log Fortran status messages
!  sk         int(mxh3,mxh2,mxh1)  Skip table (SK) contructed with skll2d()
!  ll         int(npt)       Linked list (LL) contructed with skll2d()
!
! OUTPUTS
!  Variable---Type-------------Description--------------------------------------
!  aapv       r8(nsepv,nsepp)       Counts in radial and projected separation bins
!  baapv      r8(nbts,nsepv,nsepp)  Boostrap counts in radial and projected separation bins
!
! NOTES  -----------------------------------------------------------------------
!  1. RA limits should be set to ramin=0 and ramax=360.
!  2. Remember to update the declarations in cflibfor.pyf if you add/remove 
!     in/out parameters to this exposed subroutine.

implicit none
real(kind=8)  :: dec(npt),sbound(6),ral,rau,decl,decu,dcoml,dcomu,hc1,hc2,hc3,rp2
real(kind=8)  :: x(npt),y(npt),z(npt),dc(npt),rcl,stm2,dltdec,dltra,dci,xi,yi,zi
real(kind=8)  :: sepp(nsepp+1),sepv(nsepv+1),sepp2(nsepp+1),rpmax,rpmax2,rvmax,rv,idsepv
real(kind=8)  :: aapv(nsepv,nsepp),baapv(nbts,nsepv,nsepp)
real(kind=4)  :: wbts(nbts,npt)
integer       :: sk(mxh3,mxh2,mxh1),ll(npt),mxh1,mxh2,mxh3,npt,nsepp,nsepv
integer       :: nt,nthr,nc1,nc2,nc3,jq1m,jq2m,nbts,bseed
integer       :: i,ii,j,jj,iq1,iq2,iq3,jq1,jq2,jq3,jq1min,jq2min,jq2max,jq2t
character     :: cntid*2,logf*80

!------Open log file---------------------------------
open(11,file=logf,action='write',access='append')

! Some useful print statements. Place where desired
!print *, mxh1,mxh2,mxh3
!print *, nc1,nc2,nc3
!print *, hc1, hc2, hc3
!print *, omp_get_num_threads(), omp_get_max_threads(), omp_get_num_procs()

!----------------------------------------------------
! Reset the counts, set (rpmax,rvmax) and square rp bins
aapv   = 0.0d0
baapv  = 0.0d0
rpmax  = sepp(nsepp+1)
rvmax  = sepv(nsepv+1)
sepp2  = sepp*sepp
rpmax2 = sepp2(nsepp+1)
idsepv = 1./(sepv(2)-sepv(1))  !inverse radial bin size

!----------------------------------------------------
! Unpack boundaries and get cell dimensions in (dec,ra,dcom)
ral = sbound(1) ; rau = sbound(2) ; decl = sbound(3) ; decu = sbound(4)
dcoml = sbound(5) ; dcomu = sbound(6)
nc1 = mxh1
nc2 = mxh2
nc3 = int((dcomu-dcoml)/rvmax)
if(nc3>mxh3) nc3 = mxh3
hc1 = (decu-decl)/float(nc1)     !effective nr of DEC cells
hc2 = (rau-ral)/float(nc2)       !effective nr of RA cells
hc3 = (dcomu-dcoml)/float(nc3)   !effective nr of DCOM cells

!----------------------------------------------------
! Generate bootstrap samples
call bootstrap(npt,nbts,bseed,wbts)

!----------------------------------------------------
! Set the number of threads
if(nt<=0) then
   nthr = omp_get_num_procs()
else
   nthr = nt
endif
call omp_set_num_threads(nthr)

!----------------------------------------------------
! Count pairs in SK grid
write(*,*) ' '
write(*,fmt='(a,i3,a)') '====  Counting '//cntid//' pairs in ', mxh1, ' DEC strips  ===='
!$omp parallel do reduction(+:aapv,baapv) default(shared) &
!$omp& private(iq1,iq2,iq3,jq1,jq2,jq3,rcl,stm2,dltdec,jq1m,jq1min) &
!$omp& private(dltra,jq2m,jq2min,jq2max,jq2t,j,i,xi,yi,zi,dci,rv,rp2,ii,jj) &
!$omp& schedule(guided) if(nthr>1)
do iq1=1,nc1
   !$omp critical
   write(*,fmt="(i4)",advance='no') iq1                             ! for screen
   write(11,*) cntid//' counting in DEC strip > ',iq1,'/',mxh1      ! for disk
   flush(11)
   !$omp end critical
   do iq2=1,nc2
      do iq3=1,nc3
         i = sk(iq3,iq2,iq1) ! index of ith particle 
         do while(i/=0) !---- Loop over ith particles ----
            lp_jq3: do jq3=iq3,iq3+1
               if(jq3>nc3) cycle lp_jq3
               rcl    = (jq3-1)*hc3+dcoml
               stm2   = sthmax2(rpmax,dc(i),rcl)
               dltdec = 2.0d0*asin(stm2)*rad2deg
               jq1m   = int(dltdec/hc1)+1 !dltdec/hc1+1.
               if(jq3==iq3) then
                  jq1min = iq1
               else
                  jq1min = iq1-jq1m
               end if
               lp_jq1: do jq1=jq1min,iq1+jq1m
                  if(jq1>nc1.or.jq1<1) cycle lp_jq1
                  if(jq1==iq1) then
                     dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,1)
                  else
                     dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,0)
                  end if
                  jq2m   = int(dltra/hc2)+1  !dltra/hc2+1.
                  jq2max = iq2+jq2m
                  jq2min = iq2-jq2m
                  if(jq2max-jq2min+1>nc2) jq2max = jq2min-1+nc2
                  lp_jq2: do jq2=jq2min,jq2max
                     if(jq2>nc2) then
                        jq2t = jq2-nc2
                     else if(jq2<1) then
                        jq2t = jq2+nc2
                     else
                        jq2t = jq2
                     end if
                     if(jq3==iq3.and.jq2t<iq2.and.jq1==iq1) cycle lp_jq2
                     if(jq3==iq3.and.jq2t==iq2.and.jq1==iq1) then
                        j = ll(i)
                     else
                        j = sk(jq3,jq2t,jq1)
                     endif
                     xi  = x(i)
                     yi  = y(i)
                     zi  = z(i)
                     dci = dc(i)

                     do while(j/=0)   !---- Loop over jth particles ----
                        rv = abs(dci-dc(j))
                        if(rv<=rvmax) then
                           rp2 = 4.*dci*dc(j)*((xi-x(j))**2 + (yi-y(j))**2 + (zi-z(j))**2)
                           if(rp2<=rpmax2) then
                               jj = 1 + int(rv*idsepv)       !find rv bin number
                               if(rp2>sepp2(nsepp)) then 
                                  aapv(jj,nsepp) = aapv(jj,nsepp) + 1.0d0
                                  baapv(:,jj,nsepp) = baapv(:,jj,nsepp) + wbts(:,i)*wbts(:,j)
                                  goto 70
                               endif
                               do ii=nsepp-1,1,-1
                                  if(rp2>sepp2(ii)) then
                                     aapv(jj,ii) = aapv(jj,ii) + 1.0d0
                                     baapv(:,jj,ii) = baapv(:,jj,ii) + wbts(:,i)*wbts(:,j)
                                     goto 70
                                  endif
                               enddo
                           endif !rp2<rpmax2
                        else
                           if(jq3>iq3) cycle lp_jq2  !if(rv>rvmax.and.jq3>iq3)
                        endif
                        70 j = ll(j)
                     end do
                     
                  end do lp_jq2
               end do lp_jq1
            end do lp_jq3
            i = ll(i)
         end do !---- End loop over ith particles ----
      end do !loop iq3
   end do !loop iq2
end do
!$omp end parallel do
close(11)  ! close log
write(*,*) ' '
end subroutine rppi_Ab


subroutine rppi_Ab_wg(nt,npt,dec,dc,wei,x,y,z,nsepp,sepp,nsepv,sepv,sbound, &
           mxh1,mxh2,mxh3,nbts,bseed,wfib,cntid,logf,sk,ll,aapv,baapv)
!===============================================================================
! NAME
!  rppi_Ab_wg()
!
! DESCRIPTION
!  Count weighted  pairs in projected space for ONE sample of particles. Also count
!  boostraped data pairs
!
! INPUTS
!  Variable---Type-----------Description----------------------------------------
!  nt         int            Number of thread to use (OpenMP)
!  npt        int            Number of particles
!  dec        r8(npt)        DEC of particles [deg]
!  dc         r8(npt)        DCOM of particles [Mpc/h]
!  wei        r4(npt)        WEIGHT of particles
!  x,y,z      r8(npt)        X,Y,Z coordinates of particles (see radec2xyz())
!  nsepp      int            Number of projected separation bins
!  sepp       r8(nsepp+1)    Bins in projected separation [Mpc/h]
!  nsepv      int            Number of radial separation bins
!  sepv       r8(nsepv+1)    Bins in radial separation [Mpc/h]
!  sbound     r8(6)          Survey boundaries in RA,DEC,DCOM
!                            Form is (ramin,ramax,decmin,decmax,dcmin,dcmax)
!  mxh1       int            Nr of DEC cells of skip table
!  mxh2       int            Nr of RA cells of skip table
!  mxh3       int            Nr of DCOM cells of skip table
!  nbts       int            Number of boostrap samples
!  bseed      int            Seed for RNG during boostrap resampling
!  wfib       int            If 1, also apply fiber-correction weights
!  cntid      char*2         Two-character string to identify the samples,
!                            e.g. DD, RR, DR,... (used only for ouput messages)
!  logf       char*80        String for the file to log Fortran status messages
!  sk         int(mxh3,mxh2,mxh1)  Skip table (SK) contructed with skll2d()
!  ll         int(npt)       Linked list (LL) contructed with skll2d()
!
! OUTPUTS
!  Variable---Type-------------Description--------------------------------------
!  aapv       r8(nsepv,nsepp)       Counts in radial and projected separation bins
!  baapv      r8(nbts,nsepv,nsepp)  Boostrap counts in radial and projected separation bins
!
! NOTES  -----------------------------------------------------------------------
!  1. RA limits should be set to ramin=0 and ramax=360.
!  2. Remember to update the declarations in cflibfor.pyf if you add/remove 
!     in/out parameters to this exposed subroutine.

implicit none
real(kind=8)  :: dec(npt),sbound(6),ral,rau,decl,decu,dcoml,dcomu,hc1,hc2,hc3,rp2
real(kind=4)  :: wei(npt),wpp,wi
real(kind=8)  :: x(npt),y(npt),z(npt),dc(npt),rcl,stm2,dltdec,dltra,dci,xi,yi,zi
real(kind=8)  :: sepp(nsepp+1),sepv(nsepv+1),sepp2(nsepp+1),rpmax,rpmax2,rvmax,rv,shth2,idsepv
real(kind=8)  :: aapv(nsepv,nsepp),baapv(nbts,nsepv,nsepp)
real(kind=4)  :: wbts(nbts,npt)
integer       :: sk(mxh3,mxh2,mxh1),ll(npt),mxh1,mxh2,mxh3,npt,nsepp,nsepv
integer       :: nt,nthr,nc1,nc2,nc3,jq1m,jq2m,nbts,bseed,wfib
integer       :: i,ii,j,jj,iq1,iq2,iq3,jq1,jq2,jq3,jq1min,jq2min,jq2max,jq2t
character     :: cntid*2,logf*80

!------Open log file---------------------------------
open(11,file=logf,action='write',access='append')

! Some useful print statements. Place where desired
!print *, mxh1,mxh2,mxh3
!print *, nc1,nc2,nc3
!print *, hc1, hc2, hc3
!print *, omp_get_num_threads(), omp_get_max_threads(), omp_get_num_procs()

!----------------------------------------------------
! Reset the counts, set (rpmax,rvmax) and square rp bins
aapv   = 0.0d0
baapv  = 0.0d0
rpmax  = sepp(nsepp+1)
rvmax  = sepv(nsepv+1)
sepp2  = sepp*sepp
rpmax2 = sepp2(nsepp+1)
idsepv = 1./(sepv(2)-sepv(1))  !inverse radial bin size

!----------------------------------------------------
! Unpack boundaries and get cell dimensions in (dec,ra,dcom)
ral = sbound(1) ; rau = sbound(2) ; decl = sbound(3) ; decu = sbound(4)
dcoml = sbound(5) ; dcomu = sbound(6)
nc1 = mxh1
nc2 = mxh2
nc3 = int((dcomu-dcoml)/rvmax)
if(nc3>mxh3) nc3 = mxh3
hc1 = (decu-decl)/float(nc1)     !effective nr of DEC cells
hc2 = (rau-ral)/float(nc2)       !effective nr of RA cells
hc3 = (dcomu-dcoml)/float(nc3)   !effective nr of DCOM cells

!----------------------------------------------------
! Generate bootstrap samples
call bootstrap(npt,nbts,bseed,wbts)

!----------------------------------------------------
! Set the number of threads
if(nt<=0) then
   nthr = omp_get_num_procs()
else
   nthr = nt
endif
call omp_set_num_threads(nthr)

!----------------------------------------------------
! Count pairs in SK grid
write(*,*) ' '
write(*,fmt='(a,i3,a)') '====  Counting '//cntid//' pairs in ', mxh1, ' DEC strips  ===='
!$omp parallel do reduction(+:aapv,baapv) default(shared) &
!$omp& private(iq1,iq2,iq3,jq1,jq2,jq3,rcl,stm2,dltdec,jq1m,jq1min) &
!$omp& private(dltra,jq2m,jq2min,jq2max,jq2t,j,i,xi,yi,zi,dci,rv,rp2,ii,jj,wi,wpp,shth2) &
!$omp& schedule(guided) if(nthr>1)
do iq1=1,nc1
   !$omp critical
   write(*,fmt="(i4)",advance='no') iq1                             ! for screen
   write(11,*) cntid//' counting in DEC strip > ',iq1,'/',mxh1      ! for disk
   flush(11)
   !$omp end critical
   do iq2=1,nc2
      do iq3=1,nc3
         i = sk(iq3,iq2,iq1) ! index of ith particle 
         do while(i/=0) !---- Loop over ith particles ----
            lp_jq3: do jq3=iq3,iq3+1
               if(jq3>nc3) cycle lp_jq3
               rcl    = (jq3-1)*hc3+dcoml
               stm2   = sthmax2(rpmax,dc(i),rcl)
               dltdec = 2.0d0*asin(stm2)*rad2deg
               jq1m   = int(dltdec/hc1)+1 !dltdec/hc1+1.
               if(jq3==iq3) then
                  jq1min = iq1
               else
                  jq1min = iq1-jq1m
               end if
               lp_jq1: do jq1=jq1min,iq1+jq1m
                  if(jq1>nc1.or.jq1<1) cycle lp_jq1
                  if(jq1==iq1) then
                     dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,1)
                  else
                     dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,0)
                  end if
                  jq2m   = int(dltra/hc2)+1  !dltra/hc2+1.
                  jq2max = iq2+jq2m
                  jq2min = iq2-jq2m
                  if(jq2max-jq2min+1>nc2) jq2max = jq2min-1+nc2
                  lp_jq2: do jq2=jq2min,jq2max
                     if(jq2>nc2) then
                        jq2t = jq2-nc2
                     else if(jq2<1) then
                        jq2t = jq2+nc2
                     else
                        jq2t = jq2
                     end if
                     if(jq3==iq3.and.jq2t<iq2.and.jq1==iq1) cycle lp_jq2
                     if(jq3==iq3.and.jq2t==iq2.and.jq1==iq1) then
                        j = ll(i)
                     else
                        j = sk(jq3,jq2t,jq1)
                     endif
                     xi  = x(i)
                     yi  = y(i)
                     zi  = z(i)
                     dci = dc(i)
                     wi  = wei(i)

                     do while(j/=0)   !---- Loop over jth particles ----
                        rv = abs(dci-dc(j))
                        if(rv<=rvmax) then
                           rp2 = 4.*dci*dc(j)*((xi-x(j))**2 + (yi-y(j))**2 + (zi-z(j))**2)
                           if(rp2<=rpmax2) then
                               jj  = 1 + int(rv*idsepv)      !find rv bin number
                               wpp = wi*wei(j)               !weight by input
                               if(wfib==1) then              !weight by fiber
                                   shth2 = rp2/(4.*dci*dc(j))
                                   wpp   = wpp*wfiber(shth2)
                               endif
                               if(rp2>sepp2(nsepp)) then 
                                  aapv(jj,nsepp) = aapv(jj,nsepp) + wpp
                                  baapv(:,jj,nsepp) = baapv(:,jj,nsepp) + wpp*wbts(:,i)*wbts(:,j)
                                  goto 70
                               endif
                               do ii=nsepp-1,1,-1
                                  if(rp2>sepp2(ii)) then
                                     aapv(jj,ii) = aapv(jj,ii) + wpp
                                     baapv(:,jj,ii) = baapv(:,jj,ii) + wpp*wbts(:,i)*wbts(:,j)
                                     goto 70
                                  endif
                               enddo
                           endif !rp2<rpmax2
                        else
                           if(jq3>iq3) cycle lp_jq2  !if(rv>rvmax.and.jq3>iq3)
                        endif
                        70 j = ll(j)
                     end do
                     
                  end do lp_jq2
               end do lp_jq1
            end do lp_jq3
            i = ll(i)
         end do !---- End loop over ith particles ----
      end do !loop iq3
   end do !loop iq2
end do
!$omp end parallel do
close(11)  ! close log
write(*,*) ' '
end subroutine rppi_Ab_wg


subroutine rppi_C(nt,npt,ra,dec,dc,x,y,z,npt1,dc1,x1,y1,z1, &
           nsepp,sepp,nsepv,sepv,sbound,mxh1,mxh2,mxh3,cntid,logf,sk1,ll1,cdpv)
!===============================================================================
! NAME
!  rppi_C()
!
! DESCRIPTION
!  Cross-count data pairs in projected space for TWO samples of particles
!
! INPUTS
!  Variable---Type-----------Description----------------------------------------
!  nt         int            Number of thread to use (OpenMP)
!  npt        int            [sampleC] Number of particles
!  ra         r8(npt)        [sampleC] RA of particles [deg]
!  dec        r8(npt)        [sampleC] DEC of particles [deg]
!  dc         r8(npt)        [sampleC] DCOM of particles [Mpc/h]
!  x,y,z      r8(npt)        [sampleC] X,Y,Z coordinates of particles (see radec2xyz())
!  npt1       int            [sampleD] Number of particles
!  dc1        r8(npt1)       [sampleD] DCOM of particles [Mpc/h]
!  x1,y1,z1   r8(npt1)       [sampleD] X,Y,Z coordinates of particles (see radec2xyz())
!  nsepp      int            Number of projected separation bins
!  sepp       r8(nsepp+1)    Bins in projected separation [Mpc/h]
!  nsepv      int            Number of radial separation bins
!  sepv       r8(nsepv+1)    Bins in radial separation [Mpc/h]
!  sbound     r8(6)          Survey boundaries in RA,DEC,DCOM
!                            Form is (ramin,ramax,decmin,decmax,dcmin,dcmax)
!  mxh1       int            Nr of DEC cells of skip table
!  mxh2       int            Nr of RA cells of skip table
!  mxh3       int            Nr of DCOM cells of skip table
!  cntid      char*2         Two-character string to identify the samples,
!                            e.g. DD, RR, DR,... (used only for ouput messages)
!  logf       char*80        String for the file to log Fortran status messages
!  sk1        int(mxh3,mxh2,mxh1) [sampleD] Skip table (SK) contructed with skll2d()
!  ll1        int(npt1)      [sampleD] Linked list (LL) contructed with skll2d()
!
! OUTPUTS
!  Variable---Type-------------Description--------------------------------------
!  cdpv       r8(nsepv,nsepp)  Counts in radial and projected separation bins
!
! NOTES  -----------------------------------------------------------------------
!  1. RA limits should be set to ramin=0 and ramax=360.
!  2. Remember to update the declarations in cflibfor.pyf if you add/remove 
!     in/out parameters to this exposed subroutine.

implicit none
real(kind=8)  :: ra(npt),dec(npt),sbound(6),ral,rau,decl,decu,dcoml,dcomu,hc1,hc2,hc3,rp2
real(kind=8)  :: x(npt),y(npt),z(npt),x1(npt1),y1(npt1),z1(npt1),rcl,stm2,dltdec,dltra
real(kind=8)  :: dc(npt),dc1(npt1),xi,dci
real(kind=8)  :: sepp(nsepp+1),sepv(nsepv+1),sepp2(nsepp+1),rpmax,rpmax2,rvmax,rv,idsepv
real(kind=8)  :: cdpv(nsepv,nsepp)
integer       :: sk1(mxh3,mxh2,mxh1),ll1(npt1),mxh1,mxh2,mxh3,npt,npt1,nsepp,nsepv
integer       :: nt,nthr,fracp,dpart,nadv,nc1,nc2,nc3,jq1m,jq2m
integer       :: i,ii,j,jj,iq1,iq2,iq3,jq1,jq2,jq3,jq2min,jq2max,jq2t,p1,p2
character     :: cntid*2,logf*80

!------Open log file---------------------------------
open(11,file=logf,action='write',access='append')

fracp = 0 ; dpart = 0 ; nadv = 0  !Reset progress counters

!----------------------------------------------------
! Reset the counts, set (rpmax,rvmax) and square rp bins
cdpv   = 0.d0
rpmax  = sepp(nsepp+1)
rvmax  = sepv(nsepv+1)
sepp2  = sepp*sepp
rpmax2 = sepp2(nsepp+1)
idsepv = 1./(sepv(2)-sepv(1))  !inverse radial bin size

!----------------------------------------------------
! Unpack boundaries and get cell dimensions in (dec,ra,dcom)
ral = sbound(1) ; rau = sbound(2) ; decl = sbound(3) ; decu = sbound(4)
dcoml = sbound(5) ; dcomu = sbound(6)
nc1 = mxh1
nc2 = mxh2
nc3 =int((dcomu-dcoml)/rvmax)
if(nc3>mxh3) nc3 = mxh3
hc1 = (decu-decl)/float(nc1)    !effective nr of DEC cells
hc2 = (rau-ral)/float(nc2)      !effective nr of RA cells
hc3 = (dcomu-dcoml)/float(nc3)  !effective nr of DCOM cells

!----------------------------------------------------
! Set the number of threads
if(nt<=0) then
   nthr = omp_get_num_procs()
else
   nthr = nt
endif
call omp_set_num_threads(nthr)

!----------------------------------------------------
! Some inits
dpart = int(float(npt)/float(mxh1)) !choose dpart so we get mhx1 parts

!----------------------------------------------------
! Count pairs in SK grid
write(*,*) ' '
write(*,fmt='(a,i3,a)') '====  Counting '//cntid//' pairs in ', mxh1, ' DEC strips'
!$omp parallel do reduction(+:cdpv) default(shared) &
!$omp& private(iq1,iq2,iq3,jq1,jq2,jq3,rcl,stm2,dltdec,jq1m) &
!$omp& private(dltra,jq2m,jq2min,jq2max,jq2t,j,i,xi,dci,rv,rp2,ii,jj) &
!$omp& schedule(guided) firstprivate(fracp) if(nthr>1)

do i=1,npt   !---- Loop over ith particles ----
   xi  = x(i)
   dci = dc(i)
   
   iq1 = int((dec(i)-decl)/hc1)+1
   iq2 = int((ra(i)-ral)/hc2)+1
   iq3 = int((dci-dcoml)/hc3)+1

   fracp = fracp + 1  ! accumulate particles and check when above the step size
   if(fracp>=dpart) then
       !$omp critical
       nadv = nadv + 1
       !omp flush (nadv)
       p1   = (nadv-1)*dpart + 1
       p2   = nadv*dpart
       if((nadv+1)*dpart>npt) p2=npt
       ! Note we are not really counting in mxh1 strips, just mymicking
       write(*,fmt="(i4)",advance='no') nadv                        ! for screen
       write(11,*) cntid//' counting in DEC strip > ',nadv,' (',p1,'-',p2,')' !for disk
       fracp = 0
       !$omp end critical
   endif
   
   lp_jq3: do jq3 = iq3-1,iq3+1
      if(jq3<1.or.jq3>nc3) cycle lp_jq3
      rcl    = (jq3-1)*hc3+dcoml
      stm2   = sthmax2(rpmax,dci,rcl)
      dltdec = 2.0d0*asin(stm2)*rad2deg
      jq1m   = int(dltdec/hc1)+1
      lp_jq1: do jq1=iq1-jq1m,iq1+jq1m
         if(jq1>nc1.or.jq1<1) cycle lp_jq1
         if(jq1==iq1) then
            dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,1)
         else
            dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,0)
         end if
         jq2m   = int(dltra/hc2)+1
         jq2max = iq2+jq2m
         jq2min = iq2-jq2m
         if(jq2max-jq2min+1>nc2) jq2max = jq2min-1+nc2
         lp_jq2: do jq2=jq2min,jq2max
            if(jq2>nc2) then
               jq2t = jq2-nc2
            else if(jq2<1) then
               jq2t = jq2+nc2
            else
               jq2t = jq2
            end if
            j = sk1(jq3,jq2t,jq1)
            do while(j/=0)
               !---- Loop over jth particles ----
               ! For each ith particle, index j runs through all particles in
               ! neighboring cells. Prune those outside rv and rp2 range.
               rv  = abs(dci-dc1(j))
               if(rv<=rvmax) then
                  rp2 = 4.*dci*dc1(j)*((xi-x1(j))**2 + (y(i)-y1(j))**2 + (z(i)-z1(j))**2)
                  if(rp2<=rpmax2) then
                     jj = 1 + int(rv*idsepv)       !find rv bin number
                     if(rp2>sepp2(nsepp)) then 
                        cdpv(jj,nsepp) = cdpv(jj,nsepp) + 1.0d0
                        goto 70
                     endif
                     do ii=nsepp-1,1,-1
                        if(rp2>sepp2(ii)) then
                           cdpv(jj,ii) = cdpv(jj,ii) + 1.0d0
                           goto 70
                        endif
                     enddo
                  endif
               else 
                  if(jq3>iq3) cycle lp_jq2  !if(rv>rvmax.and.jq3>iq3)
               endif
               70 j = ll1(j)
            end do
         end do lp_jq2
      end do lp_jq1
   end do lp_jq3
end do
!$omp end parallel do
close(11)  ! close log
write(*,*) ' '
end subroutine rppi_C


subroutine rppi_C_wg(nt,npt,ra,dec,dc,wei,x,y,z,npt1,dc1,wei1,x1,y1,z1, &
           nsepp,sepp,nsepv,sepv,sbound,mxh1,mxh2,mxh3,wfib,cntid,logf,sk1,ll1,cdpv)
!===============================================================================
! NAME
!  rppi_C_wg()
!
! DESCRIPTION
!  Cross-count weighted data pairs in projected space for TWO samples of particles
!
! INPUTS
!  Variable---Type-----------Description----------------------------------------
!  nt         int            Number of thread to use (OpenMP)
!  npt        int            [sampleC] Number of particles
!  ra         r8(npt)        [sampleC] RA of particles [deg]
!  dec        r8(npt)        [sampleC] DEC of particles [deg]
!  dc         r8(npt)        [sampleC] DCOM of particles [Mpc/h]
!  wei        r4(npt)        [sampleC] WEIGHT of particles
!  x,y,z      r8(npt)        [sampleC] X,Y,Z coordinates of particles (see radec2xyz())
!  npt1       int            [sampleD] Number of particles
!  dc1        r8(npt1)       [sampleD] DCOM of particles [Mpc/h]
!  wei1       r4(npt1)       [sampleD] WEIGHT of particles
!  x1,y1,z1   r8(npt1)       [sampleD] X,Y,Z coordinates of particles (see radec2xyz())
!  nsepp      int            Number of projected separation bins
!  sepp       r8(nsepp+1)    Bins in projected separation [Mpc/h]
!  nsepv      int            Number of radial separation bins
!  sepv       r8(nsepv+1)    Bins in radial separation [Mpc/h]
!  sbound     r8(6)          Survey boundaries in RA,DEC,DCOM
!                            Form is (ramin,ramax,decmin,decmax,dcmin,dcmax)
!  mxh1       int            Nr of DEC cells of skip table
!  mxh2       int            Nr of RA cells of skip table
!  mxh3       int            Nr of DCOM cells of skip table
!  wfib       int            If 1, also apply fiber-correction weights
!  cntid      char*2         Two-character string to identify the samples,
!                            e.g. DD, RR, DR,... (used only for ouput messages)
!  logf       char*80        String for the file to log Fortran status messages
!  sk1        int(mxh3,mxh2,mxh1) [sampleD] Skip table (SK) contructed with skll2d()
!  ll1        int(npt1)      [sampleD] Linked list (LL) contructed with skll2d()
!
! OUTPUTS
!  Variable---Type-------------Description--------------------------------------
!  cdpv       r8(nsepv,nsepp)  Counts in radial and projected separation bins
!
! NOTES  -----------------------------------------------------------------------
!  1. RA limits should be set to ramin=0 and ramax=360.
!  2. Remember to update the declarations in cflibfor.pyf if you add/remove 
!     in/out parameters to this exposed subroutine.

implicit none
real(kind=8)  :: ra(npt),dec(npt),sbound(6),ral,rau,decl,decu,dcoml,dcomu,hc1,hc2,hc3,rp2
real(kind=8)  :: x(npt),y(npt),z(npt),x1(npt1),y1(npt1),z1(npt1),rcl,stm2,dltdec,dltra
real(kind=8)  :: dc(npt),dc1(npt1),xi,dci
real(kind=4)  :: wei(npt),wei1(npt1),wi,wpp
real(kind=8)  :: sepp(nsepp+1),sepv(nsepv+1),sepp2(nsepp+1),rpmax,rpmax2,rvmax,rv,idsepv,shth2
real(kind=8)  :: cdpv(nsepv,nsepp)
integer       :: sk1(mxh3,mxh2,mxh1),ll1(npt1),mxh1,mxh2,mxh3,npt,npt1,nsepp,nsepv,wfib
integer       :: nt,nthr,ndp,fracp,dpart,nadv,nc1,nc2,nc3,jq1m,jq2m
integer       :: i,ii,j,jj,iq1,iq2,iq3,jq1,jq2,jq3,jq2min,jq2max,jq2t,p1,p2
character     :: cntid*2,logf*80

!------Open log file---------------------------------
open(11,file=logf,action='write',access='append')

fracp = 0 ; dpart = 0 ; nadv = 0  !Reset progress counters

!----------------------------------------------------
! Reset the counts, set (rpmax,rvmax) and square rp bins
cdpv   = 0.d0
rpmax  = sepp(nsepp+1)
rvmax  = sepv(nsepv+1)
sepp2  = sepp*sepp
rpmax2 = sepp2(nsepp+1)
idsepv = 1./(sepv(2)-sepv(1))  !inverse radial bin size

!----------------------------------------------------
! Unpack boundaries and get cell dimensions in (dec,ra,dcom)
ral = sbound(1) ; rau = sbound(2) ; decl = sbound(3) ; decu = sbound(4)
dcoml = sbound(5) ; dcomu = sbound(6)
nc1 = mxh1
nc2 = mxh2
nc3 =int((dcomu-dcoml)/rvmax)
if(nc3>mxh3) nc3 = mxh3
hc1 = (decu-decl)/float(nc1)    !effective nr of DEC cells
hc2 = (rau-ral)/float(nc2)      !effective nr of RA cells
hc3 = (dcomu-dcoml)/float(nc3)  !effective nr of DCOM cells

!----------------------------------------------------
! Set the number of threads
if(nt<=0) then
   nthr = omp_get_num_procs()
else
   nthr = nt
endif
call omp_set_num_threads(nthr)

!----------------------------------------------------
! Some inits
dpart = int(float(npt)/float(mxh1)) !choose dpart so we get mhx1 parts

!----------------------------------------------------
! Count pairs in SK grid
write(*,*) ' '
write(*,fmt='(a,i3,a)') '====  Counting '//cntid//' pairs in ', mxh1, ' DEC strips'
!$omp parallel do reduction(+:cdpv) default(shared) &
!$omp& private(iq1,iq2,iq3,jq1,jq2,jq3,rcl,stm2,dltdec,jq1m) &
!$omp& private(dltra,jq2m,jq2min,jq2max,jq2t,j,i,xi,wi,dci,rv,rp2,shth2,wpp,ii,jj) &
!$omp& schedule(guided) firstprivate(fracp) if(nthr>1)

do i=1,npt   !---- Loop over ith particles ----
   xi  = x(i)
   dci = dc(i)
   wi  = wei(i)
   
   iq1 = int((dec(i)-decl)/hc1)+1
   iq2 = int((ra(i)-ral)/hc2)+1
   iq3 = int((dci-dcoml)/hc3)+1

   fracp = fracp + 1  ! accumulate particles and check when above the step size
   if(fracp>=dpart) then
       !$omp critical
       nadv = nadv + 1
       !omp flush (nadv)
       p1   = (nadv-1)*dpart + 1
       p2   = nadv*dpart
       if((nadv+1)*dpart>npt) p2=npt
       ! Note we are not really counting in mxh1 strips, just mymicking
       write(*,fmt="(i4)",advance='no') nadv                        ! for screen
       write(11,*) cntid//' counting in DEC strip > ',nadv,' (',p1,'-',p2,')' !for disk
       fracp = 0
       !$omp end critical
   endif
   
   lp_jq3: do jq3 = iq3-1,iq3+1
      if(jq3<1.or.jq3>nc3) cycle lp_jq3
      rcl    = (jq3-1)*hc3+dcoml
      stm2   = sthmax2(rpmax,dci,rcl)
      dltdec = 2.0d0*asin(stm2)*rad2deg
      jq1m   = int(dltdec/hc1)+1
      lp_jq1: do jq1=iq1-jq1m,iq1+jq1m
         if(jq1>nc1.or.jq1<1) cycle lp_jq1
         if(jq1==iq1) then
            dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,1)
         else
            dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,0)
         end if
         jq2m   = int(dltra/hc2)+1
         jq2max = iq2+jq2m
         jq2min = iq2-jq2m
         if(jq2max-jq2min+1>nc2) jq2max = jq2min-1+nc2
         lp_jq2: do jq2=jq2min,jq2max
            if(jq2>nc2) then
               jq2t = jq2-nc2
            else if(jq2<1) then
               jq2t = jq2+nc2
            else
               jq2t = jq2
            end if
            j = sk1(jq3,jq2t,jq1)
            do while(j/=0)
               !---- Loop over jth particles ----
               ! For each ith particle, index j runs through all particles in
               ! neighboring cells. Prune those outside rv and rp2 range.
               rv  = abs(dci-dc1(j))
               if(rv<=rvmax) then
                  rp2 = 4.*dci*dc1(j)*((xi-x1(j))**2 + (y(i)-y1(j))**2 + (z(i)-z1(j))**2)
                  if(rp2<=rpmax2) then
                     jj = 1 + int(rv*idsepv)             !find rv bin number
                     wpp = wi*wei1(j)                    !weight by input
                     if(wfib==1) then                    !weight by fiber
                        shth2 = rp2/(4.*dci*dc1(j))
                        wpp   = wpp*wfiber(shth2)                        
                     endif
                     if(rp2>sepp2(nsepp)) then 
                        cdpv(jj,nsepp) = cdpv(jj,nsepp) + wpp
                        goto 70
                     endif
                     do ii=nsepp-1,1,-1
                        if(rp2>sepp2(ii)) then
                           cdpv(jj,ii) = cdpv(jj,ii) + wpp
                           goto 70
                        endif
                     enddo
                  endif
               else 
                  if(jq3>iq3) cycle lp_jq2  !if(rv>rvmax.and.jq3>iq3)
               endif
               70 j = ll1(j)
            end do
         end do lp_jq2
      end do lp_jq1
   end do lp_jq3
end do
!$omp end parallel do
close(11)  ! close log
write(*,*) ' '
end subroutine rppi_C_wg


subroutine rppi_Cb(nt,npt,ra,dec,dc,x,y,z,npt1,dc1,x1,y1,z1,nsepp,sepp,nsepv,sepv, &
           sbound,mxh1,mxh2,mxh3,nbts,bseed,cntid,logf,sk1,ll1,cdpv,bcdpv)
!===============================================================================
! NAME
!  rppi_Cb()
!
! DESCRIPTION
!  Cross-count data pairs in projected space for TWO samples of particles. Also
!  count boostraped data pairs
!
! INPUTS
!  Variable---Type-----------Description----------------------------------------
!  nt         int            Number of thread to use (OpenMP)
!  npt        int            [sampleC] Number of particles
!  ra         r8(npt)        [sampleC] RA of particles [deg]
!  dec        r8(npt)        [sampleC] DEC of particles [deg]
!  dc         r8(npt)        [sampleC] DCOM of particles [Mpc/h]
!  x,y,z      r8(npt)        [sampleC] X,Y,Z coordinates of particles (see radec2xyz())
!  npt1       int            [sampleD] Number of particles
!  dc1        r8(npt1)       [sampleD] DCOM of particles [Mpc/h]
!  x1,y1,z1   r8(npt1)       [sampleD] X,Y,Z coordinates of particles (see radec2xyz())
!  nsepp      int            Number of projected separation bins
!  sepp       r8(nsepp+1)    Bins in projected separation [Mpc/h]
!  nsepv      int            Number of radial separation bins
!  sepv       r8(nsepv+1)    Bins in radial separation [Mpc/h]
!  sbound     r8(6)          Survey boundaries in RA,DEC,DCOM
!                            Form is (ramin,ramax,decmin,decmax,dcmin,dcmax)
!  mxh1       int            Nr of DEC cells of skip table
!  mxh2       int            Nr of RA cells of skip table
!  mxh3       int            Nr of DCOM cells of skip table
!  nbts       int            Number of boostrap samples
!  bseed      int            Seed for RNG during boostrap resampling
!  cntid      char*2         Two-character string to identify the samples,
!                            e.g. DD, RR, DR,... (used only for ouput messages)
!  logf       char*80        String for the file to log Fortran status messages
!  sk1        int(mxh3,mxh2,mxh1) [sampleD] Skip table (SK) contructed with skll2d()
!  ll1        int(npt1)      [sampleD] Linked list (LL) contructed with skll2d()
!
! OUTPUTS
!  Variable---Type-------------Description--------------------------------------
!  cdpv       r8(nsepv,nsepp)      Counts in radial and projected separation bins
!  bcdpv      r8(nbts,nsepv,nsepp) Boostrap counts in radial and projected separation bins
!
! NOTES  -----------------------------------------------------------------------
!  1. RA limits should be set to ramin=0 and ramax=360.
!  2. Remember to update the declarations in cflibfor.pyf if you add/remove 
!     in/out parameters to this exposed subroutine.

implicit none
real(kind=8)  :: ra(npt),dec(npt),sbound(6),ral,rau,decl,decu,dcoml,dcomu,hc1,hc2,hc3,rp2
real(kind=8)  :: x(npt),y(npt),z(npt),x1(npt1),y1(npt1),z1(npt1),rcl,stm2,dltdec,dltra
real(kind=8)  :: dc(npt),dc1(npt1),xi,dci
real(kind=4)  :: wbts(nbts,npt),wbts1(nbts,npt1)
real(kind=8)  :: sepp(nsepp+1),sepv(nsepv+1),sepp2(nsepp+1),rpmax,rpmax2,rvmax,rv,idsepv
real(kind=8)  :: cdpv(nsepv,nsepp),bcdpv(nbts,nsepv,nsepp)
integer       :: sk1(mxh3,mxh2,mxh1),ll1(npt1),mxh1,mxh2,mxh3,npt,npt1,nsepp,nsepv,nbts,bseed
integer       :: nt,nthr,ndp,fracp,dpart,nadv,nc1,nc2,nc3,jq1m,jq2m
integer       :: i,ii,j,jj,iq1,iq2,iq3,jq1,jq2,jq3,jq2min,jq2max,jq2t,p1,p2
character     :: cntid*2,logf*80

!------Open log file---------------------------------
open(11,file=logf,action='write',access='append')

fracp = 0 ; dpart = 0 ; nadv = 0  !Reset progress counters

!----------------------------------------------------
! Reset the counts, set (rpmax,rvmax) and square rp bins
cdpv   = 0.d0
bcdpv  = 0.0
rpmax  = sepp(nsepp+1)
rvmax  = sepv(nsepv+1)
sepp2  = sepp*sepp
rpmax2 = sepp2(nsepp+1)
idsepv = 1./(sepv(2)-sepv(1))  !inverse radial bin size

!----------------------------------------------------
! Unpack boundaries and get cell dimensions in (dec,ra,dcom)
ral = sbound(1) ; rau = sbound(2) ; decl = sbound(3) ; decu = sbound(4)
dcoml = sbound(5) ; dcomu = sbound(6)
nc1 = mxh1
nc2 = mxh2
nc3 =int((dcomu-dcoml)/rvmax)
if(nc3>mxh3) nc3 = mxh3
hc1 = (decu-decl)/float(nc1)    !effective nr of DEC cells
hc2 = (rau-ral)/float(nc2)      !effective nr of RA cells
hc3 = (dcomu-dcoml)/float(nc3)  !effective nr of DCOM cells

!----------------------------------------------------
! Generate bootstrap samples
call bootstrap(npt,nbts,bseed,wbts)
call bootstrap(npt1,nbts,bseed,wbts1)

!----------------------------------------------------
! Set the number of threads
if(nt<=0) then
   nthr = omp_get_num_procs()
else
   nthr = nt
endif
call omp_set_num_threads(nthr)

!----------------------------------------------------
! Some inits
dpart = int(float(npt)/float(mxh1)) !choose dpart so we get mhx1 parts

!----------------------------------------------------
! Count pairs in SK grid
write(*,*) ' '
write(*,fmt='(a,i3,a)') '====  Counting '//cntid//' pairs in ', mxh1, ' DEC strips'
!$omp parallel do reduction(+:cdpv,bcdpv) default(shared) &
!$omp& private(iq1,iq2,iq3,jq1,jq2,jq3,rcl,stm2,dltdec,jq1m) &
!$omp& private(dltra,jq2m,jq2min,jq2max,jq2t,j,i,xi,dci,rv,rp2,ii,jj) &
!$omp& schedule(guided) firstprivate(fracp) if(nthr>1)

do i=1,npt   !---- Loop over ith particles ----
   xi  = x(i)
   dci = dc(i)
   
   iq1 = int((dec(i)-decl)/hc1)+1
   iq2 = int((ra(i)-ral)/hc2)+1
   iq3 = int((dci-dcoml)/hc3)+1

   fracp = fracp + 1  ! accumulate particles and check when above the step size
   if(fracp>=dpart) then
       !$omp critical
       nadv = nadv + 1
       !omp flush (nadv)
       p1   = (nadv-1)*dpart + 1
       p2   = nadv*dpart
       if((nadv+1)*dpart>npt) p2=npt
       ! Note we are not really counting in mxh1 strips, just mymicking
       write(*,fmt="(i4)",advance='no') nadv                        ! for screen
       write(11,*) cntid//' counting in DEC strip > ',nadv,' (',p1,'-',p2,')' !for disk
       fracp = 0
       !$omp end critical
   endif
   
   lp_jq3: do jq3 = iq3-1,iq3+1
      if(jq3<1.or.jq3>nc3) cycle lp_jq3
      rcl    = (jq3-1)*hc3+dcoml
      stm2   = sthmax2(rpmax,dci,rcl)
      dltdec = 2.0d0*asin(stm2)*rad2deg
      jq1m   = int(dltdec/hc1)+1
      lp_jq1: do jq1=iq1-jq1m,iq1+jq1m
         if(jq1>nc1.or.jq1<1) cycle lp_jq1
         if(jq1==iq1) then
            dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,1)
         else
            dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,0)
         end if
         jq2m   = int(dltra/hc2)+1
         jq2max = iq2+jq2m
         jq2min = iq2-jq2m
         if(jq2max-jq2min+1>nc2) jq2max = jq2min-1+nc2
         lp_jq2: do jq2=jq2min,jq2max
            if(jq2>nc2) then
               jq2t = jq2-nc2
            else if(jq2<1) then
               jq2t = jq2+nc2
            else
               jq2t = jq2
            end if
            j = sk1(jq3,jq2t,jq1)
            do while(j/=0)
               !---- Loop over jth particles ----
               ! For each ith particle, index j runs through all particles in
               ! neighboring cells. Prune those outside rv and rp2 range.
               rv  = abs(dci-dc1(j))
               if(rv<=rvmax) then
                  rp2 = 4.*dci*dc1(j)*((xi-x1(j))**2 + (y(i)-y1(j))**2 + (z(i)-z1(j))**2)
                  if(rp2<=rpmax2) then
                     jj = 1 + int(rv*idsepv)       !find rv bin number
                     if(rp2>sepp2(nsepp)) then 
                        cdpv(jj,nsepp) = cdpv(jj,nsepp) + 1.0d0
                        bcdpv(:,jj,nsepp) = bcdpv(:,jj,nsepp) + wbts(:,i)*wbts1(:,j)
                        goto 70
                     endif
                     do ii=nsepp-1,1,-1
                        if(rp2>sepp2(ii)) then
                           cdpv(jj,ii) = cdpv(jj,ii) + 1.0d0
                           bcdpv(:,jj,ii) = bcdpv(:,jj,ii) + wbts(:,i)*wbts1(:,j)
                           goto 70
                        endif
                     enddo
                  endif
               else 
                  if(jq3>iq3) cycle lp_jq2  !if(rv>rvmax.and.jq3>iq3)
               endif
               70 j = ll1(j)
            end do
         end do lp_jq2
      end do lp_jq1
   end do lp_jq3
end do
!$omp end parallel do
close(11)  ! close log
write(*,*) ' '
end subroutine rppi_Cb


subroutine rppi_Cb_wg(nt,npt,ra,dec,dc,wei,x,y,z,npt1,dc1,wei1,x1,y1,z1, &
           nsepp,sepp,nsepv,sepv,sbound,mxh1,mxh2,mxh3,nbts,bseed,wfib,cntid, &
           logf,sk1,ll1,cdpv,bcdpv)
!===============================================================================
! NAME
!  rppi_Cb_wg()
!
! DESCRIPTION
!  Cross-count weighted data pairs in projected space for TWO samples of particles.
!  Also count boostraped data pairs.
!
! INPUTS
!  Variable---Type-----------Description----------------------------------------
!  nt         int            Number of thread to use (OpenMP)
!  npt        int            [sampleC] Number of particles
!  ra         r8(npt)        [sampleC] RA of particles [deg]
!  dec        r8(npt)        [sampleC] DEC of particles [deg]
!  dc         r8(npt)        [sampleC] DCOM of particles [Mpc/h]
!  wei        r4(npt)        [sampleC] WEIGHT of particles
!  x,y,z      r8(npt)        [sampleC] X,Y,Z coordinates of particles (see radec2xyz())
!  npt1       int            [sampleD] Number of particles
!  dc1        r8(npt1)       [sampleD] DCOM of particles [Mpc/h]
!  wei1       r4(npt1)       [sampleD] WEIGHT of particles
!  x1,y1,z1   r8(npt1)       [sampleD] X,Y,Z coordinates of particles (see radec2xyz())
!  nsepp      int            Number of projected separation bins
!  sepp       r8(nsepp+1)    Bins in projected separation [Mpc/h]
!  nsepv      int            Number of radial separation bins
!  sepv       r8(nsepv+1)    Bins in radial separation [Mpc/h]
!  sbound     r8(6)          Survey boundaries in RA,DEC,DCOM
!                            Form is (ramin,ramax,decmin,decmax,dcmin,dcmax)
!  mxh1       int            Nr of DEC cells of skip table
!  mxh2       int            Nr of RA cells of skip table
!  mxh3       int            Nr of DCOM cells of skip table
!  nbts       int            Number of boostrap samples
!  bseed      int            Seed for RNG during boostrap resampling
!  wfib       int            If 1, also apply fiber-correction weights
!  cntid      char*2         Two-character string to identify the samples,
!                            e.g. DD, RR, DR,... (used only for ouput messages)
!  logf       char*80        String for the file to log Fortran status messages
!  sk1        int(mxh3,mxh2,mxh1) [sampleD] Skip table (SK) contructed with skll2d()
!  ll1        int(npt1)      [sampleD] Linked list (LL) contructed with skll2d()
!
! OUTPUTS
!  Variable---Type-------------Description--------------------------------------
!  cdpv       r8(nsepv,nsepp)  Counts in radial and projected separation bins
!  bcdpv      r8(nbts,nsepv,nsepp) Boostrap counts in radial and projected separation bins
!
! NOTES  -----------------------------------------------------------------------
!  1. RA limits should be set to ramin=0 and ramax=360.
!  2. Remember to update the declarations in cflibfor.pyf if you add/remove 
!     in/out parameters to this exposed subroutine.

implicit none
real(kind=8)  :: ra(npt),dec(npt),sbound(6),ral,rau,decl,decu,dcoml,dcomu,hc1,hc2,hc3,rp2
real(kind=8)  :: x(npt),y(npt),z(npt),x1(npt1),y1(npt1),z1(npt1),rcl,stm2,dltdec,dltra
real(kind=8)  :: dc(npt),dc1(npt1),xi,dci
real(kind=4)  :: wei(npt),wei1(npt1),wi,wpp
real(kind=4)  :: wbts(nbts,npt),wbts1(nbts,npt1)
real(kind=8)  :: sepp(nsepp+1),sepv(nsepv+1),sepp2(nsepp+1),rpmax,rpmax2,rvmax,rv,idsepv,shth2
real(kind=8)  :: cdpv(nsepv,nsepp),bcdpv(nbts,nsepv,nsepp)
integer       :: sk1(mxh3,mxh2,mxh1),ll1(npt1),mxh1,mxh2,mxh3,npt,npt1,nsepp,nsepv,nbts,bseed,wfib
integer       :: nt,nthr,ndp,fracp,dpart,nadv,nc1,nc2,nc3,jq1m,jq2m
integer       :: i,ii,j,jj,iq1,iq2,iq3,jq1,jq2,jq3,jq2min,jq2max,jq2t,p1,p2
character     :: cntid*2,logf*80

!------Open log file---------------------------------
open(11,file=logf,action='write',access='append')

fracp = 0 ; dpart = 0 ; nadv = 0  !Reset progress counters

!----------------------------------------------------
! Reset the counts, set (rpmax,rvmax) and square rp bins
cdpv   = 0.d0
bcdpv  = 0.0
rpmax  = sepp(nsepp+1)
rvmax  = sepv(nsepv+1)
sepp2  = sepp*sepp
rpmax2 = sepp2(nsepp+1)
idsepv = 1./(sepv(2)-sepv(1))  !inverse radial bin size

!----------------------------------------------------
! Unpack boundaries and get cell dimensions in (dec,ra,dcom)
ral = sbound(1) ; rau = sbound(2) ; decl = sbound(3) ; decu = sbound(4)
dcoml = sbound(5) ; dcomu = sbound(6)
nc1 = mxh1
nc2 = mxh2
nc3 =int((dcomu-dcoml)/rvmax)
if(nc3>mxh3) nc3 = mxh3
hc1 = (decu-decl)/float(nc1)    !effective nr of DEC cells
hc2 = (rau-ral)/float(nc2)      !effective nr of RA cells
hc3 = (dcomu-dcoml)/float(nc3)  !effective nr of DCOM cells

!----------------------------------------------------
! Generate bootstrap samples
call bootstrap(npt,nbts,bseed,wbts)
call bootstrap(npt1,nbts,bseed,wbts1)

!----------------------------------------------------
! Set the number of threads
if(nt<=0) then
   nthr = omp_get_num_procs()
else
   nthr = nt
endif
call omp_set_num_threads(nthr)

!----------------------------------------------------
! Some inits
dpart = int(float(npt)/float(mxh1)) !choose dpart so we get mhx1 parts

!----------------------------------------------------
! Count pairs in SK grid
write(11,*) 'Looping...'
write(*,*) ' '
write(*,fmt='(a,i3,a)') '====  Counting '//cntid//' pairs in ', mxh1, ' DEC strips'
!$omp parallel do reduction(+:cdpv,bcdpv) default(shared) &
!$omp& private(iq1,iq2,iq3,jq1,jq2,jq3,rcl,stm2,dltdec,jq1m) &
!$omp& private(dltra,jq2m,jq2min,jq2max,jq2t,j,i,xi,wi,dci,rv,rp2,shth2,wpp,ii,jj) &
!$omp& schedule(guided) firstprivate(fracp) if(nthr>1)

do i=1,npt   !---- Loop over ith particles ----
   xi  = x(i)
   dci = dc(i)
   wi  = wei(i)
   
   iq1 = int((dec(i)-decl)/hc1)+1
   iq2 = int((ra(i)-ral)/hc2)+1
   iq3 = int((dci-dcoml)/hc3)+1

   fracp = fracp + 1  ! accumulate particles and check when above the step size
   if(fracp>=dpart) then
       !$omp critical
       nadv = nadv + 1
       !omp flush (nadv)
       p1   = (nadv-1)*dpart + 1
       p2   = nadv*dpart
       if((nadv+1)*dpart>npt) p2=npt
       ! Note we are not really counting in mxh1 strips, just mymicking
       write(*,fmt="(i4)",advance='no') nadv                        ! for screen
       write(11,*) cntid//' counting in DEC strip > ',nadv,' (',p1,'-',p2,')' !for disk
       fracp = 0
       !$omp end critical
   endif
   
   lp_jq3: do jq3 = iq3-1,iq3+1
      if(jq3<1.or.jq3>nc3) cycle lp_jq3
      rcl    = (jq3-1)*hc3+dcoml
      stm2   = sthmax2(rpmax,dci,rcl)
      dltdec = 2.0d0*asin(stm2)*rad2deg
      jq1m   = int(dltdec/hc1)+1
      lp_jq1: do jq1=iq1-jq1m,iq1+jq1m
         if(jq1>nc1.or.jq1<1) cycle lp_jq1
         if(jq1==iq1) then
            dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,1)
         else
            dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,0)
         end if
         jq2m   = int(dltra/hc2)+1
         jq2max = iq2+jq2m
         jq2min = iq2-jq2m
         if(jq2max-jq2min+1>nc2) jq2max = jq2min-1+nc2
         lp_jq2: do jq2=jq2min,jq2max
            if(jq2>nc2) then
               jq2t = jq2-nc2
            else if(jq2<1) then
               jq2t = jq2+nc2
            else
               jq2t = jq2
            end if
            j = sk1(jq3,jq2t,jq1)
            do while(j/=0)
               !---- Loop over jth particles ----
               ! For each ith particle, index j runs through all particles in
               ! neighboring cells. Prune those outside rv and rp2 range.
               rv  = abs(dci-dc1(j))
               if(rv<=rvmax) then
                  rp2 = 4.*dci*dc1(j)*((xi-x1(j))**2 + (y(i)-y1(j))**2 + (z(i)-z1(j))**2)
                  if(rp2<=rpmax2) then
                     jj = 1 + int(rv*idsepv)             !find rv bin number
                     wpp = wi*wei1(j)                    !weight by input
                     if(wfib==1) then                    !weight by fiber
                        shth2 = rp2/(4.*dci*dc1(j))
                        wpp   = wpp*wfiber(shth2)                        
                     endif
                     if(rp2>sepp2(nsepp)) then 
                        cdpv(jj,nsepp) = cdpv(jj,nsepp) + wpp
                        bcdpv(:,jj,nsepp) = bcdpv(:,jj,nsepp) + wpp*wbts(:,i)*wbts1(:,j)
                        goto 70
                     endif
                     do ii=nsepp-1,1,-1
                        if(rp2>sepp2(ii)) then
                           cdpv(jj,ii) = cdpv(jj,ii) + wpp
                           bcdpv(:,jj,ii) = bcdpv(:,jj,ii) + wpp*wbts(:,i)*wbts1(:,j)
                           goto 70
                        endif
                     enddo
                  endif
               else 
                  if(jq3>iq3) cycle lp_jq2  !if(rv>rvmax.and.jq3>iq3)
               endif
               70 j = ll1(j)
            end do
         end do lp_jq2
      end do lp_jq1
   end do lp_jq3
end do
!$omp end parallel do
close(11)  ! close log
write(*,*) ' '
end subroutine rppi_Cb_wg


!==========================================================================================
!   COUNTING ROUTINES - 3D SPHERICAL (redshift-space)
!==========================================================================================

subroutine s_A(nt,npt,dec,dc,x,y,z,nseps,seps,sbound,mxh1,mxh2,mxh3,cntid,logf,sk,ll,aas)
!===============================================================================
! NAME
!  s_A()
!
! DESCRIPTION
!  Count data pairs in redshift space for ONE sample of particles
!
! INPUTS
!  Variable---Type-----------Description----------------------------------------
!  nt         int            Number of thread to use (OpenMP)
!  npt        int            Number of particles
!  dec        r8(npt)        DEC of particles [deg]
!  dc         r8(npt)        DCOM of particles [Mpc/h]
!  x,y,z      r8(npt)        X,Y,Z coordinates of particles (see radec2xyz())
!  nseps      int            Number of redshift-space separation bins
!  seps       r8(nseps+1)    Bins in redshift-space separation [Mpc/h]
!  sbound     r8(6)          Survey boundaries in RA,DEC,DCOM
!                            Form is (ramin,ramax,decmin,decmax,dcmin,dcmax)
!  mxh1       int            Nr of DEC cells of skip table
!  mxh2       int            Nr of RA cells of skip table
!  mxh3       int            Nr of DCOM cells of skip table
!  cntid      char*2         Two-character string to identify the samples,
!                            e.g. DD, RR, DR,... (used only for ouput messages)
!  logf       char*80        String for the file to log Fortran status messages
!  sk         int(mxh3,mxh2,mxh1)  Skip table (SK) contructed with skll2d()
!  ll         int(npt)       Linked list (LL) contructed with skll2d()
!
! OUTPUTS
!  Variable---Type-------------Description--------------------------------------
!  aas        r8(nseps)      Counts in redshift-space separation bins
!
! NOTES  -----------------------------------------------------------------------
!  1. RA limits should be set to ramin=0 and ramax=360.
!  2. Remember to update the declarations in cflibfor.pyf if you add/remove 
!     in/out parameters to this exposed subroutine.

implicit none
real(kind=8)  :: dec(npt),sbound(6),ral,rau,decl,decu,dcoml,dcomu,hc1,hc2,hc3,rp2,r2
real(kind=8)  :: x(npt),y(npt),z(npt),dc(npt),rcl,stm2,dltdec,dltra,dci,xi,yi,zi
real(kind=8)  :: seps(nseps+1),seps2(nseps+1),rsmax,rsmax2,rv
real(kind=8)  :: aas(nseps)
integer       :: sk(mxh3,mxh2,mxh1),ll(npt),mxh1,mxh2,mxh3,npt,nseps
integer       :: nt,nthr,nc1,nc2,nc3,jq1m,jq2m
integer       :: i,ii,j,jj,iq1,iq2,iq3,jq1,jq2,jq3,jq1min,jq2min,jq2max,jq2t
character     :: cntid*2,logf*80

!------Open log file---------------------------------
open(11,file=logf,action='write',access='append')

! Some useful print statements. Place where desired
!print *, mxh1,mxh2,mxh3
!print *, nc1,nc2,nc3
!print *, hc1, hc2, hc3
!print *, omp_get_num_threads(), omp_get_max_threads(), omp_get_num_procs()

!----------------------------------------------------
! Reset the counts, set rsmax and square s bins
aas    = 0.0d0
rsmax  = seps(nseps+1)
seps2  = seps*seps
rsmax2 = seps2(nseps+1)

!----------------------------------------------------
! Unpack boundaries and get cell dimensions in (dec,ra,dcom)
ral = sbound(1) ; rau = sbound(2) ; decl = sbound(3) ; decu = sbound(4)
dcoml = sbound(5) ; dcomu = sbound(6)
nc1 = mxh1
nc2 = mxh2
nc3 = int((dcomu-dcoml)/rsmax)
if(nc3>mxh3) nc3 = mxh3
hc1 = (decu-decl)/float(nc1)     !effective nr of DEC cells
hc2 = (rau-ral)/float(nc2)       !effective nr of RA cells
hc3 = (dcomu-dcoml)/float(nc3)   !effective nr of DCOM cells

!----------------------------------------------------
! Set the number of threads
if(nt<=0) then
   nthr = omp_get_num_procs()
else
   nthr = nt
endif
call omp_set_num_threads(nthr)

!----------------------------------------------------
! Count pairs in SK grid
write(*,*) ' '
write(*,fmt='(a,i3,a)') '====  Counting '//cntid//' pairs in ', mxh1, ' DEC strips  ===='
!$omp parallel do reduction(+:aas) default(shared) &
!$omp& private(iq1,iq2,iq3,jq1,jq2,jq3,rcl,stm2,dltdec,jq1m,jq1min) &
!$omp& private(dltra,jq2m,jq2min,jq2max,jq2t,j,i,xi,yi,zi,dci,rv,rp2,r2,ii) &
!$omp& schedule(guided) if(nthr>1)
do iq1=1,nc1
   !$omp critical
   write(*,fmt="(i4)",advance='no') iq1                             ! for screen
   write(11,*) cntid//' counting in DEC strip > ',iq1,'/',mxh1      ! for disk
   flush(11)
   !$omp end critical
   do iq2=1,nc2
      do iq3=1,nc3
         i = sk(iq3,iq2,iq1) ! index of ith particle 
         do while(i/=0) !---- Loop over ith particles ----
            lp_jq3: do jq3=iq3,iq3+1
               if(jq3>nc3) cycle lp_jq3
               rcl    = (jq3-1)*hc3+dcoml
               stm2   = sthmax2(rsmax,dc(i),rcl)
               dltdec = 2.0d0*asin(stm2)*rad2deg
               jq1m   = int(dltdec/hc1)+1 !dltdec/hc1+1.
               if(jq3==iq3) then
                  jq1min = iq1
               else
                  jq1min = iq1-jq1m
               end if
               lp_jq1: do jq1=jq1min,iq1+jq1m
                  if(jq1>nc1.or.jq1<1) cycle lp_jq1
                  if(jq1==iq1) then
                     dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,1)
                  else
                     dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,0)
                  end if
                  jq2m   = int(dltra/hc2)+1  !dltra/hc2+1.
                  jq2max = iq2+jq2m
                  jq2min = iq2-jq2m
                  if(jq2max-jq2min+1>nc2) jq2max = jq2min-1+nc2
                  lp_jq2: do jq2=jq2min,jq2max
                     if(jq2>nc2) then
                        jq2t = jq2-nc2
                     else if(jq2<1) then
                        jq2t = jq2+nc2
                     else
                        jq2t = jq2
                     end if
                     if(jq3==iq3.and.jq2t<iq2.and.jq1==iq1) cycle lp_jq2
                     if(jq3==iq3.and.jq2t==iq2.and.jq1==iq1) then
                        j = ll(i)
                     else
                        j = sk(jq3,jq2t,jq1)
                     endif
                     xi  = x(i)
                     yi  = y(i)
                     zi  = z(i)
                     dci = dc(i)

                     do while(j/=0)
                        rv = abs(dci-dc(j))
                        rp2 = 4.*dci*dc(j)*((xi-x(j))**2 + (yi-y(j))**2 + (zi-z(j))**2)
                        r2  = rv*rv + rp2
                        if(r2<=rsmax2) then
                           if(r2>seps2(nseps)) then
                              aas(nseps) = aas(nseps) + 1.0d0
                              goto 74
                           endif
                           if(r2>seps2(nseps-1)) then
                              aas(nseps-1) = aas(nseps-1) + 1.0d0
                              goto 74
                           endif
                           do ii=nseps-2,1,-1
                              if(r2>seps2(ii)) then
                                 aas(ii) = aas(ii) + 1.0d0
                                 goto 74
                              endif
                           enddo
                        else
                           if(rv>rsmax.and.jq3>iq3) cycle lp_jq2
                        endif
                        74 j = ll(j)
                     end do
                  end do lp_jq2
               end do lp_jq1
            end do lp_jq3
            i = ll(i)
         end do
      end do 
   end do 
end do
!$omp end parallel do
close(11)  ! close log
write(*,*) ' '
end subroutine s_A


subroutine s_A_wg(nt,npt,dec,dc,wei,x,y,z,nseps,seps,sbound,mxh1,mxh2,mxh3,& 
                  wfib,cntid,logf,sk,ll,aas)
!===============================================================================
! NAME
!  s_A_wg()
!
! DESCRIPTION
!  Count weighted data pairs in redshift space for ONE sample of particles
!
! INPUTS
!  Variable---Type-----------Description----------------------------------------
!  nt         int            Number of thread to use (OpenMP)
!  npt        int            Number of particles
!  dec        r8(npt)        DEC of particles [deg]
!  dc         r8(npt)        DCOM of particles [Mpc/h]
!  wei        r4(npt)        WEIGHT of particles
!  x,y,z      r8(npt)        X,Y,Z coordinates of particles (see radec2xyz())
!  nseps      int            Number of redshift-space separation bins
!  seps       r8(nseps+1)    Bins in redshift-space separation [Mpc/h]
!  sbound     r8(6)          Survey boundaries in RA,DEC,DCOM
!                            Form is (ramin,ramax,decmin,decmax,dcmin,dcmax)
!  mxh1       int            Nr of DEC cells of skip table
!  mxh2       int            Nr of RA cells of skip table
!  mxh3       int            Nr of DCOM cells of skip table
!  wfib       int            If 1, also apply fiber-correction weights
!  cntid      char*2         Two-character string to identify the samples,
!                            e.g. DD, RR, DR,... (used only for ouput messages)
!  logf       char*80        String for the file to log Fortran status messages
!  sk         int(mxh3,mxh2,mxh1)  Skip table (SK) contructed with skll2d()
!  ll         int(npt)       Linked list (LL) contructed with skll2d()
!
! OUTPUTS
!  Variable---Type-------------Description--------------------------------------
!  aas        r8(nseps)      Counts in redshift-space separation bins
!
! NOTES  -----------------------------------------------------------------------
!  1. RA limits should be set to ramin=0 and ramax=360.
!  2. Remember to update the declarations in cflibfor.pyf if you add/remove 
!     in/out parameters to this exposed subroutine.

implicit none
real(kind=8)  :: dec(npt),sbound(6),ral,rau,decl,decu,dcoml,dcomu,hc1,hc2,hc3,rp2,r2
real(kind=4)  :: wei(npt),wpp,wi
real(kind=8)  :: x(npt),y(npt),z(npt),dc(npt),rcl,stm2,dltdec,dltra,dci,xi,yi,zi
real(kind=8)  :: seps(nseps+1),seps2(nseps+1),rsmax,rsmax2,rv,shth2
real(kind=8)  :: aas(nseps)
integer       :: sk(mxh3,mxh2,mxh1),ll(npt),mxh1,mxh2,mxh3,npt,nseps
integer       :: nt,nthr,nc1,nc2,nc3,jq1m,jq2m
integer       :: i,ii,j,jj,iq1,iq2,iq3,jq1,jq2,jq3,jq1min,jq2min,jq2max,jq2t,wfib
character     :: cntid*2,logf*80

!------Open log file---------------------------------
open(11,file=logf,action='write',access='append')

! Some useful print statements. Place where desired
!print *, mxh1,mxh2,mxh3
!print *, nc1,nc2,nc3
!print *, hc1, hc2, hc3
!print *, omp_get_num_threads(), omp_get_max_threads(), omp_get_num_procs()

!----------------------------------------------------
! Reset the counts, set rsmax and square s bins
aas    = 0.0d0
rsmax  = seps(nseps+1)
seps2  = seps*seps
rsmax2 = seps2(nseps+1)

!----------------------------------------------------
! Unpack boundaries and get cell dimensions in (dec,ra,dcom)
ral = sbound(1) ; rau = sbound(2) ; decl = sbound(3) ; decu = sbound(4)
dcoml = sbound(5) ; dcomu = sbound(6)
nc1 = mxh1
nc2 = mxh2
nc3 = int((dcomu-dcoml)/rsmax)
if(nc3>mxh3) nc3 = mxh3
hc1 = (decu-decl)/float(nc1)     !effective nr of DEC cells
hc2 = (rau-ral)/float(nc2)       !effective nr of RA cells
hc3 = (dcomu-dcoml)/float(nc3)   !effective nr of DCOM cells

!----------------------------------------------------
! Set the number of threads
if(nt<=0) then
   nthr = omp_get_num_procs()
else
   nthr = nt
endif
call omp_set_num_threads(nthr)

!----------------------------------------------------
! Count pairs in SK grid
write(*,*) ' '
write(*,fmt='(a,i3,a)') '====  Counting '//cntid//' pairs in ', mxh1, ' DEC strips  ===='
!$omp parallel do reduction(+:aas) default(shared) &
!$omp& private(iq1,iq2,iq3,jq1,jq2,jq3,rcl,stm2,dltdec,jq1m,jq1min) &
!$omp& private(dltra,jq2m,jq2min,jq2max,jq2t,j,i,xi,yi,zi,dci,rv,rp2,r2,ii,wi,wpp,shth2) &
!$omp& schedule(guided) if(nthr>1)
do iq1=1,nc1
   !$omp critical
   write(*,fmt="(i4)",advance='no') iq1                             ! for screen
   write(11,*) cntid//' counting in DEC strip > ',iq1,'/',mxh1      ! for disk
   flush(11)
   !$omp end critical
   do iq2=1,nc2
      do iq3=1,nc3
         i = sk(iq3,iq2,iq1) ! index of ith particle 
         do while(i/=0) !---- Loop over ith particles ----
            lp_jq3: do jq3=iq3,iq3+1
               if(jq3>nc3) cycle lp_jq3
               rcl    = (jq3-1)*hc3+dcoml
               stm2   = sthmax2(rsmax,dc(i),rcl)
               dltdec = 2.0d0*asin(stm2)*rad2deg
               jq1m   = int(dltdec/hc1)+1 !dltdec/hc1+1.
               if(jq3==iq3) then
                  jq1min = iq1
               else
                  jq1min = iq1-jq1m
               end if
               lp_jq1: do jq1=jq1min,iq1+jq1m
                  if(jq1>nc1.or.jq1<1) cycle lp_jq1
                  if(jq1==iq1) then
                     dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,1)
                  else
                     dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,0)
                  end if
                  jq2m   = int(dltra/hc2)+1  !dltra/hc2+1.
                  jq2max = iq2+jq2m
                  jq2min = iq2-jq2m
                  if(jq2max-jq2min+1>nc2) jq2max = jq2min-1+nc2
                  lp_jq2: do jq2=jq2min,jq2max
                     if(jq2>nc2) then
                        jq2t = jq2-nc2
                     else if(jq2<1) then
                        jq2t = jq2+nc2
                     else
                        jq2t = jq2
                     end if
                     if(jq3==iq3.and.jq2t<iq2.and.jq1==iq1) cycle lp_jq2
                     if(jq3==iq3.and.jq2t==iq2.and.jq1==iq1) then
                        j = ll(i)
                     else
                        j = sk(jq3,jq2t,jq1)
                     endif
                     xi  = x(i)
                     yi  = y(i)
                     zi  = z(i)
                     dci = dc(i)
                     wi  = wei(i)

                     do while(j/=0)
                        rv = abs(dci-dc(j))
                        rp2 = 4.*dci*dc(j)*((xi-x(j))**2 + (yi-y(j))**2 + (zi-z(j))**2)
                        r2  = rv*rv + rp2
                        if(r2<=rsmax2) then
                           wpp = wi*wei(j)               !weight by input
                           if(wfib==1) then              !weight by fiber
                               shth2 = rp2/(4.*dci*dc(j))
                               wpp   = wpp*wfiber(shth2)
                           endif
                           if(r2>seps2(nseps)) then
                              aas(nseps) = aas(nseps) + wpp
                              goto 74
                           endif
                           if(r2>seps2(nseps-1)) then
                              aas(nseps-1) = aas(nseps-1) + wpp
                              goto 74
                           endif
                           do ii=nseps-2,1,-1
                              if(r2>seps2(ii)) then
                                 aas(ii) = aas(ii) + wpp
                                 goto 74
                              endif
                           enddo
                        else
                           if(rv>rsmax.and.jq3>iq3) cycle lp_jq2
                        endif
                        74 j = ll(j)
                     end do
                  end do lp_jq2
               end do lp_jq1
            end do lp_jq3
            i = ll(i)
         end do
      end do 
   end do 
end do
!$omp end parallel do
close(11)  ! close log
write(*,*) ' '
end subroutine s_A_wg


subroutine s_Ab(nt,npt,dec,dc,x,y,z,nseps,seps,sbound,mxh1,mxh2,mxh3,nbts,bseed, &
                cntid,logf,sk,ll,aas,baas)
!===============================================================================
! NAME
!  s_Ab()
!
! DESCRIPTION
!  Count data pairs in redshift space for ONE sample of particles. Also count
!  boostraped data pairs
!
! INPUTS
!  Variable---Type-----------Description----------------------------------------
!  nt         int            Number of thread to use (OpenMP)
!  npt        int            Number of particles
!  dec        r8(npt)        DEC of particles [deg]
!  dc         r8(npt)        DCOM of particles [Mpc/h]
!  x,y,z      r8(npt)        X,Y,Z coordinates of particles (see radec2xyz())
!  nseps      int            Number of redshift-space separation bins
!  seps       r8(nseps+1)    Bins in redshift-space separation [Mpc/h]
!  sbound     r8(6)          Survey boundaries in RA,DEC,DCOM
!                            Form is (ramin,ramax,decmin,decmax,dcmin,dcmax)
!  mxh1       int            Nr of DEC cells of skip table
!  mxh2       int            Nr of RA cells of skip table
!  mxh3       int            Nr of DCOM cells of skip table
!  nbts       int            Number of boostrap samples
!  bseed      int            Seed for RNG during boostrap resampling
!  cntid      char*2         Two-character string to identify the samples,
!                            e.g. DD, RR, DR,... (used only for ouput messages)
!  logf       char*80        String for the file to log Fortran status messages
!  sk         int(mxh3,mxh2,mxh1)  Skip table (SK) contructed with skll2d()
!  ll         int(npt)       Linked list (LL) contructed with skll2d()
!
! OUTPUTS
!  Variable---Type-------------Description--------------------------------------
!  aas        r8(nseps)      Counts in redshift-space separation bins
!  baas       r8(nbts,nseps) Booostraped counts in redshift-space separation bins
!
! NOTES  -----------------------------------------------------------------------
!  1. RA limits should be set to ramin=0 and ramax=360.
!  2. Remember to update the declarations in cflibfor.pyf if you add/remove 
!     in/out parameters to this exposed subroutine.

implicit none
real(kind=8)  :: dec(npt),sbound(6),ral,rau,decl,decu,dcoml,dcomu,hc1,hc2,hc3,rp2,r2
real(kind=8)  :: x(npt),y(npt),z(npt),dc(npt),rcl,stm2,dltdec,dltra,dci,xi,yi,zi
real(kind=8)  :: seps(nseps+1),seps2(nseps+1),rsmax,rsmax2,rv
real(kind=8)  :: aas(nseps),baas(nbts,nseps)
real(kind=4)  :: wbts(nbts,npt)
integer       :: sk(mxh3,mxh2,mxh1),ll(npt),mxh1,mxh2,mxh3,npt,nseps
integer       :: nt,nthr,nc1,nc2,nc3,jq1m,jq2m,nbts,bseed
integer       :: i,ii,j,jj,iq1,iq2,iq3,jq1,jq2,jq3,jq1min,jq2min,jq2max,jq2t
character     :: cntid*2,logf*80

!------Open log file---------------------------------
open(11,file=logf,action='write',access='append')

! Some useful print statements. Place where desired
!print *, mxh1,mxh2,mxh3
!print *, nc1,nc2,nc3
!print *, hc1, hc2, hc3
!print *, omp_get_num_threads(), omp_get_max_threads(), omp_get_num_procs()

!----------------------------------------------------
! Reset the counts, set rsmax and square s bins
aas    = 0.0d0
baas   = 0.0d0
rsmax  = seps(nseps+1)
seps2  = seps*seps
rsmax2 = seps2(nseps+1)

!----------------------------------------------------
! Unpack boundaries and get cell dimensions in (dec,ra,dcom)
ral = sbound(1) ; rau = sbound(2) ; decl = sbound(3) ; decu = sbound(4)
dcoml = sbound(5) ; dcomu = sbound(6)
nc1 = mxh1
nc2 = mxh2
nc3 = int((dcomu-dcoml)/rsmax)
if(nc3>mxh3) nc3 = mxh3
hc1 = (decu-decl)/float(nc1)     !effective nr of DEC cells
hc2 = (rau-ral)/float(nc2)       !effective nr of RA cells
hc3 = (dcomu-dcoml)/float(nc3)   !effective nr of DCOM cells

!----------------------------------------------------
! Generate bootstrap samples
call bootstrap(npt,nbts,bseed,wbts)

!----------------------------------------------------
! Set the number of threads
if(nt<=0) then
   nthr = omp_get_num_procs()
else
   nthr = nt
endif
call omp_set_num_threads(nthr)

!----------------------------------------------------
! Count pairs in SK grid
write(*,*) ' '
write(*,fmt='(a,i3,a)') '====  Counting '//cntid//' pairs in ', mxh1, ' DEC strips  ===='
!$omp parallel do reduction(+:aas,baas) default(shared) &
!$omp& private(iq1,iq2,iq3,jq1,jq2,jq3,rcl,stm2,dltdec,jq1m,jq1min) &
!$omp& private(dltra,jq2m,jq2min,jq2max,jq2t,j,i,xi,yi,zi,dci,rv,rp2,r2,ii) &
!$omp& schedule(guided) if(nthr>1)
do iq1=1,nc1
   !$omp critical
   write(*,fmt="(i4)",advance='no') iq1                             ! for screen
   write(11,*) cntid//' counting in DEC strip > ',iq1,'/',mxh1      ! for disk
   flush(11)
   !$omp end critical
   do iq2=1,nc2
      do iq3=1,nc3
         i = sk(iq3,iq2,iq1) ! index of ith particle 
         do while(i/=0) !---- Loop over ith particles ----
            lp_jq3: do jq3=iq3,iq3+1
               if(jq3>nc3) cycle lp_jq3
               rcl    = (jq3-1)*hc3+dcoml
               stm2   = sthmax2(rsmax,dc(i),rcl)
               dltdec = 2.0d0*asin(stm2)*rad2deg
               jq1m   = int(dltdec/hc1)+1 !dltdec/hc1+1.
               if(jq3==iq3) then
                  jq1min = iq1
               else
                  jq1min = iq1-jq1m
               end if
               lp_jq1: do jq1=jq1min,iq1+jq1m
                  if(jq1>nc1.or.jq1<1) cycle lp_jq1
                  if(jq1==iq1) then
                     dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,1)
                  else
                     dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,0)
                  end if
                  jq2m   = int(dltra/hc2)+1  !dltra/hc2+1.
                  jq2max = iq2+jq2m
                  jq2min = iq2-jq2m
                  if(jq2max-jq2min+1>nc2) jq2max = jq2min-1+nc2
                  lp_jq2: do jq2=jq2min,jq2max
                     if(jq2>nc2) then
                        jq2t = jq2-nc2
                     else if(jq2<1) then
                        jq2t = jq2+nc2
                     else
                        jq2t = jq2
                     end if
                     if(jq3==iq3.and.jq2t<iq2.and.jq1==iq1) cycle lp_jq2
                     if(jq3==iq3.and.jq2t==iq2.and.jq1==iq1) then
                        j = ll(i)
                     else
                        j = sk(jq3,jq2t,jq1)
                     endif
                     xi  = x(i)
                     yi  = y(i)
                     zi  = z(i)
                     dci = dc(i)

                     do while(j/=0)
                        rv = abs(dci-dc(j))
                        rp2 = 4.*dci*dc(j)*((xi-x(j))**2 + (yi-y(j))**2 + (zi-z(j))**2)
                        r2  = rv*rv + rp2
                        if(r2<=rsmax2) then
                           if(r2>seps2(nseps)) then
                              aas(nseps) = aas(nseps) + 1.0d0
                              baas(:,nseps) = baas(:,nseps) + wbts(:,i)*wbts(:,j)
                              goto 74
                           endif
                           if(r2>seps2(nseps-1)) then
                              aas(nseps-1) = aas(nseps-1) + 1.0d0
                              baas(:,nseps-1) = baas(:,nseps-1) + wbts(:,i)*wbts(:,j)
                              goto 74
                           endif
                           do ii=nseps-2,1,-1
                              if(r2>seps2(ii)) then
                                 aas(ii) = aas(ii) + 1.0d0
                                 baas(:,ii) = baas(:,ii) + wbts(:,i)*wbts(:,j)
                                 goto 74
                              endif
                           enddo
                        else
                           if(rv>rsmax.and.jq3>iq3) cycle lp_jq2
                        endif
                        74 j = ll(j)
                     end do
                  end do lp_jq2
               end do lp_jq1
            end do lp_jq3
            i = ll(i)
         end do
      end do 
   end do 
end do
!$omp end parallel do
close(11)  ! close log
write(*,*) ' '
end subroutine s_Ab


subroutine s_Ab_wg(nt,npt,dec,dc,wei,x,y,z,nseps,seps,sbound,mxh1,mxh2,mxh3,& 
                   nbts,bseed,wfib,cntid,logf,sk,ll,aas,baas)
!===============================================================================
! NAME
!  s_Ab_wg()
!
! DESCRIPTION
!  Count weighted data pairs in redshift space for ONE sample of particles. Also
!  count booostraped data pairs
!
! INPUTS
!  Variable---Type-----------Description----------------------------------------
!  nt         int            Number of thread to use (OpenMP)
!  npt        int            Number of particles
!  dec        r8(npt)        DEC of particles [deg]
!  dc         r8(npt)        DCOM of particles [Mpc/h]
!  wei        r4(npt)        WEIGHT of particles
!  x,y,z      r8(npt)        X,Y,Z coordinates of particles (see radec2xyz())
!  nseps      int            Number of redshift-space separation bins
!  seps       r8(nseps+1)    Bins in redshift-space separation [Mpc/h]
!  sbound     r8(6)          Survey boundaries in RA,DEC,DCOM
!                            Form is (ramin,ramax,decmin,decmax,dcmin,dcmax)
!  mxh1       int            Nr of DEC cells of skip table
!  mxh2       int            Nr of RA cells of skip table
!  mxh3       int            Nr of DCOM cells of skip table
!  nbts       int            Number of boostrap samples
!  bseed      int            Seed for RNG during boostrap resampling
!  wfib       int            If 1, also apply fiber-correction weights
!  cntid      char*2         Two-character string to identify the samples,
!                            e.g. DD, RR, DR,... (used only for ouput messages)
!  logf       char*80        String for the file to log Fortran status messages
!  sk         int(mxh3,mxh2,mxh1)  Skip table (SK) contructed with skll2d()
!  ll         int(npt)       Linked list (LL) contructed with skll2d()
!
! OUTPUTS
!  Variable---Type-------------Description--------------------------------------
!  aas        r8(nseps)      Counts in redshift-space separation bins
!  baas       r8(nbts,nseps) Booostraped counts in redshift-space separation bins
!
! NOTES  -----------------------------------------------------------------------
!  1. RA limits should be set to ramin=0 and ramax=360.
!  2. Remember to update the declarations in cflibfor.pyf if you add/remove 
!     in/out parameters to this exposed subroutine.

implicit none
real(kind=8)  :: dec(npt),sbound(6),ral,rau,decl,decu,dcoml,dcomu,hc1,hc2,hc3,rp2,r2
real(kind=4)  :: wei(npt),wpp,wi
real(kind=8)  :: x(npt),y(npt),z(npt),dc(npt),rcl,stm2,dltdec,dltra,dci,xi,yi,zi
real(kind=8)  :: seps(nseps+1),seps2(nseps+1),rsmax,rsmax2,rv,shth2
real(kind=8)  :: aas(nseps),baas(nbts,nseps)
real(kind=4)  :: wbts(nbts,npt)
integer       :: sk(mxh3,mxh2,mxh1),ll(npt),mxh1,mxh2,mxh3,npt,nseps
integer       :: nt,nthr,nc1,nc2,nc3,jq1m,jq2m,nbts,bseed
integer       :: i,ii,j,jj,iq1,iq2,iq3,jq1,jq2,jq3,jq1min,jq2min,jq2max,jq2t,wfib
character     :: cntid*2,logf*80

!------Open log file---------------------------------
open(11,file=logf,action='write',access='append')

! Some useful print statements. Place where desired
!print *, mxh1,mxh2,mxh3
!print *, nc1,nc2,nc3
!print *, hc1, hc2, hc3
!print *, omp_get_num_threads(), omp_get_max_threads(), omp_get_num_procs()

!----------------------------------------------------
! Reset the counts, set rsmax and square s bins
aas    = 0.0d0
baas   = 0.0d0
rsmax  = seps(nseps+1)
seps2  = seps*seps
rsmax2 = seps2(nseps+1)

!----------------------------------------------------
! Unpack boundaries and get cell dimensions in (dec,ra,dcom)
ral = sbound(1) ; rau = sbound(2) ; decl = sbound(3) ; decu = sbound(4)
dcoml = sbound(5) ; dcomu = sbound(6)
nc1 = mxh1
nc2 = mxh2
nc3 = int((dcomu-dcoml)/rsmax)
if(nc3>mxh3) nc3 = mxh3
hc1 = (decu-decl)/float(nc1)     !effective nr of DEC cells
hc2 = (rau-ral)/float(nc2)       !effective nr of RA cells
hc3 = (dcomu-dcoml)/float(nc3)   !effective nr of DCOM cells

!----------------------------------------------------
! Generate bootstrap samples
call bootstrap(npt,nbts,bseed,wbts)

!----------------------------------------------------
! Set the number of threads
if(nt<=0) then
   nthr = omp_get_num_procs()
else
   nthr = nt
endif
call omp_set_num_threads(nthr)

!----------------------------------------------------
! Count pairs in SK grid
write(*,*) ' '
write(*,fmt='(a,i3,a)') '====  Counting '//cntid//' pairs in ', mxh1, ' DEC strips  ===='
!$omp parallel do reduction(+:aas,baas) default(shared) &
!$omp& private(iq1,iq2,iq3,jq1,jq2,jq3,rcl,stm2,dltdec,jq1m,jq1min) &
!$omp& private(dltra,jq2m,jq2min,jq2max,jq2t,j,i,xi,yi,zi,dci,rv,rp2,r2,ii,wi,wpp,shth2) &
!$omp& schedule(guided) if(nthr>1)
do iq1=1,nc1
   !$omp critical
   write(*,fmt="(i4)",advance='no') iq1                             ! for screen
   write(11,*) cntid//' counting in DEC strip > ',iq1,'/',mxh1      ! for disk
   flush(11)
   !$omp end critical
   do iq2=1,nc2
      do iq3=1,nc3
         i = sk(iq3,iq2,iq1) ! index of ith particle 
         do while(i/=0) !---- Loop over ith particles ----
            lp_jq3: do jq3=iq3,iq3+1
               if(jq3>nc3) cycle lp_jq3
               rcl    = (jq3-1)*hc3+dcoml
               stm2   = sthmax2(rsmax,dc(i),rcl)
               dltdec = 2.0d0*asin(stm2)*rad2deg
               jq1m   = int(dltdec/hc1)+1 !dltdec/hc1+1.
               if(jq3==iq3) then
                  jq1min = iq1
               else
                  jq1min = iq1-jq1m
               end if
               lp_jq1: do jq1=jq1min,iq1+jq1m
                  if(jq1>nc1.or.jq1<1) cycle lp_jq1
                  if(jq1==iq1) then
                     dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,1)
                  else
                     dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,0)
                  end if
                  jq2m   = int(dltra/hc2)+1  !dltra/hc2+1.
                  jq2max = iq2+jq2m
                  jq2min = iq2-jq2m
                  if(jq2max-jq2min+1>nc2) jq2max = jq2min-1+nc2
                  lp_jq2: do jq2=jq2min,jq2max
                     if(jq2>nc2) then
                        jq2t = jq2-nc2
                     else if(jq2<1) then
                        jq2t = jq2+nc2
                     else
                        jq2t = jq2
                     end if
                     if(jq3==iq3.and.jq2t<iq2.and.jq1==iq1) cycle lp_jq2
                     if(jq3==iq3.and.jq2t==iq2.and.jq1==iq1) then
                        j = ll(i)
                     else
                        j = sk(jq3,jq2t,jq1)
                     endif
                     xi  = x(i)
                     yi  = y(i)
                     zi  = z(i)
                     dci = dc(i)
                     wi  = wei(i)

                     do while(j/=0)
                        rv = abs(dci-dc(j))
                        rp2 = 4.*dci*dc(j)*((xi-x(j))**2 + (yi-y(j))**2 + (zi-z(j))**2)
                        r2  = rv*rv + rp2
                        if(r2<=rsmax2) then
                           wpp = wi*wei(j)               !weight by input
                           if(wfib==1) then              !weight by fiber
                               shth2 = rp2/(4.*dci*dc(j))
                               wpp   = wpp*wfiber(shth2)
                           endif
                           if(r2>seps2(nseps)) then
                              aas(nseps) = aas(nseps) + wpp
                              baas(:,nseps) = baas(:,nseps) + wpp*wbts(:,i)*wbts(:,j)
                              goto 74
                           endif
                           if(r2>seps2(nseps-1)) then
                              aas(nseps-1) = aas(nseps-1) + wpp
                              baas(:,nseps-1) = baas(:,nseps-1) + wpp*wbts(:,i)*wbts(:,j)
                              goto 74
                           endif
                           do ii=nseps-2,1,-1
                              if(r2>seps2(ii)) then
                                 aas(ii) = aas(ii) + wpp
                                 baas(:,ii) = baas(:,ii) + wbts(:,i)*wbts(:,j)
                                 goto 74
                              endif
                           enddo
                        else
                           if(rv>rsmax.and.jq3>iq3) cycle lp_jq2
                        endif
                        74 j = ll(j)
                     end do
                  end do lp_jq2
               end do lp_jq1
            end do lp_jq3
            i = ll(i)
         end do
      end do 
   end do 
end do
!$omp end parallel do
close(11)  ! close log
write(*,*) ' '
end subroutine s_Ab_wg


subroutine s_C(nt,npt,ra,dec,dc,x,y,z,npt1,dc1,x1,y1,z1, &
               nseps,seps,sbound,mxh1,mxh2,mxh3,cntid,logf,sk1,ll1,cds)
!===============================================================================
! NAME
!  s_C()
!
! DESCRIPTION
!  Cross-count data pairs in redshift space for TWO samples of particles
!
! INPUTS
!  Variable---Type-----------Description----------------------------------------
!  nt         int            Number of thread to use (OpenMP)
!  npt        int            [sampleC]Number of particles
!  ra         r8(npt)        [sampleC] RA of particles [deg]
!  dec        r8(npt)        [sampleC]DEC of particles [deg]
!  dc         r8(npt)        [sampleC]DCOM of particles [Mpc/h]
!  x,y,z      r8(npt)        [sampleC]X,Y,Z coordinates of particles (see radec2xyz())
!  npt1       int            [sampleD] Number of particles
!  dc1        r8(npt1)       [sampleD] DCOM of particles [Mpc/h]
!  x1,y1,z1   r8(npt1)       [sampleD] X,Y,Z coordinates of particles (see radec2xyz())
!  nseps      int            Number of redshift-space separation bins
!  seps       r8(nseps+1)    Bins in redshift-space separation [Mpc/h]
!  sbound     r8(6)          Survey boundaries in RA,DEC,DCOM
!                            Form is (ramin,ramax,decmin,decmax,dcmin,dcmax)
!  mxh1       int            Nr of DEC cells of skip table
!  mxh2       int            Nr of RA cells of skip table
!  mxh3       int            Nr of DCOM cells of skip table
!  cntid      char*2         Two-character string to identify the samples,
!                            e.g. DD, RR, DR,... (used only for ouput messages)
!  logf       char*80        String for the file to log Fortran status messages
!  sk1        int(mxh3,mxh2,mxh1)  [sampleD] Skip table (SK) contructed with skll2d()
!  ll1        int(npt)       [sampleD] Linked list (LL) contructed with skll2d()
!
! OUTPUTS
!  Variable---Type-------------Description--------------------------------------
!  cds        r8(nseps)      Counts in redshift-space separation bins
!
! NOTES  -----------------------------------------------------------------------
!  1. RA limits should be set to ramin=0 and ramax=360.
!  2. Remember to update the declarations in cflibfor.pyf if you add/remove 
!     in/out parameters to this exposed subroutine.

implicit none
real(kind=8)  :: ra(npt),dec(npt),sbound(6),ral,rau,decl,decu,dcoml,dcomu,hc1,hc2,hc3,rp2,r2
real(kind=8)  :: x(npt),y(npt),z(npt),x1(npt1),y1(npt1),z1(npt1),rcl,stm2,dltdec,dltra
real(kind=8)  :: dc(npt),dc1(npt1),xi,dci
real(kind=8)  :: seps(nseps+1),seps2(nseps+1),rsmax,rsmax2,rv
real(kind=8)  :: cds(nseps)
integer       :: sk1(mxh3,mxh2,mxh1),ll1(npt1),mxh1,mxh2,mxh3,npt,npt1,nseps
integer       :: nt,nthr,fracp,dpart,nadv,nc1,nc2,nc3,jq1m,jq2m
integer       :: i,ii,j,jj,iq1,iq2,iq3,jq1,jq2,jq3,jq2min,jq2max,jq2t,p1,p2
character     :: cntid*2,logf*80

!------Open log file---------------------------------
open(11,file=logf,action='write',access='append')

fracp = 0 ; dpart = 0 ; nadv = 0  !Reset progress counters

!----------------------------------------------------
! Reset the counts, set (rpmax,rvmax) and square rp bins
cds    = 0.0d0
rsmax  = seps(nseps+1)
seps2  = seps*seps
rsmax2 = seps2(nseps+1)

!----------------------------------------------------
! Unpack boundaries and get cell dimensions in (dec,ra,dcom)
ral = sbound(1) ; rau = sbound(2) ; decl = sbound(3) ; decu = sbound(4)
dcoml = sbound(5) ; dcomu = sbound(6)
nc1 = mxh1
nc2 = mxh2
nc3 =int((dcomu-dcoml)/rsmax)
if(nc3>mxh3) nc3 = mxh3
hc1 = (decu-decl)/float(nc1)    !effective nr of DEC cells
hc2 = (rau-ral)/float(nc2)      !effective nr of RA cells
hc3 = (dcomu-dcoml)/float(nc3)  !effective nr of DCOM cells

!----------------------------------------------------
! Set the number of threads
if(nt<=0) then
   nthr = omp_get_num_procs()
else
   nthr = nt
endif
call omp_set_num_threads(nthr)

!----------------------------------------------------
! Some inits
dpart = int(float(npt)/float(mxh1)) !choose dpart so we get mhx1 parts

!----------------------------------------------------
! Count pairs in SK grid
write(*,*) ' '
write(*,fmt='(a,i3,a)') '====  Counting '//cntid//' pairs in ', mxh1, ' DEC strips'
!$omp parallel do reduction(+:cds) default(shared) &
!$omp& private(iq1,iq2,iq3,jq1,jq2,jq3,rcl,stm2,dltdec,jq1m) &
!$omp& private(dltra,jq2m,jq2min,jq2max,jq2t,j,i,xi,dci,rv,rp2,r2,ii) &
!$omp& schedule(guided) firstprivate(fracp) if(nthr>1)

do i=1,npt   !---- Loop over ith particles ----
   xi  = x(i)
   dci = dc(i)
   
   iq1 = int((dec(i)-decl)/hc1)+1
   iq2 = int((ra(i)-ral)/hc2)+1
   iq3 = int((dci-dcoml)/hc3)+1

   fracp = fracp + 1  ! accumulate particles and check when above the step size
   if(fracp>=dpart) then
       !$omp critical
       nadv = nadv + 1
       !omp flush (nadv)
       p1   = (nadv-1)*dpart + 1
       p2   = nadv*dpart
       if((nadv+1)*dpart>npt) p2=npt
       ! Note we are not really counting in mxh1 strips, just mymicking
       write(*,fmt="(i4)",advance='no') nadv                        ! for screen
       write(11,*) cntid//' counting in DEC strip > ',nadv,' (',p1,'-',p2,')' !for disk
       fracp = 0
       !$omp end critical
   endif
   
   lp_jq3: do jq3 = iq3-1,iq3+1
      if(jq3<1.or.jq3>nc3) cycle lp_jq3
      rcl    = (jq3-1)*hc3+dcoml
      stm2   = sthmax2(rsmax,dci,rcl)
      dltdec = 2.0d0*asin(stm2)*rad2deg
      jq1m   = int(dltdec/hc1)+1
      lp_jq1: do jq1=iq1-jq1m,iq1+jq1m
         if(jq1>nc1.or.jq1<1) cycle lp_jq1
         if(jq1==iq1) then
            dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,1)
         else
            dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,0)
         end if
         jq2m   = int(dltra/hc2)+1
         jq2max = iq2+jq2m
         jq2min = iq2-jq2m
         if(jq2max-jq2min+1>nc2) jq2max = jq2min-1+nc2
         lp_jq2: do jq2=jq2min,jq2max
            if(jq2>nc2) then
               jq2t = jq2-nc2
            else if(jq2<1) then
               jq2t = jq2+nc2
            else
               jq2t = jq2
            end if
            j = sk1(jq3,jq2t,jq1)
            do while(j/=0)
               !---- Loop over jth particles ----
               ! For each ith particle, index j runs through all particles in
               ! neighboring cells. Prune those outside rv and rp2 range.
               rv  = abs(dci-dc1(j))
               rp2 = 4.*dci*dc1(j)*((xi-x1(j))**2 + (y(i)-y1(j))**2 + (z(i)-z1(j))**2)
               r2  = rv*rv + rp2
               if(r2<=rsmax2) then
                  if(r2>seps2(nseps)) then
                     cds(nseps) = cds(nseps) + 1.0d0
                     goto 70
                  endif
                  do ii=nseps-1,1,-1
                     if(r2>seps2(ii)) then
                        cds(ii) = cds(ii) + 1.0d0
                        goto 70
                     endif
                  enddo
               else
                  if(rv>rsmax.and.jq3>iq3) cycle lp_jq2  !if(jq3>iq3) cycle lp_jq2
               endif
               70 j = ll1(j)
            end do
         end do lp_jq2
      end do lp_jq1
   end do lp_jq3
end do
!$omp end parallel do
close(11)  ! close log
write(*,*) ' '
end subroutine s_C


subroutine s_C_wg(nt,npt,ra,dec,dc,wei,x,y,z,npt1,dc1,wei1,x1,y1,z1, &
                  nseps,seps,sbound,mxh1,mxh2,mxh3,wfib,cntid,logf,sk1,ll1,cds)
!===============================================================================
! NAME
!  s_C_wg()
!
! DESCRIPTION
!  Cross-count weighted data pairs in redshift space for TWO samples of particles
!
! INPUTS
!  Variable---Type-----------Description----------------------------------------
!  nt         int            Number of thread to use (OpenMP)
!  npt        int            [sampleC]Number of particles
!  ra         r8(npt)        [sampleC] RA of particles [deg]
!  dec        r8(npt)        [sampleC]DEC of particles [deg]
!  dc         r8(npt)        [sampleC]DCOM of particles [Mpc/h]
!  wei        r4(npt)        [sampleC] WEIGHT of particles
!  x,y,z      r8(npt)        [sampleC]X,Y,Z coordinates of particles (see radec2xyz())
!  npt1       int            [sampleD] Number of particles
!  dc1        r8(npt1)       [sampleD] DCOM of particles [Mpc/h]
!  wei1       r4(npt1)       [sampleD] WEIGHT of particles
!  x1,y1,z1   r8(npt1)       [sampleD] X,Y,Z coordinates of particles (see radec2xyz())
!  nseps      int            Number of redshift-space separation bins
!  seps       r8(nseps+1)    Bins in redshift-space separation [Mpc/h]
!  sbound     r8(6)          Survey boundaries in RA,DEC,DCOM
!                            Form is (ramin,ramax,decmin,decmax,dcmin,dcmax)
!  mxh1       int            Nr of DEC cells of skip table
!  mxh2       int            Nr of RA cells of skip table
!  mxh3       int            Nr of DCOM cells of skip table
!  wfib       int            If 1, also apply fiber-correction weights
!  cntid      char*2         Two-character string to identify the samples,
!                            e.g. DD, RR, DR,... (used only for ouput messages)
!  logf       char*80        String for the file to log Fortran status messages
!  sk1        int(mxh3,mxh2,mxh1)  [sampleD] Skip table (SK) contructed with skll2d()
!  ll1        int(npt)       [sampleD] Linked list (LL) contructed with skll2d()
!
! OUTPUTS
!  Variable---Type-------------Description--------------------------------------
!  cds        r8(nseps)      Counts in redshift-space separation bins
!
! NOTES  -----------------------------------------------------------------------
!  1. RA limits should be set to ramin=0 and ramax=360.
!  2. Remember to update the declarations in cflibfor.pyf if you add/remove 
!     in/out parameters to this exposed subroutine.

implicit none
real(kind=8)  :: ra(npt),dec(npt),sbound(6),ral,rau,decl,decu,dcoml,dcomu,hc1,hc2,hc3,rp2,r2
real(kind=8)  :: x(npt),y(npt),z(npt),x1(npt1),y1(npt1),z1(npt1),rcl,stm2,dltdec,dltra
real(kind=8)  :: dc(npt),dc1(npt1),xi,dci
real(kind=4)  :: wei(npt),wei1(npt),wpp,wi
real(kind=8)  :: seps(nseps+1),seps2(nseps+1),rsmax,rsmax2,rv,shth2
real(kind=8)  :: cds(nseps)
integer       :: sk1(mxh3,mxh2,mxh1),ll1(npt1),mxh1,mxh2,mxh3,npt,npt1,nseps,wfib
integer       :: nt,nthr,fracp,dpart,nadv,nc1,nc2,nc3,jq1m,jq2m
integer       :: i,ii,j,jj,iq1,iq2,iq3,jq1,jq2,jq3,jq2min,jq2max,jq2t,p1,p2
character     :: cntid*2,logf*80

!------Open log file---------------------------------
open(11,file=logf,action='write',access='append')

fracp = 0 ; dpart = 0 ; nadv = 0  !Reset progress counters

!----------------------------------------------------
! Reset the counts, set (rpmax,rvmax) and square rp bins
cds    = 0.0d0
rsmax  = seps(nseps+1)
seps2  = seps*seps
rsmax2 = seps2(nseps+1)

!----------------------------------------------------
! Unpack boundaries and get cell dimensions in (dec,ra,dcom)
ral = sbound(1) ; rau = sbound(2) ; decl = sbound(3) ; decu = sbound(4)
dcoml = sbound(5) ; dcomu = sbound(6)
nc1 = mxh1
nc2 = mxh2
nc3 =int((dcomu-dcoml)/rsmax)
if(nc3>mxh3) nc3 = mxh3
hc1 = (decu-decl)/float(nc1)    !effective nr of DEC cells
hc2 = (rau-ral)/float(nc2)      !effective nr of RA cells
hc3 = (dcomu-dcoml)/float(nc3)  !effective nr of DCOM cells

!----------------------------------------------------
! Set the number of threads
if(nt<=0) then
   nthr = omp_get_num_procs()
else
   nthr = nt
endif
call omp_set_num_threads(nthr)

!----------------------------------------------------
! Some inits
dpart = int(float(npt)/float(mxh1)) !choose dpart so we get mhx1 parts

!----------------------------------------------------
! Count pairs in SK grid
write(*,*) ' '
write(*,fmt='(a,i3,a)') '====  Counting '//cntid//' pairs in ', mxh1, ' DEC strips'
!$omp parallel do reduction(+:cds) default(shared) &
!$omp& private(iq1,iq2,iq3,jq1,jq2,jq3,rcl,stm2,dltdec,jq1m) &
!$omp& private(dltra,jq2m,jq2min,jq2max,jq2t,j,i,xi,dci,wi,wpp,rv,rp2,r2,shth2,ii) &
!$omp& schedule(guided) firstprivate(fracp) if(nthr>1)

do i=1,npt   !---- Loop over ith particles ----
   xi  = x(i)
   dci = dc(i)
   wi  = wei(i)
   iq1 = int((dec(i)-decl)/hc1)+1
   iq2 = int((ra(i)-ral)/hc2)+1
   iq3 = int((dci-dcoml)/hc3)+1

   fracp = fracp + 1  ! accumulate particles and check when above the step size
   if(fracp>=dpart) then
       !$omp critical
       nadv = nadv + 1
       !omp flush (nadv)
       p1   = (nadv-1)*dpart + 1
       p2   = nadv*dpart
       if((nadv+1)*dpart>npt) p2=npt
       ! Note we are not really counting in mxh1 strips, just mymicking
       write(*,fmt="(i4)",advance='no') nadv                        ! for screen
       write(11,*) cntid//' counting in DEC strip > ',nadv,' (',p1,'-',p2,')' !for disk
       fracp = 0
       !$omp end critical
   endif
   
   lp_jq3: do jq3 = iq3-1,iq3+1
      if(jq3<1.or.jq3>nc3) cycle lp_jq3
      rcl    = (jq3-1)*hc3+dcoml
      stm2   = sthmax2(rsmax,dci,rcl)
      dltdec = 2.0d0*asin(stm2)*rad2deg
      jq1m   = int(dltdec/hc1)+1
      lp_jq1: do jq1=iq1-jq1m,iq1+jq1m
         if(jq1>nc1.or.jq1<1) cycle lp_jq1
         if(jq1==iq1) then
            dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,1)
         else
            dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,0)
         end if
         jq2m   = int(dltra/hc2)+1
         jq2max = iq2+jq2m
         jq2min = iq2-jq2m
         if(jq2max-jq2min+1>nc2) jq2max = jq2min-1+nc2
         lp_jq2: do jq2=jq2min,jq2max
            if(jq2>nc2) then
               jq2t = jq2-nc2
            else if(jq2<1) then
               jq2t = jq2+nc2
            else
               jq2t = jq2
            end if
            j = sk1(jq3,jq2t,jq1)
            do while(j/=0)
               !---- Loop over jth particles ----
               ! For each ith particle, index j runs through all particles in
               ! neighboring cells. Prune those outside rv and rp2 range.
               rv  = abs(dci-dc1(j))
               rp2 = 4.*dci*dc1(j)*((xi-x1(j))**2 + (y(i)-y1(j))**2 + (z(i)-z1(j))**2)
               r2  = rv*rv + rp2
               if(r2<=rsmax2) then
                  wpp = wi*wei1(j)    !weight by input
                  if(wfib==1) then    !weight by fiber if requested
                      shth2 = rp2/(4.*dci*dc1(j))
                      wpp   = wpp*wfiber(shth2)
                  endif
                  if(r2>seps2(nseps)) then
                     cds(nseps) = cds(nseps) + wpp
                     goto 70
                  endif
                  do ii=nseps-1,1,-1
                     if(r2>seps2(ii)) then
                        cds(ii) = cds(ii) + wpp
                        goto 70
                     endif
                  enddo
               else 
                  if(rv>rsmax.and.jq3>iq3) cycle lp_jq2
               endif
               70 j = ll1(j)
            end do
         end do lp_jq2
      end do lp_jq1
   end do lp_jq3
end do
!$omp end parallel do
close(11)  ! close log
write(*,*) ' '
end subroutine s_C_wg


subroutine s_Cb(nt,npt,ra,dec,dc,x,y,z,npt1,dc1,x1,y1,z1, &
                nseps,seps,sbound,mxh1,mxh2,mxh3,nbts,bseed,cntid,logf,sk1,ll1,cds,bcds)
!===============================================================================
! NAME
!  s_Cb()
!
! DESCRIPTION
!  Cross-count data pairs in redshift space for TWO samples of particles. Also
!  count boostraped data pairs.
!
! INPUTS
!  Variable---Type-----------Description----------------------------------------
!  nt         int            Number of thread to use (OpenMP)
!  npt        int            [sampleC]Number of particles
!  ra         r8(npt)        [sampleC] RA of particles [deg]
!  dec        r8(npt)        [sampleC]DEC of particles [deg]
!  dc         r8(npt)        [sampleC]DCOM of particles [Mpc/h]
!  x,y,z      r8(npt)        [sampleC]X,Y,Z coordinates of particles (see radec2xyz())
!  npt1       int            [sampleD] Number of particles
!  dc1        r8(npt1)       [sampleD] DCOM of particles [Mpc/h]
!  x1,y1,z1   r8(npt1)       [sampleD] X,Y,Z coordinates of particles (see radec2xyz())
!  nseps      int            Number of redshift-space separation bins
!  seps       r8(nseps+1)    Bins in redshift-space separation [Mpc/h]
!  sbound     r8(6)          Survey boundaries in RA,DEC,DCOM
!                            Form is (ramin,ramax,decmin,decmax,dcmin,dcmax)
!  mxh1       int            Nr of DEC cells of skip table
!  mxh2       int            Nr of RA cells of skip table
!  mxh3       int            Nr of DCOM cells of skip table
!  nbts       int            Number of boostrap samples
!  bseed      int            Seed for RNG during boostrap resampling
!  cntid      char*2         Two-character string to identify the samples,
!                            e.g. DD, RR, DR,... (used only for ouput messages)
!  logf       char*80        String for the file to log Fortran status messages
!  sk1        int(mxh3,mxh2,mxh1)  [sampleD] Skip table (SK) contructed with skll2d()
!  ll1        int(npt)       [sampleD] Linked list (LL) contructed with skll2d()
!
! OUTPUTS
!  Variable---Type-------------Description--------------------------------------
!  cds        r8(nseps)      Counts in redshift-space separation bins
!  bcds       r8(nbts,nseps) Boostrap counts in redshift-space separation bins
!
! NOTES  -----------------------------------------------------------------------
!  1. RA limits should be set to ramin=0 and ramax=360.
!  2. Remember to update the declarations in cflibfor.pyf if you add/remove 
!     in/out parameters to this exposed subroutine.

implicit none
real(kind=8)  :: ra(npt),dec(npt),sbound(6),ral,rau,decl,decu,dcoml,dcomu,hc1,hc2,hc3,rp2,r2
real(kind=8)  :: x(npt),y(npt),z(npt),x1(npt1),y1(npt1),z1(npt1),rcl,stm2,dltdec,dltra
real(kind=8)  :: dc(npt),dc1(npt1),xi,dci
real(kind=4)  :: wbts(nbts,npt),wbts1(nbts,npt1)
real(kind=8)  :: seps(nseps+1),seps2(nseps+1),rsmax,rsmax2,rv
real(kind=8)  :: cds(nseps),bcds(nbts,nseps)
integer       :: sk1(mxh3,mxh2,mxh1),ll1(npt1),mxh1,mxh2,mxh3,npt,npt1,nseps,nbts,bseed
integer       :: nt,nthr,fracp,dpart,nadv,nc1,nc2,nc3,jq1m,jq2m
integer       :: i,ii,j,jj,iq1,iq2,iq3,jq1,jq2,jq3,jq2min,jq2max,jq2t,p1,p2
character     :: cntid*2,logf*80

!------Open log file---------------------------------
open(11,file=logf,action='write',access='append')

fracp = 0 ; dpart = 0 ; nadv = 0  !Reset progress counters

!----------------------------------------------------
! Reset the counts, set (rpmax,rvmax) and square rp bins
cds    = 0.0d0
bcds   = 0.0d0
rsmax  = seps(nseps+1)
seps2  = seps*seps
rsmax2 = seps2(nseps+1)

!----------------------------------------------------
! Unpack boundaries and get cell dimensions in (dec,ra,dcom)
ral = sbound(1) ; rau = sbound(2) ; decl = sbound(3) ; decu = sbound(4)
dcoml = sbound(5) ; dcomu = sbound(6)
nc1 = mxh1
nc2 = mxh2
nc3 =int((dcomu-dcoml)/rsmax)
if(nc3>mxh3) nc3 = mxh3
hc1 = (decu-decl)/float(nc1)    !effective nr of DEC cells
hc2 = (rau-ral)/float(nc2)      !effective nr of RA cells
hc3 = (dcomu-dcoml)/float(nc3)  !effective nr of DCOM cells

!----------------------------------------------------
! Generate bootstrap samples
call bootstrap(npt,nbts,bseed,wbts)
call bootstrap(npt1,nbts,bseed,wbts1)

!----------------------------------------------------
! Set the number of threads
if(nt<=0) then
   nthr = omp_get_num_procs()
else
   nthr = nt
endif
call omp_set_num_threads(nthr)

!----------------------------------------------------
! Some inits
dpart = int(float(npt)/float(mxh1)) !choose dpart so we get mhx1 parts

!----------------------------------------------------
! Count pairs in SK grid
write(*,*) ' '
write(*,fmt='(a,i3,a)') '====  Counting '//cntid//' pairs in ', mxh1, ' DEC strips'
!$omp parallel do reduction(+:cds,bcds) default(shared) &
!$omp& private(iq1,iq2,iq3,jq1,jq2,jq3,rcl,stm2,dltdec,jq1m) &
!$omp& private(dltra,jq2m,jq2min,jq2max,jq2t,j,i,xi,dci,rv,rp2,r2,ii) &
!$omp& schedule(guided) firstprivate(fracp) if(nthr>1)

do i=1,npt   !---- Loop over ith particles ----
   xi  = x(i)
   dci = dc(i)
   iq1 = int((dec(i)-decl)/hc1)+1
   iq2 = int((ra(i)-ral)/hc2)+1
   iq3 = int((dci-dcoml)/hc3)+1

   fracp = fracp + 1  ! accumulate particles and check when above the step size
   if(fracp>=dpart) then
       !$omp critical
       nadv = nadv + 1
       !omp flush (nadv)
       p1   = (nadv-1)*dpart + 1
       p2   = nadv*dpart
       if((nadv+1)*dpart>npt) p2=npt
       ! Note we are not really counting in mxh1 strips, just mymicking
       write(*,fmt="(i4)",advance='no') nadv                        ! for screen
       write(11,*) cntid//' counting in DEC strip > ',nadv,' (',p1,'-',p2,')' !for disk
       fracp = 0
       !$omp end critical
   endif
   
   lp_jq3: do jq3 = iq3-1,iq3+1
      if(jq3<1.or.jq3>nc3) cycle lp_jq3
      rcl    = (jq3-1)*hc3+dcoml
      stm2   = sthmax2(rsmax,dci,rcl)
      dltdec = 2.0d0*asin(stm2)*rad2deg
      jq1m   = int(dltdec/hc1)+1
      lp_jq1: do jq1=iq1-jq1m,iq1+jq1m
         if(jq1>nc1.or.jq1<1) cycle lp_jq1
         if(jq1==iq1) then
            dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,1)
         else
            dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,0)
         end if
         jq2m   = int(dltra/hc2)+1
         jq2max = iq2+jq2m
         jq2min = iq2-jq2m
         if(jq2max-jq2min+1>nc2) jq2max = jq2min-1+nc2
         lp_jq2: do jq2=jq2min,jq2max
            if(jq2>nc2) then
               jq2t = jq2-nc2
            else if(jq2<1) then
               jq2t = jq2+nc2
            else
               jq2t = jq2
            end if
            j = sk1(jq3,jq2t,jq1)
            do while(j/=0)
               !---- Loop over jth particles ----
               ! For each ith particle, index j runs through all particles in
               ! neighboring cells. Prune those outside rv and rp2 range.
               rv  = abs(dci-dc1(j))
               rp2 = 4.*dci*dc1(j)*((xi-x1(j))**2 + (y(i)-y1(j))**2 + (z(i)-z1(j))**2)
               r2  = rv*rv + rp2
               if(r2<=rsmax2) then
                  if(r2>seps2(nseps)) then
                     cds(nseps) = cds(nseps) + 1.0d0
                     bcds(:,nseps) = bcds(:,nseps) + wbts(:,i)*wbts1(:,j)
                     goto 70
                  endif
                  do ii=nseps-1,1,-1
                     if(r2>seps2(ii)) then
                        cds(ii) = cds(ii) + 1.0d0
                        bcds(:,ii) = bcds(:,ii) + wbts(:,i)*wbts1(:,j)
                        goto 70
                     endif
                  enddo
               else 
                  if(rv>rsmax.and.jq3>iq3) cycle lp_jq2
               endif
               70 j = ll1(j)
            end do
         end do lp_jq2
      end do lp_jq1
   end do lp_jq3
end do
!$omp end parallel do
close(11)  ! close log
write(*,*) ' '
end subroutine s_Cb


subroutine s_Cb_wg(nt,npt,ra,dec,dc,wei,x,y,z,npt1,dc1,wei1,x1,y1,z1,nseps,seps,sbound, &
                   mxh1,mxh2,mxh3,nbts,bseed,wfib,cntid,logf,sk1,ll1,cds,bcds)
!===============================================================================
! NAME
!  s_Cb_wg()
!
! DESCRIPTION
!  Cross-count weighted data pairs in redshift space for TWO samples of particles.
!  Also count boostraped data pairs.
!
! INPUTS
!  Variable---Type-----------Description----------------------------------------
!  nt         int            Number of thread to use (OpenMP)
!  npt        int            [sampleC]Number of particles
!  ra         r8(npt)        [sampleC] RA of particles [deg]
!  dec        r8(npt)        [sampleC]DEC of particles [deg]
!  dc         r8(npt)        [sampleC]DCOM of particles [Mpc/h]
!  wei        r4(npt)        [sampleC] WEIGHT of particles
!  x,y,z      r8(npt)        [sampleC]X,Y,Z coordinates of particles (see radec2xyz())
!  npt1       int            [sampleD] Number of particles
!  dc1        r8(npt1)       [sampleD] DCOM of particles [Mpc/h]
!  wei1       r4(npt1)       [sampleD] WEIGHT of particles
!  x1,y1,z1   r8(npt1)       [sampleD] X,Y,Z coordinates of particles (see radec2xyz())
!  nseps      int            Number of redshift-space separation bins
!  seps       r8(nseps+1)    Bins in redshift-space separation [Mpc/h]
!  sbound     r8(6)          Survey boundaries in RA,DEC,DCOM
!                            Form is (ramin,ramax,decmin,decmax,dcmin,dcmax)
!  mxh1       int            Nr of DEC cells of skip table
!  mxh2       int            Nr of RA cells of skip table
!  mxh3       int            Nr of DCOM cells of skip table
!  nbts       int            Number of boostrap samples
!  bseed      int            Seed for RNG during boostrap resampling
!  wfib       int            If 1, also apply fiber-correction weights
!  cntid      char*2         Two-character string to identify the samples,
!                            e.g. DD, RR, DR,... (used only for ouput messages)
!  logf       char*80        String for the file to log Fortran status messages
!  sk1        int(mxh3,mxh2,mxh1)  [sampleD] Skip table (SK) contructed with skll2d()
!  ll1        int(npt)       [sampleD] Linked list (LL) contructed with skll2d()
!
! OUTPUTS
!  Variable---Type-------------Description--------------------------------------
!  cds        r8(nseps)      Counts in redshift-space separation bins
!  bcds       r8(nbts,nseps) Boostrap counts in redshift-space separation bins
!
! NOTES  -----------------------------------------------------------------------
!  1. RA limits should be set to ramin=0 and ramax=360.
!  2. Remember to update the declarations in cflibfor.pyf if you add/remove 
!     in/out parameters to this exposed subroutine.

implicit none
real(kind=8)  :: ra(npt),dec(npt),sbound(6),ral,rau,decl,decu,dcoml,dcomu,hc1,hc2,hc3,rp2,r2
real(kind=8)  :: x(npt),y(npt),z(npt),x1(npt1),y1(npt1),z1(npt1),rcl,stm2,dltdec,dltra
real(kind=8)  :: dc(npt),dc1(npt1),xi,dci
real(kind=4)  :: wei(npt),wei1(npt),wpp,wi
real(kind=4)  :: wbts(nbts,npt),wbts1(nbts,npt1)
real(kind=8)  :: seps(nseps+1),seps2(nseps+1),rsmax,rsmax2,rv,shth2
real(kind=8)  :: cds(nseps),bcds(nbts,nseps)
integer       :: sk1(mxh3,mxh2,mxh1),ll1(npt1),mxh1,mxh2,mxh3,npt,npt1,nseps,wfib,nbts,bseed
integer       :: nt,nthr,fracp,dpart,nadv,nc1,nc2,nc3,jq1m,jq2m
integer       :: i,ii,j,jj,iq1,iq2,iq3,jq1,jq2,jq3,jq2min,jq2max,jq2t,p1,p2
character     :: cntid*2,logf*80

!------Open log file---------------------------------
open(11,file=logf,action='write',access='append')

fracp = 0 ; dpart = 0 ; nadv = 0  !Reset progress counters

!----------------------------------------------------
! Reset the counts, set (rpmax,rvmax) and square rp bins
cds    = 0.0d0
bcds   = 0.0d0
rsmax  = seps(nseps+1)
seps2  = seps*seps
rsmax2 = seps2(nseps+1)

!----------------------------------------------------
! Unpack boundaries and get cell dimensions in (dec,ra,dcom)
ral = sbound(1) ; rau = sbound(2) ; decl = sbound(3) ; decu = sbound(4)
dcoml = sbound(5) ; dcomu = sbound(6)
nc1 = mxh1
nc2 = mxh2
nc3 =int((dcomu-dcoml)/rsmax)
if(nc3>mxh3) nc3 = mxh3
hc1 = (decu-decl)/float(nc1)    !effective nr of DEC cells
hc2 = (rau-ral)/float(nc2)      !effective nr of RA cells
hc3 = (dcomu-dcoml)/float(nc3)  !effective nr of DCOM cells

!----------------------------------------------------
! Generate bootstrap samples
call bootstrap(npt,nbts,bseed,wbts)
call bootstrap(npt1,nbts,bseed,wbts1)

!----------------------------------------------------
! Set the number of threads
if(nt<=0) then
   nthr = omp_get_num_procs()
else
   nthr = nt
endif
call omp_set_num_threads(nthr)

!----------------------------------------------------
! Some inits
dpart = int(float(npt)/float(mxh1)) !choose dpart so we get mhx1 parts

!----------------------------------------------------
! Count pairs in SK grid
write(*,*) ' '
write(*,fmt='(a,i3,a)') '====  Counting '//cntid//' pairs in ', mxh1, ' DEC strips'
!$omp parallel do reduction(+:cds,bcds) default(shared) &
!$omp& private(iq1,iq2,iq3,jq1,jq2,jq3,rcl,stm2,dltdec,jq1m) &
!$omp& private(dltra,jq2m,jq2min,jq2max,jq2t,j,i,xi,dci,wi,wpp,rv,rp2,r2,shth2,ii) &
!$omp& schedule(guided) firstprivate(fracp) if(nthr>1)

do i=1,npt   !---- Loop over ith particles ----
   xi  = x(i)
   dci = dc(i)
   wi  = wei(i)
   iq1 = int((dec(i)-decl)/hc1)+1
   iq2 = int((ra(i)-ral)/hc2)+1
   iq3 = int((dci-dcoml)/hc3)+1

   fracp = fracp + 1  ! accumulate particles and check when above the step size
   if(fracp>=dpart) then
       !$omp critical
       nadv = nadv + 1
       !omp flush (nadv)
       p1   = (nadv-1)*dpart + 1
       p2   = nadv*dpart
       if((nadv+1)*dpart>npt) p2=npt
       ! Note we are not really counting in mxh1 strips, just mymicking
       write(*,fmt="(i4)",advance='no') nadv                        ! for screen
       write(11,*) cntid//' counting in DEC strip > ',nadv,' (',p1,'-',p2,')' !for disk
       fracp = 0
       !$omp end critical
   endif
   
   lp_jq3: do jq3 = iq3-1,iq3+1
      if(jq3<1.or.jq3>nc3) cycle lp_jq3
      rcl    = (jq3-1)*hc3+dcoml
      stm2   = sthmax2(rsmax,dci,rcl)
      dltdec = 2.0d0*asin(stm2)*rad2deg
      jq1m   = int(dltdec/hc1)+1
      lp_jq1: do jq1=iq1-jq1m,iq1+jq1m
         if(jq1>nc1.or.jq1<1) cycle lp_jq1
         if(jq1==iq1) then
            dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,1)
         else
            dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,0)
         end if
         jq2m   = int(dltra/hc2)+1
         jq2max = iq2+jq2m
         jq2min = iq2-jq2m
         if(jq2max-jq2min+1>nc2) jq2max = jq2min-1+nc2
         lp_jq2: do jq2=jq2min,jq2max
            if(jq2>nc2) then
               jq2t = jq2-nc2
            else if(jq2<1) then
               jq2t = jq2+nc2
            else
               jq2t = jq2
            end if
            j = sk1(jq3,jq2t,jq1)
            do while(j/=0)
               !---- Loop over jth particles ----
               ! For each ith particle, index j runs through all particles in
               ! neighboring cells. Prune those outside rv and rp2 range.
               rv  = abs(dci-dc1(j))
               rp2 = 4.*dci*dc1(j)*((xi-x1(j))**2 + (y(i)-y1(j))**2 + (z(i)-z1(j))**2)
               r2  = rv*rv + rp2
               if(r2<=rsmax2) then
                  wpp = wi*wei1(j)    !weight by input
                  if(wfib==1) then    !weight by fiber if requested
                      shth2 = rp2/(4.*dci*dc1(j))
                      wpp   = wpp*wfiber(shth2)
                  endif
                  if(r2>seps2(nseps)) then
                     cds(nseps) = cds(nseps) + wpp
                     bcds(:,nseps) = bcds(:,nseps) + wpp*wbts(:,i)*wbts1(:,j)
                     goto 70
                  endif
                  do ii=nseps-1,1,-1
                     if(r2>seps2(ii)) then
                        cds(ii) = cds(ii) + wpp
                        bcds(:,ii) = bcds(:,ii) + wpp*wbts(:,i)*wbts1(:,j)
                        goto 70
                     endif
                  enddo
               else 
                  if(rv>rsmax.and.jq3>iq3) cycle lp_jq2
               endif
               70 j = ll1(j)
            end do
         end do lp_jq2
      end do lp_jq1
   end do lp_jq3
end do
!$omp end parallel do
close(11)  ! close log
write(*,*) ' '
end subroutine s_Cb_wg


!==========================================================================================
!   COUNTING ROUTINES - 2D SURFACE
!==========================================================================================

subroutine skll2d(mxh1,mxh2,npt,ra,dec,sbound,sk,ll)
!===============================================================================
! NAME
!  skll2d()
!
! PURPOSE
!  Construct a skip table (SK) and linked list (LL) for a set of particles with 
!  (ra,dec) coordinates. 
!  The SK table of size (mxh2,mxh1) is constructed within the area 
!  specified by sbound=(ramin,ramax,decmin,decmax)
!
! INPUTS
!  Variable----Type------------Description--------------------------------------
!  mxh1        int             Nr. of DEC cells in SK table
!  mxh2        int             Nr. of RA cells in SK table
!  npt          int             Number of particles
!  ra          real8(npt)       RA of particles
!  dec         real8(npt)       DEC of particles
!  sbound      real8(8)        Survey boundaries in RA,DEC
!
! OUTPUTS
!  Variable---Type---------------Description------------------------------------
!  sk          int(mxh1,mxh2)    SK table
!  ll          int(npt)           Linked list for particles
!
! NOTES  -----------------------------------------------------------------------
!  1. Angles are in degrees
!  2. The boundaries of RA should be 0 - 360

implicit none
integer      :: mxh1,mxh2,npt,q(2)
real(kind=8) :: ra(npt),dec(npt),sbound(4),hc1,hc2
integer      :: sk(mxh2,mxh1),ll(npt),i,i0

sk = 0
!-------------------------------------
! Get cell dimensions in (dec,ra)
hc1 = (sbound(4)-sbound(3))/float(mxh1)
hc2 = (sbound(2)-sbound(1))/float(mxh2)

!------------------------------
! Build SK and LL
do i=1, npt
   ! Find cell coordinates of the i-th particle
   q(1) = int((dec(i)-sbound(3))/hc1)+1
   q(2) = int((ra(i)-sbound(1))/hc2)+1
   
   ! Check coordinates are within range
   if(q(1)>mxh1.or.q(1)<1) cycle
   if(q(2)>mxh2) then
      q(2)=q(2)-mxh2
   else if(q(2)<1) then
      q(2)=q(2)+mxh2
   end if

   ! Assign LL and SK entry
   ll(i) = sk(q(2),q(1))
   sk(q(2),q(1)) = i

   ! Sort in increasing th(i) so that the sk element has the smallest th
   if(ll(i)==0) cycle
   if(dec(i)<=dec(ll(i))) cycle
   sk(q(2),q(1)) = ll(i)
   i0 = ll(i)
   11 if(ll(i0)/=0)then
      if(dec(i)<=dec(ll(i0)))then
         ll(i)  = ll(i0)
         ll(i0) = i
         cycle
      else
         i0 = ll(i0)
         goto 11
      endif
   else
      ll(i0) = i
      ll(i)  = 0
   endif
end do

end subroutine skll2d


subroutine th_A(nt,npt,dec,x,y,z,nsep,sep,sbound,mxh1,mxh2,cntid,logf,sk,ll,aa)
!===============================================================================
! NAME
!  th_A()
!
! DESCRIPTION
!  Count data pairs in angular space for ONE sample of particles
!
! INPUTS
!  Variable---Type-------------Description--------------------------------------
!  nt         int              Number of thread to use (OpenMP)
!  npt        int              Number of particles
!  dec        r8(npt)          DEC of particles [deg]
!  x,y,z      r8(npt)          X,Y,Z coordinates of particles (see radec2xyz())
!  nsep       int              Number of angular separation bins
!  sep        r8(nsep+1)       Bins in angular separation [deg]
!  sbound     r8(4)            Survey boundaries in RA,DEC. Form is (ramin,ramax,decmin,decmax)
!  mxh1       int              Nr of DEC cells of skip table
!  mxh2       int              Nr of RA cells of skip table
!  cntid      char*2           Two-character string to identify the samples,
!                              e.g. DD, RR, DR,... (used only for ouput messages)
!  logf       char*80          String for the file to log Fortran status messages
!  sk         int(mxh2,mxh1)   Skip table (SK) contructed with skll2d()
!  ll         int(npt)         Linked list (LL) contructed with skll2d()
!
! OUTPUTS
!  Variable---Type-------------Description--------------------------------------
!  aa        r8(nsep)          Counts in angular separation bins
!
! NOTES  -----------------------------------------------------------------------
!  1. RA limits should be set to ramin=0 and ramax=360.
!  2. Remember to update the declarations in cflibfor.pyf if you add/remove 
!     in/out parameters to this exposed subroutine.

implicit none
real(kind=8) :: dec(npt),ral,rau,decl,decu,hc1,hc2,shth2,stm2,dltra,dltdec
real(kind=8) :: x(npt),y(npt),z(npt),sep(nsep+1),sep2(nsep+1),sep2max,sbound(4),xi,yi,zi
real(kind=8) :: aa(nsep)
integer      :: nt,nthr,npt,mxh1,mxh2,nsep,sk(mxh2,mxh1),ll(npt),nc1,nc2
integer      :: jq1m,iq1,iq2,jq1,jq2,i,j,jq2max,jq2min,jq2m,jq2t,ii
character    :: cntid*2,logf*80

!------Open log file---------------------------------
open(11,file=logf,action='write',access='append')

! Some useful print statements. Place where desired
!print *, mxh1,mxh2
!print *, nc1,nc2
!print *, hc1, hc2
!print *, omp_get_num_threads(), omp_get_max_threads(), omp_get_num_procs()

!------------------
! reset the counts, set max. ang. distance and square bins
aa      = 0.d0
dltdec  = sep(nsep+1)
sep2    = (sin(0.5*sep*deg2rad))**2
sep2max = sep2(nsep+1)
 
!-----------------------------------
! Unpack boundaries and get cell dimensions in (dec,ra)
ral = sbound(1) ; rau = sbound(2) ; decl = sbound(3) ; decu = sbound(4)
nc1 = mxh1
nc2 = mxh2
hc1 = (decu-decl)/float(nc1)   !effective nr of DEC cells
hc2 = (rau-ral)/float(nc2)     !effective nr of RA cells

!----------------------------------------------------
! Set the number of threads
if(nt<=0) then
   nthr = omp_get_num_procs()
else
   nthr = nt
endif
call omp_set_num_threads(nthr)


stm2 = sin(dltdec*0.5*deg2rad)
jq1m = int(dltdec/hc1)+1

!----------------------------------------------------
! Count pairs in SK grid
write(*,*) ' '
write(*,fmt='(a,i3,a)') '====  Counting '//cntid//' pairs in ', mxh1, ' DEC strips  ===='
!$omp parallel do reduction(+:aa) default(shared) &
!$omp& private(iq1,iq2,jq1,jq2,dltra,jq2m,jq2min,jq2max,jq2t,j,i) &
!$omp& private(xi,yi,zi,shth2,ii) &
!$omp& schedule(guided) if(nthr>1)
do iq1=1,nc1
   !$omp critical
   write(*,fmt="(i4)",advance='no') iq1                             ! for screen
   write(11,*) cntid//' counting in DEC strip > ',iq1,'/',mxh1      ! for disk
   flush(11)
   !$omp end critical
   do iq2=1,nc2
      i = sk(iq2,iq1) ! index of ith particle 
      do while(i/=0)  !---- Loop over ith particles ----
         lp_jq1: do jq1=iq1,iq1+jq1m
            if(jq1>nc1.or.jq1<1) cycle lp_jq1
            if(jq1==iq1) then
               dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,1)
            else
               dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,0)
            end if
            jq2m   = int(dltra/hc2)+1
            jq2max = iq2+jq2m
            jq2min = iq2-jq2m
            if(jq2max-jq2min+1>nc2) jq2max=jq2min-1+nc2
            lp_jq2: do jq2=jq2min,jq2max
               if(jq2>nc2) then
                  jq2t = jq2-nc2
               else if(jq2<1) then
                  jq2t = jq2+nc2
               else
                  jq2t = jq2
               end if
               if(jq2t<iq2.and.jq1==iq1) cycle lp_jq2
               if(jq2t==iq2.and.jq1==iq1) then
                  j = ll(i)
               else
                  j = sk(jq2t,jq1)
               endif
               xi = x(i)
               yi = y(i)
               zi = z(i)
               
               do while(j/=0)
                  shth2 = (xi-x(j))**2 + (yi-y(j))**2 + (zi-z(j))**2
                  if(shth2<=sep2max) then
                      ! Now count the pair by finding its (ii) bin in vector of ang-space bins
                      if(shth2>sep2(nsep)) then
                          aa(nsep) = aa(nsep) + 1.0d0
                          goto 78
                      endif
                      if(shth2>sep2(nsep-1)) then
                          aa(nsep-1) = aa(nsep-1) + 1.0d0
                          goto 78
                      endif
                      if(shth2>sep2(nsep-2)) then
                          aa(nsep-2) = aa(nsep-2) + 1.0d0
                          goto 78
                      endif
                      if(shth2>sep2(nsep-3)) then
                          aa(nsep-3) = aa(nsep-3) + 1.0d0
                          goto 78
                      endif
                      do ii=nsep-4,1,-1
                         if(shth2>sep2(ii)) then
                            aa(ii) = aa(ii) + 1.0d0
                            goto 78
                         endif
                      enddo
                  endif
                  78 j = ll(j)
               end do
               
            end do lp_jq2
         end do lp_jq1
         i = ll(i)
      end do
   end do !iq2 loop
end do !iq1 loop
!$omp end parallel do
close(11)  ! close log
write(*,*) ' '
end subroutine th_A


subroutine th_A_wg(nt,npt,dec,wei,x,y,z,nsep,sep,sbound,mxh1,mxh2,wfib,cntid,logf,sk,ll,aa)
!===============================================================================
! NAME
!  th_A_wg()
!
! DESCRIPTION
!  Count weighted data pairs in angular space for ONE sample of particles
!
! INPUTS
!  Variable---Type-------------Description--------------------------------------
!  nt         int              Number of thread to use (OpenMP)
!  npt        int              Number of particles
!  dec        r8(npt)          DEC of particles [deg]
!  wei        r4(npt)          WEIGHT of particles
!  x,y,z      r8(npt)          X,Y,Z coordinates of particles (see radec2xyz())
!  nsep       int              Number of angular separation bins
!  sep        r8(nsep+1)       Bins in angular separation [deg]
!  sbound     r8(4)            Survey boundaries in RA,DEC. Form is (ramin,ramax,decmin,decmax)
!  mxh1       int              Nr of DEC cells of skip table
!  mxh2       int              Nr of RA cells of skip table
!  wfib       int              If 1, also apply fiber-correction weights
!  cntid      char*2           Two-character string to identify the samples,
!                              e.g. DD, RR, DR,... (used only for ouput messages)
!  logf       char*80          String for the file to log Fortran status messages
!  sk         int(mxh2,mxh1)   Skip table (SK) contructed with skll2d()
!  ll         int(npt)         Linked list (LL) contructed with skll2d()
!
! OUTPUTS
!  Variable---Type-------------Description--------------------------------------
!  aa        r8(nsep)          Counts in angular separation bins
!
! NOTES  -----------------------------------------------------------------------
!  1. RA limits should be set to ramin=0 and ramax=360.
!  2. Remember to update the declarations in cflibfor.pyf if you add/remove 
!     in/out parameters to this exposed subroutine.

implicit none
real(kind=8) :: dec(npt),ral,rau,decl,decu,hc1,hc2,shth2,stm2,dltra,dltdec
real(kind=8) :: x(npt),y(npt),z(npt),sep(nsep+1),sep2(nsep+1),sep2max,sbound(4),xi,yi,zi
real(kind=4) :: wei(npt),wi,wpp
real(kind=8) :: aa(nsep)
integer      :: nt,nthr,npt,mxh1,mxh2,nsep,sk(mxh2,mxh1),ll(npt),nc1,nc2
integer      :: jq1m,iq1,iq2,jq1,jq2,i,j,jq2max,jq2min,jq2m,jq2t,ii,wfib
character    :: cntid*2,logf*80

!------Open log file---------------------------------
open(11,file=logf,action='write',access='append')

! Some useful print statements. Place where desired
!print *, mxh1,mxh2
!print *, nc1,nc2
!print *, hc1, hc2
!print *, omp_get_num_threads(), omp_get_max_threads(), omp_get_num_procs()

!------------------
! reset the counts, set max. ang. distance and square bins
aa      = 0.d0
dltdec  = sep(nsep+1)
sep2    = (sin(0.5*sep*deg2rad))**2
sep2max = sep2(nsep+1)
 
!-----------------------------------
! Unpack boundaries and get cell dimensions in (dec,ra)
ral = sbound(1) ; rau = sbound(2) ; decl = sbound(3) ; decu = sbound(4)
nc1 = mxh1
nc2 = mxh2
hc1 = (decu-decl)/float(nc1)   !effective nr of DEC cells
hc2 = (rau-ral)/float(nc2)     !effective nr of RA cells

!----------------------------------------------------
! Set the number of threads
if(nt<=0) then
   nthr = omp_get_num_procs()
else
   nthr = nt
endif
call omp_set_num_threads(nthr)


stm2 = sin(dltdec*0.5*deg2rad)
jq1m = int(dltdec/hc1)+1

!----------------------------------------------------
! Count pairs in SK grid
write(*,*) ' '
write(*,fmt='(a,i3,a)') '====  Counting '//cntid//' pairs in ', mxh1, ' DEC strips  ===='
!$omp parallel do reduction(+:aa) default(shared) &
!$omp& private(iq1,iq2,jq1,jq2,dltra,jq2m,jq2min,jq2max,jq2t,j,i) &
!$omp& private(xi,yi,zi,wi,wpp,shth2,ii) &
!$omp& schedule(guided) if(nthr>1)
do iq1=1,nc1
   !$omp critical
   write(*,fmt="(i4)",advance='no') iq1                             ! for screen
   write(11,*) cntid//' counting in DEC strip > ',iq1,'/',mxh1      ! for disk
   flush(11)
   !$omp end critical
   do iq2=1,nc2
      i = sk(iq2,iq1) ! index of ith particle 
      do while(i/=0)  !---- Loop over ith particles ----
         lp_jq1: do jq1=iq1,iq1+jq1m
            if(jq1>nc1.or.jq1<1) cycle lp_jq1
            if(jq1==iq1) then
               dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,1)
            else
               dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,0)
            end if
            jq2m   = int(dltra/hc2)+1
            jq2max = iq2+jq2m
            jq2min = iq2-jq2m
            if(jq2max-jq2min+1>nc2) jq2max=jq2min-1+nc2
            lp_jq2: do jq2=jq2min,jq2max
               if(jq2>nc2) then
                  jq2t = jq2-nc2
               else if(jq2<1) then
                  jq2t = jq2+nc2
               else
                  jq2t = jq2
               end if
               if(jq2t<iq2.and.jq1==iq1) cycle lp_jq2
               if(jq2t==iq2.and.jq1==iq1) then
                  j = ll(i)
               else
                  j = sk(jq2t,jq1)
               endif
               xi = x(i)
               yi = y(i)
               zi = z(i)
               wi = wei(i)
               
               do while(j/=0)
                  shth2 = (xi-x(j))**2 + (yi-y(j))**2 + (zi-z(j))**2
                  if(shth2<=sep2max) then
                      wpp = wi*wei(j)          !weight by input
                      if(wfib==1) then         !weight by fiber
                          wpp = wpp*wfiber(shth2)  !TODO ???
                      endif
                      ! Now count the pair by finding its (ii) bin in vector of ang-space bins
                      if(shth2>sep2(nsep)) then
                          aa(nsep) = aa(nsep) + wpp
                          goto 78
                      endif
                      if(shth2>sep2(nsep-1)) then
                          aa(nsep-1) = aa(nsep-1) + wpp
                          goto 78
                      endif
                      if(shth2>sep2(nsep-2)) then
                          aa(nsep-2) = aa(nsep-2) + wpp
                          goto 78
                      endif
                      if(shth2>sep2(nsep-3)) then
                          aa(nsep-3) = aa(nsep-3) + wpp
                          goto 78
                      endif
                      do ii=nsep-4,1,-1
                         if(shth2>sep2(ii)) then
                            aa(ii) = aa(ii) + wpp
                            goto 78
                         endif
                      enddo
                  endif
                  78 j = ll(j)
               end do
               
            end do lp_jq2
         end do lp_jq1
         i = ll(i)
      end do
   end do !iq2 loop
end do !iq1 loop
!$omp end parallel do
close(11)  ! close log
write(*,*) ' '
end subroutine th_A_wg


subroutine th_Ab(nt,npt,dec,x,y,z,nsep,sep,sbound,mxh1,mxh2,nbts,bseed,cntid,logf, &
           sk,ll,aa,baa)
!===============================================================================
! NAME
!  th_Ab()
!
! DESCRIPTION
!  Count data pairs in angular space for ONE sample of particles. Also count
!  booostraped data pairs.
!
! INPUTS
!  Variable---Type-------------Description--------------------------------------
!  nt         int              Number of thread to use (OpenMP)
!  npt        int              Number of particles
!  dec        r8(npt)          DEC of particles [deg]
!  x,y,z      r8(npt)          X,Y,Z coordinates of particles (see radec2xyz())
!  nsep       int              Number of angular separation bins
!  sep        r8(nsep+1)       Bins in angular separation [deg]
!  sbound     r8(4)            Survey boundaries in RA,DEC. Form is (ramin,ramax,decmin,decmax)
!  mxh1       int              Nr of DEC cells of skip table
!  mxh2       int              Nr of RA cells of skip table
!  nbts       int              Number of boostrap samples
!  bseed      int              Seed for RNG during boostrap resampling
!  cntid      char*2           Two-character string to identify the samples,
!                              e.g. DD, RR, DR,... (used only for ouput messages)
!  logf       char*80          String for the file to log Fortran status messages
!  sk         int(mxh2,mxh1)   Skip table (SK) contructed with skll2d()
!  ll         int(npt)         Linked list (LL) contructed with skll2d()
!
! OUTPUTS
!  Variable---Type-------------Description--------------------------------------
!  aa        r8(nsep)          Counts in angular separation bins
!  baa       r8(nbts,nsep)     Boostrap counts in angular separation bins
!
! NOTES  -----------------------------------------------------------------------
!  1. RA limits should be set to ramin=0 and ramax=360.
!  2. Remember to update the declarations in cflibfor.pyf if you add/remove 
!     in/out parameters to this exposed subroutine.

implicit none
real(kind=8) :: dec(npt),ral,rau,decl,decu,hc1,hc2,shth2,stm2,dltra,dltdec
real(kind=8) :: x(npt),y(npt),z(npt),sep(nsep+1),sep2(nsep+1),sep2max,sbound(4),xi,yi,zi
real(kind=4) :: wbts(nbts,npt)
real(kind=8) :: aa(nsep),baa(nbts,nsep)
integer      :: nt,nthr,npt,mxh1,mxh2,nsep,sk(mxh2,mxh1),ll(npt),nc1,nc2,nbts,bseed
integer      :: jq1m,iq1,iq2,jq1,jq2,i,j,jq2max,jq2min,jq2m,jq2t,ii
character    :: cntid*2,logf*80

!------Open log file---------------------------------
open(11,file=logf,action='write',access='append')

! Some useful print statements. Place where desired
!print *, mxh1,mxh2
!print *, nc1,nc2
!print *, hc1, hc2
!print *, omp_get_num_threads(), omp_get_max_threads(), omp_get_num_procs()

!------------------
! reset the counts, set max. ang. distance and square bins
aa      = 0.d0
baa     = 0.d0
dltdec  = sep(nsep+1)
sep2    = (sin(0.5*sep*deg2rad))**2
sep2max = sep2(nsep+1)
 
!-----------------------------------
! Unpack boundaries and get cell dimensions in (dec,ra)
ral = sbound(1) ; rau = sbound(2) ; decl = sbound(3) ; decu = sbound(4)
nc1 = mxh1
nc2 = mxh2
hc1 = (decu-decl)/float(nc1)   !effective nr of DEC cells
hc2 = (rau-ral)/float(nc2)     !effective nr of RA cells

!----------------------------------------------------
! Generate bootstrap samples
call bootstrap(npt,nbts,bseed,wbts)

!----------------------------------------------------
! Set the number of threads
if(nt<=0) then
   nthr = omp_get_num_procs()
else
   nthr = nt
endif
call omp_set_num_threads(nthr)


stm2 = sin(dltdec*0.5*deg2rad)
jq1m = int(dltdec/hc1)+1

!----------------------------------------------------
! Count pairs in SK grid
write(*,*) ' '
write(*,fmt='(a,i3,a)') '====  Counting '//cntid//' pairs in ', mxh1, ' DEC strips  ===='
!$omp parallel do reduction(+:aa,baa) default(shared) &
!$omp& private(iq1,iq2,jq1,jq2,dltra,jq2m,jq2min,jq2max,jq2t,j,i) &
!$omp& private(xi,yi,zi,shth2,ii) &
!$omp& schedule(guided) if(nthr>1)
do iq1=1,nc1
   !$omp critical
   write(*,fmt="(i4)",advance='no') iq1                             ! for screen
   write(11,*) cntid//' counting in DEC strip > ',iq1,'/',mxh1      ! for disk
   flush(11)
   !$omp end critical
   do iq2=1,nc2
      i = sk(iq2,iq1) ! index of ith particle 
      do while(i/=0)  !---- Loop over ith particles ----
         lp_jq1: do jq1=iq1,iq1+jq1m
            if(jq1>nc1.or.jq1<1) cycle lp_jq1
            if(jq1==iq1) then
               dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,1)
            else
               dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,0)
            end if
            jq2m   = int(dltra/hc2)+1
            jq2max = iq2+jq2m
            jq2min = iq2-jq2m
            if(jq2max-jq2min+1>nc2) jq2max=jq2min-1+nc2
            lp_jq2: do jq2=jq2min,jq2max
               if(jq2>nc2) then
                  jq2t = jq2-nc2
               else if(jq2<1) then
                  jq2t = jq2+nc2
               else
                  jq2t = jq2
               end if
               if(jq2t<iq2.and.jq1==iq1) cycle lp_jq2
               if(jq2t==iq2.and.jq1==iq1) then
                  j = ll(i)
               else
                  j = sk(jq2t,jq1)
               endif
               xi = x(i)
               yi = y(i)
               zi = z(i)
               do while(j/=0)
                  shth2 = (xi-x(j))**2 + (yi-y(j))**2 + (zi-z(j))**2
                  if(shth2<=sep2max) then
                      ! Now count the pair by finding its (ii) bin in vector of ang-space bins
                      if(shth2>sep2(nsep)) then
                          aa(nsep) = aa(nsep) + 1.0d0
                          baa(:,nsep) = baa(:,nsep) + wbts(:,i)*wbts(:,j)
                          goto 78
                      endif
                      if(shth2>sep2(nsep-1)) then
                          aa(nsep-1) = aa(nsep-1) + 1.0d0
                          baa(:,nsep-1) = baa(:,nsep-1) + wbts(:,i)*wbts(:,j)
                          goto 78
                      endif
                      if(shth2>sep2(nsep-2)) then
                          aa(nsep-2) = aa(nsep-2) + 1.0d0
                          baa(:,nsep-2) = baa(:,nsep-2) + wbts(:,i)*wbts(:,j)
                          goto 78
                      endif
                      if(shth2>sep2(nsep-3)) then
                          aa(nsep-3) = aa(nsep-3) + 1.0d0
                          baa(:,nsep-3) = baa(:,nsep-3) + wbts(:,i)*wbts(:,j)
                          goto 78
                      endif
                      do ii=nsep-4,1,-1
                         if(shth2>sep2(ii)) then
                            aa(ii) = aa(ii) + 1.0d0
                            baa(:,ii) = baa(:,ii) + wbts(:,i)*wbts(:,j)
                            goto 78
                         endif
                      enddo
                  endif
                  78 j = ll(j)
               end do
            end do lp_jq2
         end do lp_jq1
         i = ll(i)
      end do
   end do !iq2 loop
end do !iq1 loop
!$omp end parallel do
close(11)  ! close log
write(*,*) ' '
end subroutine th_Ab


subroutine th_Ab_wg(nt,npt,dec,wei,x,y,z,nsep,sep,sbound,mxh1,mxh2,nbts,bseed,wfib, &
                    cntid,logf,sk,ll,aa,baa)
!===============================================================================
! NAME
!  th_Ab_wg()
!
! DESCRIPTION
!  Count weighted data pairs in angular space for ONE sample of particles. Also
!  count booostraped data pairs.
!
! INPUTS
!  Variable---Type-------------Description--------------------------------------
!  nt         int              Number of thread to use (OpenMP)
!  npt        int              Number of particles
!  dec        r8(npt)          DEC of particles [deg]
!  wei        r4(npt)          WEIGHT of particles
!  x,y,z      r8(npt)          X,Y,Z coordinates of particles (see radec2xyz())
!  nsep       int              Number of angular separation bins
!  sep        r8(nsep+1)       Bins in angular separation [deg]
!  sbound     r8(4)            Survey boundaries in RA,DEC. Form is (ramin,ramax,decmin,decmax)
!  mxh1       int              Nr of DEC cells of skip table
!  mxh2       int              Nr of RA cells of skip table
!  nbts       int              Number of boostrap samples
!  bseed      int              Seed for RNG during boostrap resampling
!  wfib       int              If 1, also apply fiber-correction weights
!  cntid      char*2           Two-character string to identify the samples,
!                              e.g. DD, RR, DR,... (used only for ouput messages)
!  logf       char*80          String for the file to log Fortran status messages
!  sk         int(mxh2,mxh1)   Skip table (SK) contructed with skll2d()
!  ll         int(npt)         Linked list (LL) contructed with skll2d()
!
! OUTPUTS
!  Variable---Type-------------Description--------------------------------------
!  aa        r8(nsep)          Counts in angular separation bins
!  baa       r8(nbts,nsep)     Boostrap counts in angular separation bins
!
! NOTES  -----------------------------------------------------------------------
!  1. RA limits should be set to ramin=0 and ramax=360.
!  2. Remember to update the declarations in cflibfor.pyf if you add/remove 
!     in/out parameters to this exposed subroutine.

implicit none
real(kind=8) :: dec(npt),ral,rau,decl,decu,hc1,hc2,shth2,stm2,dltra,dltdec
real(kind=8) :: x(npt),y(npt),z(npt),sep(nsep+1),sep2(nsep+1),sep2max,sbound(4),xi,yi,zi
real(kind=4) :: wei(npt),wi,wpp
real(kind=8) :: aa(nsep),baa(nbts,nsep)
real(kind=4) :: wbts(nbts,npt)
integer      :: nt,nthr,npt,mxh1,mxh2,nsep,sk(mxh2,mxh1),ll(npt),nc1,nc2,nbts,bseed,wfib
integer      :: jq1m,iq1,iq2,jq1,jq2,i,j,jq2max,jq2min,jq2m,jq2t,ii
character    :: cntid*2,logf*80

!------Open log file---------------------------------
open(11,file=logf,action='write',access='append')

! Some useful print statements. Place where desired
!print *, mxh1,mxh2
!print *, nc1,nc2
!print *, hc1, hc2
!print *, omp_get_num_threads(), omp_get_max_threads(), omp_get_num_procs()

!------------------
! reset the counts, set max. ang. distance and square bins
aa      = 0.d0
baa     = 0.d0
dltdec  = sep(nsep+1)
sep2    = (sin(0.5*sep*deg2rad))**2
sep2max = sep2(nsep+1)
 
!-----------------------------------
! Unpack boundaries and get cell dimensions in (dec,ra)
ral = sbound(1) ; rau = sbound(2) ; decl = sbound(3) ; decu = sbound(4)
nc1 = mxh1
nc2 = mxh2
hc1 = (decu-decl)/float(nc1)   !effective nr of DEC cells
hc2 = (rau-ral)/float(nc2)     !effective nr of RA cells

!----------------------------------------------------
! Generate bootstrap samples
call bootstrap(npt,nbts,bseed,wbts)

!----------------------------------------------------
! Set the number of threads
if(nt<=0) then
   nthr = omp_get_num_procs()
else
   nthr = nt
endif
call omp_set_num_threads(nthr)


stm2 = sin(dltdec*0.5*deg2rad)
jq1m = int(dltdec/hc1)+1

!----------------------------------------------------
! Count pairs in SK grid
write(*,*) ' '
write(*,fmt='(a,i3,a)') '====  Counting '//cntid//' pairs in ', mxh1, ' DEC strips  ===='
!$omp parallel do reduction(+:aa,baa) default(shared) &
!$omp& private(iq1,iq2,jq1,jq2,dltra,jq2m,jq2min,jq2max,jq2t,j,i) &
!$omp& private(xi,yi,zi,wi,wpp,shth2,ii) &
!$omp& schedule(guided) if(nthr>1)
do iq1=1,nc1
   !$omp critical
   write(*,fmt="(i4)",advance='no') iq1                             ! for screen
   write(11,*) cntid//' counting in DEC strip > ',iq1,'/',mxh1      ! for disk
   flush(11)
   !$omp end critical
   do iq2=1,nc2
      i = sk(iq2,iq1) ! index of ith particle 
      do while(i/=0)  !---- Loop over ith particles ----
         lp_jq1: do jq1=iq1,iq1+jq1m
            if(jq1>nc1.or.jq1<1) cycle lp_jq1
            if(jq1==iq1) then
               dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,1)
            else
               dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,0)
            end if
            jq2m   = int(dltra/hc2)+1
            jq2max = iq2+jq2m
            jq2min = iq2-jq2m
            if(jq2max-jq2min+1>nc2) jq2max=jq2min-1+nc2
            lp_jq2: do jq2=jq2min,jq2max
               if(jq2>nc2) then
                  jq2t = jq2-nc2
               else if(jq2<1) then
                  jq2t = jq2+nc2
               else
                  jq2t = jq2
               end if
               if(jq2t<iq2.and.jq1==iq1) cycle lp_jq2
               if(jq2t==iq2.and.jq1==iq1) then
                  j = ll(i)
               else
                  j = sk(jq2t,jq1)
               endif
               xi = x(i)
               yi = y(i)
               zi = z(i)
               wi = wei(i)
               
               do while(j/=0)
                  shth2 = (xi-x(j))**2 + (yi-y(j))**2 + (zi-z(j))**2
                  if(shth2<=sep2max) then
                      wpp = wi*wei(j)          !weight by input
                      if(wfib==1) then         !weight by fiber
                          wpp = wpp*wfiber(shth2)  !TODO ???
                      endif
                      ! Now count the pair by finding its (ii) bin in vector of ang-space bins
                      if(shth2>sep2(nsep)) then
                          aa(nsep) = aa(nsep) + wpp
                          baa(:,nsep) = baa(:,nsep) + wpp*wbts(:,i)*wbts(:,j)
                          goto 78
                      endif
                      if(shth2>sep2(nsep-1)) then
                          aa(nsep-1) = aa(nsep-1) + wpp
                          baa(:,nsep-1) = baa(:,nsep-1) + wpp*wbts(:,i)*wbts(:,j)
                          goto 78
                      endif
                      if(shth2>sep2(nsep-2)) then
                          aa(nsep-2) = aa(nsep-2) + wpp
                          baa(:,nsep-2) = baa(:,nsep-2) + wpp*wbts(:,i)*wbts(:,j)
                          goto 78
                      endif
                      if(shth2>sep2(nsep-3)) then
                          aa(nsep-3) = aa(nsep-3) + wpp
                          baa(:,nsep-3) = baa(:,nsep-3) + wpp*wbts(:,i)*wbts(:,j)
                          goto 78
                      endif
                      do ii=nsep-4,1,-1
                         if(shth2>sep2(ii)) then
                            aa(ii) = aa(ii) + wpp
                            baa(:,ii) = baa(:,ii) + wpp*wbts(:,i)*wbts(:,j)
                            goto 78
                         endif
                      enddo
                  endif
                  78 j = ll(j)
               end do
               
            end do lp_jq2
         end do lp_jq1
         i = ll(i)
      end do
   end do !iq2 loop
end do !iq1 loop
!$omp end parallel do
close(11)  ! close log
write(*,*) ' '
end subroutine th_Ab_wg


subroutine th_C(nt,npt,ra,dec,x,y,z,npt1,x1,y1,z1, &
                nsep,sep,sbound,mxh1,mxh2,cntid,logf,sk1,ll1,cdth)
!===============================================================================
! NAME
!  th_C()
!
! DESCRIPTION
!  Cross-count data pairs in angular space for TWO samples of particles
!
! INPUTS
!  Variable---Type-------------Description--------------------------------------
!  nt         int              Number of thread to use (OpenMP)
!  npt        int              Number of particles
!  ra         r8(npt)          [sampleC] RA of particles [deg]
!  dec        r8(npt)          [sampleC] DEC of particles [deg]
!  x,y,z      r8(npt)          [sampleC] X,Y,Z coordinates of particles (see radec2xyz())
!  npt1       int              [sampleD] Number of particles
!  x1,y1,z1   r8(npt1)         [sampleD] X,Y,Z coordinates of particles (see radec2xyz())
!  nsep       int              Number of angular separation bins
!  sep        r8(nsep+1)       Bins in angular separation [deg]
!  sbound     r8(4)            Survey boundaries in RA,DEC. Form is (ramin,ramax,decmin,decmax)
!  mxh1       int              Nr of DEC cells of skip table
!  mxh2       int              Nr of RA cells of skip table
!  cntid      char*2           Two-character string to identify the samples,
!                              e.g. DD, RR, DR,... (used only for ouput messages)
!  logf       char*80          String for the file to log Fortran status messages
!  sk1        int(mxh2,mxh1)   [sampleD] Skip table (SK) contructed with skll2d()
!  ll1        int(npt1)        [sampleD] Linked list (LL) contructed with skll2d()
!
! OUTPUTS
!  Variable---Type-------------Description--------------------------------------
!  cdth       r8(nsep)         Counts in angular separation bins
!
! NOTES  -----------------------------------------------------------------------
!  1. RA limits should be set to ramin=0 and ramax=360.
!  2. Remember to update the declarations in cflibfor.pyf if you add/remove 
!     in/out parameters to this exposed subroutine.

implicit none
real(kind=8) :: ra(npt),dec(npt),ral,rau,decl,decu,hc1,hc2,shth2,stm2,dltra,dltdec
real(kind=8) :: x(npt),y(npt),z(npt),x1(npt1),y1(npt1),z1(npt1) 
real(kind=8) :: sep(nsep+1),sep2(nsep+1),sbound(4),sep2max,xi,yi,zi
real(kind=8) :: cdth(nsep)
integer      :: nt,nthr,sk1(mxh2,mxh1),ll1(npt1),mxh1,mxh2,npt,npt1,nsep,ndp
integer      :: fracp,dpart,nadv,nc1,nc2,i,ii,j,iq1,iq2,p1,p2,jq1,jq2,jq1m
integer      :: jq2m,jq2max,jq2min,jq2t
character    :: cntid*2,logf*80

!------Open log file---------------------------------
open(11,file=logf,action='write',access='append')

fracp = 0 ; dpart = 0 ; nadv = 0  !Reset progress counters

!----------------------------------------------------
! reset the counts, set max. ang. distance and square bins
cdth    = 0.0d0
dltdec  = sep(nsep+1)
sep2    = (sin(0.5*sep*deg2rad))**2
sep2max = sep2(nsep+1)

!-----------------------------------
! Unpack boundaries and get cell dimensions in (dec,ra)
ral = sbound(1) ; rau = sbound(2) ; decl = sbound(3) ; decu = sbound(4)
nc1 = mxh1
nc2 = mxh2
hc1 = (decu-decl)/float(nc1)   !effective nr of DEC cells
hc2 = (rau-ral)/float(nc2)     !effective nr of RA cells

!----------------------------------------------------
! Set the number of threads
if(nt<=0) then
   nthr = omp_get_num_procs()
else
   nthr = nt
endif
call omp_set_num_threads(nthr)

!----------------------------------------------------
! Some inits
dpart = int(float(npt)/float(mxh1)) !choose dpart so we get mhx1 parts
stm2  = sin(dltdec*0.5*deg2rad)
jq1m  = int(dltdec/hc1)+1

!----------------------------------------------------
! Count pairs in SK grid
write(*,*) ' '
write(*,fmt='(a,i3,a)') '====  Counting '//cntid//' pairs in ', mxh1, ' DEC strips'
!$omp parallel do reduction(+:cdth) default(shared) &
!$omp& private(iq1,iq2,jq1,jq2) &
!$omp& private(dltra,jq2m,jq2min,jq2max,jq2t,j,i,xi,yi,zi,shth2,ii,fracp,p1,p2) &
!$omp& schedule(guided) if(nthr>1)

do i=1,npt
   iq1 = int((dec(i)-decl)/hc1)+1
   iq2 = int((ra(i)-ral)/hc2)+1
   fracp = fracp + 1  ! accumulate particles and check when above the step size
   if(fracp>=dpart) then
       !$omp critical
       nadv = nadv + 1
       !omp flush (nadv)
       p1   = (nadv-1)*dpart + 1
       p2   = nadv*dpart
       if((nadv+1)*dpart>npt) p2=npt
       !note we are not really counting in mxh1 strips, just mymicking
       write(*,fmt="(i4)",advance='no') nadv                        ! for screen
       write(11,*) cntid//' counting in DEC strip > ',nadv,' (',p1,'-',p2,')' !for disk
       flush(11)
       fracp = 0
       !$omp end critical
   endif
   lp_jq1: do jq1=iq1-jq1m,iq1+jq1m
      if(jq1>nc1.or.jq1<1) cycle lp_jq1
      if(jq1==iq1) then
         dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,1)
      else
         dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,0)
      end if
      jq2m   = int(dltra/hc2)+1
      jq2max = iq2+jq2m
      jq2min = iq2-jq2m
      if(jq2max-jq2min+1>nc2) jq2max=jq2min-1+nc2
      lp_jq2: do jq2=jq2min,jq2max
         if(jq2>nc2) then
            jq2t = jq2-nc2
         else if(jq2<1) then
            jq2t = jq2+nc2
         else
            jq2t = jq2
         end if
         j = sk1(jq2t,jq1)
         xi = x(i)
         yi = y(i)
         zi = z(i)

         do while(j/=0)
            shth2 = (xi-x1(j))**2 + (yi-y1(j))**2 + (zi-z1(j))**2
            !th = 2.*sqrt(shth2)*180./3.1415
            if(shth2<=sep2max) then
                ! Now count the pair by finding its (ii) bin in vector of ang-space bins
                if(shth2>sep2(nsep)) then
                    cdth(nsep) = cdth(nsep) + 1.0d0
                    goto 78
                endif
                if(shth2>sep2(nsep-1)) then
                    cdth(nsep-1) = cdth(nsep-1) + 1.0d0
                    goto 78
                endif
                if(shth2>sep2(nsep-2)) then
                    cdth(nsep-2) = cdth(nsep-2) + 1.0d0
                    goto 78
                endif
                if(shth2>sep2(nsep-3)) then
                    cdth(nsep-3) = cdth(nsep-3) + 1.0d0
                    goto 78
                endif
                do ii=nsep-4,1,-1
                   if(shth2>sep2(ii)) then
                      cdth(ii) = cdth(ii) + 1.0d0
                      goto 78
                   endif
                enddo
            endif
            78 j = ll1(j)
         end do
         
      end do lp_jq2
   end do lp_jq1
end do
!$omp end parallel do
close(11)  ! close log
write(*,*) ' '
end subroutine th_C


subroutine th_C_wg(nt,npt,ra,dec,wei,x,y,z,npt1,wei1,x1,y1,z1, &
                   nsep,sep,sbound,mxh1,mxh2,wfib,cntid,logf,sk1,ll1,cdth)
!===============================================================================
! NAME
!  th_C_wg()
!
! DESCRIPTION
!  Cross-count weighted data pairs in angular space for TWO samples of particles
!
! INPUTS
!  Variable---Type-------------Description--------------------------------------
!  nt         int              Number of thread to use (OpenMP)
!  npt        int              Number of particles
!  ra         r8(npt)          [sampleC] RA of particles [deg]
!  dec        r8(npt)          [sampleC] DEC of particles [deg]
!  wei        r4(npt)          [sampleC] WEIGHT of particles
!  x,y,z      r8(npt)          [sampleC] X,Y,Z coordinates of particles (see radec2xyz())
!  npt1       int              [sampleD] Number of particles
!  wei1       r4(npt1)         [sampleD] WEIGHT of particles
!  x1,y1,z1   r8(npt1)         [sampleD] X,Y,Z coordinates of particles (see radec2xyz())
!  nsep       int              Number of angular separation bins
!  sep        r8(nsep+1)       Bins in angular separation [deg]
!  sbound     r8(4)            Survey boundaries in RA,DEC. Form is (ramin,ramax,decmin,decmax)
!  mxh1       int              Nr of DEC cells of skip table
!  mxh2       int              Nr of RA cells of skip table
!  wfib       int              If 1, also apply fiber-correction weights
!  cntid      char*2           Two-character string to identify the samples,
!                              e.g. DD, RR, DR,... (used only for ouput messages)
!  logf       char*80          String for the file to log Fortran status messages
!  sk1        int(mxh2,mxh1)   [sampleD] Skip table (SK) contructed with skll2d()
!  ll1        int(npt1)        [sampleD] Linked list (LL) contructed with skll2d()
!
! OUTPUTS
!  Variable---Type-------------Description--------------------------------------
!  cdth       r8(nsep)         Counts in angular separation bins
!
! NOTES  -----------------------------------------------------------------------
!  1. RA limits should be set to ramin=0 and ramax=360.
!  2. Remember to update the declarations in cflibfor.pyf if you add/remove 
!     in/out parameters to this exposed subroutine.

implicit none
real(kind=8) :: ra(npt),dec(npt),ral,rau,decl,decu,hc1,hc2,shth2,stm2,dltra,dltdec
real(kind=8) :: x(npt),y(npt),z(npt),x1(npt1),y1(npt1),z1(npt1) 
real(kind=8) :: sep(nsep+1),sep2(nsep+1),sbound(4),sep2max,xi,yi,zi
real(kind=4) :: wei(npt),wei1(npt1),wi,wpp
real(kind=8) :: cdth(nsep)
integer      :: nt,nthr,sk1(mxh2,mxh1),ll1(npt1),mxh1,mxh2,npt,npt1,nsep,ndp
integer      :: fracp,dpart,nadv,nc1,nc2,i,ii,j,iq1,iq2,p1,p2,jq1,jq2,jq1m,wfib
integer      :: jq2m,jq2max,jq2min,jq2t
character    :: cntid*2,logf*80

!------Open log file---------------------------------
open(11,file=logf,action='write',access='append')

fracp = 0 ; dpart = 0 ; nadv = 0  !Reset progress counters

!----------------------------------------------------
! reset the counts, set max. ang. distance and square bins
cdth    = 0.0d0
dltdec  = sep(nsep+1)
sep2    = (sin(0.5*sep*deg2rad))**2
sep2max = sep2(nsep+1)

!-----------------------------------
! Unpack boundaries and get cell dimensions in (dec,ra)
ral = sbound(1) ; rau = sbound(2) ; decl = sbound(3) ; decu = sbound(4)
nc1 = mxh1
nc2 = mxh2
hc1 = (decu-decl)/float(nc1)   !effective nr of DEC cells
hc2 = (rau-ral)/float(nc2)     !effective nr of RA cells

!----------------------------------------------------
! Set the number of threads
if(nt<=0) then
   nthr = omp_get_num_procs()
else
   nthr = nt
endif
call omp_set_num_threads(nthr)

!----------------------------------------------------
! Some inits
dpart = int(float(npt)/float(mxh1)) !choose dpart so we get mhx1 parts
stm2  = sin(dltdec*0.5*deg2rad)
jq1m  = int(dltdec/hc1)+1

!----------------------------------------------------
! Count pairs in SK grid
write(*,*) ' '
write(*,fmt='(a,i3,a)') '====  Counting '//cntid//' pairs in ', mxh1, ' DEC strips'
!$omp parallel do reduction(+:cdth) default(shared) &
!$omp& private(iq1,iq2,jq1,jq2,dltra,jq2m,jq2min,jq2max) &
!$omp& private(jq2t,j,i,xi,yi,zi,wi,wpp,shth2,ii,p1,p2) &
!$omp& schedule(guided) firstprivate(fracp) if(nthr>1)

do i=1,npt
   iq1 = int((dec(i)-decl)/hc1)+1
   iq2 = int((ra(i)-ral)/hc2)+1
   fracp = fracp + 1  ! accumulate particles and check when above the step size
   
   if(fracp>=dpart) then
       !$omp critical
       nadv = nadv + 1
       !omp flush (nadv)
       p1   = (nadv-1)*dpart + 1
       p2   = nadv*dpart
       if((nadv+1)*dpart>npt) p2=npt
       !note we are not really counting in mxh1 strips, just mymicking
       write(*,fmt="(i4)",advance='no') nadv                        ! for screen
       write(11,*) cntid//' counting in DEC strip > ',nadv,' (',p1,'-',p2,')' !for disk
       fracp = 0
       !$omp end critical
   endif
   lp_jq1: do jq1=iq1-jq1m,iq1+jq1m
      if(jq1>nc1.or.jq1<1) cycle lp_jq1
      if(jq1==iq1) then
         dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,1)
      else
         dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,0)
      end if
      jq2m   = int(dltra/hc2)+1
      jq2max = iq2+jq2m
      jq2min = iq2-jq2m
      if(jq2max-jq2min+1>nc2) jq2max=jq2min-1+nc2
      lp_jq2: do jq2=jq2min,jq2max
         if(jq2>nc2) then
            jq2t = jq2-nc2
         else if(jq2<1) then
            jq2t = jq2+nc2
         else
            jq2t = jq2
         end if
         j = sk1(jq2t,jq1)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         wi = wei(i)

         do while(j/=0)
            shth2 = (xi-x1(j))**2 + (yi-y1(j))**2 + (zi-z1(j))**2
            !th = 2.*sqrt(shth2)*180./3.1415
            if(shth2<=sep2max) then
                wpp = wi*wei1(j)                    !weight by input
                if(wfib==1) then                    !weight by fiber
                    wpp = wpp*wfiber(shth2)  !TODO ???
                endif
                ! Now count the pair by finding its (ii) bin in vector of ang-space bins
                if(shth2>sep2(nsep)) then
                    cdth(nsep) = cdth(nsep) + wpp
                    goto 78
                endif
                if(shth2>sep2(nsep-1)) then
                    cdth(nsep-1) = cdth(nsep-1) + wpp
                    goto 78
                endif
                if(shth2>sep2(nsep-2)) then
                    cdth(nsep-2) = cdth(nsep-2) + wpp
                    goto 78
                endif
                if(shth2>sep2(nsep-3)) then
                    cdth(nsep-3) = cdth(nsep-3) + wpp
                    goto 78
                endif
                do ii=nsep-4,1,-1
                   if(shth2>sep2(ii)) then
                      cdth(ii) = cdth(ii) + wpp
                      goto 78
                   endif
                enddo
            endif
            78 j = ll1(j)
         end do
         
      end do lp_jq2
   end do lp_jq1
end do
!$omp end parallel do
close(11)  ! close log
write(*,*) ' '
end subroutine th_C_wg


subroutine th_Cb(nt,npt,ra,dec,x,y,z,npt1,x1,y1,z1, &
                 nsep,sep,sbound,mxh1,mxh2,nbts,bseed,cntid,logf,sk1,ll1,cdth,bcdth)
!===============================================================================
! NAME
!  th_Cb()
!
! DESCRIPTION
!  Cross-count data pairs in angular space for TWO samples of particles. Also
!  count booostraped data pairs
!
! INPUTS
!  Variable---Type-------------Description--------------------------------------
!  nt         int              Number of thread to use (OpenMP)
!  npt        int              Number of particles
!  ra         r8(npt)          [sampleC] RA of particles [deg]
!  dec        r8(npt)          [sampleC] DEC of particles [deg]
!  x,y,z      r8(npt)          [sampleC] X,Y,Z coordinates of particles (see radec2xyz())
!  npt1       int              [sampleD] Number of particles
!  x1,y1,z1   r8(npt1)         [sampleD] X,Y,Z coordinates of particles (see radec2xyz())
!  nsep       int              Number of angular separation bins
!  sep        r8(nsep+1)       Bins in angular separation [deg]
!  sbound     r8(4)            Survey boundaries in RA,DEC. Form is (ramin,ramax,decmin,decmax)
!  mxh1       int              Nr of DEC cells of skip table
!  mxh2       int              Nr of RA cells of skip table
!  nbts       int              Number of boostrap samples
!  bseed      int              Seed for RNG during boostrap resampling
!  cntid      char*2           Two-character string to identify the samples,
!                              e.g. DD, RR, DR,... (used only for ouput messages)
!  logf       char*80          String for the file to log Fortran status messages
!  sk1        int(mxh2,mxh1)   [sampleD] Skip table (SK) contructed with skll2d()
!  ll1        int(npt1)        [sampleD] Linked list (LL) contructed with skll2d()
!
! OUTPUTS
!  Variable---Type-------------Description--------------------------------------
!  cdth       r8(nsep)         Counts in angular separation bins
!  bcdth      r8(nbts,nsep)    Boostrap counts in angular separation bins
!
! NOTES  -----------------------------------------------------------------------
!  1. RA limits should be set to ramin=0 and ramax=360.
!  2. Remember to update the declarations in cflibfor.pyf if you add/remove 
!     in/out parameters to this exposed subroutine.

implicit none
real(kind=8) :: ra(npt),dec(npt),ral,rau,decl,decu,hc1,hc2,shth2,stm2,dltra,dltdec
real(kind=8) :: x(npt),y(npt),z(npt),x1(npt1),y1(npt1),z1(npt1) 
real(kind=8) :: sep(nsep+1),sep2(nsep+1),sbound(4),sep2max,xi,yi,zi
real(kind=4) :: wbts(nbts,npt),wbts1(nbts,npt1)
real(kind=8) :: cdth(nsep),bcdth(nbts,nsep)
integer      :: nt,nthr,sk1(mxh2,mxh1),ll1(npt1),mxh1,mxh2,npt,npt1,nsep,ndp
integer      :: fracp,dpart,nadv,nc1,nc2,i,ii,j,iq1,iq2,p1,p2,jq1,jq2,jq1m,nbts,bseed
integer      :: jq2m,jq2max,jq2min,jq2t
character    :: cntid*2,logf*80

!------Open log file---------------------------------
open(11,file=logf,action='write',access='append')

fracp = 0 ; dpart = 0 ; nadv = 0  !Reset progress counters

!----------------------------------------------------
! reset the counts, set max. ang. distance and square bins
cdth    = 0.0d0
dltdec  = sep(nsep+1)
sep2    = (sin(0.5*sep*deg2rad))**2
sep2max = sep2(nsep+1)

!-----------------------------------
! Unpack boundaries and get cell dimensions in (dec,ra)
ral = sbound(1) ; rau = sbound(2) ; decl = sbound(3) ; decu = sbound(4)
nc1 = mxh1
nc2 = mxh2
hc1 = (decu-decl)/float(nc1)   !effective nr of DEC cells
hc2 = (rau-ral)/float(nc2)     !effective nr of RA cells

!----------------------------------------------------
! Generate bootstrap samples
call bootstrap(npt,nbts,bseed,wbts)
call bootstrap(npt1,nbts,bseed,wbts1)

!----------------------------------------------------
! Set the number of threads
if(nt<=0) then
   nthr = omp_get_num_procs()
else
   nthr = nt
endif
call omp_set_num_threads(nthr)

!----------------------------------------------------
! Some inits
dpart = int(float(npt)/float(mxh1)) !choose dpart so we get mhx1 parts
stm2  = sin(dltdec*0.5*deg2rad)
jq1m  = int(dltdec/hc1)+1

!----------------------------------------------------
! Count pairs in SK grid
write(*,*) ' '
write(*,fmt='(a,i3,a)') '====  Counting '//cntid//' pairs in ', mxh1, ' DEC strips'
!$omp parallel do reduction(+:cdth,bcdth) default(shared) &
!$omp& private(iq1,iq2,jq1,jq2) &
!$omp& private(dltra,jq2m,jq2min,jq2max,jq2t,j,i,xi,yi,zi,shth2,ii,p1,p2) &
!$omp& schedule(guided) firstprivate(fracp) if(nthr>1)

do i=1,npt
   iq1 = int((dec(i)-decl)/hc1)+1
   iq2 = int((ra(i)-ral)/hc2)+1
   fracp = fracp + 1  ! accumulate particles and check when above the step size
   if(fracp>=dpart) then
       !$omp critical
       nadv = nadv + 1
       !omp flush (nadv)
       p1   = (nadv-1)*dpart + 1
       p2   = nadv*dpart
       if((nadv+1)*dpart>npt) p2=npt
       !note we are not really counting in mxh1 strips, just mymicking
       write(*,fmt="(i4)",advance='no') nadv                        ! for screen
       write(11,*) cntid//' counting in DEC strip > ',nadv,' (',p1,'-',p2,')' !for disk
       flush(11)
       fracp = 0
       !$omp end critical
   endif
   
   lp_jq1: do jq1=iq1-jq1m,iq1+jq1m
      if(jq1>nc1.or.jq1<1) cycle lp_jq1
      if(jq1==iq1) then
         dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,1)
      else
         dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,0)
      end if
      jq2m   = int(dltra/hc2)+1
      jq2max = iq2+jq2m
      jq2min = iq2-jq2m
      if(jq2max-jq2min+1>nc2) jq2max=jq2min-1+nc2
      lp_jq2: do jq2=jq2min,jq2max
         if(jq2>nc2) then
            jq2t = jq2-nc2
         else if(jq2<1) then
            jq2t = jq2+nc2
         else
            jq2t = jq2
         end if
         j = sk1(jq2t,jq1)
         xi = x(i)
         yi = y(i)
         zi = z(i)

         do while(j/=0)
            shth2 = (xi-x1(j))**2 + (yi-y1(j))**2 + (zi-z1(j))**2
            !th = 2.*sqrt(shth2)*180./3.1415
            if(shth2<=sep2max) then
                ! Now count the pair by finding its (ii) bin in vector of ang-space bins
                if(shth2>sep2(nsep)) then
                    cdth(nsep) = cdth(nsep) + 1.0d0
                    bcdth(:,nsep) = bcdth(:,nsep) + wbts(:,i)*wbts1(:,j)
                    goto 78
                endif
                if(shth2>sep2(nsep-1)) then
                    cdth(nsep-1) = cdth(nsep-1) + 1.0d0
                    bcdth(:,nsep-1) = bcdth(:,nsep-1) + wbts(:,i)*wbts1(:,j)
                    goto 78
                endif
                if(shth2>sep2(nsep-2)) then
                    cdth(nsep-2) = cdth(nsep-2) + 1.0d0
                    bcdth(:,nsep-2) = bcdth(:,nsep-2) + wbts(:,i)*wbts1(:,j)
                    goto 78
                endif
                if(shth2>sep2(nsep-3)) then
                    cdth(nsep-3) = cdth(nsep-3) + 1.0d0
                    bcdth(:,nsep-3) = bcdth(:,nsep-3) + wbts(:,i)*wbts1(:,j)
                    goto 78
                endif
                do ii=nsep-4,1,-1
                   if(shth2>sep2(ii)) then
                      cdth(ii) = cdth(ii) + 1.0d0
                      bcdth(:,ii) = bcdth(:,ii) + wbts(:,i)*wbts1(:,j)
                      goto 78
                   endif
                enddo
            endif
            78 j = ll1(j)
         end do
         
      end do lp_jq2
   end do lp_jq1
end do
!$omp end parallel do
close(11)  ! close log
write(*,*) ' '
end subroutine th_Cb


subroutine th_Cb_wg(nt,npt,ra,dec,wei,x,y,z,npt1,wei1,x1,y1,z1, &
                    nsep,sep,sbound,mxh1,mxh2,nbts,bseed,wfib,cntid,logf,sk1,ll1,cdth,bcdth)
!===============================================================================
! NAME
!  th_Cb_wg()
!
! DESCRIPTION
!  Cross-count weighted data pairs in angular space for TWO samples of particles.
!  Also count boostraped data pairs.
!
! INPUTS
!  Variable---Type-------------Description--------------------------------------
!  nt         int              Number of thread to use (OpenMP)
!  npt        int              Number of particles
!  ra         r8(npt)          [sampleC] RA of particles [deg]
!  dec        r8(npt)          [sampleC] DEC of particles [deg]
!  wei        r4(npt)          [sampleC] WEIGHT of particles
!  x,y,z      r8(npt)          [sampleC] X,Y,Z coordinates of particles (see radec2xyz())
!  npt1       int              [sampleD] Number of particles
!  wei1       r4(npt1)         [sampleD] WEIGHT of particles
!  x1,y1,z1   r8(npt1)         [sampleD] X,Y,Z coordinates of particles (see radec2xyz())
!  nsep       int              Number of angular separation bins
!  sep        r8(nsep+1)       Bins in angular separation [deg]
!  sbound     r8(4)            Survey boundaries in RA,DEC. Form is (ramin,ramax,decmin,decmax)
!  mxh1       int              Nr of DEC cells of skip table
!  mxh2       int              Nr of RA cells of skip table
!  nbts       int              Number of boostrap samples
!  bseed      int              Seed for RNG during boostrap resampling
!  wfib       int              If 1, also apply fiber-correction weights
!  cntid      char*2           Two-character string to identify the samples,
!                              e.g. DD, RR, DR,... (used only for ouput messages)
!  logf       char*80          String for the file to log Fortran status messages
!  sk1        int(mxh2,mxh1)   [sampleD] Skip table (SK) contructed with skll2d()
!  ll1        int(npt1)        [sampleD] Linked list (LL) contructed with skll2d()
!
! OUTPUTS
!  Variable---Type-------------Description--------------------------------------
!  cdth       r8(nsep)         Counts in angular separation bins
!  bcdth      r8(nbts,nsep)    Boostrap counts in angular separation bins
!
! NOTES  -----------------------------------------------------------------------
!  1. RA limits should be set to ramin=0 and ramax=360.
!  2. Remember to update the declarations in cflibfor.pyf if you add/remove 
!     in/out parameters to this exposed subroutine.

implicit none
real(kind=8) :: ra(npt),dec(npt),ral,rau,decl,decu,hc1,hc2,shth2,stm2,dltra,dltdec
real(kind=8) :: x(npt),y(npt),z(npt),x1(npt1),y1(npt1),z1(npt1) 
real(kind=8) :: sep(nsep+1),sep2(nsep+1),sbound(4),sep2max,xi,yi,zi
real(kind=4) :: wbts(nbts,npt),wbts1(nbts,npt1)
real(kind=4) :: wei(npt),wei1(npt1),wi,wpp
real(kind=8) :: cdth(nsep),bcdth(nbts,nsep)
integer      :: nt,nthr,sk1(mxh2,mxh1),ll1(npt1),mxh1,mxh2,npt,npt1,nsep,ndp
integer      :: fracp,dpart,nadv,nc1,nc2,i,ii,j,iq1,iq2,p1,p2,jq1,jq2,jq1m,nbts,bseed,wfib
integer      :: jq2m,jq2max,jq2min,jq2t
character    :: cntid*2,logf*80

!------Open log file---------------------------------
open(11,file=logf,action='write',access='append')

fracp = 0 ; dpart = 0 ; nadv = 0  !Reset progress counters

!----------------------------------------------------
! reset the counts, set max. ang. distance and square bins
cdth    = 0.0d0
bcdth   = 0.0d0
dltdec  = sep(nsep+1)
sep2    = (sin(0.5*sep*deg2rad))**2
sep2max = sep2(nsep+1)

!-----------------------------------
! Unpack boundaries and get cell dimensions in (dec,ra)
ral = sbound(1) ; rau = sbound(2) ; decl = sbound(3) ; decu = sbound(4)
nc1 = mxh1
nc2 = mxh2
hc1 = (decu-decl)/float(nc1)   !effective nr of DEC cells
hc2 = (rau-ral)/float(nc2)     !effective nr of RA cells

!----------------------------------------------------
! Generate bootstrap samples
call bootstrap(npt,nbts,bseed,wbts)
call bootstrap(npt1,nbts,bseed,wbts1)

!----------------------------------------------------
! Set the number of threads
if(nt<=0) then
   nthr = omp_get_num_procs()
else
   nthr = nt
endif
call omp_set_num_threads(nthr)

!----------------------------------------------------
! Some inits
dpart = int(float(npt)/float(mxh1)) !choose dpart so we get mhx1 parts
stm2  = sin(dltdec*0.5*deg2rad)
jq1m  = int(dltdec/hc1)+1

!----------------------------------------------------
! Count pairs in SK grid
write(*,*) ' '
write(*,fmt='(a,i3,a)') '====  Counting '//cntid//' pairs in ', mxh1, ' DEC strips'
!$omp parallel do reduction(+:cdth,bcdth) default(shared) &
!$omp& private(iq1,iq2,jq1,jq2) &
!$omp& private(dltra,jq2m,jq2min,jq2max,jq2t,j,i,xi,yi,zi,wi,wpp,shth2,ii,p1,p2) &
!$omp& schedule(guided) firstprivate(fracp) if(nthr>1)

do i=1,npt
   iq1 = int((dec(i)-decl)/hc1)+1
   iq2 = int((ra(i)-ral)/hc2)+1
   fracp = fracp + 1  ! accumulate particles and check when above the step size
   if(fracp>=dpart) then
       !$omp critical
       nadv = nadv + 1
       !omp flush (nadv)
       p1   = (nadv-1)*dpart + 1
       p2   = nadv*dpart
       if((nadv+1)*dpart>npt) p2=npt
       !note we are not really counting in mxh1 strips, just mymicking
       write(*,fmt="(i4)",advance='no') nadv                        ! for screen
       write(11,*) cntid//' counting in DEC strip > ',nadv,' (',p1,'-',p2,')' !for disk
       flush(11)
       fracp = 0
       !$omp end critical
   endif
   
   lp_jq1: do jq1=iq1-jq1m,iq1+jq1m
      if(jq1>nc1.or.jq1<1) cycle lp_jq1
      if(jq1==iq1) then
         dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,1)
      else
         dltra = dalp(stm2,dec(i),dec(i)*deg2rad,jq1,hc1,decl,0)
      end if
      jq2m   = int(dltra/hc2)+1
      jq2max = iq2+jq2m
      jq2min = iq2-jq2m
      if(jq2max-jq2min+1>nc2) jq2max=jq2min-1+nc2
      lp_jq2: do jq2=jq2min,jq2max
         if(jq2>nc2) then
            jq2t = jq2-nc2
         else if(jq2<1) then
            jq2t = jq2+nc2
         else
            jq2t = jq2
         end if
         j = sk1(jq2t,jq1)
         xi = x(i)
         yi = y(i)
         zi = z(i)
         wi = wei(i)

         do while(j/=0)
            shth2 = (xi-x1(j))**2 + (yi-y1(j))**2 + (zi-z1(j))**2
            !th = 2.*sqrt(shth2)*180./3.1415
            if(shth2<=sep2max) then
                wpp = wi*wei1(j)                    !weight by input
                if(wfib==1) then                    !weight by fiber
                    wpp = wpp*wfiber(shth2)  !TODO ???
                endif
                ! Now count the pair by finding its (ii) bin in vector of ang-space bins
                if(shth2>sep2(nsep)) then
                    cdth(nsep) = cdth(nsep) + wpp
                    bcdth(:,nsep) = bcdth(:,nsep) + wpp*wbts(:,i)*wbts1(:,j)
                    goto 78
                endif
                if(shth2>sep2(nsep-1)) then
                    cdth(nsep-1) = cdth(nsep-1) + wpp
                    bcdth(:,nsep-1) = bcdth(:,nsep-1) + wpp*wbts(:,i)*wbts1(:,j)
                    goto 78
                endif
                if(shth2>sep2(nsep-2)) then
                    cdth(nsep-2) = cdth(nsep-2) + wpp
                    bcdth(:,nsep-2) = bcdth(:,nsep-2) + wpp*wbts(:,i)*wbts1(:,j)
                    goto 78
                endif
                if(shth2>sep2(nsep-3)) then
                    cdth(nsep-3) = cdth(nsep-3) + wpp
                    bcdth(:,nsep-3) = bcdth(:,nsep-3) + wpp*wbts(:,i)*wbts1(:,j)
                    goto 78
                endif
                do ii=nsep-4,1,-1
                   if(shth2>sep2(ii)) then
                      cdth(ii) = cdth(ii) + wpp
                      bcdth(:,ii) = bcdth(:,ii) + wpp*wbts(:,i)*wbts1(:,j)
                      goto 78
                   endif
                enddo
            endif
            78 j = ll1(j)
         end do
         
      end do lp_jq2
   end do lp_jq1
end do
!$omp end parallel do
close(11)  ! close log
write(*,*) ' '
end subroutine th_Cb_wg


end module mod

