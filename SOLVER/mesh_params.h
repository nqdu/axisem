! Proc   2: Header for mesh information to run static solver
! created by the mesher on 10/18/2026, at 05h 22min
 
!:::::::::::::::::::: Input parameters :::::::::::::::::::::::::::
!   Background model     :            external
!   Dominant period [s]  :   20.0000
!   Elements/wavelength  :    1.5000
!   Courant number       :    0.6000
!   Coarsening levels    :         3
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 
 integer, parameter ::         npol =         4  !            polynomial order
 integer, parameter ::        nelem =     13311  !                   proc. els
 integer, parameter ::       npoint =    332775  !               proc. all pts
 integer, parameter ::    nel_solid =     11151  !             proc. solid els
 integer, parameter ::    nel_fluid =      2160  !             proc. fluid els
 integer, parameter :: npoint_solid =    278775  !             proc. solid pts
 integer, parameter :: npoint_fluid =     54000  !             proc. fluid pts
 integer, parameter ::  nglob_fluid =     35649  !            proc. flocal pts
 integer, parameter ::     nel_bdry =       162  ! proc. solid-fluid bndry els
 integer, parameter ::        ndisc =         3  !   # disconts in bkgrd model
 integer, parameter ::   nproc_mesh =         2  !        number of processors
 integer, parameter :: lfbkgrdmodel =         8  !   length of bkgrdmodel name
 
!:::::::::::::::::::: Output parameters ::::::::::::::::::::::::::
!   Time step [s]        :    0.1274
!   Min(h/vp),dt/courant :    1.2662    0.8494
!   max(h/vs),T0/wvlngth :   13.2086   13.3333
!   Inner core r_min [km]:  799.5618
!   Max(h) r/ns(icb) [km]:   35.4011
!   Max(h) precalc.  [km]:   35.6947
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 
