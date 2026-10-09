!
!    Copyright 2013, Tarje Nissen-Meyer, Alexandre Fournier, Martin van Driel
!                    Simon Stähler, Kasra Hosseini, Stefanie Hempel
!
!    This file is part of AxiSEM.
!    It is distributed from the webpage <http://www.axisem.info>
!
!    AxiSEM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    AxiSEM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with AxiSEM.  If not, see <http://www.gnu.org/licenses/>.
!

!=========================================================================================
module utlity

  use global_parameters
  implicit none
  
  public :: compute_coordinates, scoord, zcoord, rcoord, thetacoord
  public :: dblreldiff_small, reldiff_small
  public :: dblereldiff, reldiff
  public :: dbleabsreldiff, absreldiff
  public :: to_lower

  !nqdu
  public :: inside_element,lagrange_interpol_2D_td
  private

contains

!-----------------------------------------------------------------------------------------
pure logical function dblreldiff_small(x1,x2)

  real(kind=dp), intent(in) :: x1,x2

  dblreldiff_small = .false.

  if (x1 /= zero) then 
     if (abs((x1-x2)/x1) <= smallval_dble) dblreldiff_small = .true.
  elseif (x2 /=zero) then
     if (abs((x1-x2)/x2) <= smallval_dble) dblreldiff_small = .true.
  else
     dblreldiff_small = .true.
  endif

end function dblreldiff_small
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure logical function reldiff_small(x1,x2)

  real(kind=realkind), intent(in) :: x1,x2
  real(kind=realkind)             ::  smallval1

  if (realkind==sp) smallval1 = smallval_sngl
  if (realkind==dp) smallval1 = smallval_dble

  reldiff_small = .false.

  if (x1 /= zero) then 
     if (abs((x1-x2)/x1) <= smallval1) reldiff_small = .true.
  elseif (x2 /=zero) then
     if (abs((x1-x2)/x2) <= smallval1) reldiff_small = .true.
  else
     reldiff_small = .true.
  endif

end function reldiff_small
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=realkind) function reldiff(x1,x2)

  real(kind=realkind), intent(in) :: x1,x2

  if (x1/=zero) then
     reldiff=(x1-x2)/x1
  elseif (x2/=zero) then
     reldiff=(x1-x2)/x2
  else
     reldiff=zero
  endif

end function reldiff
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp) function dblereldiff(x1,x2)

  real(kind=dp), intent(in) :: x1,x2

  if (x1/=zero) then
     dblereldiff=(x1-x2)/x1
  elseif (x2/=zero) then 
     dblereldiff=(x1-x2)/x2
  else
     dblereldiff=zero
  endif

end function dblereldiff
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=realkind) function absreldiff(x1,x2)

  real(kind=realkind), intent(in) :: x1,x2

  if (x1/=zero) then
     absreldiff=abs((x1-x2)/x1)
  elseif (x2/=zero) then 
     absreldiff=abs((x1-x2)/x2)
  else
     absreldiff=zero
  endif

end function absreldiff
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp) function dbleabsreldiff(x1,x2)

  real(kind=dp), intent(in) :: x1,x2

  if (x1/=zero) then
     dbleabsreldiff=abs((x1-x2)/x1)
  elseif (x2/=zero) then 
     dbleabsreldiff=abs((x1-x2)/x2)
  else
     dbleabsreldiff=zero
  endif

end function dbleabsreldiff
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure subroutine compute_coordinates(s,z,r,theta,ielem,ipol,jpol)
!< Given the elemental grid point index, outputs s,z,r,theta coordinate [m,rad].
!! These coordinates are by default ALWAYS global (no solid or fluid domains).
  
  use data_mesh,            only: min_distance_dim
  use data_mesh,            only: lnods, crd_nodes, axis
  use data_spec,            only: xi_k, eta
  use analytic_mapping,     only: mapping
  
  real(kind=dp), intent(out)    :: s,z,r,theta
  integer, intent(in)           :: ielem,ipol,jpol
  integer                       :: ipt,inode
  real(kind=dp)                 :: nodes_crd(8,2)

  do inode = 1, 8
     ipt = lnods(ielem,inode)
     nodes_crd(inode,1) = crd_nodes(ipt,1)
     nodes_crd(inode,2) = crd_nodes(ipt,2)
  end do

  ! Fill global coordinate array
  if ( axis(ielem) ) then 

     s= mapping( xi_k(ipol),eta(jpol),nodes_crd,1,ielem)
     z= mapping( xi_k(ipol),eta(jpol),nodes_crd,2,ielem)

  else 
     s= mapping(eta(ipol),eta(jpol),nodes_crd,1,ielem)
     z= mapping(eta(ipol),eta(jpol),nodes_crd,2,ielem)
  end if

  ! Eliminate roundoff errors
  if (abs(s) < min_distance_dim) s=zero
  if (abs(z) < min_distance_dim) z=zero

  r = dsqrt(s**2+z**2)
  theta = datan(s/(z+epsi))
  if ( zero > theta ) theta = pi + theta
  if (theta == zero .and. z < 0) theta = pi

end subroutine compute_coordinates
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp) function scoord(ipol,jpol,ielem)
!< Given the elemental grid point index, outputs the s coordinate [m].
!! These coordinates are by default ALWAYS global (no solid or fluid domains).
  
  use data_mesh,            only: min_distance_dim
  use data_mesh,            only: lnods, crd_nodes, axis
  use data_spec,            only: xi_k, eta
  use analytic_mapping,     only: mapping
  
  integer, intent(in)  :: ielem, ipol, jpol
  integer              :: ipt, inode
  real(kind=dp)        :: nodes_crd(8,2)

  do inode = 1, 8
     ipt = lnods(ielem,inode)
     nodes_crd(inode,1) = crd_nodes(ipt,1)
     nodes_crd(inode,2) = crd_nodes(ipt,2)
  end do

  ! Fill global coordinate array
  if ( axis(ielem) ) then 
     scoord = mapping( xi_k(ipol),eta(jpol),nodes_crd,1,ielem)
  else 
     scoord = mapping(eta(ipol),eta(jpol),nodes_crd,1,ielem)
  end if

  ! Eliminate roundoff errors
  if (abs(scoord) < min_distance_dim) scoord=zero

end function scoord
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp)    function zcoord(ipol,jpol,ielem)
!< Given the elemental grid point index, outputs the z coordinate [m].
!! These coordinates are by default ALWAYS global (no solid or fluid domains).
  
  use data_mesh,            only: min_distance_dim
  use data_mesh,            only: lnods, crd_nodes, axis
  use data_spec,            only: xi_k, eta
  use analytic_mapping,     only: mapping
  
  integer, intent(in)  :: ielem, ipol, jpol
  integer              :: ipt, inode
  real(kind=dp)        :: nodes_crd(8,2)

  do inode = 1, 8
     ipt = lnods(ielem,inode)
     nodes_crd(inode,1) = crd_nodes(ipt,1)
     nodes_crd(inode,2) = crd_nodes(ipt,2)
  end do

  ! Fill global coordinate array
  if ( axis(ielem) ) then 
     zcoord = mapping( xi_k(ipol),eta(jpol),nodes_crd,2,ielem)
  else 
     zcoord = mapping(eta(ipol),eta(jpol),nodes_crd,2,ielem)
  end if

  ! Eliminate roundoff errors
  if (abs(zcoord) < min_distance_dim) zcoord=zero

end function zcoord
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp)    function rcoord(ipol,jpol,ielem)
!< Given the elemental grid point index, outputs the radius coordinate [m].
!! These coordinates are by default ALWAYS global (no solid or fluid domains).
  
  use data_mesh,            only: min_distance_dim
  use data_mesh,            only: lnods, crd_nodes, axis
  use data_spec,            only: xi_k, eta
  use analytic_mapping,     only: mapping
  
  integer, intent(in)  :: ielem, ipol, jpol
  integer              :: ipt, inode
  real(kind=dp)        :: nodes_crd(8,2),s,z

  do inode = 1, 8
     ipt = lnods(ielem,inode)
     nodes_crd(inode,1) = crd_nodes(ipt,1)
     nodes_crd(inode,2) = crd_nodes(ipt,2)
  end do

  ! Fill global coordinate array
  if ( axis(ielem) ) then 
     s = mapping( xi_k(ipol),eta(jpol),nodes_crd,1,ielem)
     z = mapping( xi_k(ipol),eta(jpol),nodes_crd,2,ielem)
  else 
     s = mapping(eta(ipol),eta(jpol),nodes_crd,1,ielem)
     z = mapping(eta(ipol),eta(jpol),nodes_crd,2,ielem)
  end if

  rcoord = sqrt(s**2 + z**2)
  ! Eliminate roundoff errors
  if (abs(rcoord) < min_distance_dim) rcoord=zero

end function rcoord
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
pure real(kind=dp) function thetacoord(ipol,jpol,ielem)
!< Given the elemental grid point index, outputs the theta coordinate [rad].
!! These coordinates are by default ALWAYS global (no solid or fluid domains).
  
  use data_mesh,            only: min_distance_dim
  use data_mesh,            only: lnods, crd_nodes,axis
  use data_spec,            only: xi_k, eta
  use analytic_mapping,     only: mapping
  
  integer, intent(in)  :: ielem, ipol, jpol
  integer              :: ipt, inode
  real(kind=dp)        :: nodes_crd(8,2),s,z

  do inode = 1, 8
     ipt = lnods(ielem,inode)
     nodes_crd(inode,1) = crd_nodes(ipt,1)
     nodes_crd(inode,2) = crd_nodes(ipt,2)
  end do

  ! Fill global coordinate array
  if ( axis(ielem) ) then 
     s = mapping(xi_k(ipol), eta(jpol), nodes_crd, 1, ielem)
     z = mapping(xi_k(ipol), eta(jpol), nodes_crd, 2, ielem)
  else 
     s = mapping(eta(ipol), eta(jpol), nodes_crd, 1, ielem)
     z = mapping(eta(ipol), eta(jpol), nodes_crd, 2, ielem)
  end if

  thetacoord = datan(s/(z+epsi))
  if ( zero > thetacoord ) thetacoord = pi + thetacoord
  if (thetacoord == zero .and. z < 0) thetacoord = pi

end function thetacoord
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
function to_lower(strIn) result(strOut)
!< Converts string to lowercase, adapted from http://www.star.le.ac.uk/~cgp/fortran.html
    implicit none

    character(len=*), intent(in) :: strIn
    character(len=len(strIn))    :: strOut
    integer                      :: i,j

    do i = 1, len(strIn)
        j = iachar(strIn(i:i))
        if (j>= iachar("A") .and. j<=iachar("Z") ) then
            strOut(i:i) = achar(iachar(strIn(i:i))+32)
        else
            strOut(i:i) = strIn(i:i)
        end if
    end do

end function to_lower
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------

function lagrange_interpol_2D_td(points1, points2, coefficients, x1, x2) result(s)
!> computes the Lagrangian interpolation polynomial of a function defined by its values at
!  a set of collocation points in 2D, where the points are a tensorproduct of two sets of
!  points in 1D, for time dependent coefficients
   implicit none

   real(kind=dp), intent(in)  :: points1(0:), points2(0:)
   real(kind=dp), intent(in)  :: coefficients(0:,0:)
   real(kind=dp), intent(in)  :: x1, x2
   real(kind=dp)              :: s
   real(kind=dp)              :: l_i(0:size(points1)-1), l_j(0:size(points2)-1)

   integer               :: i, j, m1, m2, n1, n2

   n1 = size(points1) - 1
   n2 = size(points2) - 1

   do i=0, n1
      l_i(i) = 1
      do m1=0, n1
         if (m1 == i) cycle
         l_i(i) = l_i(i) * (x1 - points1(m1)) / (points1(i) - points1(m1))
      enddo
   enddo

   do j=0, n2
      l_j(j) = 1
      do m2=0, n2
         if (m2 == j) cycle
         l_j(j) = l_j(j) * (x2 - points2(m2)) / (points2(j) - points2(m2))
      enddo
   enddo

   s = 0.0_dp

   do i=0, n1
      do j=0, n2
         s = s + coefficients(i,j) * l_i(i) * l_j(j)
      enddo
   enddo

end function lagrange_interpol_2D_td

!-----------------------------------------------------------------------------------------
subroutine inside_element(s,z,iel,xi,eta,sloc,zloc,in_it)
! check if (s,z) is in element (iel), if true, also compute it's local coordinates (xi,eta)
   use data_mesh,only : crd_nodes, lnods
   use analytic_mapping, only : mapping,compute_partial_derivatives
   implicit none

   real(kind=dp),intent(in)  :: s,z
   real(kind=dp),intent(out) :: xi,eta,sloc,zloc
   integer, intent(in)       :: iel
   logical,intent(out)       :: in_it

   !local variables
   integer                   :: i
   integer,parameter         :: maxiter = 20
   real(kind=dp),parameter   :: tol = 1.0e-3
   real(kind=dp)             :: nodes(8,2),ds,dz,dist,jaco_det,jaco(2,2),inv_jaco(2,2)

   in_it = .false.
   ! compute control coordinates in this element
   nodes(:,1) = crd_nodes(lnods(iel,:),1)
   nodes(:,2) = crd_nodes(lnods(iel,:),2)

   ! start value
   xi = 0.0_dp
   eta = 0.0_dp

   do i=1,maxiter
     sloc = mapping(xi,eta,nodes,1,iel)
     zloc = mapping(xi,eta,nodes,2,iel)
     ds = s - sloc
     dz = z - zloc

     ! check convergence
     dist = hypot(ds,dz)
     if (dist < 1.0e-7_dp * max(1.0_dp, hypot(s,z))) then
       exit
     endif

     ! update
     call compute_partial_derivatives(jaco(1,1),jaco(2,1),jaco(1,2),jaco(2,2),xi,eta,nodes,iel)
     jaco_det = jaco(1,1) * jaco(2,2) - jaco(1,2) * jaco(2,1)
     if (abs(jaco_det) <= tiny(jaco_det)) exit
     inv_jaco(1,1) = jaco(2,2) / jaco_det
     inv_jaco(2,1) = -jaco(2,1) / jaco_det
     inv_jaco(1,2) = -jaco(1,2) / jaco_det
     inv_jaco(2,2) = jaco(1,1) / jaco_det
     xi  =  xi + inv_jaco(1,1) * ds + inv_jaco(1,2) * dz
     eta = eta + inv_jaco(2,1) * ds + inv_jaco(2,2) * dz

     if (abs(xi) > 1.0 + tol) xi = xi / abs(xi) * 1.1
     if (abs(eta) > 1.0 + tol) eta = sign(1.1_dp, eta)
   enddo
   sloc = mapping(xi,eta,nodes,1,iel)
   zloc = mapping(xi,eta,nodes,2,iel)

   ! check inside this element
   in_it = (hypot(s-sloc,z-zloc) <= 1.0e-7_dp * max(1.0_dp,hypot(s,z)) .and. &
           xi >= -1 - tol .and. &
           xi <= 1 + tol .and. &
           eta >= -1 - tol .and. &
           eta <= 1 + tol)

end subroutine inside_element

end module utlity
!=========================================================================================
