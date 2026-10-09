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
module data_mesh

  ! Arrays here pertain to some sort of mesh peculiarities and mainly serve 
  ! as information or parameters for many "if"-decisions such as
  ! axis, north, element type, solid-fluid boundary mapping, coarsening,
  ! and specifically also related to the background model such as
  ! solid-fluid boundary mapping, discontinuities, and element arrays to 
  ! map solid/fluid to global element domains.
  ! These quantities are in active memory throughout the simulation.
  ! 
  ! Any global arrays containing properties inside elements are defined in data_matr.
  
  use global_parameters
  use kdtree2_module, only: kdkind, kdtree2, kdtree2_result, &
                            kdtree2_create, kdtree2_destroy, kdtree2_n_nearest
  implicit none
  
  public 

  ! Very basic mesh parameters, have been in mesh_params.h before
  integer, protected ::         npol !<            polynomial order
  integer, protected ::        nelem !<                   proc. els
  integer, protected ::       npoint !<               proc. all pts
  integer, protected ::    nel_solid !<             proc. solid els
  integer, protected ::    nel_fluid !<             proc. fluid els
  integer, protected :: npoint_solid !<             proc. solid pts
  integer, protected :: npoint_fluid !<             proc. fluid pts
  integer, protected ::  nglob_solid !<            proc. slocal pts
  integer, protected ::  nglob_fluid !<            proc. flocal pts
  integer, protected ::     nel_bdry !< proc. solid-fluid bndry els
  integer, protected ::        ndisc !<   # disconts in bkgrd model
  integer, protected ::   nproc_mesh !<        number of processors
  integer, protected :: lfbkgrdmodel !<   length of bkgrdmodel name

  ! global number in solid varies across procs due to central cube domain decomposition
  integer                                       :: nglob
  ! global numbering array for the solid and fluid assembly
  integer, protected, allocatable, dimension(:) :: igloc_solid ! (npoint_solid)
  integer, protected, allocatable, dimension(:) :: igloc_fluid ! (npoint_fluid)

  ! Element centers and search trees used while locating receivers and faces.
  integer, parameter :: nearest_element_count = 10
  real(kind=kdkind), pointer :: element_midpoint(:,:) => null()
  real(kind=kdkind), pointer :: solid_midpoint(:,:) => null()
  real(kind=kdkind), pointer :: fluid_midpoint(:,:) => null()
  type(kdtree2), pointer :: element_tree => null()
  type(kdtree2), pointer :: solid_tree => null()
  type(kdtree2), pointer :: fluid_tree => null()
  integer, allocatable :: solid_element_index(:)
  integer, allocatable :: fluid_element_index(:)

  ! Misc definitions
  integer                           :: nsize
  logical                           :: do_mesh_tests

  real(kind=realkind), allocatable  :: gvec_solid(:,:) 
  real(kind=realkind), allocatable  :: gvec_fluid(:)

  ! Deprecated elemental mesh (radius & colatitude of elemental midpoint)
  ! This is used for blow up localization. Might want to remove this when 
  ! everything is running smoothly in all eternities...
  real(kind=realkind), allocatable  :: mean_rad_colat_solid(:,:)
  real(kind=realkind), allocatable  :: mean_rad_colat_fluid(:,:)
  
  ! Global mesh informations
  real(kind=dp)                     :: router ! Outer radius (surface)
  character(len=100)                :: model_name_ext_model ! name of external model

  ! critical mesh parameters (spacing/velocity, characteristic lead time etc)
  real(kind=dp)                     :: pts_wavelngth
  real(kind=dp)                     :: hmin_glob, hmax_glob
  real(kind=dp)                     :: min_distance_dim, min_distance_nondim
  real(kind=dp)                     :: char_time_max
  integer                           :: char_time_max_globel
  real(kind=dp)                     :: char_time_max_rad, char_time_max_theta
  real(kind=dp)                     :: char_time_min
  integer                           :: char_time_min_globel
  real(kind=dp)                     :: char_time_min_rad, char_time_min_theta
  real(kind=dp)                     :: vpmin, vsmin, vpmax, vsmax
  real(kind=dp)                     :: vpminr, vsminr, vpmaxr, vsmaxr
  integer, dimension(3)             :: vpminloc, vsminloc, vpmaxloc, vsmaxloc

  !----------------------------------------------------------------------
  ! Axial elements
  integer                           :: naxel, naxel_solid, naxel_fluid
  integer, protected, allocatable   :: ax_el(:), ax_el_solid(:), ax_el_fluid(:)
  logical,            allocatable   :: axis_solid(:)
  logical,            allocatable   :: axis_fluid(:)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Background Model related
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  ! Solid-Fluid boundary----------------------------------------------------

  ! mapping nel_bdry elements into the solid/fluid element numbers
  integer, protected, allocatable :: bdry_solid_el(:), bdry_fluid_el(:)

  ! mapping the z coordinate of the boundary for each element 
  ! (depending on north/south, above/below)
  integer, protected, allocatable :: bdry_jpol_solid(:), bdry_jpol_fluid(:)

  ! Boolean to determine whether proc has solid-fluid boundary elements
  logical :: have_bdry_elem 

  !Not used anywhere
  !! integer array of size nel_bdry containing the "global" element number 
  !! for 1:nel_bdry 
  !integer, dimension(nel_bdry) :: ibdryel
  

  ! Background model--------------------------------------------------------
  character(len=100)          :: bkgrdmodel
  character(len=100)          :: meshname
  logical                     :: have_fluid
  real(kind=dp), allocatable  :: discont(:)
  logical, allocatable        :: solid_domain(:)
  integer, allocatable        :: idom_fluid(:)
  real(kind=dp)               :: rmin, minh_ic, maxh_ic, maxh_icb
  logical                     :: anel_true ! anelastic model?
  !--------------------------------------------------------------------------

  ! Receiver locations
  integer                      :: maxind_glob, maxind, ind_first, ind_last
  integer                      :: num_rec, num_rec_tot
  integer, allocatable         :: surfelem(:), jsurfel(:)
  real(kind=sp), allocatable   :: surfcoord(:)
  integer                      :: ielepi, ielantipode, ielequ
  integer, allocatable         :: recfile_el(:,:), loc2globrec(:)
  ! Tensor-product interpolation weights for each receiver's containing element.
  real(kind=dp), allocatable   :: recfile_weights(:,:,:)
  logical                      :: have_epi, have_equ, have_antipode
  real                         :: dtheta_rec
  
  ! CMB receivers (same as receivers, just above CMB instead)
  integer                      :: num_cmb
  integer, allocatable         :: cmbfile_el(:,:), loc2globcmb(:)
  !--------------------------------------------------------------------------

  ! for xdmf plotting
  integer                      :: nelem_plot, npoint_plot
  logical, allocatable         :: plotting_mask(:,:,:)
  integer, allocatable         :: mapping_ijel_iplot(:,:,:)

  ! for kernel wavefields in displ_only mode
  integer                      :: nelem_kwf_global, nelem_kwf
  integer                      :: npoint_kwf_global, npoint_kwf
  integer                      :: npoint_solid_kwf, npoint_fluid_kwf
  logical, allocatable         :: kwf_mask(:,:,:)
  integer, allocatable         :: mapping_ijel_ikwf(:,:,:)
  integer, allocatable         :: midpoint_mesh_kwf(:), eltype_kwf(:), axis_kwf(:)
  integer, allocatable         :: fem_mesh_kwf(:,:)
  integer, allocatable         :: sem_mesh_kwf(:,:,:)

  ! Only needed before the simulation and later deallocated
  ! Global mesh informations
  integer, dimension(:,:), allocatable          :: lnods ! (nelem,1:8)
  character(len=6), dimension(:), allocatable   :: eltype ! (nelem)
  integer                                       :: npoin
  real(kind=dp), dimension(:,:), allocatable    :: crd_nodes ! (npoin,2)
  logical, dimension(:), allocatable            :: coarsing,north ! (nelem)
  integer                                       :: num_spher_radii
  real(kind=dp)   , dimension(:), allocatable   :: spher_radii
  ! Axial elements
  logical, dimension(:), allocatable            :: axis ! (nelem)

  ! Mapping between solid/fluid elements:
  ! integer array of size nel_fluid containing the glocal (global per-proc)
  ! element number for 1:nel_solid/fluid
  integer, dimension(:), allocatable            :: ielsolid ! (nel_solid)
  integer, dimension(:), allocatable            :: ielfluid ! (nel_fluid)

  ! Mask for the free surface boundary in fluids
  real(kind=realkind), dimension(:,:,:), allocatable :: fluid_free_surface_mask

  ! damping factor for absorbing boundaries
  logical                                              :: have_absorbing_bc
  real(kind=realkind), dimension(:,:,:), allocatable   :: fluid_absorbing_gamma
  real(kind=realkind), dimension(:,:,:), allocatable   :: solid_absorbing_gamma

contains

!-----------------------------------------------------------------------------------------
!> Cache physical element centers and build phase-specific local search trees.
subroutine build_element_trees(midpoints)
   real(kind=dp), intent(in) :: midpoints(2,nelem)
   integer :: i

   if (associated(element_midpoint)) call destroy_element_trees
   allocate(element_midpoint(2,nelem),solid_element_index(nelem),fluid_element_index(nelem))
   element_midpoint = real(midpoints,kind=kdkind)
   solid_element_index = 0
   fluid_element_index = 0
   do i=1,nel_solid
      solid_element_index(ielsolid(i)) = i
   enddo
   do i=1,nel_fluid
      fluid_element_index(ielfluid(i)) = i
   enddo

   if (nelem > 1) element_tree => kdtree2_create(element_midpoint,sort=.true.)
   if (nel_solid > 0) then
      allocate(solid_midpoint(2,nel_solid))
      solid_midpoint = element_midpoint(:,ielsolid)
      if (nel_solid > 1) solid_tree => kdtree2_create(solid_midpoint,sort=.true.)
   endif
   if (nel_fluid > 0) then
      allocate(fluid_midpoint(2,nel_fluid))
      fluid_midpoint = element_midpoint(:,ielfluid)
      if (nel_fluid > 1) fluid_tree => kdtree2_create(fluid_midpoint,sort=.true.)
   endif
end subroutine build_element_trees

!-----------------------------------------------------------------------------------------
!> Return up to ten nearest local element IDs (1-based global mesh order).
!! phase: 0=all, 1=solid, 2=fluid.
subroutine nearest_mesh_elements(s,z,phase,elements,nfound)
   real(kind=dp), intent(in) :: s,z
   integer, intent(in) :: phase
   integer, intent(out) :: elements(nearest_element_count),nfound
   real(kind=kdkind), target :: query(2)
   type(kdtree2_result) :: results(nearest_element_count)
   integer :: i,navailable

   elements = 0
   select case(phase)
   case(0)
      navailable = nelem
   case(1)
      navailable = nel_solid
   case(2)
      navailable = nel_fluid
   case default
      error stop 'Invalid element search phase'
   end select
   nfound = min(nearest_element_count,navailable)
   if (nfound == 0) return
   if (nfound == 1) then
      select case(phase)
      case(0)
         elements(1) = 1
      case(1)
         elements(1) = ielsolid(1)
      case(2)
         elements(1) = ielfluid(1)
      end select
      return
   endif

   query = [real(s,kind=kdkind),real(z,kind=kdkind)]
   select case(phase)
   case(0)
      call kdtree2_n_nearest(element_tree,query,nfound,results)
      elements(1:nfound) = [(results(i)%idx,i=1,nfound)]
   case(1)
      call kdtree2_n_nearest(solid_tree,query,nfound,results)
      elements(1:nfound) = [(ielsolid(results(i)%idx),i=1,nfound)]
   case(2)
      call kdtree2_n_nearest(fluid_tree,query,nfound,results)
      elements(1:nfound) = [(ielfluid(results(i)%idx),i=1,nfound)]
   end select
end subroutine nearest_mesh_elements

!-----------------------------------------------------------------------------------------
!> Release the search data after station and face preparation.
subroutine destroy_element_trees
   if (associated(element_tree)) then
      call kdtree2_destroy(element_tree)
      nullify(element_tree)
   endif
   if (associated(solid_tree)) then
      call kdtree2_destroy(solid_tree)
      nullify(solid_tree)
   endif
   if (associated(fluid_tree)) then
      call kdtree2_destroy(fluid_tree)
      nullify(fluid_tree)
   endif
   if (associated(element_midpoint)) deallocate(element_midpoint)
   if (associated(solid_midpoint)) deallocate(solid_midpoint)
   if (associated(fluid_midpoint)) deallocate(fluid_midpoint)
   if (allocated(solid_element_index)) deallocate(solid_element_index)
   if (allocated(fluid_element_index)) deallocate(fluid_element_index)
end subroutine destroy_element_trees

!-----------------------------------------------------------------------------------------
!> Read parameters formerly in mesh_params.h 
!! It is slightly dirty to have this routine in a data module
!! but it allows to define the variables as 'protected', i.e.
!! fixed outside of this module.
subroutine read_mesh_basics(iounit)

   integer, intent(in)   :: iounit

   read(iounit) nproc_mesh
   read(iounit) npol
   read(iounit) nelem
   read(iounit) npoint
   read(iounit) nel_solid
   read(iounit) nel_fluid
   read(iounit) npoint_solid
   read(iounit) npoint_fluid
   read(iounit) nglob_solid
   read(iounit) nglob_fluid
   read(iounit) nel_bdry
   read(iounit) ndisc
   read(iounit) lfbkgrdmodel

end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine read_mesh_advanced(iounit)
   use data_io, only     : verbose 
   use data_spec
   integer, intent(in)  :: iounit
   integer              :: iptcp, iel, inode

   allocate(eta(0:npol))
   allocate(dxi(0:npol))
   allocate(wt(0:npol))
   allocate(xi_k(0:npol))
   allocate(wt_axial_k(0:npol))
   allocate(G1(0:npol,0:npol))
   allocate(G1T(0:npol,0:npol))
   allocate(G2(0:npol,0:npol))
   allocate(G2T(0:npol,0:npol))
   allocate(G0(0:npol))

   ! spectral stuff
   read(iounit) xi_k        
   read(iounit) eta 
   read(iounit) dxi       
   read(iounit) wt        
   read(iounit) wt_axial_k
   read(iounit) G0
   read(iounit) G1
   read(iounit) G1T
   read(iounit) G2
   read(iounit) G2T

   read(iounit) npoin
   
   if (verbose > 1) then
      write(69,*) 'reading coordinates/control points...'
      write(69,*) 'global number of control points:',npoin
   endif
   allocate(crd_nodes(1:npoin,1:2))
   
   read(iounit) crd_nodes(:,1)
   read(iounit) crd_nodes(:,2)
   do iptcp = 1, npoin 
      if(abs(crd_nodes(iptcp,2)) < 1.e-8) crd_nodes(iptcp,2) = zero
   end do

   allocate(lnods(1:nelem,1:8))
   do iel = 1, nelem
      read(iounit) (lnods(iel,inode), inode=1,8)
   end do


   ! Number of global distinct points (slightly differs for each processor!)
   read(iounit) nglob
   if (verbose > 1) write(69,*) '  global number:', nglob
 
   ! Element type
   allocate(eltype(1:nelem), coarsing(1:nelem))
   read(iounit) eltype
   read(iounit) coarsing
 
   !!!!!!!!!!! SOLID/FLUID !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
   ! mapping from sol/flu (1:nel_fluid) to global element numbers (1:neltot) 
   if (verbose > 1) write(69,*) 'reading solid/fluid domain info...'
   allocate(ielsolid(1:nel_solid))
   allocate(ielfluid(1:nel_fluid))
   read(iounit) ielsolid
   read(iounit) ielfluid
 
   ! slocal numbering 
   allocate(igloc_solid(npoint_solid))
   read(iounit) igloc_solid(1:npoint_solid)
 
   ! flocal numbering 
   allocate(igloc_fluid(npoint_fluid))
   read(iounit) igloc_fluid(1:npoint_fluid)
 
   ! Solid-Fluid boundary
   if (verbose > 1) write(69,*) 'reading solid/fluid boundary info...'
   read(iounit) have_bdry_elem
 
   if (have_bdry_elem) then
      allocate(bdry_solid_el(1:nel_bdry))
      allocate(bdry_fluid_el(1:nel_bdry))
      allocate(bdry_jpol_solid(1:nel_bdry))
      allocate(bdry_jpol_fluid(1:nel_bdry))
      read(iounit) bdry_solid_el(1:nel_bdry)
      read(iounit) bdry_fluid_el(1:nel_bdry)
      read(iounit) bdry_jpol_solid(1:nel_bdry)
      read(iounit) bdry_jpol_fluid(1:nel_bdry)
   endif
end subroutine
!-----------------------------------------------------------------------------------------

!-----------------------------------------------------------------------------------------
subroutine read_mesh_axel(iounit)
   integer, intent(in) :: iounit

   allocate(ax_el(naxel))
   allocate(ax_el_solid(1:naxel_solid))
   allocate(ax_el_fluid(1:naxel_fluid))
   allocate(axis_solid(nel_solid))
   allocate(axis_fluid(nel_fluid))

   read(iounit) ax_el(1:naxel)
   read(iounit) ax_el_solid(1:naxel_solid)
   read(iounit) ax_el_fluid(1:naxel_fluid)
end subroutine
!-----------------------------------------------------------------------------------------

end module data_mesh
!=========================================================================================
