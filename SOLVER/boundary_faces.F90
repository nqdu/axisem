! GLL wavefields on elements selected by ordered 5 x 5 boundary faces.
module boundary_faces
#ifdef enable_netcdf
  use netcdf
#endif
  use global_parameters
  use data_io, only: datapath, lfdata, save_bdry_faces
  use data_mesh, only: nelem, nel_solid, nel_fluid, npol, lnods, crd_nodes, axis, &
                       router, ielsolid, ielfluid, eltype, nearest_element_count, &
                       nearest_mesh_elements, solid_element_index, fluid_element_index
  use data_spec, only: xi_k, eta, G0, G1, G2
  use data_source, only: rot_src, src_type
  use data_proc, only: mynum, nproc
  use data_time, only: niter, deltat, strain_it
  use data_io, only: nstrain, strain_t0
  use commun, only: pmin, psum_int, pcheck, barrier, comm_elem_number
  use utlity, only: inside_element, compute_coordinates
  use nc_helpers, only: check
  implicit none
  private
  public :: prepare_boundary_faces, cache_boundary_materials, create_boundary_output, &
            sample_boundary_faces, finish_boundary_output

  integer :: nfaces = 0, npoints = 0, nlocal = 0, buffered = 0, first_sample = 1
  integer :: nsolid_local = 0, nfluid_local = 0, nsolid_global = 0, nfluid_global = 0
  integer :: solid_first = 1, fluid_first = 1, mesh_first = 1
  integer, allocatable :: face_phase(:), point_id(:), element_id(:), global_element_id(:)
  integer, allocatable :: solid_elements(:), fluid_elements(:), unique_elements(:)
  integer, allocatable :: point_slot(:), point_mesh_index(:)
  real(kind=dp), allocatable :: longitude(:), latitude(:), depth_km(:), &
                                phi_src(:), xi_loc(:), eta_loc(:)
  real(kind=dp), allocatable :: material(:,:,:,:)
  real(kind=realkind), allocatable :: solid_samples(:,:,:,:,:), fluid_samples(:,:,:,:)
  character(len=*), parameter :: output_name = '/boundary_wavefields.nc4'

contains

  subroutine prepare_boundary_faces
    integer :: io, i, iface, owner, phase, local_phase, rank, j
    integer :: row, line_count, local_elem, local_count
    integer :: global_nelem, global_first, global_last, first, last, total
    real(kind=dp) :: s, z, phi, local_xi, local_eta, estimated_s, estimated_z
    real(kind=dp), allocatable :: location_error_m(:)
    logical, allocatable :: selected(:)
    character(len=512) :: line, errmsg
    character(len=32) :: input_file
    logical :: file_exists

    if (.not. save_bdry_faces) return
    input_file = 'boundary_faces.dat'
    inquire(file=trim(input_file),exist=file_exists)
    if (.not. file_exists) then
       input_file = 'boundary_surfaces.dat'
       inquire(file=trim(input_file),exist=file_exists)
       if (.not. file_exists) input_file = 'boundar_surfaces.dat'
    endif
    open(newunit=io, file=trim(input_file), status='old', action='read', iostat=local_count)
    call pcheck(local_count /= 0, 'Cannot open boundary_faces.dat, boundary_surfaces.dat, or boundar_surfaces.dat')
    line_count = 0
    do
       read(io,'(a)',iostat=local_count) line
       if (local_count < 0) exit
       call pcheck(local_count /= 0, 'Cannot read '//trim(input_file))
       if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
       line_count = line_count + 1
    enddo
    call pcheck(line_count == 0 .or. mod(line_count,25) /= 0, &
                trim(input_file)//' must have 25 rows per face')
    npoints = line_count
    nfaces = npoints / 25
    allocate(longitude(npoints), latitude(npoints), depth_km(npoints), &
             phi_src(npoints), face_phase(nfaces), xi_loc(npoints), eta_loc(npoints))
    rewind(io)
    row = 0
    do
       read(io,'(a)',iostat=local_count) line
       if (local_count < 0) exit
       if (len_trim(line) == 0 .or. line(1:1) == '#') cycle
       row = row + 1
       read(line,*,iostat=local_count) longitude(row), latitude(row), depth_km(row)
       write(errmsg,'(a,a,a,i0)') 'Invalid ',trim(input_file),' row ',row
       call pcheck(local_count /= 0, errmsg)
       call pcheck(abs(latitude(row)) > 90.0_dp .or. depth_km(row) < 0.0_dp .or. &
                   depth_km(row)*1000.0_dp > router, errmsg)
    enddo
    close(io)

    ! Determine the phase from the center (row 13) of each face.
    do iface=1,nfaces
       i = (iface-1)*25 + 13
       call coordinates(i,s,z,phi)
       call locate_point(s,z,0,local_elem,local_xi,local_eta,estimated_s,estimated_z)
       owner = nproc
       if (local_elem > 0) owner = mynum
       owner = int(pmin(real(owner,dp)))
       if (owner == nproc) then
          call locate_point(s,z,0,local_elem,local_xi,local_eta,estimated_s,estimated_z,.true.)
          if (local_elem > 0) owner = mynum
          owner = int(pmin(real(owner,dp)))
       endif
       write(errmsg,'(a,i0)') 'Cannot locate center of boundary face ', iface
       call pcheck(owner == nproc, errmsg)
       local_phase = 0
       if (mynum == owner) then
          if (fluid_element_index(local_elem) > 0) local_phase = 1
       endif
       face_phase(iface) = psum_int(local_phase)
    enddo

    call comm_elem_number(nelem,global_nelem,global_first,global_last)
    allocate(selected(nelem))
    selected = .false.
    allocate(point_id(npoints), element_id(npoints), global_element_id(npoints), &
             point_slot(npoints), point_mesh_index(npoints), location_error_m(npoints))
    nlocal = 0
    do i=1,npoints
       iface = (i-1)/25 + 1
       call coordinates(i,s,z,phi)
       phi_src(i) = phi
       phase = face_phase(iface)
       call locate_point(s,z,phase+1,local_elem,local_xi,local_eta,estimated_s,estimated_z)
       owner = nproc
       if (local_elem > 0) owner = mynum
       owner = int(pmin(real(owner,dp)))
       if (owner == nproc) then
          call locate_point(s,z,phase+1,local_elem,local_xi,local_eta,estimated_s,estimated_z,.true.)
          if (local_elem > 0) owner = mynum
          owner = int(pmin(real(owner,dp)))
       endif
       write(errmsg,'(a,i0)') 'Cannot locate boundary point in its center phase: ', i
       call pcheck(owner == nproc, errmsg)
       if (owner /= mynum) cycle
       nlocal = nlocal + 1
       point_id(nlocal) = i
       selected(local_elem) = .true.
       global_element_id(nlocal) = global_first + local_elem - 2
       if (phase == 0) then
          element_id(nlocal) = solid_element_index(local_elem)
       else
          element_id(nlocal) = fluid_element_index(local_elem)
       endif
       xi_loc(nlocal) = local_xi
       eta_loc(nlocal) = local_eta
       location_error_m(nlocal) = hypot(s-estimated_s,z-estimated_z)
    enddo
    nsolid_local = count(selected(ielsolid))
    nfluid_local = count(selected(ielfluid))
    allocate(solid_elements(nsolid_local), fluid_elements(nfluid_local), &
             unique_elements(nsolid_local+nfluid_local))
    j = 0
    do i=1,nel_solid
       if (.not. selected(ielsolid(i))) cycle
       j = j+1
       solid_elements(j) = i
       unique_elements(j) = ielsolid(i)
    enddo
    j = 0
    do i=1,nel_fluid
       if (.not. selected(ielfluid(i))) cycle
       j = j+1
       fluid_elements(j) = i
       unique_elements(nsolid_local+j) = ielfluid(i)
    enddo
    call comm_elem_number(nsolid_local,nsolid_global,solid_first,last)
    call comm_elem_number(nfluid_local,nfluid_global,fluid_first,last)
    call comm_elem_number(nsolid_local+nfluid_local,total,mesh_first,last)
    do i=1,nlocal
       if (face_phase((point_id(i)-1)/25+1) == 0) then
          do j=1,nsolid_local
             if (solid_elements(j) /= element_id(i)) cycle
             point_slot(i) = solid_first+j-2
             point_mesh_index(i) = mesh_first+j-2
             exit
          enddo
       else
          do j=1,nfluid_local
             if (fluid_elements(j) /= element_id(i)) cycle
             point_slot(i) = fluid_first+j-2
             point_mesh_index(i) = mesh_first+nsolid_local+j-2
             exit
          enddo
       endif
    enddo
    deallocate(selected)
    if (mynum == 0) write(6,'(a,i0,a,i0)') &
         'Boundary faces: ',nfaces,' faces, points: ',npoints
    do rank=0,nproc-1
       call barrier
       if (mynum /= rank) cycle
       do i=1,nlocal
          write(6,'(a,i0,a,i0,a,i0,a,es14.5)') &
               '  boundary point=',point_id(i), ' proc=',mynum, &
               ' global element=',global_element_id(i), &
               ' distance to real location [m]=',location_error_m(i)
       enddo
       call flush(6)
    enddo
    call barrier
    deallocate(location_error_m)
    call pcheck(nstrain < 1, 'No boundary samples: adjust DUMP_T0 or simulation length')
    allocate(solid_samples(0:npol,0:npol,3,nsolid_local,min(nstrain,max(1,nc_buffer_size()))))
    allocate(fluid_samples(0:npol,0:npol,nfluid_local,min(nstrain,max(1,nc_buffer_size()))))
    solid_samples = 0.0_realkind
    fluid_samples = 0.0_realkind
    allocate(material(0:npol,0:npol,8,size(unique_elements)))
    material = 0.0_dp
  end subroutine prepare_boundary_faces

  subroutine cache_boundary_materials(rho,lambda,mu,xi_ani,phi_ani,eta_ani)
    real(kind=dp), intent(in) :: rho(0:,0:,:),lambda(0:,0:,:),mu(0:,0:,:)
    real(kind=dp), intent(in) :: xi_ani(0:,0:,:),phi_ani(0:,0:,:),eta_ani(0:,0:,:)
    integer :: i,elem
    if (.not. save_bdry_faces) return
    do i=1,size(unique_elements)
       elem = unique_elements(i)
       material(:,:,1,i) = rho(:,:,elem)
       material(:,:,2,i) = lambda(:,:,elem)
       material(:,:,3,i) = mu(:,:,elem)
       material(:,:,4,i) = sqrt((lambda(:,:,elem)+2.0_dp*mu(:,:,elem))/rho(:,:,elem))
       material(:,:,5,i) = sqrt(mu(:,:,elem)/rho(:,:,elem))
       material(:,:,6,i) = xi_ani(:,:,elem)
       material(:,:,7,i) = phi_ani(:,:,elem)
       material(:,:,8,i) = eta_ani(:,:,elem)
    enddo
  end subroutine cache_boundary_materials

  subroutine coordinates(i,s,z,phi)
    use data_io, only: trans_rot_mat
    integer, intent(in) :: i
    real(kind=dp), intent(out) :: s,z,phi
    real(kind=dp) :: theta, lambda, r, xyz(3), rotated(3)
    theta = (90.0_dp-latitude(i))*pi/180.0_dp
    lambda = longitude(i)*pi/180.0_dp
    r = router - depth_km(i)*1000.0_dp
    xyz = [r*sin(theta)*cos(lambda), r*sin(theta)*sin(lambda), r*cos(theta)]
    rotated = xyz
    if (rot_src) rotated = matmul(trans_rot_mat,xyz)
    s = hypot(rotated(1),rotated(2))
    z = rotated(3)
    phi = atan2(rotated(2),rotated(1))
  end subroutine coordinates

  subroutine locate_point(s,z,wanted_phase,elem,xil,etal,estimated_s,estimated_z,full_search)
    integer, intent(in) :: wanted_phase ! 0: either, 1: solid, 2: fluid
    real(kind=dp), intent(in) :: s,z
    integer, intent(out) :: elem
    real(kind=dp), intent(out) :: xil,etal,estimated_s,estimated_z
    logical, intent(in), optional :: full_search
    integer :: ie, i, nfound, candidates(nearest_element_count)
    real(kind=dp) :: x,y,sloc,zloc,err,best,margin
    logical :: inside, scan_all
    elem = 0
    best = huge(1.0_dp)
    margin = max(1.0e-3_dp,1.0e-8_dp*router)
    scan_all = .false.
    if (present(full_search)) scan_all = full_search
    if (scan_all) then
       nfound = nelem
    else
       call nearest_mesh_elements(s,z,wanted_phase,candidates,nfound)
    endif
    do i=1,nfound
       ie = i
       if (.not. scan_all) ie = candidates(i)
       if (wanted_phase == 1 .and. solid_element_index(ie) == 0) cycle
       if (wanted_phase == 2 .and. fluid_element_index(ie) == 0) cycle
       if (.not. scan_all) then
          if (s < minval(crd_nodes(lnods(ie,:),1))-margin .or. &
              s > maxval(crd_nodes(lnods(ie,:),1))+margin .or. &
              z < minval(crd_nodes(lnods(ie,:),2))-margin .or. &
              z > maxval(crd_nodes(lnods(ie,:),2))+margin) cycle
       endif
       call inside_element(s,z,ie,x,y,sloc,zloc,inside)
       if (.not. inside) cycle
       err = hypot(s-sloc,z-zloc)
       if (err >= best) cycle
       best = err
       elem = ie
       xil = max(-1.0_dp,min(1.0_dp,x))
       etal = max(-1.0_dp,min(1.0_dp,y))
       estimated_s = sloc
       estimated_z = zloc
    enddo
  end subroutine locate_point

  integer function nc_buffer_size()
    use nc_routines, only: nc_dumpbuffersize
    nc_buffer_size = nc_dumpbuffersize
  end function nc_buffer_size

  subroutine create_boundary_output
#ifdef enable_netcdf
    integer :: ncid, meshid, dimpoint, dimtime, dimface, dimelem, dimpol, dimcontrol
    integer :: dimsolid, dimfluid, i, it, ivar, rank, pointids(13), meshvar
    integer :: face_index(npoints), point_in_face(npoints)
    real(kind=dp) :: times(nstrain)
    character(len=12), parameter :: material_names(8) = &
         ['mesh_rho    ','mesh_lambda ','mesh_mu     ','mesh_vp     ', &
          'mesh_vs     ','mesh_xi     ','mesh_phi    ','mesh_eta    ']
    if (.not. save_bdry_faces) return
    if (mynum == 0) then
       call check(nf90_create(datapath(1:lfdata)//output_name, &
                              ior(NF90_CLOBBER,NF90_NETCDF4),ncid))
       call check(nf90_def_dim(ncid,'point',npoints,dimpoint))
       call check(nf90_def_dim(ncid,'time',nstrain,dimtime))
       call check(nf90_def_dim(ncid,'face',nfaces,dimface))
       call check(nf90_def_dim(ncid,'npol',npol+1,dimpol))
       if (nsolid_global > 0) call check(nf90_def_dim(ncid,'solid_element',nsolid_global,dimsolid))
       if (nfluid_global > 0) call check(nf90_def_dim(ncid,'fluid_element',nfluid_global,dimfluid))
       call check(nf90_def_var(ncid,'longitude',NF90_DOUBLE,[dimpoint],pointids(1)))
       call check(nf90_def_var(ncid,'latitude',NF90_DOUBLE,[dimpoint],pointids(2)))
       call check(nf90_def_var(ncid,'depth_km',NF90_DOUBLE,[dimpoint],pointids(3)))
       call check(nf90_def_var(ncid,'source_phi',NF90_DOUBLE,[dimpoint],pointids(4)))
       call check(nf90_def_var(ncid,'face_phase',NF90_INT,[dimface],pointids(5)))
       call check(nf90_def_var(ncid,'sample_time',NF90_DOUBLE,[dimtime],pointids(6)))
       call check(nf90_def_var(ncid,'face_index',NF90_INT,[dimpoint],pointids(7)))
       call check(nf90_put_att(ncid,pointids(7),'description','0-based boundary face index'))
       call check(nf90_def_var(ncid,'point_in_face',NF90_INT,[dimpoint],pointids(8)))
       call check(nf90_def_var(ncid,'element_id',NF90_INT,[dimpoint],pointids(9)))
       call check(nf90_def_var(ncid,'xi',NF90_DOUBLE,[dimpoint],pointids(10)))
       call check(nf90_def_var(ncid,'eta',NF90_DOUBLE,[dimpoint],pointids(11)))
       call check(nf90_def_var(ncid,'point_element_slot',NF90_INT,[dimpoint],pointids(12)))
       call check(nf90_def_var(ncid,'point_mesh_index',NF90_INT,[dimpoint],pointids(13)))
       call check(nf90_put_att(ncid,pointids(9),'description', &
                               '0-based global solver element index'))
       call check(nf90_put_att(ncid,pointids(12),'description', &
                               '0-based index into solid_element or fluid_element, according to face_phase'))
       call check(nf90_put_att(ncid,pointids(13),'description', &
                               '0-based index into Mesh/elements'))
       call check(nf90_def_grp(ncid,'Mesh',meshid))
       call check(nf90_def_dim(meshid,'elements',nsolid_global+nfluid_global,dimelem))
       call check(nf90_def_dim(meshid,'control_points',4,dimcontrol))
       call check(nf90_def_var(meshid,'axis',NF90_INT,[dimelem],meshvar))
       call check(nf90_def_var(meshid,'eltype',NF90_INT,[dimelem],meshvar))
       call check(nf90_def_var(meshid,'midpoint_mesh',NF90_INT,[dimelem],meshvar))
       call check(nf90_def_var(meshid,'fem_mesh',NF90_INT,[dimcontrol,dimelem],meshvar))
       call check(nf90_def_var(meshid,'sem_mesh',NF90_INT,[dimpol,dimpol,dimelem],meshvar))
       call check(nf90_def_var(meshid,'mesh_S',NF90_DOUBLE,[dimpol,dimpol,dimelem],meshvar))
       call check(nf90_def_var(meshid,'mesh_Z',NF90_DOUBLE,[dimpol,dimpol,dimelem],meshvar))
       do i=1,8
          call check(nf90_def_var(meshid,trim(material_names(i)),NF90_FLOAT, &
                                  [dimpol,dimpol,dimelem],meshvar))
       enddo
       call check(nf90_def_var(meshid,'gll',NF90_DOUBLE,[dimpol],meshvar))
       call check(nf90_def_var(meshid,'glj',NF90_DOUBLE,[dimpol],meshvar))
       call check(nf90_def_var(meshid,'G0',NF90_DOUBLE,[dimpol],meshvar))
       call check(nf90_def_var(meshid,'G1',NF90_DOUBLE,[dimpol,dimpol],meshvar))
       call check(nf90_def_var(meshid,'G2',NF90_DOUBLE,[dimpol,dimpol],meshvar))
       if (nsolid_global > 0) then
          call check(nf90_def_var(ncid,'solid_element_to_mesh',NF90_INT,[dimsolid],ivar))
          call check(nf90_def_var(ncid,'disp_s',NF90_FLOAT,[dimpol,dimpol,dimsolid,dimtime],ivar))
          call check(nf90_def_var(ncid,'disp_z',NF90_FLOAT,[dimpol,dimpol,dimsolid,dimtime],ivar))
          if (src_type(1) /= 'monopole') &
             call check(nf90_def_var(ncid,'disp_p',NF90_FLOAT,[dimpol,dimpol,dimsolid,dimtime],ivar))
       endif
       if (nfluid_global > 0) then
          call check(nf90_def_var(ncid,'fluid_element_to_mesh',NF90_INT,[dimfluid],ivar))
          call check(nf90_def_var(ncid,'chi',NF90_FLOAT,[dimpol,dimpol,dimfluid,dimtime],ivar))
       endif
       call check(nf90_put_att(ncid,NF90_GLOBAL,'points_per_face',25))
       call check(nf90_put_att(ncid,NF90_GLOBAL,'face_phase_codes','0=solid, 1=fluid; from point 13'))
       call check(nf90_put_att(ncid,NF90_GLOBAL,'components', &
                               'source-centered cylindrical modal components at element GLL nodes'))
       call check(nf90_enddef(ncid))
       call check(nf90_put_var(ncid,pointids(1),longitude))
       call check(nf90_put_var(ncid,pointids(2),latitude))
       call check(nf90_put_var(ncid,pointids(3),depth_km))
       call check(nf90_put_var(ncid,pointids(4),phi_src))
       call check(nf90_put_var(ncid,pointids(5),face_phase))
       face_index = [((i-1)/25,i=1,npoints)]
       point_in_face = [(mod(i-1,25)+1,i=1,npoints)]
       call check(nf90_put_var(ncid,pointids(7),face_index))
       call check(nf90_put_var(ncid,pointids(8),point_in_face))
       it = 0
       do i=1,niter
          if (mod(i,strain_it) /= 0 .or. i*deltat <= strain_t0) cycle
          it = it+1
          if (it <= nstrain) times(it) = i*deltat
       enddo
       call check(nf90_put_var(ncid,pointids(6),times))
       call check(nf90_inq_varid(meshid,'gll',meshvar))
       call check(nf90_put_var(meshid,meshvar,eta))
       call check(nf90_inq_varid(meshid,'glj',meshvar))
       call check(nf90_put_var(meshid,meshvar,xi_k))
       call check(nf90_inq_varid(meshid,'G0',meshvar))
       call check(nf90_put_var(meshid,meshvar,real(G0,dp)))
       call check(nf90_inq_varid(meshid,'G1',meshvar))
       call check(nf90_put_var(meshid,meshvar,real(G1,dp)))
       call check(nf90_inq_varid(meshid,'G2',meshvar))
       call check(nf90_put_var(meshid,meshvar,real(G2,dp)))
       call check(nf90_close(ncid))
    endif
    do rank=0,nproc-1
       call barrier
       if (mynum /= rank .or. nlocal == 0) cycle
       call check(nf90_open(datapath(1:lfdata)//output_name,NF90_WRITE,ncid))
       call check(nf90_inq_varid(ncid,'element_id',pointids(9)))
       call check(nf90_inq_varid(ncid,'xi',pointids(10)))
       call check(nf90_inq_varid(ncid,'eta',pointids(11)))
       call check(nf90_inq_varid(ncid,'point_element_slot',pointids(12)))
       call check(nf90_inq_varid(ncid,'point_mesh_index',pointids(13)))
       do i=1,nlocal
          call check(nf90_put_var(ncid,pointids(9),global_element_id(i),start=[point_id(i)]))
          call check(nf90_put_var(ncid,pointids(10),xi_loc(i),start=[point_id(i)]))
          call check(nf90_put_var(ncid,pointids(11),eta_loc(i),start=[point_id(i)]))
          call check(nf90_put_var(ncid,pointids(12),point_slot(i),start=[point_id(i)]))
          call check(nf90_put_var(ncid,pointids(13),point_mesh_index(i),start=[point_id(i)]))
       enddo
       if (nsolid_local > 0) then
          call check(nf90_inq_varid(ncid,'solid_element_to_mesh',ivar))
          do i=1,nsolid_local
             call check(nf90_put_var(ncid,ivar,mesh_first+i-2,start=[solid_first+i-1]))
          enddo
       endif
       if (nfluid_local > 0) then
          call check(nf90_inq_varid(ncid,'fluid_element_to_mesh',ivar))
          do i=1,nfluid_local
             call check(nf90_put_var(ncid,ivar,mesh_first+nsolid_local+i-2,start=[fluid_first+i-1]))
          enddo
       endif
       call check(nf90_inq_ncid(ncid,'Mesh',meshid))
       call write_boundary_mesh(meshid)
       call check(nf90_close(ncid))
    enddo
    call barrier
    deallocate(material)
#endif
  end subroutine create_boundary_output

  subroutine write_boundary_mesh(meshid)
#ifdef enable_netcdf
    integer, intent(in) :: meshid
    integer :: ids(7), material_ids(8), q, elem, pos, i, j, code, n
    integer :: sem(0:npol,0:npol), fem(4)
    real(kind=dp) :: scoord(0:npol,0:npol), zcoord(0:npol,0:npol), r, theta
    character(len=12), parameter :: material_names(8) = &
         ['mesh_rho    ','mesh_lambda ','mesh_mu     ','mesh_vp     ', &
          'mesh_vs     ','mesh_xi     ','mesh_phi    ','mesh_eta    ']
    n = npol+1
    call check(nf90_inq_varid(meshid,'axis',ids(1)))
    call check(nf90_inq_varid(meshid,'eltype',ids(2)))
    call check(nf90_inq_varid(meshid,'midpoint_mesh',ids(3)))
    call check(nf90_inq_varid(meshid,'fem_mesh',ids(4)))
    call check(nf90_inq_varid(meshid,'sem_mesh',ids(5)))
    call check(nf90_inq_varid(meshid,'mesh_S',ids(6)))
    call check(nf90_inq_varid(meshid,'mesh_Z',ids(7)))
    do i=1,8
       call check(nf90_inq_varid(meshid,trim(material_names(i)),material_ids(i)))
    enddo
    do q=1,size(unique_elements)
       elem = unique_elements(q)
       pos = mesh_first+q-1
       do j=0,npol
          do i=0,npol
             call compute_coordinates(scoord(i,j),zcoord(i,j),r,theta,elem,i,j)
             sem(i,j) = (pos-1)*n*n+j*n+i
          enddo
       enddo
       fem = [sem(0,0),sem(npol,0),sem(npol,npol),sem(0,npol)]
       code = -1
       select case (trim(eltype(elem)))
       case ('curved'); code = 0
       case ('linear'); code = 1
       case ('semino'); code = 2
       case ('semiso'); code = 3
       end select
       call check(nf90_put_var(meshid,ids(1),merge(1,0,axis(elem)),start=[pos]))
       call check(nf90_put_var(meshid,ids(2),code,start=[pos]))
       call check(nf90_put_var(meshid,ids(3),sem(npol/2,npol/2),start=[pos]))
       call check(nf90_put_var(meshid,ids(4),reshape(fem,[4,1]),start=[1,pos],count=[4,1]))
       call check(nf90_put_var(meshid,ids(5),reshape(sem,[n,n,1]), &
                               start=[1,1,pos],count=[n,n,1]))
       call check(nf90_put_var(meshid,ids(6),reshape(scoord,[n,n,1]), &
                               start=[1,1,pos],count=[n,n,1]))
       call check(nf90_put_var(meshid,ids(7),reshape(zcoord,[n,n,1]), &
                               start=[1,1,pos],count=[n,n,1]))
       do i=1,8
          call check(nf90_put_var(meshid,material_ids(i), &
                                  reshape(real(material(:,:,i,q),realkind),[n,n,1]), &
                                  start=[1,1,pos],count=[n,n,1]))
       enddo
    enddo
#endif
  end subroutine write_boundary_mesh

  subroutine sample_boundary_faces(u,chi,isample)
    real(kind=realkind), intent(in) :: u(0:,0:,:,:),chi(0:,0:,:)
    integer, intent(in) :: isample
    integer :: i,elem
    if (.not. save_bdry_faces) return
    buffered = buffered+1
    do i=1,nsolid_local
       elem = solid_elements(i)
       if (src_type(1) == 'dipole') then
          solid_samples(:,:,1,i,buffered) = u(:,:,elem,1)+u(:,:,elem,2)
          solid_samples(:,:,3,i,buffered) = u(:,:,elem,1)-u(:,:,elem,2)
       else
          solid_samples(:,:,1,i,buffered) = u(:,:,elem,1)
          if (src_type(1) /= 'monopole') solid_samples(:,:,3,i,buffered) = u(:,:,elem,2)
       endif
       solid_samples(:,:,2,i,buffered) = u(:,:,elem,3)
    enddo
    do i=1,nfluid_local
       fluid_samples(:,:,i,buffered) = chi(:,:,fluid_elements(i))
    enddo
    if (buffered == size(solid_samples,5) .or. isample == nstrain) call flush_boundary_output
  end subroutine sample_boundary_faces

  subroutine flush_boundary_output
#ifdef enable_netcdf
    integer :: rank, ncid, ids(4), i,j,n
    n = npol+1
    do rank=0,nproc-1
       call barrier
       if (mynum /= rank .or. nsolid_local+nfluid_local == 0) cycle
       call check(nf90_open(datapath(1:lfdata)//output_name,NF90_WRITE,ncid))
       if (nsolid_local > 0) then
          call check(nf90_inq_varid(ncid,'disp_s',ids(1)))
          call check(nf90_inq_varid(ncid,'disp_z',ids(2)))
          if (src_type(1) /= 'monopole') call check(nf90_inq_varid(ncid,'disp_p',ids(3)))
          do i=1,nsolid_local
             do j=1,3
                if (j == 3 .and. src_type(1) == 'monopole') cycle
                call check(nf90_put_var(ncid,ids(j), &
                     reshape(solid_samples(:,:,j,i,1:buffered),[n,n,1,buffered]), &
                     start=[1,1,solid_first+i-1,first_sample],count=[n,n,1,buffered]))
             enddo
          enddo
       endif
       if (nfluid_local > 0) then
          call check(nf90_inq_varid(ncid,'chi',ids(4)))
          do i=1,nfluid_local
             call check(nf90_put_var(ncid,ids(4), &
                  reshape(fluid_samples(:,:,i,1:buffered),[n,n,1,buffered]), &
                  start=[1,1,fluid_first+i-1,first_sample],count=[n,n,1,buffered]))
          enddo
       endif
       call check(nf90_close(ncid))
    enddo
    call barrier
#endif
    first_sample = first_sample+buffered
    buffered = 0
  end subroutine flush_boundary_output

  subroutine finish_boundary_output
    if (.not. save_bdry_faces) return
    if (buffered > 0) call flush_boundary_output
  end subroutine finish_boundary_output
end module boundary_faces
