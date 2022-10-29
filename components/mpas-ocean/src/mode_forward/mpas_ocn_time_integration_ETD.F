! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integration_ETD
!
!> \brief MPAS ocean ETD Time integration scheme for tracers
!> \author Sara Calandrini
!> \date   October 2020
!> \details
!>  This module contains the ETD time integration routine for tracers. 
!
!-----------------------------------------------------------------------

module ocn_time_integration_ETD

   use mpas_derived_types
   use mpas_pool_routines
   use mpas_constants
   use mpas_dmpar
   use mpas_threading
   use mpas_vector_reconstruction
   use mpas_spline_interpolation
   use mpas_timer

   use ocn_constants
   use ocn_tendency
   use ocn_diagnostics_variables
   use ocn_mesh
   use ocn_gm

   use ocn_equation_of_state
   use ocn_vmix_cvmix
   use ocn_vmix
   use ocn_time_average_coupled
   use ocn_wetting_drying

   use ocn_tracer_nonlocalflux
   use ocn_surface_land_ice_fluxes
   use ocn_effective_density_in_land_ice

   use mpas_tracer_advection_helpers

   use ocn_tracer_advection_std
   use ocn_tracer_hmix_del2
   use ocn_tracer_hmix_Redi
   use ocn_tracer_advection_vert
   use ocn_tracer_surface_flux_to_tend
   use ocn_tracer_ecosys
   use ocn_tracer_DMS
   use ocn_tracer_MacroMolecules

   implicit none
   private
   save

   ! private module variables
   real (kind=RKIND) :: &
      coef3rdOrder       !< coefficient for blending high-order terms

   integer :: vertOrder  !< choice of order for vertical advection
   integer, parameter :: &! enumerator for supported vertical adv order
      vertOrder2=2,      &!< 2nd order
      vertOrder3=3,      &!< 3rd order
      vertOrder4=4        !< 4th order

   !logical :: del2On
   !logical :: tracerVmixOn
   real (kind=RKIND) :: eddyDiff2

   !--------------------------------------------------------------------
   !
   ! Public parameters
   !
   !--------------------------------------------------------------------

   !--------------------------------------------------------------------
   !
   ! Public member functions
   !
   !--------------------------------------------------------------------

   public :: ocn_time_integrator_ETD

   contains

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_time_integrator_ETD
!
!> \brief MPAS ocean ETD Time integration scheme for tracers
!> \author Sara Calandrini
!> \date   October 2020
!> \details
!>  This routine integrates one timestep (dt) using an ETD time integrator.
!
!-----------------------------------------------------------------------

   subroutine ocn_time_integrator_ETD(dt, vertAleTransportAvg, layerThicknessOld, &
      layerThicknessNew, tracersGroupCur, tracersGroupNew, tracersGroupTend)
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   ! Advance model state forward in time by the specified time step using
   ! a 2nd order exponential integrator 
   !
   ! Input: domain - current model state in time level 1 (e.g.,
   ! time_levs(1)state%h(:,:)) plus mesh meta-data
   ! Output: domain - upon exit, time level 2 (e.g., time_levs(2)%state%h(:,:))
   ! contains  model state advanced forward in time by dt seconds
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      implicit none

      !-----------------------------------------------------------------
      ! Input variables
      !-----------------------------------------------------------------

      real (kind=RKIND), intent(in) :: &
         dt              !< [in] time step (sec) to move forward

      real(KIND=RKIND), dimension(:,:), intent(in) :: &
         layerThicknessOld, &
         layerThicknessNew

      real(KIND=RKIND), dimension(:,:), intent(in) :: &
         vertAleTransportAvg

      !-----------------------------------------------------------------
      ! Input/output variables
      !-----------------------------------------------------------------

      real(KIND=RKIND), dimension(:,:,:), intent(inout) :: &
         tracersGroupTend, &
         tracersGroupCur, &
         tracersGroupNew

      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:), allocatable :: &
         JacZ            ! local Jacobian stored as a vector

      real (kind=RKIND), dimension(:,:,:), allocatable :: &
         phi1JTot        ! phi_1 matrix for all cells 

      integer :: k, k1, k2, err1, err2, err, iCell, cell1, cell2, iEdge, iTracer, CFL_pow, NLayers, &
                 num_tracers, num_tracersActive, num_tracersEco, num_tracersDMS, num_tracersMacro
      real (kind=RKIND), dimension(:,:,:), allocatable :: TracerGroupTendHor, TracerGroupTendVert
!      real (kind=RKIND), dimension(:,:,:), allocatable :: TracerGroupTendHorActive, TracerGroupTendVertActive
!      real (kind=RKIND), dimension(:,:,:), allocatable :: TracerGroupTendHorEcosysTracers, &
!                                                          TracerGroupTendHorDMSTracers, &
!                                                          TracerGroupTendHorMacroMoleculesTracers
      real (kind=RKIND), dimension(:), allocatable :: rhs, tracer_cur
      real (kind=RKIND), dimension(:,:), allocatable :: A
      real (kind=RKIND), dimension(:), allocatable :: vec

      logical :: computeBudgets ! flag to compute active tracer budget

      CFL_pow = 2 !5  
      vertOrder = 2
      coef3rdOrder = 0.25_RKIND !0.25
      !eddyDiff2 = 10.0 
      eddyDiff2 = 0.0001_RKIND

      call mpas_timer_start('time-stepping time for each iteration')

      !--------------------
      ! FIRST STAGE
      !--------------------

!      block => domain % blocklist

!      call mpas_pool_get_subpool(block % structs, 'state', statePool)
!      call mpas_pool_get_subpool(block % structs, 'mesh', meshPool)
!      call mpas_pool_get_subpool(block % structs, 'scratch', scratchPool)
!      call mpas_pool_get_subpool(block % structs, 'tend', tendPool)
!      call mpas_pool_get_subpool(block % structs, 'diagnostics', diagnosticsPool)
!      call mpas_pool_get_subpool(block % structs, 'forcing', forcingPool)
!      call mpas_pool_get_subpool(statePool, 'tracers', tracersPool)
!      call mpas_pool_get_subpool(tendPool, 'tracersTend', tracersTendPool)
!      call mpas_pool_get_subpool(forcingPool, 'tracersSurfaceFlux', tracersSurfaceFluxPool)
!      call mpas_pool_get_subpool(block % structs, 'verticalMesh', verticalMeshPool)

!      call mpas_pool_get_array(diagnosticsPool, 'normalTransportVelocity', normalTransportVelocity)
!      call mpas_pool_get_array(diagnosticsPool, 'layerThicknessEdge', layerThicknessEdge)
!      call mpas_pool_get_array(diagnosticsPool, 'vertAleTransportTop', vertAleTransportTop)
!      call mpas_pool_get_array(diagnosticsPool, 'vertNonLocalFlux', vertNonLocalFlux)
!      call mpas_pool_get_array(diagnosticsPool, 'vertDiffTopOfCell', vertDiffTopOfCell)
!      call mpas_pool_get_array(diagnosticsPool, 'vertViscTopOfEdge', vertViscTopOfEdge)
!      call mpas_pool_get_array(diagnosticsPool, 'vertViscTopOfCell', vertViscTopOfCell)
!      call mpas_pool_get_array(diagnosticsPool, 'velocityCell', velocityCell)
!      call mpas_pool_get_array(diagnosticsPool, 'zMid', zMid)
!      call mpas_pool_get_array(diagnosticsPool, 'slopeTriadUp', slopeTriadUp)
!      call mpas_pool_get_array(diagnosticsPool, 'slopeTriadDown', slopeTriadDown)

!      call mpas_pool_get_array(forcingPool, 'fractionAbsorbed',fractionAbsorbed)
!      call mpas_pool_get_array(forcingPool, 'fractionAbsorbedRunoff', fractionAbsorbedRunoff)

!      call mpas_pool_get_array(statePool, 'ssh', sshCur, 1) !TODO non sono sicura ci vada l'1, guardare nella chiamata delle tendencies
!      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 1)

!      call mpas_pool_get_array(tracersSurfaceFluxPool, 'debugTracersSurfaceFlux', tracerGroupSurfaceFlux)
!      call mpas_pool_get_array(tracersSurfaceFluxPool, 'ecosysTracersSurfaceFlux', tracerGroupEcosysTracersSurfaceFlux)
!      call mpas_pool_get_array(tracersSurfaceFluxPool, 'debugTracersSurfaceFlux', tracerGroupSurfaceFluxRunoff)
!      call mpas_pool_get_array(tracersSurfaceFluxPool, 'ecosysTracersSurfaceFlux', tracerGroupEcosysTracersSurfaceFluxRunoff)

!      call mpas_pool_get_array(tracersPool, 'debugTracers', tracerGroup, 1)
!      call mpas_pool_get_array(tracersPool, 'debugTracers', tracerGroupNew, 2)
!      call mpas_pool_get_array(tracersTendPool, 'debugTracersTend', TracerGroupTend)
!      call mpas_pool_get_array(tracersPool, 'activeTracers', tracerGroupActive, 1)
!      call mpas_pool_get_array(tracersPool, 'activeTracers', tracerGroupNewActive, 2)
!      call mpas_pool_get_array(tracersTendPool, 'activeTracersTend', TracerGroupTendActive)
!      call mpas_pool_get_array(tracersPool, 'ecosysTracers', tracerGroupEcosysTracers, 1)
!      call mpas_pool_get_array(tracersPool, 'ecosysTracers', tracerGroupNewEcosysTracers, 2)
!      call mpas_pool_get_array(tracersTendPool, 'ecosysTracersTend', TracerGroupTendEcosysTracers)

      num_tracers = size(tracersGroupCur,dim=1)
      !num_tracersActive  = size(tracerGroupActive,dim=1)
      !num_tracersEco = size(tracerGroupEcosysTracers,dim=1)
      allocate(TracerGroupTendHor(num_tracers, nVertLevels, nCellsAll))
      allocate(TracerGroupTendVert(num_tracers, nVertLevels, nCellsAll))
      !allocate(TracerGroupTendHorActive(num_tracersActive, nVertLevels, nCellsAll))
      !allocate(TracerGroupTendVertActive(num_tracersActive, nVertLevels, nCellsAll))
      !allocate(TracerGroupTendHorEcosysTracers(num_tracersEco, nVertLevels, nCellsAll))
      allocate(phi1JTot(nVertLevels, nVertLevels, nCellsAll))
      allocate(JacZ(nVertLevels*nVertLevels))
      allocate(rhs(nVertLevels))
      allocate(tracer_cur(nVertLevels))
      allocate(A(nVertLevels, nVertLevels))
      allocate(vec(nVertLevels))

      !call mpas_pool_initialize_time_levels(tracersPool)

!      do iCell = 1, nCellsAll
!         TracerGroupTendActive(:,:,iCell) = 0.0_RKIND
!         TracerGroupTendHorActive(:,:,iCell) = 0.0_RKIND
!         TracerGroupTendVertActive(:,:,iCell) = 0.0_RKIND
!      end do

      !Halo update for pv_edge
      !call mpas_dmpar_field_halo_exch(domain, 'normalizedRelativeVorticityEdge')

      ! determine whether active tracer budgets should be computed
      !computeBudgets = (budgetDiagsOn .and. tracerGroupName == 'activeTracers')
      computeBudgets = .false.

      call mpas_timer_start("1st stage-tendency computations")

      !print*, 'ETD ', num_tracers, num_tracersEco
      if (config_use_debugTracers) then
         !Copying hor tendencies
         do iCell = 1, nCellsAll
            TracerGroupTendHor(:,:,iCell)=tracersGroupTend(:,:,iCell)
         end do

!         !Computation of the vertical advection TODO metti la subroutine nuova
!         call ocn_tracer_vert_advection_std(tracersGroupCur, vertAleTransportAvg, &
!                                           layerThicknessOld, dt, tracersGroupTend)
      end if

      call mpas_timer_stop("1st stage-tendency computations")

      call mpas_timer_start('building jacobian and phi function for all the elements')
      do iCell = 1, nCellsOwned
         call mpas_timer_start('jacobian for each element')
         call jacobianVert(JacZ, layerAvg, vertAleTransportAvg, iCell)
         call mpas_timer_stop('jacobian for each element')
         Nlayers = maxLevelCell(iCell)
         call mpas_timer_start('phi function for each element')
         !call ocn_phi_function(JacZ, CFL_pow, phi1JTot(:,:,iCell), NLayers, dt)
         call ocn_phi_function_taylor(JacZ, phi1JTot(:,:,iCell), NLayers, dt)
         call mpas_timer_stop('phi function for each element')
      end do
      call mpas_timer_stop('building jacobian and phi function for all the elements')

      call mpas_timer_start("ETD step 1st-stage")
      if (config_use_debugTracers) then
         do iTracer = 1, num_tracers
            do iCell = 1, nCellsOwned
               Nlayers = maxLevelCell(iCell)
               tracer_cur(:) = 0.0_RKIND
               rhs(:) = 0.0_RKIND
               do k=1,Nlayers
                  tracer_cur(k) = tracersGroupCur(iTracer, k, iCell) &
                                * layerThicknessOld(k,iCell)
                  rhs(k) = tracersGroupTend(iTracer, k, iCell)
               end do

               CALL DGEMM('N','N',nVertLevels,1,nVertLevels,dt,phi1JTot(:,:,iCell), &
                          nVertLevels,rhs,nVertLevels,1.0_RKIND,tracer_cur,nVertLevels)
               do k = 1, Nlayers
                  if (layerThicknessNew(k,iCell)==0.0) then
                     tracersGroupNew(iTracer, k, iCell) = -1e+34
                  else
                     tracersGroupNew(iTracer, k, iCell) = tracer_cur(k) &
                                                    / layerThicknessNew(k,iCell)
                  end if
               end do
               !tracerGroup(iTracer, :, iCell) = tracer_cur(:) /
               !layerThickness(:,iCell)  
            end do
         end do
      end if
      call mpas_timer_stop("ETD step 1st-stage")

!      call mpas_timer_start("1st stage-prognostic halo update")
!      !Halo update tracer values
!      if (config_use_debugTracers) then
!         call mpas_dmpar_field_halo_exch(domain, 'debugTracers', timeLevel=2)
!      end if
!      if (config_use_ecosysTracers) then
!         call mpas_dmpar_field_halo_exch(domain, 'ecosysTracers', timeLevel=1)
!      end if
!      if (config_use_DMSTracers) then
!         call mpas_dmpar_field_halo_exch(domain, 'DMSTracers', timeLevel=1)
!      end if
!      if (config_use_MacroMoleculesTracers) then
!         call mpas_dmpar_field_halo_exch(domain, 'MacroMoleculesTracers', timeLevel=1)
!      end if
!      call mpas_timer_stop("1st stage-prognostic halo update")

!      do iTracer = 1, num_tracers
!         do iCell = 1, nCellsOwned
!            do k=1,maxLevelCell(iCell)
!               if(abs(tracerGroupNew(iTracer, k, iCell))>1.5_RKIND) then
!                  print*, '1st stage Old ', tracerGroup(iTracer, k, iCell), iTracer, k, maxLevelCell(iCell), iCell
!                  print*, '1st stage New ', tracerGroupNew(iTracer, k, iCell), iTracer, k, maxLevelCell(iCell), iCell
!               end if
!            end do
!         end do
!      end do

      !deallocate(TracerGroupTendHorActive, TracerGroupTendVertActive)
      !deallocate(TracerGroupTendHorEcosysTracers)
      deallocate(phi1JTot, JacZ, rhs, tracer_cur, &
                 TracerGroupTendHor, TracerGroupTendVert, &
                 A, vec)

      !print*,'second stage complete'

      call mpas_timer_stop('time-stepping time for each iteration')

   end subroutine ocn_time_integrator_ETD

   subroutine ocn_tracer_vert_advection_std(tracers, w, layerThickness, &
                                           dt, tend)

      !-----------------------------------------------------------------
      ! Input/Output variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(inout) :: &
         tend            !< [inout] Accumulated tracer tendency

      !-----------------------------------------------------------------
      ! Input variables
      !-----------------------------------------------------------------

      real (kind=RKIND), dimension(:,:,:), intent(in) :: &
         tracers               !< [in] current tracer values

      real (kind=RKIND), dimension(:,:), intent(in) :: &
         w,                   &!< [in] Vertical velocity
         layerThickness        !< [in] Thickness

      real (kind=RKIND), intent(in) :: &
         dt                    !< [in] Timestep

      !-----------------------------------------------------------------
      ! Local variables
      !-----------------------------------------------------------------

      integer :: &
         i, k,         &! loop indices for neighbors, vertical
         kmin, kmax,   &! min, max active vertical layers
         iCell, iEdge, &! loop indices for cells, edges
         iTracer,      &! tracer index
         cell1, cell2, &! neighbor cell indices across edge
         numTracers

      real (kind=RKIND) :: &
         tracerWgt,        &! local temporary
         invAreaCell1,     &! inverse cell area
         vertDivFactor      ! vertical divergence factor

      real (kind=RKIND), dimension(:,:), allocatable :: &
         tracerCur,        &! reordered tracer at current time
         highOrderFlxHorz, &! high-order flux in horizontal
         highOrderFlxVert   ! high-order flux in vertical

      real (kind=RKIND), parameter :: &
         eps = 1.e-10_RKIND  ! small value to avoid div by zero

      ! End preamble
      !-----------------------------------------------------------------
      ! Begin code

      ! Get dimensions
      numTracers = size(tracers,dim=1)

      ! Allocate some arrays
      allocate(tracerCur       (nVertLevels  , nCellsAll), &
               highOrderFlxHorz(nVertLevels  , nEdgesAll), &
               highOrderFlxVert(nVertLevels+1, nCellsAll))

      !$acc enter data &
      !$acc    create(tracerCur, highOrderFlxHorz, highOrderFlxVert)

      ! Loop over tracers. One tracer is advected at a time. It is
      ! copied into a temporary array in order to improve locality
      do iTracer = 1, numTracers

        ! Initialize variables for use in this iTracer iteration
#ifdef MPAS_OPENACC
        !$acc parallel loop collapse(2) &
        !$acc    present(tracerCur, tracers)
#else
        !$omp parallel
        !$omp do schedule(runtime) private(k)
#endif
        do iCell = 1, nCellsAll
        do k=1,nVertLevels
           tracerCur(k,iCell) = tracers(iTracer,k,iCell)
        end do
        end do
#ifndef MPAS_OPENACC
        !$omp end do
#endif

#ifdef MPAS_OPENACC
        !$acc parallel loop collapse(2) &
        !$acc    present(highOrderFlxHorz)
#else
        !$omp do schedule(runtime) private(k)
#endif
        do iEdge = 1, nEdgesAll
        do k=1,nVertLevels
           highOrderFlxHorz(k, iEdge) = 0.0_RKIND
        end do
        end do
#ifndef MPAS_OPENACC
        !$omp end do
#endif

        ! Compute the high order vertical flux.

        call ocn_tracer_advection_vert_flx(tracerCur, w, &
                           layerThickness, highOrderFlxVert)

        ! Accumulate the scaled high order vertical tendencies.
        vertDivFactor = 1.0_RKIND

#ifdef MPAS_OPENACC
        !$acc parallel loop &
        !$acc    present(minLevelCell, maxLevelCell, tend, &
        !$acc            highOrderFlxVert) &
        !$acc    private(k, kmin, kmax)
#else
        !$omp do schedule(runtime) private(k, kmin, kmax)
#endif
        do iCell = 1, nCellsOwned
           kmin = minLevelCell(iCell)
           kmax = maxLevelCell(iCell)

           do k = kmin,kmax
              tend(iTracer, k, iCell) = tend(iTracer, k, iCell) + &
                         vertDivFactor*(highOrderFlxVert(k+1, iCell) &
                                      - highOrderFlxVert(k, iCell))
           end do ! k loop
        end do ! iCell loop
#ifndef MPAS_OPENACC
        !$omp end do
        !$omp end parallel
#endif
      end do ! iTracer loop

      !$acc exit data &
      !$acc    delete(tracerCur, highOrderFlxHorz, highOrderFlxVert)

      deallocate(tracerCur, highOrderFlxHorz, highOrderFlxVert)

   end subroutine ocn_tracer_vert_advection_std!}}}

   subroutine jacobianVert(JacZ, layerThickness, w, iCell)

      real (kind=RKIND), dimension(:), intent(out) :: &
         JacZ            ! Jacobian at each point stored as vector

      real (kind=RKIND), dimension(:,:), intent(in) :: layerThickness    !< Input: thickness at cell center
      real (kind=RKIND), dimension(:,:), intent(in) :: w                 !< Input: Vertical velocity 
      integer, intent(in) :: iCell

      real (kind=RKIND), dimension(:,:), allocatable :: JacVertAdv
      integer :: i, j, k, N

      N = maxLevelCell(iCell)
      allocate(JacVertAdv(N,N))

      JacZ(:) = 0.0_RKIND
      JacVertAdv(:,:) = 0.0_RKIND

      select case (vertOrder)
      case (vertOrder2)

         JacVertAdv(1,1) = w(2,iCell) * (layerThickness(2, iCell) / (layerThickness(2, iCell) + layerThickness(1, iCell))) / layerThickness(1, iCell)
         JacVertAdv(1,2) = w(2,iCell) * (layerThickness(1, iCell) / (layerThickness(2, iCell) + layerThickness(1, iCell))) / layerThickness(2, iCell)

         do k = 2, N-1
            JacVertAdv(k,k-1) = - w(k,iCell) * (layerThickness(k, iCell) / (layerThickness(k, iCell) + layerThickness(k-1, iCell))) / layerThickness(k-1, iCell)
            JacVertAdv(k,k) = w(k+1,iCell) * (layerThickness(k+1, iCell) / (layerThickness(k, iCell) + layerThickness(k+1, iCell))) / layerThickness(k, iCell) - &
                              w(k,iCell) * (layerThickness(k-1, iCell) / (layerThickness(k, iCell) + layerThickness(k-1, iCell))) / layerThickness(k, iCell)
            JacVertAdv(k,k+1) = w(k+1,iCell) * (layerThickness(k, iCell) / (layerThickness(k, iCell) + layerThickness(k+1, iCell))) / layerThickness(k+1, iCell)
         end do

         JacVertAdv(N,N-1) = - w(N,iCell) * (layerThickness(N, iCell) / (layerThickness(N, iCell) + layerThickness(N-1, iCell))) / layerThickness(N-1, iCell)
         JacVertAdv(N,N) = - w(N,iCell) * (layerThickness(N-1, iCell) / (layerThickness(N, iCell) + layerThickness(N-1, iCell))) / layerThickness(N, iCell)

      end select ! vertOrder

      !do i = 1, N
      !   do j = 1, N
      !      JacVert(i,j) = JacVertAdv(i,j) + JacVertDiff(i,j)
      !   end do
      !end do

      !Assemble the jacobian by colum using JacVert
      do i = 1, N
         do j = 1, N
            JacZ((j-1)*N+i) = JacVertAdv(i,j)
         end do
      end do

      deallocate(JacVertAdv)

   end subroutine jacobianVert

   subroutine ocn_phi_function(Jac, CFL_pow, phi1J, NLayers, dt)!{{{

   !***********************************************************************
   !
   !  routine ocn_phi_function
   !
   !> \brief   Computes phi1(A) for the ETD time-stepping scheme
   !> \author  Sara Calandrini
   !> \date    July 2020
   !> \details
   !>  This routine computes the phi1 function of the linear operator for
   !>  the tracers based on the scaling and squaring algorithm.
   !
   !-----------------------------------------------------------------------

      integer, intent(in) :: CFL_pow, NLayers

      real (kind=RKIND), intent(in), dimension(:) :: &
         Jac             ! Jacobian at each point stored as vector

      real (kind=RKIND), intent(inout), dimension(:,:) :: phi1J
      real (kind=RKIND), intent(in) :: dt

      !real (kind=RKIND), dimension(16) :: JacProva
      real (kind=RKIND), dimension(NLayers,NLayers) :: A
      real (kind=RKIND), dimension(NLayers,NLayers) :: Temp
      real (kind=RKIND), dimension(NLayers,NLayers) :: phi1A
      real (kind=RKIND), dimension(NLayers,NLayers) :: Id

      real (kind=RKIND) :: power = 0.5_RKIND
      integer :: k1, k2

      !if (NLayers==4) then
      !   JacProva =
      !   (/-0.18486511111111,0.18486511111111,0.0,0.0,4.e-06,-0.24648948148148,0.24648548148148,0.0,0.0,4e-06,-0.18486911111111,0.18486511111111,0.0,0.0,4e-06,-4e-06/)
      !else 
      !   if (iCell==10 .or. iCell==90) then 
      !      print*, Jac(1:50)
      !   end if 
      !end if 

      phi1J(:,:)=0.0_RKIND
      phi1A(:,:)=0.0_RKIND
      Id(:,:)=0.0_RKIND

      do k1 = 1, NLayers
         phi1A(k1,k1)=1.0_RKIND
         Id(k1,k1)=1.0_RKIND
         do k2 = 1, NLayers
           A(k1,k2) = (dt/(2**CFL_pow))*Jac((k2-1)*NLayers+k1)
         end do
      end do

      CALL DGEMM('N','N',NLayers,NLayers,NLayers,1.0_RKIND,A,Nlayers,A,NLayers,0.0_RKIND,Temp,Nlayers) !A^2
      CALL DGEMM('N','N',NLayers,NLayers,NLayers,0.5_RKIND,A,Nlayers,Id,NLayers,1.0_RKIND,phi1A,Nlayers) !I+1/2A
      CALL DGEMM('N','N',NLayers,NLayers,NLayers,1./6._RKIND,Temp,Nlayers,Id,NLayers,1.0_RKIND,phi1A,Nlayers) !I+1/2A+1/6A^2
      CALL DGEMM('N','N',NLayers,NLayers,NLayers,1./24._RKIND,Temp,Nlayers,A,NLayers,1.0_RKIND,phi1A,Nlayers) !I+1/2A+1/6A^2+1/24A^3

      do k1 = 1, CFL_pow
         Temp(:,:)=0.0_RKIND
         CALL DGEMM('N','N',NLayers,NLayers,NLayers,1.0_RKIND,phi1A,Nlayers,phi1A,NLayers,0.0_RKIND,Temp,Nlayers) !phi1A^2
         if (k1==1) then
            power = 0.5_RKIND
         elseif (k1==2) then
            power = 1.0_RKIND
         elseif (k1>2) then
           power = 2**(k1-2)
         end if
         CALL DGEMM('N','N',NLayers,NLayers,NLayers,power,A,Nlayers,Temp,NLayers,1.0_RKIND,phi1A,Nlayers) !phi1A + 2^(k1-1)*A*phi1A^2
      end do

      do k1 = 1, NLayers
         do k2 = 1, NLayers
            phi1J(k1,k2) = phi1A(k1,k2)
         end do
      end do

   end subroutine ocn_phi_function!}}

   subroutine ocn_phi_function_taylor(Jac, phi1J, NLayers, dt)!{{{

   !***********************************************************************
   !
   !  routine ocn_phi_function_taylor
   !
   !> \brief   Computes phi1(A) for the ETD time-stepping scheme
   !> \author  Sara Calandrini
   !> \date    June 20201
   !> \details
   !>  This routine computes the phi1 function of the linear operator for
   !>  the tracers using just the taylor expansion of exp.
   !
   !-----------------------------------------------------------------------

      integer, intent(in) :: NLayers
      real (kind=RKIND), dimension(:), intent(in) :: &
         Jac             ! Jacobian at each point stored as vector
      real (kind=RKIND), intent(inout), dimension(:,:) :: phi1J
      real (kind=RKIND), intent(in) :: dt

      real (kind=RKIND), dimension(NLayers,NLayers) :: A
      real (kind=RKIND), dimension(NLayers,NLayers) :: phi1A
      real (kind=RKIND), dimension(NLayers,NLayers) :: Temp
      real (kind=RKIND), dimension(NLayers,NLayers) :: Id

      integer :: k1, k2

      phi1J(:,:)=0.0_RKIND
      phi1A(:,:)=0.0_RKIND
      Id(:,:)=0.0_RKIND

      do k1 = 1, NLayers
         phi1A(k1,k1)=1.0_RKIND
         Id(k1,k1)=1.0_RKIND
         do k2 = 1, NLayers
           A(k1,k2) = dt*Jac((k2-1)*NLayers+k1)
         end do
      end do

      CALL DGEMM('N','N',NLayers,NLayers,NLayers,1.0_RKIND,A,Nlayers,A,NLayers,0.0_RKIND,Temp,Nlayers) !A^2
      CALL DGEMM('N','N',NLayers,NLayers,NLayers,0.5_RKIND,A,Nlayers,Id,NLayers,1.0_RKIND,phi1A,Nlayers) !I+1/2A
      CALL DGEMM('N','N',NLayers,NLayers,NLayers,1./6._RKIND,Temp,Nlayers,Id,NLayers,1.0_RKIND,phi1A,Nlayers) !I+1/2A+1/6A^2
      CALL DGEMM('N','N',NLayers,NLayers,NLayers,1./24._RKIND,Temp,Nlayers,A,NLayers,1.0_RKIND,phi1A,Nlayers) !I+1/2A+1/6A^2+1/24A^3
      !CALL DGEMM('N','N',NLayers,NLayers,NLayers,1./120._RKIND,Temp,Nlayers,Temp,NLayers,1.0_RKIND,phi1A,Nlayers) !I+1/2A+1/6A^2+1/24A^3+1/120A^4

      do k1 = 1, NLayers
         do k2 = 1, NLayers
            phi1J(k1,k2) = phi1A(k1,k2)
         end do
      end do

      !print*, '----------------------'

   end subroutine ocn_phi_function_taylor!}}

end module ocn_time_integration_ETD
