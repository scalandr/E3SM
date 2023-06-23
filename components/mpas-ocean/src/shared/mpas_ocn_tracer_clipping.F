! Copyright (c) 2013,  Los Alamos National Security, LLC (LANS)
! and the University Corporation for Atmospheric Research (UCAR).
!
! Unless noted otherwise source code is licensed under the BSD license.
! Additional copyright and license information can be found in the LICENSE file
! distributed with this code, or at http://mpas-dev.github.com/license.html
!
!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
!
!  ocn_tracer_clipping
!
!> \brief MPAS ocean driver to clip the tracers values after the 
!> time-stepping to maintaing monotonicity
!> \author Sara Calandrini
!> \date   April 2022
!> \details
!>  This module contains the main driver routine for clipping tracer
!>  values after thte time-stepping.
!
!-----------------------------------------------------------------------

module ocn_tracer_clipping

   use mpas_timer
   use mpas_derived_types
   use mpas_pool_routines
   use mpas_threading

   use ocn_constants
   use ocn_config
   use ocn_mesh
   use ocn_diagnostics_variables

   implicit none
   private
   save

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

   public :: ocn_tracer_caas, &
             ocn_tracer_caas_init

   !--------------------------------------------------------------------
   !
   ! Private module variables
   !
   !--------------------------------------------------------------------

   logical :: caasOn


!***********************************************************************

contains

!***********************************************************************
!
!  routine ocn_tracer_caas
!
!> \brief   Clippes passive tracer values after the time-stepping 
!> \author  Sara Calandrini
!> \date    April 2022
!> \details
!>  This routine clippes the passive tracers values after the 
!>  time-stepping based on the CAAS algorithm.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_caas(statePool, tracersPool, err)!{{{
      !-----------------------------------------------------------------
      !
      ! input variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(in) :: statePool !< Input: State information

      !-----------------------------------------------------------------
      !
      ! input/output variables
      !
      !-----------------------------------------------------------------

      type (mpas_pool_type), intent(inout) :: tracersPool !< Input: Tracers information

      !-----------------------------------------------------------------
      !
      ! output variables
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      !-----------------------------------------------------------------
      !
      ! local variables
      !
      !-----------------------------------------------------------------

      integer :: iCell, cell, i, k, iTracer, num_tracers, num_tracers_eco

      real (kind=RKIND), dimension(:,:,:), pointer :: tracerGroup
      real (kind=RKIND), dimension(:,:,:), pointer :: tracerGroupOld
      real (kind=RKIND), dimension(:,:,:), pointer :: tracerGroupEcosysTracersOld
      real (kind=RKIND), dimension(:,:,:), pointer :: tracerGroupEcosysTracers
      real (kind=RKIND), dimension(:,:), pointer :: layerThickness

      real (kind=RKIND), dimension(:), allocatable :: q_g, q_bar_star, &
                                               x_bar_mod, b, m, v_mod, &
                                               neighbors
      real (kind=RKIND), dimension(:,:,:), allocatable :: x_bar, v, &
                                   upper, lower, upperEco, lowerEco, &
                                   boundUpper, boundLower

      err = 0 

      if (.not.caasOn) return

      call mpas_timer_start("tracer caas")

      call mpas_pool_get_array(tracersPool, 'debugTracers', tracerGroupOld, 1)
      call mpas_pool_get_array(tracersPool, 'debugTracers', tracerGroup, 2)
      call mpas_pool_get_array(tracersPool, 'ecosysTracers', tracerGroupEcosysTracersOld, 1)
      call mpas_pool_get_array(tracersPool, 'ecosysTracers', tracerGroupEcosysTracers, 2)

      call mpas_pool_get_array(statePool, 'layerThickness', layerThickness, 2)

      !num_tracers = size(tracerGroup, dim=1)
      num_tracers = 3
      !num_tracers_eco = size(tracerGroupEcosysTracers, dim=1)
      num_tracers_eco = 3

      allocate(q_g(num_tracers_eco))
      allocate(q_bar_star(num_tracers_eco))
      allocate(x_bar_mod(num_tracers_eco))
      allocate(b(num_tracers_eco))
      allocate(m(num_tracers_eco))
      allocate(v_mod(num_tracers_eco))
      allocate(x_bar(num_tracers_eco, nVertLevels, nCellsAll))
      allocate(v(num_tracers_eco, nVertLevels, nCellsAll))
      allocate(upper(num_tracers_eco, nVertLevels, nCellsAll))
      allocate(lower(num_tracers_eco, nVertLevels, nCellsAll))

      if (config_use_ecosysTracers) then
         allocate(upperEco(num_tracers_eco, nVertLevels, nCellsAll))
         allocate(lowerEco(num_tracers_eco, nVertLevels, nCellsAll))
         allocate(boundUpper(num_tracers_eco, nVertLevels, nCellsAll))
         allocate(boundLower(num_tracers_eco, nVertLevels, nCellsAll))
         allocate(neighbors(maxval(nEdgesOnCell)+2))
      end if

      !print*, 'clipping ', num_tracers, num_tracers_eco, nVertLevels, nCells

      if (config_use_debugTracers) then

         do iTracer = 1, 2
            do iCell = 1, nCellsHalo(1)
               !upper(iTracer,:,iCell) = 1.000001_RKIND - tracerGroup(iTracer, :, iCell)
               !lower(iTracer,:,iCell) = 0.999999_RKIND - tracerGroup(iTracer, :, iCell)
               upper(iTracer,:,iCell) = 1.000001_RKIND - tracerGroup(iTracer, :, iCell)
               lower(iTracer,:,iCell) = -0.000001_RKIND - tracerGroup(iTracer, :, iCell)
               !upper(iTracer,:,iCell) = 1.001_RKIND - tracerGroup(iTracer, :, iCell) 1.5
               !lower(iTracer,:,iCell) = -0.001_RKIND - tracerGroup(iTracer, :, iCell) -0.5
            end do
         end do

         do iCell = 1, nCellsHalo(1)
            upper(3,:,iCell) = 1.000001_RKIND - tracerGroup(3, :, iCell)
            lower(3,:,iCell) = 0.999999_RKIND - tracerGroup(3, :, iCell)
         end do
         
         do iTracer = 1, num_tracers
            v(iTracer,:,:) = 0.0_RKIND
            x_bar(iTracer,:,:) = 0.0_RKIND
            x_bar_mod(iTracer) = 0.0_RKIND
            q_g(iTracer) = 0.0_RKIND
            q_bar_star(iTracer) = 0.0_RKIND
            b(iTracer) = 0.0_RKIND
            m(iTracer) = 0.0_RKIND
            do iCell = 1, nCellsHalo(1)
               do k=1, maxLevelCell(iCell)
                  x_bar(iTracer,k,iCell) = max(lower(iTracer,k,iCell), min(upper(iTracer,k,iCell),0.0_RKIND))
                  x_bar_mod(iTracer) = x_bar_mod(iTracer) + x_bar(iTracer,k,iCell)
               end do
            end do
            do iCell = 1, nCellsHalo(1)
               do k=1, maxLevelCell(iCell)
                  !q_g(iTracer) = q_g(iTracer) + tracerGroup(iTracer, k, iCell) &
                  !   * layerThickness(k,iCell) * areaCell(iCell)
                  !if(abs(tracerGroupOld(iTracer, k, iCell))>2.0_RKIND) print*, &
                  !   'old ', tracerGroupOld(iTracer, k, iCell), iTracer, k, iCell
                  q_g(iTracer) = q_g(iTracer) + tracerGroupOld(iTracer, k, iCell)
                  !if(abs(tracerGroup(iTracer, k, iCell))>2.0_RKIND) print*, &
                  !   'new ', tracerGroup(iTracer, k, iCell), iTracer, k, iCell
                  q_bar_star(iTracer) = q_bar_star(iTracer) + &
                     tracerGroup(iTracer, k, iCell)
               end do
            end do
            b(iTracer) = q_g(iTracer) - q_bar_star(iTracer)
            m(iTracer) = b(iTracer) - x_bar_mod(iTracer)
            !if(iTracer==1)  print*, 'b ', b(iTracer)
            !print*, 'm ', m(iTracer)
            if (m(iTracer) > 1.0e-12) then
               v_mod(iTracer) = 0.0_RKIND
               do iCell = 1, nCellsHalo(1)
                  do k=1, maxLevelCell(iCell)
                     v(iTracer,k,iCell) = upper(iTracer,k,iCell) - x_bar(iTracer,k,iCell)
                     v_mod(iTracer) = v_mod(iTracer) + v(iTracer,k,iCell)
                  end do
               end do
!               print*, 'DEBUG', iTracer, x_bar_mod(iTracer), b(iTracer), m(iTracer), v_mod(iTracer)
               do iCell = 1, nCellsHalo(1)
                  v(iTracer,:,iCell) = v(iTracer,:,iCell)/v_mod(iTracer)
               end do
            elseif (m(iTracer) < -1.0e-12) then
               v_mod(iTracer) = 0.0_RKIND
               do iCell = 1, nCellsHalo(1)
                  do k=1, maxLevelCell(iCell)
                     v(iTracer,k,iCell) = x_bar(iTracer,k,iCell) - lower(iTracer,k,iCell)
                     v_mod(iTracer) = v_mod(iTracer) + v(iTracer,k,iCell)
                  end do
               end do
!               print*, 'DEBUG', iTracer, x_bar_mod(iTracer), b(iTracer), m(iTracer), v_mod(iTracer)
               do iCell = 1, nCellsHalo(1)
                  v(iTracer,:,iCell) = v(iTracer,:,iCell)/v_mod(iTracer)
               end do
            elseif (m(iTracer)<1.0e-12 .or. m(iTracer)>-1.0e-12) then
               v(iTracer,:,:) = 0.0_RKIND
!               print*, 'DEBUG', iTracer, x_bar_mod(iTracer), b(iTracer), m(iTracer), v_mod(iTracer)
            end if
            do iCell = 1, nCellsHalo(1)
               do k=1, maxLevelCell(iCell)
                  tracerGroup(iTracer,k,iCell) = tracerGroup(iTracer,k,iCell) + &
                     x_bar(iTracer,k,iCell) + m(iTracer)*v(iTracer,k,iCell)
                  if (abs(tracerGroup(iTracer,k,iCell))>2.0_RKIND) print*, &
                     tracerGroup(iTracer,k,iCell), iTracer, k, iCell, m(iTracer)
               end do
            end do
         end do
 
      end if

      if (config_use_ecosysTracers) then

         do iTracer = 1, num_tracers_eco
            do iCell = 1, nCellsAll
               neighbors(:) = 0.0_RKIND
               do k=1, maxLevelCell(iCell)
                  do i = 1, nEdgesOnCell(iCell)
                     cell = cellsOnCell(i,iCell)
                     neighbors(i) = tracerGroupEcosysTracers(iTracer, k, cell)
                  end do
                  if (k==1) then
                     neighbors(nEdgesOnCell(iCell)+1) = tracerGroupEcosysTracers(iTracer, k+1, iCell)
                     boundUpper(iTracer,k,iCell) = maxval(neighbors(1:nEdgesOnCell(iCell)+1))
                     boundLower(iTracer,k,iCell) = max(1e-15, minval(neighbors(1:nEdgesOnCell(iCell)+1)))
                  elseif (k==maxLevelCell(iCell)) then
                     neighbors(nEdgesOnCell(iCell)+1) = tracerGroupEcosysTracers(iTracer, k-1, iCell)
                     boundUpper(iTracer,k,iCell) = maxval(neighbors(1:nEdgesOnCell(iCell)+1))
                     boundLower(iTracer,k,iCell) = max(1e-15, minval(neighbors(1:nEdgesOnCell(iCell)+1)))
                  elseif (k>1 .and. k<maxLevelCell(iCell)) then
                     neighbors(nEdgesOnCell(iCell)+1) = tracerGroupEcosysTracers(iTracer, k-1, iCell)
                     neighbors(nEdgesOnCell(iCell)+2) = tracerGroupEcosysTracers(iTracer, k+1, iCell)
                     boundUpper(iTracer,k,iCell) = maxval(neighbors(1:nEdgesOnCell(iCell)+2))
                     boundLower(iTracer,k,iCell) = max(1e-15, minval(neighbors(1:nEdgesOnCell(iCell)+2)))
                  end if
                  !print*, boundLower(iTracer,k,iCell), boundUpper(iTracer,k,iCell)
               end do
               upperEco(iTracer,:,iCell) = boundUpper(iTracer,:,iCell) - tracerGroupEcosysTracers(iTracer, :, iCell)
               lowerEco(iTracer,:,iCell) = boundLower(iTracer,:,iCell) - tracerGroupEcosysTracers(iTracer, :, iCell)
            end do
         end do
         do k=1, maxLevelCell(1)
            !print*, boundUpper(24,k,1), boundLower(24,k,1)
         end do

         do iTracer = 1, num_tracers_eco
            v(iTracer,:,:) = 0.0_RKIND
            x_bar(iTracer,:,:) = 0.0_RKIND
            x_bar_mod(iTracer) = 0.0_RKIND
            q_g(iTracer) = 0.0_RKIND
            q_bar_star(iTracer) = 0.0_RKIND
            b(iTracer) = 0.0_RKIND
            m(iTracer) = 0.0_RKIND
            do iCell = 1, nCellsAll
               do k=1, maxLevelCell(iCell)
                  x_bar(iTracer,k,iCell) = max(lowerEco(iTracer,k,iCell), min(upperEco(iTracer,k,iCell),0.0_RKIND))
                  x_bar_mod(iTracer) = x_bar_mod(iTracer) + x_bar(iTracer,k,iCell)
               end do
            end do
            do iCell = 1, nCellsAll
               do k=1, maxLevelCell(iCell)
                  !q_g(iTracer) = q_g(iTracer) + tracerGroup(iTracer, k, iCell) &
                  !   * layerThickness(k,iCell) * areaCell(iCell)
                  !if(abs(tracerGroupOld(iTracer, k, iCell))>2.0_RKIND) print*, &
                  !   'old ', tracerGroupOld(iTracer, k, iCell), iTracer, k, iCell
                  q_g(iTracer) = q_g(iTracer) + tracerGroupEcosysTracersOld(iTracer, k, iCell)
                  !if(abs(tracerGroup(iTracer, k, iCell))>2.0_RKIND) print*, &
                  !   'new ', tracerGroup(iTracer, k, iCell), iTracer, k, iCell
                  q_bar_star(iTracer) = q_bar_star(iTracer) + &
                     tracerGroupEcosysTracers(iTracer, k, iCell)
               end do
            end do
            b(iTracer) = q_g(iTracer) - q_bar_star(iTracer)
            m(iTracer) = b(iTracer) - x_bar_mod(iTracer)
            !if(iTracer==1)  print*, 'b ', b(iTracer)
            !print*, 'm ', m(iTracer)
            if (m(iTracer) > 1.0e-12) then
               v_mod(iTracer) = 0.0_RKIND
               do iCell = 1, nCellsAll
                  do k=1, maxLevelCell(iCell)
                     v(iTracer,k,iCell) = upperEco(iTracer,k,iCell) - x_bar(iTracer,k,iCell)
                     v_mod(iTracer) = v_mod(iTracer) + v(iTracer,k,iCell)
                  end do
               end do
!               print*, iTracer, x_bar_mod(iTracer), b(iTracer), m(iTracer), v_mod(iTracer)
               do iCell = 1, nCellsAll
                  v(iTracer,:,iCell) = v(iTracer,:,iCell)/v_mod(iTracer)
               end do
            elseif (m(iTracer) < -1.0e-12) then
               v_mod(iTracer) = 0.0_RKIND
               do iCell = 1, nCellsAll
                  do k=1, maxLevelCell(iCell)
                     v(iTracer,k,iCell) = x_bar(iTracer,k,iCell) - lowerEco(iTracer,k,iCell)
                     v_mod(iTracer) = v_mod(iTracer) + v(iTracer,k,iCell)
                  end do
               end do
!               print*, iTracer, x_bar_mod(iTracer), b(iTracer), m(iTracer), v_mod(iTracer)
               do iCell = 1, nCellsAll
                  v(iTracer,:,iCell) = v(iTracer,:,iCell)/v_mod(iTracer)
               end do
            elseif (m(iTracer)<1.0e-12 .or. m(iTracer)>-1.0e-12) then
               v(iTracer,:,:) = 0.0_RKIND
            end if
            do iCell = 1, nCellsAll
               do k=1, maxLevelCell(iCell)
                  tracerGroupEcosysTracers(iTracer,k,iCell) = tracerGroupEcosysTracers(iTracer,k,iCell) + &
                     x_bar(iTracer,k,iCell) + m(iTracer)*v(iTracer,k,iCell)
                  !if (abs(tracerGroupEcosysTracers(iTracer,k,iCell))>2.0_RKIND) print*, 'ECO', &
                  !   tracerGroupEcosysTracers(iTracer,k,iCell), iTracer, k, iCell, m(iTracer)
               end do
            end do
         end do

      end if

     deallocate(q_g,q_bar_star,x_bar_mod,b,m,v_mod)
     deallocate(x_bar,v,upper,lower)
     if (config_use_ecosysTracers) then
        deallocate(upperEco,lowerEco,boundUpper,boundLower,neighbors)
     end if

      call mpas_timer_stop("tracer caas")

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_caas!}}}

!***********************************************************************
!
!  routine ocn_tracer_caas_init
!
!> \brief   Initializes ocean tracer horizontal mixing quantities
!> \author  Sara Calandrini
!> \date    April 2022
!> \details
!>  This routine initializes a variety of quantities related to
!>  Laplacian horizontal velocity mixing in the ocean.
!
!-----------------------------------------------------------------------

   subroutine ocn_tracer_caas_init(err)!{{{

   !--------------------------------------------------------------------

      !-----------------------------------------------------------------
      !
      ! call individual init routines for each parameterization
      !
      !-----------------------------------------------------------------

      integer, intent(out) :: err !< Output: error flag

      err = 0

      caasOn = .false.

      if ( config_use_tracer_clipping ) then
            caasOn = .true.
      endif

   !--------------------------------------------------------------------

   end subroutine ocn_tracer_caas_init!}}}

!***********************************************************************

end module ocn_tracer_clipping

!|||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
! vim: foldmethod=marker
