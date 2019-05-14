!/*****************************************************************************/
! *
! *  Elmer/Ice, a glaciological add-on to Elmer
! *  http://elmerice.elmerfem.org
! *
! * 
! *  This program is free software; you can redistribute it and/or
! *  modify it under the terms of the GNU General Public License
! *  as published by the Free Software Foundation; either version 2
! *  of the License, or (at your option) any later version.
! * 
! *  This program is distributed in the hope that it will be useful,
! *  but WITHOUT ANY WARRANTY; without even the implied warranty of
! *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! *  GNU General Public License for more details.
! *
! *  You should have received a copy of the GNU General Public License
! *  along with this program (in file fem/GPL-2); if not, write to the 
! *  Free Software Foundation, Inc., 51 Franklin Street, Fifth Floor, 
! *  Boston, MA 02110-1301, USA.
! *
! *****************************************************************************/
! ******************************************************************************
! *
! *  Authors: Denis Cohen, Thomas Zwinger, and Peter Raback
! *  Email:   denis.cohen@gmail.com
! *  Web:     http://elmerice.elmerfem.org
! *
! *  Original Date: 
! *   2016/02/16. Denis Cohen
! *   2016/08/05. Last update. Denis Cohen
! *****************************************************************************
! Contains three functions to model permafrost in rock beneath ice sheet:
! - PermafrostEnthalpy
! - PermafrostDensity
! - PermafrostConductivity
!
! Permafrost is modeled as a three-component mixture (rock + ice + water) with
! an effective heat capacity that depends on water content to mimick phase
! change near the freezing point.
!
! Rock density, conductivity and heat capacity can be functions of depth below 
! ground level
!
! Two models for water content as a function of temperature are available:
! - power-law (e.g. Cutler et al, 2000, A numerical investigation of 
!   ice-lobe-permafrost interaction around the southern Laurentide ice sheet, 
!   J. Glaciol., 46, 311-325)
! - exponential (e.g. Willeit and Ganopolski, 2015, Coupled Northern Hemisphere 
!   permafrost-ice-sheet evolution over the last glacial cycle, Clim. Past, 11,
!   1165-1180)
! *****************************************************************************

!==============================================================================
FUNCTION PermafrostEnthalpy(Model, Node, Temp) RESULT(enthalpy)
!==============================================================================

  USE DefUtils
  USE types
  USE CoordinateSystems
  USE SolverUtils
  USE ElementDescription

  IMPLICIT None

  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: Temp, enthalpy

  ! Local variables
  TYPE(Element_t),POINTER :: Element
  TYPE(ValueList_t), POINTER :: Material
  REAL(KIND=dp), ALLOCATABLE :: Porosity(:)
  REAL(KIND=dp), ALLOCATABLE :: Capacity(:)
  REAL(KIND=dp), ALLOCATABLE :: Density(:)

  TYPE(Variable_t), POINTER :: DepthVar, DepthVar2
  REAL(KIND=dp) :: Depth, Depth2

  REAL(KIND=dp) :: por                 ! Porosity
  REAL(KIND=dp) :: porscale            ! Porosity length scale
  REAL(KIND=dp) :: pordepth            ! Porosity as a function of depth
  ! r = rock, w = water, i = ice
  REAL(KIND=dp) :: phir, phiw          ! Volume fractions
  REAL(KIND=dp) :: rhor, rhow, rhoi    ! Densities 
  REAL(KIND=dp) :: Cr, Cw, Ci          ! Heat capacities
  REAL(KIND=dp) :: L                   ! Latent heat
  REAL(KIND=dp) :: ClausiusClapeyron   ! Clausius-Clapeyron constant
  REAL(KIND=dp) :: Tpmp                ! Melting point temperature of ice
  REAL(KIND=dp) :: a, b, Tstar, dT     ! Params for powerlaw/exponential model
  REAL(KIND=dp) :: fw                  ! Params for exponential model
  REAL(KIND=dp) :: iceDepth, pice, prock, press ! For computing pressures
  REAL(KIND=dp), PARAMETER :: factor=1.0d-06
  
  INTEGER :: N, istat, i, k

  CHARACTER(LEN=MAX_NAME_LEN) :: PermafrostModel, DepthVarName

  LOGICAL :: FirstTime = .TRUE.
  LOGICAL :: Found, UnfoundFatal, ScaleSystem

  SAVE porscale
  SAVE rhow, rhoi
  SAVE Cw, Ci, L
  SAVE a, b, dT
  SAVE FirstTime, PermafrostModel
  SAVE Porosity, Capacity, Density, ClausiusClapeyron

  Element => Model % CurrentElement
  Material => GetMaterial(Element)
  IF (.NOT.ASSOCIATED(Material)) THEN
    CALL FATAL('Permafrost', 'No Material found')
  END IF

  IF (FirstTime) THEN
    FirstTime = .FALSE.

    N = Model % MaxElementNodes
    ALLOCATE(Porosity(N),     & 
         Capacity(N),     & 
         Density(N),      & 
         STAT=istat)
    IF (istat /= 0) THEN
      CALL FATAL(  'USF_Permafrost', 'Enthalpy memory allocation error' )
    ELSE
      WRITE(Message,'(a)') 'Enthalpy memory allocation done'
      CALL INFO("Permafrost",Message,Level=4)
    END IF

    !-----------------------------------------------
    ! Read parameters from sif file Material section
    !-----------------------------------------------
    PermafrostModel = GetString( Material, 'Permafrost Model', Found )
    IF (.NOT. Found) THEN
      CALL FATAL('Permafrost', 'Could not find Permafrost Model')
    ENDIF

    ! Keep as is for the moment
    porscale = GetCReal( Material, 'Permafrost Porosity Depth Scale', Found )
    IF (.NOT. Found) THEN
      porscale = -1.0_dp ! Negative value means no depth dependence for por
    ENDIF

    L = GetCReal( Material, 'Latent Heat', Found )
    IF (.NOT. Found) THEN
      CALL FATAL('Permafrost', 'Could not find Latent Heat')
    ENDIF

    ClausiusClapeyron = GetCReal( Material, 'Clausius Clapeyron Constant', Found )
    IF (.NOT. Found) THEN
      CALL INFO('Permafrost', 'Could not find Latent Heat, using default SI value ', Level=3)
      ClausiusClapeyron = 9.8d-08
    ENDIF

    !--- Rock parameters ---
    !rhor = GetCReal( Material, 'Permafrost Density Rock', Found )
    !IF (.NOT. Found) THEN
    !   CALL FATAL('Permafrost', 'Could not find Permafrost Density Rock')
    !ENDIF

    !Cr = GetCReal( Material, 'Permafrost Heat Capacity Rock', Found )
    !IF (.NOT. Found) THEN
    !   CALL FATAL('Permafrost', 'Could not find Permafrost Heat Capacity Rock')
    !ENDIF

    !--- Water parameters ---
    rhow = GetCReal( Material, 'Permafrost Density Water', Found )
    IF (.NOT. Found) THEN
      CALL FATAL('Permafrost', 'Could not find Permafrost Density Water')
    ENDIF
    Cw = GetCReal( Material, 'Permafrost Heat Capacity Water', Found )
    IF (.NOT. Found) THEN
      CALL FATAL('Permafrost', 'Could not find Permafrost Heat Capacity Water')
    ENDIF

    !--- Ice parameters ---
    rhoi = GetCReal( Material, 'Permafrost Density Ice', Found )
    IF (.NOT. Found) THEN
      CALL FATAL('Permafrost', 'Could not find Permafrost Density Ice')
    ENDIF
    Ci = GetCReal( Material, 'Permafrost Heat Capacity Ice', Found )
    IF (.NOT. Found) THEN
      CALL FATAL('Permafrost', 'Could not find Permafrost Heat Capacity Ice')
    ENDIF

    !--- Power Law model parameters ---
    IF (TRIM(PermafrostModel) .EQ. "power law") THEN
      a = GetCReal( Material, 'Permafrost Power Law Factor', Found )
      IF (.NOT. Found) THEN
        CALL FATAL('Permafrost', 'Could not find Permafrost Power Law Factor')
      ENDIF
      b = GetCReal( Material, 'Permafrost Power Law Exponent', Found )
      IF (.NOT. Found) THEN
        CALL FATAL('Permafrost', 'Could not find Permafrost Power Law Exponent')
      ENDIF
      dT = GetCReal( Material, 'Permafrost Power law Temperature Offset', Found )
      IF (.NOT. Found) THEN
        CALL FATAL('Permafrost', 'Could not find Permafrost Power law Temperature Offset')
      ENDIF
      !--- Exponential model parameters ---
    ELSE IF (TRIM(PermafrostModel) .EQ. "exponential") THEN
      a = GetCReal( Material, 'Permafrost Exponential Temperature Interval', Found )
      IF (.NOT. Found) THEN
        CALL FATAL('Permafrost', 'Could not find Permafrost Exponential Temperature Interval')
      ENDIF
    ELSE
      CALL FATAL('Permafrost', 'Unknown Permafrost Model')
    ENDIF
  ENDIF ! End of FirstTime

  !-----------------------------------------------
  ! Get the depth of lower layer (rock layer)
  !-----------------------------------------------
  DepthVarName = GetString( Model % Solver % Values , 'Lower Depth Name', Found ) 
  IF (.NOT.Found) THEN
    WRITE(DepthVarName,'(A)') 'lower depth'
  END IF
  DepthVar => VariableGet(Model % Mesh % Variables, TRIM(DepthVarName))
  IF ( ASSOCIATED(DepthVar) ) THEN
    Depth = DepthVar % Values ( DepthVar % Perm(Node) )
  ELSE 
    Depth = 0.0_dp
  END IF

  !-----------------------------------------------
  ! Get the total depth below all layers
  !-----------------------------------------------
  DepthVarName = GetString( Model % Solver % Values , 'Total Depth Name', Found )
  IF (.NOT.Found) THEN
    WRITE(DepthVarName,'(A)') 'depth'
  END IF
  
  DepthVar2 => VariableGet(Model % Mesh % Variables, TRIM(DepthVarName))
  
  IF ( ASSOCIATED(DepthVar2) ) THEN
    Depth2 = DepthVar2 % Values ( DepthVar2 % Perm(Node) )
  ELSE
    Depth2 = 0.0_dp
  END IF

  ! In case there is no lower layer, use the only layer
  if (Depth == 0.0_dp) Depth = Depth2

  ! Get porosity
  N = GetElementNOFNodes(Element)
  Porosity(1:N) = ListGetReal ( Material, 'Permafrost Porosity', &
       N, Element % NodeIndexes, Found, & 
       UnfoundFatal=UnfoundFatal )
  IF (.NOT. Found) THEN
    CALL FATAL('Permafrost', 'Could not find Permafrost Porosity')
  ENDIF
  Capacity(1:N) = ListGetReal ( Material, 'Permafrost Heat Capacity Rock', &
       N, Element % NodeIndexes, Found, & 
       UnfoundFatal=UnfoundFatal )
  IF (.NOT. Found) THEN
    CALL FATAL('Permafrost', 'Could not find Permafrost Heat Capacity')
  ENDIF
  Density(1:N) = ListGetReal ( Material, 'Permafrost Density Rock', &
       N, Element % NodeIndexes, Found, & 
       UnfoundFatal=UnfoundFatal )
  IF (.NOT. Found) THEN
    CALL FATAL('Permafrost', 'Could not find Permafrost Density')
  ENDIF

  DO i = 1, N
    k = Element % NodeIndexes(i)
    If (Node .EQ. k) THEN
      por = Porosity(i)
      Cr = Capacity(i)
      rhor = Density(i)  ! Solid density
      EXIT
    END IF
  END DO

  !WRITE(MESSAGE,*) Node, Temp, Depth, Depth2, por
  !CALL Info('Permafrost Enthalpy', MESSAGE, level=11)

  !-------------------
  ! Start calculations
  !-------------------
  IF (porscale .LE. 0.0_dp) THEN
    pordepth = por
  ELSE
    pordepth = por * EXP(-Depth/porscale)
  ENDIF
  phir = 1.0_dp - pordepth

  !----------------------
  ! Compute Tpmp at depth
  !----------------------
  iceDepth = Depth2 - Depth
  pice = iceDepth * rhoi * 9.81_dp
  prock = Depth * rhow * 9.81_dp   ! This is for the fluid so use hydrostatic pressure
  press = pice + prock
  Tpmp = 273.15_dp - ClausiusClapeyron*press

  IF (TRIM(PermafrostModel) .EQ. "power law") THEN
    Tstar = Tpmp - Temp
    IF (Tstar <= dT) THEN
      phiw =  pordepth
    ELSE 
      phiw = (rhor * (1.0_dp - pordepth)/ rhow) * a * (Tstar)**b
    ENDIF
  ELSE IF (TRIM(PermafrostModel) .EQ. "exponential") THEN
    IF (Temp > Tpmp) then 
      fw = 1.0_dp
    ELSE
      fw = EXP(-((Temp - Tpmp)/a)**2)
    ENDIF
    phiw =  pordepth * fw
  ELSE
    CALL FATAL('Permafrost', 'Unknown Permafrost Model')
  ENDIF

  ! Check that phiw does not exceed porosity
  IF (phiw > pordepth) THEN
    phiw = pordepth
  ENDIF

  enthalpy = (phir*rhor*Cr + (pordepth-phiw)*rhoi*Ci + phiw*rhow*Cw)*(Temp) + phiw*rhow*L
  
  ! scale with remaining scaling factor of rhoi*Ci or rhow*L
  ScaleSystem = ListGetLogical( Material , 'Scale System', Found )
  IF (.NOT.Found) ScaleSystem=.FALSE.
  IF (ScaleSystem)  enthalpy = enthalpy * factor

END FUNCTION PermafrostEnthalpy

!==============================================================================
FUNCTION PermafrostCapacity(Model, Node, Temp) RESULT(effectivecapacity)
!==============================================================================

  USE DefUtils
  USE types
  USE CoordinateSystems
  USE SolverUtils
  USE ElementDescription

  IMPLICIT None

  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: Temp, effectivecapacity

  ! Local variables
  TYPE(Element_t),POINTER :: Element
  TYPE(ValueList_t), POINTER :: Material
  REAL(KIND=dp), ALLOCATABLE :: Porosity(:)
  REAL(KIND=dp), ALLOCATABLE :: Capacity(:)
  REAL(KIND=dp), ALLOCATABLE :: Density(:)

  TYPE(Variable_t), POINTER :: DepthVar, DepthVar2
  REAL(KIND=dp) :: Depth, Depth2

  REAL(KIND=dp) :: por                 ! Porosity
  REAL(KIND=dp) :: porscale            ! Porosity length scale
  REAL(KIND=dp) :: pordepth            ! Porosity as a function of depth
  ! r = rock, w = water, i = ice
  REAL(KIND=dp) :: phir, phiw          ! Volume fractions
  REAL(KIND=dp) :: rhor, rhow, rhoi    ! Densities 
  REAL(KIND=dp) :: Cr, Cw, Ci          ! Heat capacities
  REAL(KIND=dp) :: L                   ! Latent heat
  REAL(KIND=dp) :: ClausiusClapeyron   ! Clausius-Clapeyron constant
  REAL(KIND=dp) :: Tpmp                ! Melting point temperature of ice
  REAL(KIND=dp) :: a, b, Tstar, dT     ! Params for powerlaw/exponential model
  REAL(KIND=dp) :: fw                  ! Params for exponential model
  REAL(KIND=dp) :: iceDepth, pice, prock, press ! For computing pressures
  REAL(KIND=dp), PARAMETER :: factor=(31556926.0_dp)**(2.0_dp)

  INTEGER :: N, istat, i, k

  CHARACTER(LEN=MAX_NAME_LEN) :: PermafrostModel, DepthVarName

  LOGICAL :: FirstTime = .TRUE.
  LOGICAL :: Found, UnfoundFatal, ScaleSystem

  SAVE porscale
  SAVE rhow, rhoi
  SAVE Cw, Ci, L
  SAVE a, b, dT
  SAVE FirstTime, PermafrostModel
  SAVE Porosity, Capacity, Density

  Element => Model % CurrentElement
  Material => GetMaterial(Element)
  IF (.NOT.ASSOCIATED(Material)) THEN
    CALL FATAL('PermafrostCapacity', 'No Material found')
  END IF

  IF (FirstTime) THEN
    FirstTime = .FALSE.

    N = Model % MaxElementNodes
    ALLOCATE(Porosity(N),     & 
         Capacity(N),     & 
         Density(N),      & 
         STAT=istat)
    IF (istat /= 0) THEN
      CALL FATAL(  'USF_Permafrost', 'Effectivecapacity memory allocation error' )
    ELSE
      WRITE(Message,'(a)') 'Effectivecapacity memory allocation done'
      CALL INFO("Permafrost",Message,Level=4)
    END IF

    !-----------------------------------------------
    ! Read parameters from sif file Material section
    !-----------------------------------------------
    PermafrostModel = GetString( Material, 'Permafrost Model', Found )
    IF (.NOT. Found) THEN
      CALL FATAL('Permafrost', 'Could not find Permafrost Model')
    ENDIF

    ! Keep as is for the moment
    porscale = GetCReal( Material, 'Permafrost Porosity Depth Scale', Found )
    IF (.NOT. Found) THEN
      porscale = -1.0 ! Negative value means no depth dependence for por
    ENDIF

    L = GetCReal( Material, 'Latent Heat', Found )
    IF (.NOT. Found) THEN
      CALL FATAL('Permafrost', 'Could not find Latent Heat')
    ENDIF

    ClausiusClapeyron = GetCReal( Material, 'Clausius Clapeyron Constant', Found )
    IF (.NOT. Found) THEN
      CALL INFO('Permafrost', 'Could not find Latent Heat, using default SI value ', Level=3)
      ClausiusClapeyron = 9.8d-08 ! ### 7.42E-08 for air free water
    ENDIF

    !--- Rock parameters ---
    !rhor = GetCReal( Material, 'Permafrost Density Rock', Found )
    !IF (.NOT. Found) THEN
    !   CALL FATAL('Permafrost', 'Could not find Permafrost Density Rock')
    !ENDIF
    !Cr = GetCReal( Material, 'Permafrost Heat Capacity Rock', Found )
    !IF (.NOT. Found) THEN
    !   CALL FATAL('Permafrost', 'Could not find Permafrost Heat Capacity Rock')
    !ENDIF

    !--- Water parameters ---
    rhow = GetCReal( Material, 'Permafrost Density Water', Found )
    IF (.NOT. Found) THEN
      CALL FATAL('Permafrost', 'Could not find Permafrost Density Water')
    ENDIF
    Cw = GetCReal( Material, 'Permafrost Heat Capacity Water', Found )
    IF (.NOT. Found) THEN
      CALL FATAL('Permafrost', 'Could not find Permafrost Heat Capacity Water')
    ENDIF

    !--- Ice parameters ---
    rhoi = GetCReal( Material, 'Permafrost Density Ice', Found )
    IF (.NOT. Found) THEN
      CALL FATAL('Permafrost', 'Could not find Permafrost Density Ice')
    ENDIF
    Ci = GetCReal( Material, 'Permafrost Heat Capacity Ice', Found )
    IF (.NOT. Found) THEN
      CALL FATAL('Permafrost', 'Could not find Permafrost Heat Capacity Ice')
    ENDIF

    !--- Power Law model parameters ---
    IF (TRIM(PermafrostModel) .EQ. "power law") THEN
      a = GetCReal( Material, 'Permafrost Power Law Factor', Found )
      IF (.NOT. Found) THEN
        CALL FATAL('Permafrost', 'Could not find Permafrost Power Law Factor')
      ENDIF
      b = GetCReal( Material, 'Permafrost Power Law Exponent', Found )
      IF (.NOT. Found) THEN
        CALL FATAL('Permafrost', 'Could not find Permafrost Power Law Exponent')
      ENDIF
      dT = GetCReal( Material, 'Permafrost Power law Temperature Offset', Found )
      IF (.NOT. Found) THEN
        CALL FATAL('Permafrost', 'Could not find Permafrost Power law Temperature Offset')
      ENDIF
      !--- Exponential model parameters ---
    ELSE IF (TRIM(PermafrostModel) .EQ. "exponential") THEN
      a = GetCReal( Material, 'Permafrost Exponential Temperature Interval', Found )
      IF (.NOT. Found) THEN
        CALL FATAL('Permafrost', 'Could not find Permafrost Exponential Temperature Interval')
      ENDIF
    ELSE
      CALL FATAL('Permafrost', 'Unknown Permafrost Model')
    ENDIF
  ENDIF ! End of FirstTime

  !-----------------------------------------------
  ! Get the depth of lower layer (rock layer)
  !-----------------------------------------------
  DepthVarName = GetString( Model % Solver % Values , 'Lower Depth Name', Found ) 
  IF (.NOT.Found) THEN
    WRITE(DepthVarName,'(A)') 'lower depth'
  END IF
  DepthVar => VariableGet(Model % Mesh % Variables, TRIM(DepthVarName))
  !DepthVar => VariableGet(Model % Mesh % Variables, "lower depth")
  IF ( ASSOCIATED(DepthVar) ) THEN
    Depth = DepthVar % Values ( DepthVar % Perm(Node) )
  ELSE 
    Depth = 0.0_dp
  END IF

  !-----------------------------------------------
  ! Get the total depth below all layers
  !-----------------------------------------------
  DepthVarName = GetString( Model % Solver % Values , 'Total Depth Name', Found ) 
  IF (.NOT.Found) THEN
    WRITE(DepthVarName,'(A)') 'depth'
  END IF
  DepthVar => VariableGet(Model % Mesh % Variables, TRIM(DepthVarName))
  DepthVar2 => VariableGet(Model % Mesh % Variables, "depth")
  IF ( ASSOCIATED(DepthVar2) ) THEN
    Depth2 = DepthVar2 % Values ( DepthVar2 % Perm(Node) )
  ELSE
    Depth2 = 0.0_dp
  END IF

  ! In case there is no lower layer, use the only layer
  if (Depth == 0.0_dp) Depth = Depth2

  ! Get porosity
  N = GetElementNOFNodes(Element)
  Porosity(1:N) = ListGetReal ( Material, 'Permafrost Porosity', &
       N, Element % NodeIndexes, Found, & 
       UnfoundFatal=UnfoundFatal )
  IF (.NOT. Found) THEN
    CALL FATAL('Permafrost', 'Could not find Permafrost Porosity')
  ENDIF
  Capacity(1:N) = ListGetReal ( Material, 'Permafrost Heat Capacity Rock', &
       N, Element % NodeIndexes, Found, & 
       UnfoundFatal=UnfoundFatal )
  IF (.NOT. Found) THEN
    CALL FATAL('Permafrost', 'Could not find Permafrost Heat Capacity')
  ENDIF
  Density(1:N) = ListGetReal ( Material, 'Permafrost Density Rock', &
       N, Element % NodeIndexes, Found, & 
       UnfoundFatal=UnfoundFatal )
  IF (.NOT. Found) THEN
    CALL FATAL('Permafrost', 'Could not find Permafrost Density')
  ENDIF

  DO i = 1, N
    k = Element % NodeIndexes(i)
    If (Node .EQ. k) THEN
      por = Porosity(i)
      Cr = Capacity(i)
      rhor = Density(i)
      EXIT
    END IF
  END DO

  !WRITE(MESSAGE,*) Node, Temp, Depth, Depth2, por
  !CALL Info('Permafrost Enthalpy', MESSAGE, level=11)

  !-------------------
  ! Start calculations
  !-------------------
  IF (porscale .LE. 0.0_dp) THEN
    pordepth = por
  ELSE
    pordepth = por * EXP(-Depth/porscale)
  ENDIF

  phir = 1 - pordepth

  !----------------------
  ! Compute Tpmp at depth
  !----------------------
  iceDepth = Depth2 - Depth
  pice = iceDepth * rhoi * 9.81_dp
  prock = Depth * rhow * 9.81_dp   ! This is for the fluid so use hydrostatic pressure
  press = pice + prock
  Tpmp = 273.15_dp - ClausiusClapeyron*press 

  IF (TRIM(PermafrostModel) .EQ. "power law") THEN
    Tstar = Tpmp - Temp
    IF (Tstar <= dT) THEN
      phiw =  pordepth
    ELSE 
      phiw = (rhor * (1.0_dp - pordepth)/ rhow) * a * (Tstar)**b
    ENDIF
  ELSE IF (TRIM(PermafrostModel) .EQ. "exponential") THEN
    IF (Temp > Tpmp) then 
      fw = 1.0_dp
    ELSE
      fw = EXP(-((Temp - Tpmp)/a)**2)
    ENDIF
    phiw =  pordepth * fw
  ELSE
    CALL FATAL('Permafrost', 'Unknown Permafrost Model')
  ENDIF

  ! Check that phiw does not exceed porosity
  IF (phiw > pordepth) THEN
    phiw = pordepth
  ENDIF

  effectivecapacity = phir*Cr + (pordepth-phiw)*Ci + phiw*Cw

  ScaleSystem = ListGetLogical( Material , 'Scale System', Found )
  IF (.NOT.Found) ScaleSystem=.FALSE.
  IF (ScaleSystem)  effectivecapacity = effectivecapacity * factor
END FUNCTION PermafrostCapacity

!==============================================================================
FUNCTION PermafrostDensity(Model, Node, Temp) RESULT(Dens)
!==============================================================================

  USE DefUtils
  USE types
  USE CoordinateSystems
  USE SolverUtils
  USE ElementDescription

  IMPLICIT None

  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: Temp, Dens

  ! Local variables
  TYPE(Element_t),POINTER :: Element
  TYPE(ValueList_t), POINTER :: Material
  REAL(KIND=dp), ALLOCATABLE :: Porosity(:)
  REAL(KIND=dp), ALLOCATABLE :: Density(:)

  TYPE(Variable_t), POINTER :: DepthVar, DepthVar2
  REAL(KIND=dp) :: Depth, Depth2

  REAL(KIND=dp) :: por                 ! Porosity
  REAL(KIND=dp) :: porscale            ! Porosity length scale
  REAL(KIND=dp) :: pordepth            ! Porosity as a function of depth
  REAL(KIND=dp) :: phir, phiw          ! Volume fractions
  ! r = rock, w = water, i = ice
  REAL(KIND=dp) :: rhor, rhow, rhoi    ! Densities 
  REAL(KIND=dp) :: Tpmp                ! Melting point temperature of ice
  REAL(KIND=dp) :: ClausiusClapeyron   ! Clausius-Clapeyron constant
  REAL(KIND=dp) :: a, b, Tstar, dT     ! Params for powerlaw/exponential model
  REAL(KIND=dp) :: fw                  ! Params for exponential model
  REAL(KIND=dp) :: iceDepth, pice, prock, press
  REAL(KIND=dp), PARAMETER :: factor=1.0_dp/(1.0d06 * 31557600.0_dp**2.0_dp)
  
  INTEGER :: N, istat, i, k

  CHARACTER(LEN=MAX_NAME_LEN) :: PermafrostModel

  LOGICAL :: FirstTime = .TRUE.
  LOGICAL :: Found, UnfoundFatal, ScaleSystem

  SAVE porscale
  SAVE rhow, rhoi
  SAVE a, b, dT
  SAVE FirstTime, PermafrostModel
  SAVE Porosity, Density

  Element => Model % CurrentElement
  Material => GetMaterial(Element)

  IF (.NOT.ASSOCIATED(Material)) THEN
    CALL FATAL('Permafrost', 'No Material found')
  END IF

  IF (FirstTime) THEN
    FirstTime = .FALSE.

    !DEALLOCATE(Porosity)
    N = Model % MaxElementNodes
    ALLOCATE(Porosity(N),     & 
         Density(N),      & 
         STAT=istat)
    IF (istat /= 0) THEN
      CALL FATAL(  'USF_Permafrost', 'Density memory allocation error' )
    ELSE
      WRITE(Message,'(a)') 'Density memory allocation done'
      CALL INFO("Permafrost",Message,Level=4)
    END IF

    !-----------------------------------------------
    ! Read parameters from sif file Material section
    !-----------------------------------------------
    PermafrostModel = GetString( Material, 'Permafrost Model', Found )
    IF (.NOT. Found) THEN
      CALL FATAL('Permafrost', 'Could not find Permafrost Model')
    ENDIF

    ClausiusClapeyron = GetCReal( Material, 'Clausius Clapeyron Constant', Found )
    IF (.NOT. Found) THEN
      CALL INFO('Permafrost', 'Could not find Latent Heat, using default SI value ', Level=3)
      ClausiusClapeyron = 9.8d-08 ! ### 7.42E-08 for air free water
    ENDIF

    porscale = GetCReal( Material, 'Permafrost Porosity Depth Scale', Found )
    IF (.NOT. Found) THEN
      WRITE(MESSAGE,*) 'No depth scale for porosity. Set to zero.'
      CALL Info('Permafrost', MESSAGE, level=3)
      porscale = -1.0 ! Negative value means no depth depedence for por
    ENDIF

    !rhor = GetCReal( Material, 'Permafrost Density Rock', Found )
    !IF (.NOT. Found) THEN
    !  CALL FATAL('Permafrost', 'Could not find Permafrost Density Rock')
    !ENDIF

    rhow = GetCReal( Material, 'Permafrost Density Water', Found )
    IF (.NOT. Found) THEN
      CALL FATAL('Permafrost', 'Could not find Permafrost Density Water')
    ENDIF

    rhoi = GetCReal( Material, 'Permafrost Density Ice', Found )
    IF (.NOT. Found) THEN
      CALL FATAL('Permafrost', 'Could not find Permafrost Density Ice')
    ENDIF

    !--- Power Law model parameters ---
    IF (TRIM(PermafrostModel) .EQ. "power law") THEN
      a = GetCReal( Material, 'Permafrost Power Law Factor', Found )
      IF (.NOT. Found) THEN
        CALL FATAL('Permafrost', 'Could not find Permafrost Power Law Factor')
      ENDIF
      b = GetCReal( Material, 'Permafrost Power Law Exponent', Found )
      IF (.NOT. Found) THEN
        CALL FATAL('Permafrost', 'Could not find Permafrost Power Law Exponent')
      ENDIF
      dT = GetCReal( Material, 'Permafrost Power law Temperature Offset', Found )
      IF (.NOT. Found) THEN
        CALL FATAL('Permafrost', 'Could not find Permafrost Power law Temperature Offset')
      ENDIF
      !--- Exponential model parameters ---
    ELSE IF (TRIM(PermafrostModel) .EQ. "exponential") THEN
      a = GetCReal( Material, 'Permafrost Exponential Temperature Interval', Found )
      IF (.NOT. Found) THEN
        CALL FATAL('Permafrost', 'Could not find Permafrost Exponential Temperature Interval')
      ENDIF
    ELSE
      CALL FATAL('Permafrost', 'Unknown Permafrost Model')
    ENDIF
  ENDIF


  !-----------------------------------------------
  ! Get the depth of lower layer
  !-----------------------------------------------
  DepthVar => VariableGet(Model % Mesh % Variables, "lower depth")
  IF ( ASSOCIATED(DepthVar) ) THEN
    Depth = DepthVar % Values ( DepthVar % Perm(Node) )
  ELSE 
    Depth = 0.0_dp
  END IF

  !-----------------------------------------------
  ! Get the total depth below all layers
  !-----------------------------------------------
  DepthVar2 => VariableGet(Model % Mesh % Variables, "depth")
  IF ( ASSOCIATED(DepthVar2) ) THEN
    Depth2 = DepthVar2 % Values ( DepthVar2 % Perm(Node) )
  ELSE
    Depth2 = 0.0_dp
  END IF

  ! In case there is no lower layer, use the only layer
  if (Depth == 0.0) Depth = Depth2

  ! Get porosity
  N = GetElementNOFNodes(Element)
  Porosity(1:N) = ListGetReal ( Material, 'Permafrost Porosity', &
       N, Element % NodeIndexes, Found, & 
       UnfoundFatal=UnfoundFatal )
  IF (.NOT. Found) THEN
    CALL FATAL('Permafrost', 'Could not find Permafrost Porosity')
  ENDIF
  Density(1:N) = ListGetReal ( Material, 'Permafrost Density Rock', &
       N, Element % NodeIndexes, Found, & 
       UnfoundFatal=UnfoundFatal )
  IF (.NOT. Found) THEN
    CALL FATAL('Permafrost', 'Could not find Permafrost Density')
  ENDIF


  DO i = 1, N
    k = Element % NodeIndexes(i)
    If (Node .EQ. k) THEN
      por = Porosity(i)
      rhor = Density(i)
      EXIT
    END IF
  END DO

  !WRITE(MESSAGE,*) Node, Temp, Depth, Depth2, por
  !CALL Info('Permafrost Density', MESSAGE, level=11)

  !-------------------
  ! Start calculations
  !-------------------
  IF (porscale .LE. 0.0_dp) THEN
    pordepth = por
  ELSE
    pordepth = por * EXP(-Depth/porscale)
  ENDIF

  phir = 1.0_dp - pordepth

  !----------------------
  ! Compute Tpmp at depth
  !----------------------
  iceDepth = Depth2 - Depth
  pice = iceDepth * rhoi * 9.81_dp
  prock = Depth * rhow * 9.81_dp   ! This is for the fluid so use hydrostatic pressure
  press = pice + prock
  Tpmp = 273.15_dp - ClausiusClapeyron*press ! ### 7.42E-08 for air free water - lets move to read in value with 9.8E-08 as default

  IF (TRIM(PermafrostModel) .EQ. "power law") THEN
    Tstar = Tpmp - Temp
    IF (Tstar <= dT) THEN
      phiw =  pordepth
    ELSE ! equation (6) in Cutler 2000
      phiw = (rhor * (1.0_dp - pordepth)/ rhow) * a * (Tstar)**b
    ENDIF

  ELSE IF (TRIM(PermafrostModel) .EQ. "exponential") THEN
    IF (Temp > Tpmp) then
      fw = 1.0_dp
    ELSE
      fw = EXP(-((Temp - Tpmp)/a)**2.0_dp)
    ENDIF
    phiw =  pordepth * fw
  ELSE
    CALL FATAL('Permafrost', 'Unknown Permafrost Model')

  ENDIF

  ! Check that phiw does not exceed porosity
  IF (phiw > pordepth) THEN
    phiw = pordepth
  ENDIF
  
  Dens = phir*rhor + (pordepth-phiw)*rhoi + phiw*rhow

    Material => GetMaterial(Model % CurrentElement)
  IF (.NOT.ASSOCIATED(Material)) THEN
    CALL FATAL('Permafrost', 'No Material found')
  END IF
  ScaleSystem = ListGetLogical( Material , 'Scale System', Found )
  IF (.NOT.Found) ScaleSystem=.FALSE.
  IF (ScaleSystem) Dens = Dens * factor

END FUNCTION PermafrostDensity

!==============================================================================
FUNCTION PermafrostConductivity(Model, Node, Temp) RESULT(Cond)
!==============================================================================

  USE DefUtils
  USE types
  USE CoordinateSystems
  USE SolverUtils
  USE ElementDescription

  IMPLICIT None

  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: Temp, Cond

  ! Local variables
  TYPE(Element_t),POINTER :: Element
  TYPE(ValueList_t), POINTER :: Material
  REAL(KIND=dp), ALLOCATABLE :: Porosity(:)
  REAL(KIND=dp), ALLOCATABLE :: Conductivity(:)
  REAL(KIND=dp), ALLOCATABLE :: Density(:)

  TYPE(Variable_t), POINTER :: DepthVar, DepthVar2
  REAL(KIND=dp) :: Depth, Depth2

  REAL(KIND=dp) :: por                 ! Porosity
  REAL(KIND=dp) :: porscale            ! Porosity length scale
  REAL(KIND=dp) :: pordepth            ! Porosity as a function of depth
  ! r = rock, w = water, i = ice
  REAL(KIND=dp) :: phir, phiw          ! Volume fractions
  REAL(KIND=dp) :: rhoi, rhor, rhow    ! Densities 
  REAL(KIND=dp) :: Kr, Kw, Ki          ! Heat conductivities
  REAL(KIND=dp) :: L                   ! Latent heat
  REAL(KIND=dp) :: ClausiusClapeyron   ! Clausius-Clapeyron constant
  REAL(KIND=dp) :: Tpmp                ! Melting point temperature of ice
  REAL(KIND=dp) :: a, b, Tstar, dT     ! Params for powerlaw/exponential model
  REAL(KIND=dp) :: fw                  ! Params for exponential model
  REAL(KIND=dp) :: iceDepth, pice, prock, press
  REAL(KIND=dp), PARAMETER ::  factor = 31.5576000_dp
  
  INTEGER :: N, istat, i, k

  CHARACTER(LEN=MAX_NAME_LEN) :: PermafrostModel

  LOGICAL :: FirstTime = .TRUE.
  LOGICAL :: Found, UnfoundFatal, ScaleSystem

  SAVE porscale
  SAVE rhoi, rhow
  SAVE Kw, Ki, L
  SAVE a, b, dT
  SAVE FirstTime, PermafrostModel
  SAVE Porosity, Conductivity, Density

  Element => Model % CurrentElement
  Material => GetMaterial(Element)
  IF (.NOT.ASSOCIATED(Material)) THEN
    CALL FATAL('Permafrost', 'No Material found')
  END IF

  IF (FirstTime) THEN
    FirstTime = .FALSE.

    N = Model % MaxElementNodes
    ALLOCATE(Porosity(N),     & 
         Conductivity(N), & 
         Density(N),      & 
         STAT=istat)
    IF (istat /= 0) THEN
      CALL FATAL(  'USF_Permafrost', 'Conductivity memory allocation error' )
    ELSE
      WRITE(Message,'(a)') 'Conductivity memory allocation done'
      CALL INFO("Permafrost", Message, Level=4)
    END IF

    !-----------------------------------------------
    ! Read parameters from sif file Material section
    !-----------------------------------------------
    PermafrostModel = GetString( Material, 'Permafrost Model', Found )
    IF (.NOT. Found) THEN
      CALL FATAL('Permafrost', 'Could not find Permafrost Model')
    ENDIF

    ClausiusClapeyron = GetCReal( Material, 'Clausius Clapeyron Constant', Found )
    IF (.NOT. Found) THEN
      CALL INFO('Permafrost', 'Could not find Latent Heat, using default SI value ', Level=3)
      ClausiusClapeyron = 9.8d-08 ! ### 7.42E-08 for air free water
    ENDIF

    porscale = GetCReal( Material, 'Permafrost Porosity Depth Scale', Found )
    IF (.NOT. Found) THEN
      porscale = -1.0 ! Negative value means no depth depedence for por
    ENDIF

    rhow = GetCReal( Material, 'Permafrost Density Water', Found )
    IF (.NOT. Found) THEN
      CALL FATAL('Permafrost', 'Could not find Permafrost Density Water')
    ENDIF

    rhoi = GetCReal( Material, 'Permafrost Density Ice', Found )
    IF (.NOT. Found) THEN
      CALL FATAL('Permafrost', 'Could not find Permafrost Density Ice')
    ENDIF

    Kw = GetCReal( Material, 'Permafrost Heat Conductivity Water', Found )
    IF (.NOT. Found) THEN
      CALL FATAL('Permafrost', 'Could not find Permafrost Heat Conductivity Water')
    ENDIF

    Ki = GetCReal( Material, 'Permafrost Heat Conductivity Ice', Found )
    IF (.NOT. Found) THEN
      CALL FATAL('Permafrost', 'Could not find Permafrost Heat Conductivity Ice')
    ENDIF

    !--- Power Law model parameters ---
    IF (TRIM(PermafrostModel) .EQ. "power law") THEN
      a = GetCReal( Material, 'Permafrost Power Law Factor', Found )
      IF (.NOT. Found) THEN
        CALL FATAL('Permafrost', 'Could not find Permafrost Power Law Factor')
      ENDIF
      b = GetCReal( Material, 'Permafrost Power Law Exponent', Found )
      IF (.NOT. Found) THEN
        CALL FATAL('Permafrost', 'Could not find Permafrost Power Law Exponent')
      ENDIF
      dT = GetCReal( Material, 'Permafrost Power law Temperature Offset', Found )
      IF (.NOT. Found) THEN
        CALL FATAL('Permafrost', 'Could not find Permafrost Power law Temperature Offset')
      ENDIF
      !--- Exponential model parameters ---
    ELSE IF (TRIM(PermafrostModel) .EQ. "exponential") THEN
      a = GetCReal( Material, 'Permafrost Exponential Temperature Interval', Found )
      IF (.NOT. Found) THEN
        CALL FATAL('Permafrost', 'Could not find Permafrost Exponential Temperature Interval')
      ENDIF
    ELSE
      CALL FATAL('Permafrost', 'Unknown Permafrost Model')
    ENDIF
  ENDIF

  !-----------------------------------------------
  ! Get the depth of lower layer
  !-----------------------------------------------
  DepthVar => VariableGet(Model % Mesh % Variables, "lower depth")
  IF ( ASSOCIATED(DepthVar) ) THEN
    Depth = DepthVar % Values ( DepthVar % Perm(Node) )
  ELSE 
    Depth = 0.0_dp
  END IF

  !-----------------------------------------------
  ! Get the total depth below all layers
  !-----------------------------------------------
  DepthVar2 => VariableGet(Model % Mesh % Variables, "depth")
  IF ( ASSOCIATED(DepthVar2) ) THEN
    Depth2 = DepthVar2 % Values ( DepthVar2 % Perm(Node) )
  ELSE
    Depth2 = 0.0_dp
  END IF

  !-----------------------------------------------
  ! In case there is no lower layer, use the only layer
  !-----------------------------------------------
  if (Depth == 0.0_dp) Depth = Depth2

  !----------------------------------------------------------------------------
  ! Get the porosity
  !----------------------------------------------------------------------------
  N = GetElementNOFNodes(Element)
  Porosity(1:N) = ListGetReal ( Material, 'Permafrost Porosity', &
       N, Element % NodeIndexes, Found, &
       UnfoundFatal=UnfoundFatal )
  IF (.NOT. Found) THEN
    CALL FATAL('Permafrost', 'Could not find Permafrost Porosity')
  ENDIF
  Conductivity(1:N) = ListGetReal ( Material, 'Permafrost Heat Conductivity Rock', &
       N, Element % NodeIndexes, Found, & 
       UnfoundFatal=UnfoundFatal )
  IF (.NOT. Found) THEN
    CALL FATAL('Permafrost', 'Could not find Permafrost Heat Conductivity')
  ENDIF
  Density(1:N) = ListGetReal ( Material, 'Permafrost Density Rock', &
       N, Element % NodeIndexes, Found, & 
       UnfoundFatal=UnfoundFatal )
  IF (.NOT. Found) THEN
    CALL FATAL('Permafrost', 'Could not find Permafrost Density')
  ENDIF


  DO i = 1, N
    k = Element % NodeIndexes(i)
    If (Node .EQ. k) THEN
      por = Porosity(i)
      Kr = Conductivity(i)
      rhor = Density(i)
      EXIT
    END IF
  END DO

  !-------------------
  ! Start calculations
  !-------------------
  IF (porscale .LE. 0.0_dp) THEN
    pordepth = por
  ELSE
    pordepth = por * EXP(-Depth/porscale)
  ENDIF

  phir = 1 - pordepth

  !----------------------
  ! Compute Tpmp at depth
  !----------------------
  iceDepth = Depth2 - Depth
  pice = iceDepth * rhoi * 9.81_dp
  prock = Depth * rhow * 9.81_dp   ! This is for the fluid so use hydrostatic pressure
  press = pice + prock
  Tpmp = 273.15_dp - ClausiusClapeyron*press

  IF (TRIM(PermafrostModel) .EQ. "power law") THEN
    Tstar = Tpmp - Temp
    IF (Tstar <= dT) THEN
      phiw =  pordepth
    ELSE 
      phiw = (rhor * (1.0_dp - pordepth)/ rhow) * a * (Tstar)**b
    ENDIF
  ELSE IF (TRIM(PermafrostModel) .EQ. "exponential") THEN
    IF (Temp > Tpmp) then
      fw = 1.0_dp
    ELSE
      fw = EXP(-((Temp - Tpmp)/a)**2.0_dp)
    ENDIF
    phiw =  pordepth * fw
  ELSE
    CALL FATAL('Permafrost', 'Unknown Permafrost Model')
  ENDIF

  ! Check that phiw does not exceed porosity
  IF (phiw > pordepth) THEN
    phiw = pordepth
  ENDIF

  Cond = Kr**phir * Ki**(pordepth-phiw) * Kw**phiw
  
  Material => GetMaterial(Model % CurrentElement)
  IF (.NOT.ASSOCIATED(Material)) THEN
    CALL FATAL('Permafrost', 'No Material found')
  END IF
  ScaleSystem = ListGetLogical( Material , 'Scale System', Found )
  IF (.NOT.Found) ScaleSystem=.FALSE.
  IF (ScaleSystem) THEN
    CALL INFO('PermafrostConductivity','Applying MPa-m-a scaling',Level=9)
    Cond = Cond * factor
  END IF
  !write (*,*) Depth, Kr, Kw, phir, phiw, pordepth, Kr**phir * Ki**(pordepth-phiw) * Kw**phiw

END FUNCTION PermafrostConductivity

!==============================================================================
FUNCTION PermafrostPressure(Model, Node, dumm) RESULT(pressure)
!==============================================================================

  USE DefUtils
  USE types
  USE CoordinateSystems
  USE SolverUtils
  USE ElementDescription

  IMPLICIT None

  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: dumm, pressure

  ! Local variables
  TYPE(Element_t),POINTER :: Element
  TYPE(ValueList_t), POINTER :: Material

  TYPE(Variable_t), POINTER :: DepthVar, DepthVar2
  REAL(KIND=dp) :: Depth, Depth2
  REAL(KIND=dp), PARAMETER :: factor = 1.0d-06

  REAL(KIND=dp) :: rhor, rhow, rhoi     ! Densities 
  REAL(KIND=dp) :: iceDepth, pice, prock, press ! For computing pressures

  LOGICAL :: FirstTime = .TRUE.
  LOGICAL :: Found, ScaleSystem

  SAVE rhor, rhow, rhoi
  SAVE FirstTime

  Element => Model % CurrentElement
  Material => GetMaterial(Element)
  IF (FirstTime) THEN
    FirstTime = .FALSE.

    !--- Rock parameters ---
    !rhor = GetCReal( Material, 'Permafrost Density Rock', Found )
    !IF (.NOT. Found) THEN
    !  CALL FATAL('Permafrost', 'Could not find Permafrost Density Rock')
    !ENDIF

    !--- Water parameters ---
    rhow = GetCReal( Material, 'Permafrost Density Water', Found )
    IF (.NOT. Found) THEN
       CALL FATAL('Permafrost', 'Could not find Permafrost Density Water')
    ENDIF

    !--- Ice parameters ---
    rhoi = GetCReal( Material, 'Permafrost Density Ice', Found )
    IF (.NOT. Found) THEN
      CALL FATAL('Permafrost', 'Could not find Permafrost Density Ice')
    ENDIF

  ENDIF ! End of FirstTime

  !-----------------------------------------------
  ! Get the depth of lower layer (rock layer)
  !-----------------------------------------------
  DepthVar => VariableGet(Model % Mesh % Variables, "lower depth")
  IF ( ASSOCIATED(DepthVar) ) THEN
    Depth = DepthVar % Values ( DepthVar % Perm(Node) )
  ELSE 
    Depth = 0.0_dp
  END IF

  !-----------------------------------------------
  ! Get the total depth below all layers
  !-----------------------------------------------
  DepthVar2 => VariableGet(Model % Mesh % Variables, "depth")
  IF ( ASSOCIATED(DepthVar2) ) THEN
    Depth2 = DepthVar2 % Values ( DepthVar2 % Perm(Node) )
  ELSE
    Depth2 = 0.0_dp
  END IF

  ! In case there is no lower layer, use the only layer
  if (Depth == 0.0_dp) Depth = Depth2

  !----------------------
  ! Compute fluid pressure at depth
  !----------------------
  iceDepth = Depth2 - Depth
  pice = iceDepth * rhoi * 9.81_dp
  prock = Depth * rhow * 9.81_dp   ! This is for the fluid so use hydrostatic pressure
  pressure = pice + prock
  
  Material => GetMaterial(Model % CurrentElement)
  IF (.NOT.ASSOCIATED(Material)) THEN
    CALL FATAL('Permafrost', 'No Material found')
  END IF
  ScaleSystem = ListGetLogical( Material , 'Scale System', Found )
  IF (.NOT.Found) ScaleSystem=.FALSE.
  IF (ScaleSystem) THEN
    CALL INFO('PermafrostPressure','Applying Mpa-m-a scaling',Level=9)
    pressure = pressure * factor
  END IF

END FUNCTION PermafrostPressure

!==============================================================================
FUNCTION PermafrostIceConductivity(Model, Node, temp) RESULT(cond)
!==============================================================================

  USE DefUtils
  USE IceProperties

  IMPLICIT None

  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: temp, cond

  ! Local variables
  TYPE(Variable_t), POINTER :: HeightVar, HeightVar2, TotalHeightVar
  TYPE(ValueList_t), POINTER :: Material
  REAL(KIND=dp) :: Height, TotalHeight, MinimumIceDepth, factor 
  REAL(KIND=dp), PARAMETER :: DHeight = 50.0000000_dp 
  REAL(KIND=dp) :: c0, c1, sigmoid

  CHARACTER(LEN=MAX_NAME_LEN) :: HeightVarName, Height2VarName
  CHARACTER(LEN=MAX_NAME_LEN) :: TotalHeightVarName

  LOGICAL :: Found, Found2, ScaleSystem=.FALSE.

  !---------------------------------------------------------------------------
  ! Get the height of ice layer 
  ! Default is either upper layer or the only layer
  ! Optional is to give 'Lower Depth Name" and "Total Depth Name"
  !---------------------------------------------------------------------------
  Height = 100.0_dp ! Assume there is ice
  HeightVarName = GetString( Model % Solver % Values , 'Lower Depth Name', Found )
  IF (.NOT.Found) THEN
    WRITE(HeightVarName,'(A)') 'max upper depth'
    WRITE(Height2VarName,'(A)') 'max depth'
  ELSE
    WRITE(Height2VarName,'(A)') 'Depth'
  END IF

  !HeightVar => VariableGet(Model % Mesh % Variables, "max upper depth")
  HeightVar => VariableGet(Model % Mesh % Variables, TRIM(HeightVarName))
  IF ( ASSOCIATED(HeightVar) ) THEN
    Height = HeightVar % Values ( HeightVar % Perm(Node) )
  ELSE
    !HeightVar2 => VariableGet(Model % Mesh % Variables, "max depth")
    HeightVar2 => VariableGet(Model % Mesh % Variables, TRIM(Height2VarName))
    IF ( ASSOCIATED(HeightVar2) ) THEN
      Height = HeightVar2 % Values ( HeightVar2 % Perm(Node) )
    ELSE
      CALL FATAL(' PermafrostIceConductivity', 'Could not find depth or upper depth')
    END IF
  END IF

  IF (Found) THEN
    TotalHeightVarName = GetString( Model % Solver % Values , 'Total Depth Name', Found2 )
    IF (.NOT.Found2) &
         CALL FATAL('PermafrostIceConductivity', 'Could not find Total Depth Name')
    TotalHeightVar => VariableGet(Model % Mesh % Variables, TRIM(TotalHeightVarName))
    TotalHeight = TotalHeightVar % Values ( TotalHeightVar % Perm(Node) )
    Height = TotalHeight - Height
  ENDIF

  Material => GetMaterial(Model % CurrentElement)
  IF (.NOT.ASSOCIATED(Material)) THEN
    CALL FATAL('Permafrost', 'No Material found')
  END IF

  MinimumIceDepth = GetConstReal(Material, 'Minimum Height', Found)
  IF (.NOT.Found) THEN
    CALL FATAL('PermafrostIceConductivity','No variable >Minimum Height< found in Material')
  END IF
  
  IF (Height > MinimumIceDepth + DHeight) THEN  ! This is ice + DHeight
    cond = GetIceConductivity(temp)!IceConductivity(Model,temp)
  ELSE                       ! A very conductive layer
    c0 = GetIceConductivity(temp) !IceConductivity(Model,temp)
    c1 = c0*5
    sigmoid = 1.0/(1.0 + EXP( (-Height + (MinimumIceDepth+DHeight)/2.0)/3.0 ) )
    cond = c0 + (1.0 - sigmoid) * (c1-c0)
  ENDIF
  ScaleSystem = ListGetLogical( Material , 'Scale System', Found )
  IF (ScaleSystem) THEN
    factor = 31.5576000_dp
    CALL INFO('PermafrostIceConductivity','Applying Mpa-m-a scaling',Level=9)
  ELSE
    factor = GetConstReal(Material,"Heat Capacity Scaling Factor",ScaleSystem)
  ENDIF
  IF (ScaleSystem) THEN    
    cond = cond * factor
  END IF
  
END FUNCTION  PermafrostIceConductivity


!==============================================================================
FUNCTION  PermafrostIceCapacity(Model, Node, temp) RESULT(capac)
!==============================================================================

  USE DefUtils
  USE IceProperties


  IMPLICIT None

  TYPE(Model_t) :: Model
  INTEGER :: Node
  REAL(KIND=dp) :: temp, capac

  ! Local variables
  TYPE(Variable_t), POINTER :: HeightVar, HeightVar2, TotalHeightVar
  TYPE(ValueList_t), POINTER :: Material
  REAL(KIND=dp) :: Height, TotalHeight, MinimumIceDepth, factor 
  REAL(KIND=dp), PARAMETER :: DHeight = 50.0000000_dp 
  REAL(KIND=dp) :: c0, c1, sigmoid

  CHARACTER(LEN=MAX_NAME_LEN) :: HeightVarName, Height2VarName
  CHARACTER(LEN=MAX_NAME_LEN) :: TotalHeightVarName

  LOGICAL :: Found, Found2, ScaleSystem=.FALSE.

  !---------------------------------------------------------------------------
  ! Get the height of ice layer 
  ! Default is either upper layer or the only layer
  ! Optional is to give 'Lower Depth Name" and "Total Depth Name"
  !---------------------------------------------------------------------------
  Height = 100 ! Assume there is ice
  HeightVarName = GetString( Model % Solver % Values , 'Lower Depth Name', Found )
  IF (.NOT.Found) THEN
    WRITE(HeightVarName,'(A)') 'max upper depth'
    WRITE(Height2VarName,'(A)') 'max depth'
  ELSE
    WRITE(Height2VarName,'(A)') 'Depth'
  END IF

  !HeightVar => VariableGet(Model % Mesh % Variables, "max upper depth")
  HeightVar => VariableGet(Model % Mesh % Variables, TRIM(HeightVarName))
  IF ( ASSOCIATED(HeightVar) ) THEN
    Height = HeightVar % Values ( HeightVar % Perm(Node) )
  ELSE
    !HeightVar2 => VariableGet(Model % Mesh % Variables, "max depth")
    HeightVar2 => VariableGet(Model % Mesh % Variables, TRIM(Height2VarName))
    IF ( ASSOCIATED(HeightVar2) ) THEN
      Height = HeightVar2 % Values ( HeightVar2 % Perm(Node) )
    ELSE
      CALL FATAL('PermafrostIceCapacity', 'Could not find depth or upper depth')
    END IF
  END IF

  IF (Found) THEN
    TotalHeightVarName = GetString( Model % Solver % Values , 'Total Depth Name', Found2 )
    IF (.NOT. Found2) &
         CALL FATAL('PermafrostIceCapacity', 'Could not find Total Depth Name')
    TotalHeightVar => VariableGet(Model % Mesh % Variables, TRIM(TotalHeightVarName))
    TotalHeight = TotalHeightVar % Values ( TotalHeightVar % Perm(Node) )
    Height = TotalHeight - Height
  ENDIF

  Material => GetMaterial(Model % CurrentElement)
  IF (.NOT.ASSOCIATED(Material)) THEN
    CALL FATAL('Permafrost', 'No Material found')
  END IF
  
  MinimumIceDepth = GetConstReal(Material, 'Minimum Height', Found)
  IF (.NOT.Found) THEN
    CALL FATAL('PermafrostIceConductivity','No variable >Minimum Height< found in Material')
  END IF
  If (Height > MinimumIceDepth + DHeight) THEN ! This is ice + DHeight
    capac = GetIceCapacity(temp)!IceCapacity(Model,temp)
  ELSE                      ! A low heat capacity layer
    !capac = 200.0_dp
    c0 = GetIceCapacity(temp)
    c1 = c0*5
    sigmoid = 1.0/(1.0 + EXP( (-Height + (MinimumIceDepth+DHeight)/2.0)/3.0 ) )
    capac = c0 + (1.0 - sigmoid) * (c1-c0)
  ENDIF
 
  ScaleSystem = ListGetLogical( Material , 'Scale System', Found )
  IF (ScaleSystem) THEN
    factor = 31557600.0_dp**2.0_dp
    CALL INFO('PermafrostIceCapacity','Applying Mpa-m-a scaling',Level=9)
  ELSE
    factor = GetConstReal(Material,"Heat Capacity Scaling Factor",ScaleSystem)
  ENDIF
  IF (ScaleSystem) THEN
    capac = capac * factor
  END IF
END FUNCTION PermafrostIceCapacity




