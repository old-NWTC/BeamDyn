MODULE BeamDyn

   USE BeamDyn_Types
   USE NWTC_Library

   IMPLICIT NONE

   PRIVATE

   TYPE(ProgDesc), PARAMETER:: BeamDyn_Ver = ProgDesc('BeamDyn_RK4', 'v1.00.00','12-March-2014')

   ! ..... Public Subroutines....................................................................

   PUBLIC :: BeamDyn_Init                           ! Initialization routine
   PUBLIC :: BeamDyn_End                            ! Ending routine (includes clean up)

   PUBLIC :: BeamDyn_UpdateStates                   ! Loose coupling routine for solving for constraint states, integrating
                                                 !   continuous states, and updating discrete states
   PUBLIC :: BeamDyn_CalcOutput                     ! Routine for computing outputs

   PUBLIC :: BeamDyn_CalcConstrStateResidual        ! Tight coupling routine for returning the constraint state residual
   PUBLIC :: BeamDyn_CalcContStateDeriv             ! Tight coupling routine for computing derivatives of continuous states
   PUBLIC :: BeamDyn_UpdateDiscState                ! Tight coupling routine for updating discrete states

CONTAINS

INCLUDE 'NodeLoc.f90'
INCLUDE 'Tilde.f90'
INCLUDE 'CrvMatrixR.f90'
INCLUDE 'CrvMatrixH.f90'
INCLUDE 'CrvCompose.f90'
INCLUDE 'ElemNodalDispGL.f90'
INCLUDE 'ElemNodalStifGL.f90'
INCLUDE 'ElemNodalMassGL.f90'
INCLUDE 'NodalRelRotGL.f90'
INCLUDE 'BldGaussPointWeight.f90'
INCLUDE 'diffmtc.f90'
INCLUDE 'BldComputeJacobianLSGL.f90'
INCLUDE 'BldGaussPointDataAt0.f90'
INCLUDE 'BldGaussPointData.f90'
INCLUDE 'ElasticForce.f90'
INCLUDE 'BldGaussPointDataMass.f90'
INCLUDE 'MassMatrix.f90'
INCLUDE 'GyroForce.f90'
INCLUDE 'ElementMatrix.f90'
INCLUDE 'AssembleStiffKGL.f90'
INCLUDE 'AssembleRHSGL.f90'
INCLUDE 'GenerateDynamicElement.f90'
INCLUDE 'ludcmp.f90'
INCLUDE 'lubksb.f90'
INCLUDE 'AppliedNodalLoad.f90'
INCLUDE 'PrescribedRootMotion.f90'
INCLUDE 'DynamicSolution.f90'
INCLUDE 'BeamDyn_RK4.f90'
INCLUDE 'CrvMatrixHinv.f90'
INCLUDE 'ComputeUDN.f90'
INCLUDE 'BeamDyn_CalcContStateDeriv.f90'

INCLUDE 'ReadPrimaryFile.f90'
INCLUDE 'ReadBladeFile.f90'
INCLUDE 'BeamDyn_ReadInput.f90'
INCLUDE 'MemberArcLength.f90'
INCLUDE 'BldComputeMemberLength.f90'
INCLUDE 'ComputeIniNodalPosition.f90'
INCLUDE 'Norm.f90'
INCLUDE 'CrossProduct.f90'
INCLUDE 'CrvExtractCrv.f90'
INCLUDE 'ComputeIniNodalCrv.f90'


   SUBROUTINE BeamDyn_Init( InitInp, u, p, x, xd, z, OtherState, y, Interval, InitOut, ErrStat, ErrMsg )
!
! This routine is called at the start of the simulation to perform initialization steps.
! The parameters are set here and not changed during the simulation.
! The initial states and initial guess for the input are defined.
!.....................................................................................................................

   TYPE(BD_InitInputType),       INTENT(IN   )  :: InitInp     ! Input data for initialization routine
   TYPE(BD_InputType),           INTENT(  OUT)  :: u           ! An initial guess for the input; input mesh must be defined
   TYPE(BD_ParameterType),       INTENT(  OUT)  :: p           ! Parameters
   TYPE(BD_ContinuousStateType), INTENT(  OUT)  :: x           ! Initial continuous states
   TYPE(BD_DiscreteStateType),   INTENT(  OUT)  :: xd          ! Initial discrete states
   TYPE(BD_ConstraintStateType), INTENT(  OUT)  :: z           ! Initial guess of the constraint states
   TYPE(BD_OtherStateType),      INTENT(  OUT)  :: OtherState  ! Initial other/optimization states
   TYPE(BD_OutputType),          INTENT(  OUT)  :: y           ! Initial system outputs (outputs are not calculated;
                                                                    !    only the output mesh is initialized)
   REAL(DbKi),                        INTENT(INOUT)  :: Interval    ! Coupling interval in seconds: the rate that
                                                                    !   (1) Mod1_UpdateStates() is called in loose coupling &
                                                                    !   (2) Mod1_UpdateDiscState() is called in tight coupling.
                                                                    !   Input is the suggested time from the glue code;
                                                                    !   Output is the actual coupling interval that will be used
                                                                    !   by the glue code.
   TYPE(BD_InitOutputType),           INTENT(  OUT)  :: InitOut     ! Output for initialization routine
   INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   ! local variables
   TYPE(BD_InputFile)      :: InputFileData    ! Data stored in the module's input file
   INTEGER(IntKi)          :: i                ! do-loop counter
   INTEGER(IntKi)          :: j                ! do-loop counter
   INTEGER(IntKi)          :: k                ! do-loop counter
   INTEGER(IntKi)          :: temp_int
   INTEGER(IntKi)          :: temp_id
   INTEGER(IntKi)          :: temp_id2
   REAL(ReKi)              :: temp_EP1(3)
   REAL(ReKi)              :: temp_EP2(3)
   REAL(ReKi)              :: temp_MID(3)
   REAL(ReKi)              :: temp_twist(3)
   REAL(ReKi)              :: temp_phi
   REAL(ReKi)              :: temp_POS(3)
   REAL(ReKi)              :: temp_CRV(3)
   REAL(ReKi),ALLOCATABLE  :: temp_GLL(:)
   REAL(ReKi),ALLOCATABLE  :: temp_w(:)
   REAL(ReKi),ALLOCATABLE  :: temp_ratio(:,:)
   INTEGER(IntKi)               :: ErrStat2                     ! Temporary Error status
   CHARACTER(LEN(ErrMsg))       :: ErrMsg2                      ! Temporary Error message

   Real(ReKi)              :: xl               ! left most point
   Real(ReKi)              :: xr               ! right most point
   REAL(ReKi)              :: blength          !beam length: xr - xl
   REAL(ReKi)              :: elem_length          !beam length: xr - xl
   REAL(ReKi),ALLOCATABLE  :: dloc(:)

   REAL(ReKi)             :: TmpPos(3)

  ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = "" 

   xl = 0.0D0
   xr = 10.0D0  ! QW
   blength = xr - xl

   ! Initialize the NWTC Subroutine Library

   CALL NWTC_Init( )

   ! Display the module information

   CALL DispNVD( BeamDyn_Ver )

   CALL BeamDyn_ReadInput(InitInp%InputFile,InputFileData,ErrStat,ErrMsg)

   CALL AllocAry(p%member_length,InputFileData%member_total,2,'member length array',ErrStat2,ErrMsg2)
   p%member_length(:,:) = 0.0D0
   p%blade_length = 0.0D0

   CALL BldComputeMemberLength(InputFileData%member_total,InputFileData%kp_coordinate,p%member_length,p%blade_length)
   p%elem_total = InputFileData%member_total
   p%node_elem  = InputFileData%order_elem + 1       ! node per element
   p%ngp        = p%node_elem - 1
   p%dof_node   = 6
   temp_int     = p%node_elem * p%dof_node

   CALL AllocAry(p%uuN0,temp_int,p%elem_total,'uuN0 (initial position) array',ErrStat2,ErrMsg2)
   p%uuN0(:,:) = 0.0D0

   CALL AllocAry(temp_GLL,p%node_elem,'GLL points array',ErrStat2,ErrMsg2)
   temp_GLL(:) = 0.0D0
   CALL AllocAry(temp_w,p%node_elem,'GLL weight array',ErrStat2,ErrMsg2)
   temp_w(:) = 0.0D0
   CALL BD_gen_gll_LSGL(p%node_elem-1,temp_GLL,temp_w)
   DO i=1,InputFileData%member_total
       temp_id = (i-1)*2
       temp_EP1(1:3) = InputFileData%kp_coordinate(temp_id+1,1:3)
       temp_MID(1:3) = InputFileData%kp_coordinate(temp_id+2,1:3)
       temp_EP2(1:3) = InputFileData%kp_coordinate(temp_id+3,1:3)
       temp_twist(1) = InputFileData%initial_twist(temp_id+1)
       temp_twist(2) = InputFileData%initial_twist(temp_id+2)
       temp_twist(3) = InputFileData%initial_twist(temp_id+3)
       DO j=1,p%node_elem
           CALL ComputeIniNodalPosition(temp_EP1,temp_EP2,temp_MID,temp_GLL(j),temp_POS)
           temp_phi = temp_twist(1) + (temp_twist(3)-temp_twist(1))*(temp_GLL(j)+1.0D0)/2.0D0
           CALL ComputeIniNodalCrv(temp_EP1,temp_EP2,temp_MID,temp_phi,temp_GLL(j),temp_CRV)
           temp_id2 = (j-1)*p%dof_node 
           p%uuN0(temp_id2+1,i) = temp_POS(1)
           p%uuN0(temp_id2+2,i) = temp_POS(2)
           p%uuN0(temp_id2+3,i) = temp_POS(3)
           p%uuN0(temp_id2+4,i) = temp_CRV(1)
           p%uuN0(temp_id2+5,i) = temp_CRV(2)
           p%uuN0(temp_id2+6,i) = temp_CRV(3)
       ENDDO
   ENDDO
   DEALLOCATE(temp_GLL)
   DEALLOCATE(temp_w)

   CALL AllocAry(temp_ratio,p%ngp,p%elem_total,'temp_ratio',ErrStat2,ErrMsg2) 
   temp_ratio(:,:) = 0.0D0
   CALL AllocAry(temp_GLL,p%ngp,'temp_GL',ErrStat2,ErrMsg2) 
   temp_GLL(:) = 0.0D0
   CALL AllocAry(temp_w,p%ngp,'temp_weight_GL',ErrStat2,ErrMsg2) 
   temp_w(:) = 0.0D0

   CALL BldGaussPointWeight(p%ngp,temp_GLL,temp_w)

   DO i=1,p%ngp
       temp_GLL(i) = (temp_GLL(i) + 1.0D0)/2.0D0
   ENDDO

   DO i=1,p%elem_total
       IF(i .EQ. 1) THEN
           DO j=1,p%ngp
               temp_ratio(j,i) = temp_GLL(j)*p%member_length(i,2)
           ENDDO
       ELSE
           DO j=1,i-1
               temp_ratio(:,i) = temp_ratio(:,i) + p%member_length(j,2)
           ENDDO
           DO j=1,p%ngp
               temp_ratio(j,i) = temp_ratio(j,i) + temp_GLL(j)*p%member_length(i,2)
           ENDDO
       ENDIF
   ENDDO

   CALL AllocAry(p%Stif0_GL,6,6,p%ngp*p%elem_total,'Stif0_GL',ErrStat2,ErrMsg2) 
   p%Stif0_GL(:,:,:) = 0.0D0
   CALL AllocAry(p%Mass0_GL,6,6,p%ngp*p%elem_total,'Mass0_GL',ErrStat2,ErrMsg2) 
   p%Mass0_GL(:,:,:) = 0.0D0
   DO i=1,p%elem_total
       DO j=1,p%node_elem-1
           temp_id = (i-1)*p%ngp+j
           DO k=1,InputFileData%InpBl%station_total
               IF(temp_ratio(j,i) <= InputFileData%InpBl%station_eta(k)) THEN
                   IF(temp_ratio(j,i) == InputFileData%InpBl%station_eta(k)) THEN
                       p%Stif0_GL(1:6,1:6,temp_id) = InputFileData%InpBl%stiff0(1:6,1:6,k)
                       p%Mass0_GL(1:6,1:6,temp_id) = InputFileData%InpBl%mass0(1:6,1:6,k)
                   ELSE 
                       p%Stif0_GL(1:6,1:6,temp_id) = 0.5D0*(InputFileData%InpBl%stiff0(1:6,1:6,k-1)+&
                                                            InputFileData%InpBl%stiff0(1:6,1:6,k))
                       p%Mass0_GL(1:6,1:6,temp_id) = 0.5D0*(InputFileData%InpBl%mass0(1:6,1:6,k-1)+&
                                                            InputFileData%InpBl%mass0(1:6,1:6,k))
                   ENDIF
                   EXIT
               ENDIF
           ENDDO
       ENDDO
   ENDDO

   DEALLOCATE(temp_GLL)
   DEALLOCATE(temp_w)
   DEALLOCATE(temp_ratio)

   WRITE(*,*) "Finished Read Input"
   WRITE(*,*) "member_total = ", InputFileData%member_total
   WRITE(*,*) "temp_GL: ", temp_GLL(:)
   DO i=1,InputFileData%member_total*2+1
       WRITE(*,*) "kp_coordinate:", InputFileData%kp_coordinate(i,:)
       WRITE(*,*) "initial_twist:", InputFileData%initial_twist(i)
   ENDDO
   DO i=1,InputFiledata%member_total
       WRITE(*,*) "ith_member_length",i,p%member_length(i,:)
!       WRITE(*,*) "temp_ratio: ", temp_ratio(:,i)
       DO j=1,p%node_elem
           WRITE(*,*) "Nodal Position:",j
           WRITE(*,*) p%uuN0((j-1)*6+1,i),p%uuN0((j-1)*6+2,i),p%uuN0((j-1)*6+3,i)
           WRITE(*,*) p%uuN0((j-1)*6+4,i),p%uuN0((j-1)*6+5,i),p%uuN0((j-1)*6+6,i)
       ENDDO
   ENDDO
   WRITE(*,*) "Blade Length: ", p%blade_length
   WRITE(*,*) "node_elem: ", p%node_elem
!   WRITE(*,*) "Stiff0: ", InputFileData%InpBl%stiff0(4,:,1)
!   WRITE(*,*) "Stiff0: ", InputFileData%InpBl%stiff0(4,:,2)
!   WRITE(*,*) "Stiff0: ", InputFileData%InpBl%stiff0(4,:,3)

!   WRITE(*,*) "Stiff0_GL: ", p%Stif0_GL(4,:,1)
!   WRITE(*,*) "Stiff0_GL: ", p%Stif0_GL(4,:,2)
!   WRITE(*,*) "Stiff0_GL: ", p%Stif0_GL(4,:,3)
!   WRITE(*,*) "Stiff0_GL: ", p%Stif0_GL(4,:,4)
!   STOP
   ! Define parameters here:

   p%node_total  = p%elem_total*(p%node_elem-1) + 1         ! total number of node  
   p%dof_total   = p%node_total*p%dof_node   ! total number of dof
   p%dt = Interval

   ! Allocate OtherState if using multi-step method; initialize n


   ! Allocate continuous states and define initial system states here:

   CALL AllocAry(x%q,p%dof_total,'x%q',ErrStat2,ErrMsg2)
   x%q = 0.0D0
   CALL AllocAry(x%dqdt,p%dof_total,'x%dqdt',ErrStat2,ErrMsg2)
   x%dqdt = 0.0D0

   ! Define system output initializations (set up mesh) here:
   CALL MeshCreate( BlankMesh        = u%PointMesh            &
                   ,IOS              = COMPONENT_INPUT        &
                   ,NNodes           = 1                      &
                   , TranslationDisp = .TRUE. &
                   , TranslationVel  = .TRUE. &
                   , TranslationAcc  = .TRUE. &
                   , Orientation     = .TRUE. &
                   , RotationVel     = .TRUE. &
                   , RotationAcc     = .TRUE. &
                   ,nScalars        = 0                      &
                   ,ErrStat         = ErrStat               &
                   ,ErrMess         = ErrMsg                )

   CALL MeshConstructElement ( Mesh = u%PointMesh            &
                             , Xelement = ELEMENT_POINT      &
                             , P1       = 1                  &
                             , ErrStat  = ErrStat            &
                             , ErrMess  = ErrMsg             )

   ! place single node at origin; position affects mapping/coupling with other modules
   TmpPos(1) = 0.
   TmpPos(2) = 0.
   TmpPos(3) = 0.

   CALL MeshPositionNode ( Mesh = u%PointMesh          &
                         , INode = 1                &
                         , Pos = TmpPos             &
                         , ErrStat   = ErrStat      &
                         , ErrMess   = ErrMsg       )

   CALL MeshCommit ( Mesh    = u%PointMesh     &
                    ,ErrStat = ErrStat         &
                    ,ErrMess = ErrMsg          )

   CALL MeshCopy ( SrcMesh  = u%PointMesh      &
                 , DestMesh = y%PointMesh      &
                 , CtrlCode = MESH_SIBLING     &
                 , Force           = .TRUE.    &
                 , Moment          = .TRUE.    &
                 , ErrStat  = ErrStat          &
                 , ErrMess  = ErrMsg           )

   ! Define initialization-routine input here:

   u%PointMesh%TranslationDisp(1,1) = 0.
   u%PointMesh%TranslationDisp(2,1) = 0.
   u%PointMesh%TranslationDisp(3,1) = 0.

   u%PointMesh%TranslationVel(1,1) = 0.
   u%PointMesh%TranslationVel(2,1) = 0.
   u%PointMesh%TranslationVel(3,1) = 0.

   u%PointMesh%TranslationAcc(1,1) = 0.
   u%PointMesh%TranslationAcc(2,1) = 0.
   u%PointMesh%TranslationAcc(3,1) = 0.

   u%PointMesh%Orientation = 0.
   u%PointMesh%Orientation(1,1,1) = 1.
   u%PointMesh%Orientation(2,2,1) = 1.
   u%PointMesh%Orientation(3,3,1) = 1.

   u%PointMesh%RotationVel(1,1) = 0.
   u%PointMesh%RotationVel(2,1) = 0.
   u%PointMesh%RotationVel(3,1) = 0.

   u%PointMesh%RotationAcc(1,1) = 0.
   u%PointMesh%RotationAcc(2,1) = 0.
   u%PointMesh%RotationAcc(3,1) = 0.

   ! Define initial guess for the system outputs here:

   y%PointMesh%Force(1,1)   = 0.
   y%PointMesh%Force(2,1)   = 0.
   y%PointMesh%Force(3,1)   = 0.

   y%PointMesh%Moment(1,1)   = 0.
   y%PointMesh%Moment(2,1)   = 0.
   y%PointMesh%Moment(3,1)   = 0.


   ! set remap flags to true
   y%PointMesh%RemapFlag = .True.
   u%PointMesh%RemapFlag = .True.


   END SUBROUTINE BeamDyn_Init

   !----------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE BeamDyn_End( u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
   !
   ! This routine is called at the end of the simulation.
   !..................................................................................................................................

   TYPE(BD_InputType),           INTENT(INOUT)  :: u           ! System inputs
   TYPE(BD_ParameterType),       INTENT(INOUT)  :: p           ! Parameters
   TYPE(BD_ContinuousStateType), INTENT(INOUT)  :: x           ! Continuous states
   TYPE(BD_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Discrete states
   TYPE(BD_ConstraintStateType), INTENT(INOUT)  :: z           ! Constraint states
   TYPE(BD_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
   TYPE(BD_OutputType),          INTENT(INOUT)  :: y           ! System outputs
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = "" 

   ! Place any last minute operations or calculations here:

   ! Close files here:

   ! Destroy the input data:

   CALL BD_DestroyInput( u, ErrStat, ErrMsg )

   ! Destroy the parameter data:

   CALL BD_DestroyParam( p, ErrStat, ErrMsg )

   ! Destroy the state data:

   CALL BD_DestroyContState(   x,           ErrStat, ErrMsg )
   CALL BD_DestroyDiscState(   xd,          ErrStat, ErrMsg )
   CALL BD_DestroyConstrState( z,           ErrStat, ErrMsg )
   CALL BD_DestroyOtherState(  OtherState,  ErrStat, ErrMsg )

   ! Destroy the output data:

   CALL BD_DestroyOutput( y, ErrStat, ErrMsg )


   END SUBROUTINE BeamDyn_End

   SUBROUTINE BeamDyn_UpdateStates( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )
!
! Routine for solving for constraint states, integrating continuous states, and updating discrete states
! Constraint states are solved for input t; Continuous and discrete states are updated for t + p%dt
! (stepsize dt assumed to be in ModName parameter)
!..................................................................................................................................

   REAL(DbKi),                      INTENT(IN   ) :: t          ! Current simulation time in seconds
   INTEGER(IntKi),                  INTENT(IN   ) :: n          ! Current simulation time step n = 0,1,...
   TYPE(BD_InputType),            INTENT(INOUT) :: u(:)       ! Inputs at utimes
   REAL(DbKi),                      INTENT(IN   ) :: utimes(:)  ! Times associated with u(:), in seconds
   TYPE(BD_ParameterType),        INTENT(IN   ) :: p          ! Parameters
   TYPE(BD_ContinuousStateType),  INTENT(INOUT) :: x          ! Input: Continuous states at t;
                                                                      !   Output: Continuous states at t + Interval
   TYPE(BD_DiscreteStateType),    INTENT(INOUT) :: xd         ! Input: Discrete states at t;
                                                                      !   Output: Discrete states at t  + Interval
   TYPE(BD_ConstraintStateType),  INTENT(INOUT) :: z          ! Input: Initial guess of constraint states at t+dt;
                                                                      !   Output: Constraint states at t+dt
   TYPE(BD_OtherStateType),       INTENT(INOUT) :: OtherState ! Other/optimization states
   INTEGER(IntKi),                  INTENT(  OUT) :: ErrStat    ! Error status of the operation
   CHARACTER(*),                    INTENT(  OUT) :: ErrMsg     ! Error message if ErrStat /= ErrID_None

   ! local variables
   INTEGER(IntKi):: i
   INTEGER(IntKi):: j
   INTEGER(IntKi):: temp_id
   REAL(ReKi):: temp
   REAL(ReKi):: temp_pp(3)
   REAL(ReKi):: temp_qq(3)
   REAL(ReKi):: temp_rr(3)
   ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = "" 

   CALL BeamDyn_RK4( t, n, u, utimes, p, x, xd, z, OtherState, ErrStat, ErrMsg )

!   DO i=1,3
!       x%q(i) = u(1)%PointMesh%TranslationDisp(i,1)
!       x%q(i+3) = 0.0D0
!       x%dqdt(i) = u(1)%PointMesh%TranslationVel(i,1)
!       x%dqdt(i+3) = 0.0D0
!   ENDDO

   temp_pp(:) = 0.0D0
   temp_qq(:) = 0.0D0
   temp_rr(:) = 0.0D0
   x%q(5) = -4.0D0*TAN((3.1415926D0*t/3.0D0)/4.0D0)
   DO i=1,3
       temp_pp(i) = x%q(i+3)
   ENDDO
   CALL CrvCompose(temp_rr,temp_pp,temp_qq,0)
   DO i=1,3
       x%q(i+3) = temp_rr(i)
   ENDDO
!   IF(ABS(x%q(5)) .GT. 4.0D0) THEN
!       x%q(5) = -4.0D0*TAN((3.1415926D0*t/3.0D0+2.0D0*3.1415926D0)/4.0D0)
!   ENDIF
   x%dqdt(5) = -3.1415926D0/3.0D0

   DO i=2,p%node_total
       temp_id = (i-1)*6
       temp_pp(:) = 0.0D0
       temp_qq(:) = 0.0D0
       temp_rr(:) = 0.0D0
       DO j=1,3
           temp_pp(j) = x%q(temp_id+3+j)
       ENDDO
       CALL CrvCompose(temp_rr,temp_pp,temp_qq,0)
       DO j=1,3
           x%q(temp_id+3+j) = temp_rr(j)
       ENDDO
!       IF(ABS(x%q(temp_id+5)) .GT. 4.0D0) THEN
!           temp = -4.0D0*ATAN(x%q(temp_id+5)/4.0D0)
!           x%q(temp_id+5) = -4.0D0*TAN((temp+2.0D0*3.1415926D0)/4.0D0)
!       ENDIF
   ENDDO

   END SUBROUTINE BeamDyn_UpdateStates

   !----------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE BeamDyn_CalcOutput( t, u, p, x, xd, z, OtherState, y, ErrStat, ErrMsg )
   !
   ! Routine for computing outputs, used in both loose and tight coupling.
   !..................................................................................................................................

   REAL(DbKi),                        INTENT(IN   )  :: t           ! Current simulation time in seconds
   TYPE(BD_InputType),           INTENT(IN   )  :: u           ! Inputs at t
   TYPE(BD_ParameterType),       INTENT(IN   )  :: p           ! Parameters
   TYPE(BD_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
   TYPE(BD_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t
   TYPE(BD_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t
   TYPE(BD_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
   TYPE(BD_OutputType),          INTENT(INOUT)  :: y           ! Outputs computed at t (Input only so that mesh con-
                                                                    !   nectivity information does not have to be recalculated)
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = "" 

   ! see Eqs. (12), (13)  in Gasmi et al. (2013)
   !y%q    = x%q
   !y%dqdt = x%dqdt

   ! assume system is aligned with x-axis

   y%PointMesh%Force(1,1)   = x%q(1)
   y%PointMesh%Force(2,1)   = x%q(2)
   y%PointMesh%Force(3,1)   = x%q(3)

   y%PointMesh%Moment(1,1)   = x%q(4)
   y%PointMesh%Moment(2,1)   = x%q(5)
   y%PointMesh%Moment(3,1)   = x%q(6)

!   WRITE(67,*) t,  y%PointMesh%TranslationAcc(1,1)

   END SUBROUTINE BeamDyn_CalcOutput

   !----------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE BeamDyn_UpdateDiscState( t, n, u, p, x, xd, z, OtherState, ErrStat, ErrMsg )
   !
   ! Routine for updating discrete states
   !..................................................................................................................................

   REAL(DbKi),                        INTENT(IN   )  :: t           ! Current simulation time in seconds
   INTEGER(IntKi),                    INTENT(IN   )  :: n           ! Current step of the simulation: t = n*Interval
   TYPE(BD_InputType),           INTENT(IN   )  :: u           ! Inputs at t
   TYPE(BD_ParameterType),       INTENT(IN   )  :: p           ! Parameters
   TYPE(BD_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
   TYPE(BD_DiscreteStateType),   INTENT(INOUT)  :: xd          ! Input: Discrete states at t;
                                                                    !   Output: Discrete states at t + Interval
   TYPE(BD_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t
   TYPE(BD_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
   INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

      ! Initialize ErrStat

      ErrStat = ErrID_None
      ErrMsg  = "" 

      ! Update discrete states here:

!      xd%DummyDiscState = 0.0

END SUBROUTINE BeamDyn_UpdateDiscState
   !----------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE BeamDyn_CalcConstrStateResidual( t, u, p, x, xd, z, OtherState, Z_residual, ErrStat, ErrMsg )
   !
   ! Routine for solving for the residual of the constraint state equations
   !..................................................................................................................................

   REAL(DbKi),                        INTENT(IN   )  :: t           ! Current simulation time in seconds
   TYPE(BD_InputType),           INTENT(IN   )  :: u           ! Inputs at t
   TYPE(BD_ParameterType),       INTENT(IN   )  :: p           ! Parameters
   TYPE(BD_ContinuousStateType), INTENT(IN   )  :: x           ! Continuous states at t
   TYPE(BD_DiscreteStateType),   INTENT(IN   )  :: xd          ! Discrete states at t
   TYPE(BD_ConstraintStateType), INTENT(IN   )  :: z           ! Constraint states at t (possibly a guess)
   TYPE(BD_OtherStateType),      INTENT(INOUT)  :: OtherState  ! Other/optimization states
   TYPE(BD_ConstraintStateType), INTENT(  OUT)  :: Z_residual  ! Residual of the constraint state equations using
                                                                    !     the input values described above
   INTEGER(IntKi),                 INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                   INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None


   ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = "" 


   ! Solve for the constraint states here:

   Z_residual%DummyConstrState = 0

   END SUBROUTINE BeamDyn_CalcConstrStateResidual

   !----------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE BD_gen_gll_LSGL(N, x, w)
   !
   ! This subroutine determines the (N+1) Gauss-Lobatto-Legendre points x and weights w
   !
   ! For details, see
   ! @book{Deville-etal:2002,
   !  author =    {M. O. Deville and P. F. Fischer and E. H. Mund},
   !  title =     {High-Order Methods for Incompressible Fluid Flow},
   !  publisher = {Cambridge University Press},
   !  address = {Cambridge},
   !  year =      2002
   !}
   !
   !..................................................................................................................................

   ! input variables

   INTEGER(IntKi),                 INTENT(IN   )  :: N           ! Order of spectral element
   REAL(ReKi),                     INTENT(  OUT)  :: x(N+1)      ! location of GLL nodes
   REAL(ReKi),                     INTENT(  OUT)  :: w(N+1)      ! quadrature weights at GLL nodes


   ! local variables  

   REAL(ReKi)          :: tol       ! tolerance for newton-raphson solve
   INTEGER(IntKi)      :: maxit     ! maximum allowable iterations in newton-raphson solve
   REAL(ReKi)          :: x_it      ! current NR-iteration value
   REAL(ReKi)          :: x_old     ! last NR-iteration value

   REAL(ReKi)          :: dleg(N+1)   ! legendre polynomial

   INTEGER(IntKi)      :: N1        ! N+1

   INTEGER(IntKi)      :: i         ! do-loop counter
   INTEGER(IntKi)      :: j         ! do-loop counter
   INTEGER(IntKi)      :: k         ! do-loop counter


   tol = 1e-15

   N1 = N+1

   maxit = 1e3  

   ! enter known endpoints  [-1.0, 1.0]
   x(1) = -1.
   x(N1) = 1.

   pi = ACOS(-1.)  ! perhaps use NWTC library value, but does not matter here; just used to guess at solution

   DO i = 1, N1

      x_it = -COS(pi * FLOAT(i-1) / N) ! initial guess - chebyshev points

      DO j = 1, maxit
         x_old = x_it
         dleg(1) = 1.0
         dleg(2) = x_it
         DO k = 2,N
            dleg(k+1) = (  (2.0*DFLOAT(k) - 1.0) * dleg(k) * x_it &
                            - (DFLOAT(k)-1.0)*dleg(k-1) ) / DFLOAT(k)
         ENDDO

         x_it = x_it - ( x_it * dleg(N1) - dleg(N) ) / &
                       (DFLOAT(N1) * dleg(N1) )

         IF (ABS(x_it - x_old) .lt. tol) THEN
            EXIT
         ENDIF
      ENDDO

      x(i) = x_it
      w(i) = 2.0 / (DFLOAT(N * N1) * dleg(N1)**2 )

   ENDDO

   END SUBROUTINE BD_gen_gll_LSGL


END MODULE BeamDyn