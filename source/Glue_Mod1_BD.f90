!..................................................................................................................................
! LICENSING
! Copyright (C) 2013  National Renewable Energy Laboratory
!
!    Glue is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as
!    published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
!
!    This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty
!    of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License along with Module2.
!    If not, see <http://www.gnu.org/licenses/>.
!
!**********************************************************************************************************************************
!    This file is the "glue code" example for the new FAST modularization.
!
!    This code embodies the formulation and numerical methods for "Predictor-Corrector Loose Coupling."   The methods
!    employed here share similarities with those described in Gasmi et al. (2013).   However, the method used here is a
!    "symmetric" predictor-corrector approach (order of module UpdateStates does not matter).    Further, the method reduces
!    to an explicit method when only the prediction step is taken (pc_max = 1 below).   However, for all modules, inputs
!    and outputs are stored for up to three steps, allowing up to quadratic interpolation or exptrapolation of input and
!    output data.
!
!    The test problem is a simple two-degree-of-freedom damped oscillator, where each "mass" is treated by a module; see
!    Gasmi et al. (2013) for details.
!
!    Three fourth-order explicit numerical time integrators are included: Runge-Kutta (RK4), Adams-Bashforth (AB4), and
!    Adams-Bashforth-Moulton (ABM4).    RK4 and ABM4 have an implcit dependence on other-module data.
!
!    Numerical experiments have shown that, if quadratic interpolation of inputs and outpus is employed, order of accuracy of
!    the methods with pc_max predictor-corrector iterations are as follows (with Mod1 & Mod2 using same integrator):
!
!    RK4, PC1: third order
!    RK4, PC2: third order (but more accurate than PC1)
!
!    AB4, PC1: fourth order
!    AB4, PC2: fourth order (should be identical to PC1; no implicit dependence on other-module data)
!
!    ABM4, PC1: third order
!    ABM4, PC2: fourth order
!
!    NOTE: These convergence results can be obtained only when the multi-step methods have their first three steps initialized
!          with the exact benchmark solution.
!
!    References:
!
!    Gasmi, A., M. A. Sprague, J. M. Jonkman, and W. B. Jones, Numerical stability and accuracy of temporally coupled
!    multi-physics modules in wind turbine CAE tools. In proceedings of the 32nd ASME Wind Energy Symposium, 51st AIAA
!    Aerospace Sciences Meeting including the New Horizons Forum and Aerospace Exposition, Grapevine, TX, January 7-10,
!    2013.   Also published as NREL Report No. CP-2C00-57298.   Available in pdf format at:
!    http://www.nrel.gov/docs/fy13osti/57298.pdf
!
!**********************************************************************************************************************************
MODULE Mod1_BD_MappingModule

   USE Module1
   USE Module1_Types

   USE BeamDyn
   USE BeamDyn_Subs
   USE BeamDyn_Types

   USE NWTC_Library
   USE NWTC_LAPACK

   implicit none

CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE Mod1_BD_InputOutputSolve(time, &
                   Mod1_Input, Mod1_Parameter, Mod1_ContinuousState, Mod1_DiscreteState, &
                   Mod1_ConstraintState, Mod1_OtherState, Mod1_Output, &
                   BD_Input, BD_Parameter, BD_ContinuousState, BD_DiscreteState, &
                   BD_ConstraintState, BD_OtherState, BD_Output,  &
!                   Map_Mod1_P_BD_P, Map_BD_P_Mod1_P, &
                   ErrStat, ErrMsg)
!
! Solve input-output relations for Module 1 coupled to BeamDyn; this section of code corresponds to Eq. (35) in
! Gasmi et al. (2013). This code will be specific to the underlying modules
!...................................................................................................................................


   ! Module1 Derived-types variables; see Registry_Module1.txt for details

   TYPE(Mod1_InputType),           INTENT(INOUT) :: Mod1_Input
   TYPE(Mod1_ParameterType),       INTENT(IN   ) :: Mod1_Parameter
   TYPE(Mod1_ContinuousStateType), INTENT(IN   ) :: Mod1_ContinuousState
   TYPE(Mod1_DiscreteStateType),   INTENT(IN   ) :: Mod1_DiscreteState
   TYPE(Mod1_ConstraintStateType), INTENT(INOUT) :: Mod1_ConstraintState
   TYPE(Mod1_OtherStateType),      INTENT(INOUT) :: Mod1_OtherState
   TYPE(Mod1_OutputType),          INTENT(INOUT) :: Mod1_Output

   ! BeamDyn Derived-types variables; see Registry_BeamDyn.txt

   TYPE(BD_InputType),           INTENT(INOUT) :: BD_Input
   TYPE(BD_ParameterType),       INTENT(IN   ) :: BD_Parameter
   TYPE(BD_ContinuousStateType), INTENT(INOUT) :: BD_ContinuousState
   TYPE(BD_DiscreteStateType),   INTENT(IN   ) :: BD_DiscreteState
   TYPE(BD_ConstraintStateType), INTENT(INOUT) :: BD_ConstraintState
   TYPE(BD_OtherStateType),      INTENT(INOUT) :: BD_OtherState
   TYPE(BD_OutputType),          INTENT(INOUT) :: BD_Output

   ! mapping stuff

   INTEGER(IntKi),               INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                 INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None
   REAL(DbKi),                   INTENT(IN   )  :: time        ! Current simulation time in seconds

   REAL(ReKi)                                   :: BD_Force
   REAL(ReKi)                                   :: BD_RootAcc
   REAL(ReKi)                                   :: Mod1_Force
   REAL(ReKi)                                   :: Mod1_Acc
   REAL(ReKi)                                   :: RHS(2)
   REAL(ReKi)                                   :: Coef(2,2)
   REAL(ReKi)                                   :: eps
   REAL(ReKi)                                   :: d
   REAL(ReKi)                                   :: uinc(2)
   REAL(ReKi),                        PARAMETER :: TOLF = 1.0D-10
   INTEGER(IntKi)                               :: indx(2)
   INTEGER(IntKi)                               :: i
   INTEGER(IntKi),                    PARAMETER :: iter_max = 10
   TYPE(BD_OutputType)                          :: OT_tmp
   TYPE(Mod1_OutputType)                          :: Mod1OT_tmp
   INTEGER(IntKi)          :: ErrStat2                     ! Temporary Error status
   CHARACTER(ErrMsgLen)    :: ErrMsg2                      ! Temporary Error message
   CHARACTER(*), PARAMETER :: RoutineName = 'BD_InputOutputSolve'

   ErrStat = ErrID_None
   ErrMsg  = ""


   ! eps is the perturbation magnitude in the approximation of the Jacobian
   ! mas: I think the system is linear here, and should thus be insensitive to changes in eps; Qi verified that results
   ! are indeed insensitive.
    
   eps = 0.01D+00

   ! These BD inputs come directly from the Mod1 states; no solve required
   BD_Input%RootMotion%TranslationDisp(3,1) = Mod1_Output%PointMesh%TranslationDisp(1,1)
   BD_Input%RootMotion%TranslationVel(3,1) = Mod1_Output%PointMesh%TranslationVel(1,1)
!   BD_Input%RootMotion%TranslationDisp(1,1) = 0.0D0
!   BD_Input%RootMotion%TranslationDisp(2,1) = 0.0D0
!   BD_Input%RootMotion%TranslationVel(1,1) = 0.0D0
!   BD_Input%RootMotion%TranslationVel(2,1) = 0.0D0

   DO i=1,iter_max
!WRITE(*,*) 'i=',i
       CALL Mod1_CalcOutput( time, Mod1_Input, Mod1_Parameter, Mod1_ContinuousState, Mod1_DiscreteState, &
                    Mod1_ConstraintState, Mod1_OtherState, Mod1_Output, ErrStat, ErrMsg )

!       CALL BD_InputClean(BD_Input)
       CALL BD_CalcOutput( time, BD_Input, BD_Parameter, BD_ContinuousState, BD_DiscreteState, &
                    BD_ConstraintState, BD_OtherState, BD_Output, ErrStat, ErrMsg )
!IF(i .EQ. 2) EXIT

       BD_Force = BD_Output%ReactionForce%Force(3,1)
       BD_RootAcc = BD_Input%RootMotion%TranslationAcc(3,1)


       ! calculate RHS, which is the current iteration residual
       RHS(:) = 0.0D0
       RHS(1) = -(Mod1_Input%PointMesh%Force(1,1) - BD_Output%ReactionForce%Force(3,1))
       RHS(2) = -(BD_Input%RootMotion%TranslationAcc(3,1) - Mod1_Output%PointMesh%TranslationAcc(1,1))

       IF(TwoNorm(RHS) .LE. TOLF) THEN
           RETURN
       ENDIF
    
       Coef(:,:) = 0.0D0
       Coef(1,1) = 1.0D0

       BD_Input%RootMotion%TranslationAcc(3,1) = BD_Input%RootMotion%TranslationAcc(3,1) + eps

!       CALL BD_InputClean(BD_Input)
       CALL BD_CopyOutput(BD_Output, OT_tmp, MESH_NEWCOPY, ErrStat, ErrMsg)

       CALL BD_CalcOutput( time, BD_Input, BD_Parameter, BD_ContinuousState, BD_DiscreteState, &
                    BD_ConstraintState, BD_OtherState, OT_tmp, ErrStat, ErrMsg )

       Coef(1,2) = -((OT_tmp%ReactionForce%Force(3,1)-BD_Force)/eps)

       ! replace BD_Input with unperturbed value
       BD_Input%RootMotion%TranslationAcc(3,1) = BD_RootAcc

       Coef(2,2) = 1.0D0

       Mod1_Force = Mod1_Input%PointMesh%Force(1,1)
       Mod1_Acc = Mod1_Output%PointMesh%TranslationAcc(1,1)

       Mod1_Input%PointMesh%Force(1,1) = Mod1_Input%PointMesh%Force(1,1) + eps
       CALL Mod1_CopyOutput(Mod1_Output, Mod1OT_tmp, MESH_NEWCOPY, ErrStat, ErrMsg)
       CALL Mod1_CalcOutput( time, Mod1_Input, Mod1_Parameter, Mod1_ContinuousState, Mod1_DiscreteState, &
                             Mod1_ConstraintState, Mod1_OtherState, Mod1OT_tmp, ErrStat, ErrMsg )
       Coef(2,1) = -((Mod1OT_tmp%PointMesh%TranslationAcc(1,1) - Mod1_Acc)/eps)
       Mod1_Input%PointMesh%Force(1,1) = Mod1_Force
 
       CALL LAPACK_getrf( 2, 2, Coef,indx, ErrStat2, ErrMsg2) 
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )
       CALL LAPACK_getrs( 'N',2, Coef,indx,RHS, ErrStat2, ErrMsg2) 
          CALL SetErrStat( ErrStat2, ErrMsg2, ErrStat, ErrMsg, RoutineName )

       uinc(:) = RHS(:)

!WRITE(*,*) 'uinc'
!WRITE(*,*) uinc
       IF(TwoNorm(uinc) .LE. TOLF) RETURN

       Mod1_Input%PointMesh%Force(1,1) = Mod1_Input%PointMesh%Force(1,1) + uinc(1)
       BD_Input%RootMotion%TranslationAcc(3,1) = BD_Input%RootMotion%TranslationAcc(3,1) + uinc(2)
      
       IF(i .EQ. iter_max) THEN
           WRITE(*,*) "InputOutputSolve does not converge after the maximum number of iterations"
!           RETURN
           STOP 
       ENDIF

   ENDDO


END SUBROUTINE Mod1_BD_InputOutputSolve
!----------------------------------------------------------------------------------------------------------------------------------
   SUBROUTINE BD_InputClean(BD_Input)

   TYPE(BD_InputType),           INTENT(INOUT) :: BD_Input

   BD_Input%RootMotion%TranslationDisp(2:3,1) = 0.0D0
   BD_Input%RootMotion%TranslationVel(2:3,1)  = 0.0D0
!   BD_Input%RootMotion%TranslationAcc(2:3,1)  = 0.0D0
   BD_Input%RootMotion%Orientation(:,:,:) = 0.0D0
   BD_Input%RootMotion%Orientation(1,1,1) = 1.0D0
   BD_Input%RootMotion%Orientation(2,2,1) = 1.0D0
   BD_Input%RootMotion%Orientation(3,3,1) = 1.0D0
   BD_Input%RootMotion%RotationVel(:,:)   = 0.0D0
   BD_Input%RootMotion%RotationAcc(:,:)   = 0.0D0

   END SUBROUTINE BD_InputClean

!----------------------------------------------------------------------------------------------------------------------------------
END MODULE Mod1_BD_MappingModule
!----------------------------------------------------------------------------------------------------------------------------------
PROGRAM MAIN

   use Mod1_BD_MappingModule

   USE Module1
   USE Module1_Types

   USE BeamDyn
   USE BeamDyn_Subs  ! for crv extract routines
   USE BeamDyn_Types

   USE NWTC_Library

   IMPLICIT NONE

   ! global glue-code-specific variables

   INTEGER(IntKi)                     :: ErrStat          ! Error status of the operation
   CHARACTER(1024)                    :: ErrMsg           ! Error message if ErrStat /= ErrID_None

   REAL(DbKi)                         :: dt_global        ! fixed/constant global time step
   REAL(DbKi)                         :: t_initial        ! time at initialization
   REAL(DbKi)                         :: t_final          ! time at simulation end
   REAL(DbKi)                         :: t_global         ! global-loop time marker

   INTEGER(IntKi)                     :: n_t_final        ! total number of time steps
   INTEGER(IntKi)                     :: n_t_global       ! global-loop time counter

   INTEGER(IntKi)                     :: pc_max           ! 1:explicit loose; 2:pc loose
   INTEGER(IntKi)                     :: pc               ! counter for pc iterations

   INTEGER(IntKi)                     :: Mod1_interp_order     ! order of interpolation/extrapolation
   INTEGER(IntKi)                     :: BD_interp_order     ! order of interpolation/extrapolation

   INTEGER(IntKi)                     :: MaxPtsInMap      ! the maximum number of points in a mapping
   
   ! Module1 Derived-types variables; see Registry_Module1.txt for details

   TYPE(Mod1_InitInputType)           :: Mod1_InitInput
   TYPE(Mod1_ParameterType)           :: Mod1_Parameter
   TYPE(Mod1_ContinuousStateType)     :: Mod1_ContinuousState
   TYPE(Mod1_ContinuousStateType)     :: Mod1_ContinuousStateDeriv
   TYPE(Mod1_InitOutputType)          :: Mod1_InitOutput
   TYPE(Mod1_DiscreteStateType)       :: Mod1_DiscreteState
   TYPE(Mod1_ConstraintStateType)     :: Mod1_ConstraintState
   TYPE(Mod1_OtherStateType)          :: Mod1_OtherState

   TYPE(Mod1_InputType),Dimension(:),ALLOCATABLE   :: Mod1_Input
   REAL(DbKi) , DIMENSION(:), ALLOCATABLE          :: Mod1_InputTimes

   TYPE(Mod1_OutputType),Dimension(:),ALLOCATABLE  :: Mod1_Output
   REAL(DbKi) , DIMENSION(:), ALLOCATABLE          :: Mod1_OutputTimes

   TYPE(Mod1_InputType)   :: u1    ! local variable for extrapolated inputs
   TYPE(Mod1_OutputType)  :: y1    ! local variable for extrapolated outputs

   ! Module 1 deived data typed needed in pc-coupling; predicted states

   TYPE(Mod1_ContinuousStateType)     :: Mod1_ContinuousState_pred
   TYPE(Mod1_DiscreteStateType)       :: Mod1_DiscreteState_pred
   TYPE(Mod1_ConstraintStateType)     :: Mod1_ConstraintState_pred

   ! Module2 Derived-types variables; see Registry_Module2.txt

   TYPE(BD_InitInputType)           :: BD_InitInput
   TYPE(BD_ParameterType)           :: BD_Parameter
   TYPE(BD_ContinuousStateType)     :: BD_ContinuousState
   TYPE(BD_ContinuousStateType)     :: BD_ContinuousStateDeriv
   TYPE(BD_InitOutputType)          :: BD_InitOutput
   TYPE(BD_DiscreteStateType)       :: BD_DiscreteState
   TYPE(BD_ConstraintStateType)     :: BD_ConstraintState
   TYPE(BD_OtherStateType)          :: BD_OtherState

! Module 2 deived data typed needed in pc-coupling; predicted states

   TYPE(BD_ContinuousStateType)     :: BD_ContinuousState_pred
   TYPE(BD_DiscreteStateType)       :: BD_DiscreteState_pred
   TYPE(BD_ConstraintStateType)     :: BD_ConstraintState_pred
   TYPE(BD_OtherStateType)          :: BD_OtherState_pred

   TYPE(BD_InputType),DIMENSION(:),ALLOCATABLE   :: BD_Input
   REAL(DbKi) , DIMENSION(:), ALLOCATABLE        :: BD_InputTimes
   TYPE(BD_OutputType),DIMENSION(:),ALLOCATABLE  :: BD_Output
   REAL(DbKi) , DIMENSION(:), ALLOCATABLE        :: BD_OutputTimes

   TYPE(BD_InputType)   :: u2    ! local variable for extrapolated inputs
   TYPE(BD_OutputType)  :: y2    ! local variable for extrapolated outputs

   ! local variables
   Integer(IntKi)                     :: i               ! counter for various loops
   REAL(ReKi):: temp_R(3,3)
   REAL(ReKi):: temp_vec(3)
   REAL(ReKi):: temp_cc(3)
   INTEGER(IntKi),PARAMETER:: QiDisUnit = 20
   INTEGER(IntKi),PARAMETER:: BDForce = 30
   INTEGER(IntKi),PARAMETER:: Mod1Disp = 40
   INTEGER(IntKi),PARAMETER:: Mod1Vel = 50
   INTEGER(IntKi),PARAMETER:: Mod1Acc = 60

REAL(R8Ki):: start, finish
   ! -------------------------------------------------------------------------
   ! MAPPING STUFF; Likely needs to be added to ModMesh
   ! -------------------------------------------------------------------------

!   TYPE(MeshMapType) :: Map_Mod2_P_Mod1_P
!   TYPE(MeshMapType) :: Map_Mod1_P_Mod2_P

   OPEN(unit = QiDisUnit, file = 'QiDisp_Tip_Mod1_BD.out', status = 'REPLACE',ACTION = 'WRITE')
!   OPEN(unit = BDForce, file = 'Qi_Force_Mod1_BD.out', status = 'REPLACE',ACTION = 'WRITE')
   OPEN(unit = Mod1Disp, file = 'Qi_Mod1Disp.out', status = 'REPLACE',ACTION = 'WRITE')
   OPEN(unit = Mod1Vel, file = 'Qi_Mod1Vel.out', status = 'REPLACE',ACTION = 'WRITE')
!   OPEN(unit = Mod1Acc, file = 'Qi_Mod1Acc.out', status = 'REPLACE',ACTION = 'WRITE')
   ! -------------------------------------------------------------------------
   ! Initialization of glue-code time-step variables
   ! -------------------------------------------------------------------------

   t_initial = 0.d0
   t_final   = 5.0D+00

   pc_max = 1  ! Number of predictor-corrector iterations; 1 corresponds to an explicit calculation where UpdateStates
               ! is called only once  per time step for each module; inputs and outputs are extrapolated in time and
               ! are available to modules that have an implicit dependence on other-module data

   ! specify time increment; currently, all modules will be time integrated with this increment size
   ! some notes on dt_global/pc_max choice for stable solutions:
   ! -- pc_max = 1 => dt_global <= 2e-5
   ! -- pc_max = 2 => dt_global <= 5e-5
   ! -- pc_max = 3 => dt_global <= 7e-4
   ! -- pc_max = 4 => dt_global <= 1e-3
   dt_global = 2.0D-05!*0.5

   n_t_final = ((t_final - t_initial) / dt_global )

   t_global = t_initial

   ! define polynomial-order for ModName_Input_ExtrapInterp and ModName_Output_ExtrapInterp
   ! Must be 0, 1, or 2
   Mod1_interp_order = 2
   BD_interp_order   = 2

   !Module1: allocate Input and Output arrays; used for interpolation and extrapolation
   ALLOCATE(Mod1_Input(Mod1_interp_order + 1))
   ALLOCATE(Mod1_InputTimes(Mod1_interp_order + 1)) 
   ALLOCATE(Mod1_Output(Mod1_interp_order + 1))
   ALLOCATE(Mod1_OutputTimes(Mod1_interp_order + 1))

   ! Module2: allocate Input and Output arrays; used for interpolation and extrapolation
   ALLOCATE(BD_Input(BD_interp_order + 1))
   ALLOCATE(BD_InputTimes(BD_interp_order + 1))
   ALLOCATE(BD_Output(BD_interp_order + 1))
   ALLOCATE(BD_OutputTimes(BD_interp_order + 1))

   ! -------------------------------------------------------------------------
   ! Initialization of Modules
   !  note that in this example, we are assuming that dt_global is not changed
   !  in the modules, i.e., that both modules are called at the same glue-code
   !  defined coupling interval.
   ! -------------------------------------------------------------------------

   CALL Mod1_Init( Mod1_InitInput, Mod1_Input(1), Mod1_Parameter, Mod1_ContinuousState, Mod1_DiscreteState, &
                   Mod1_ConstraintState, Mod1_OtherState, Mod1_Output(1), dt_global, Mod1_InitOutput, ErrStat, ErrMsg )

   call Mod1_CopyInput(  Mod1_Input(1), u1, MESH_NEWCOPY, ErrStat, ErrMsg )
   call Mod1_CopyOutput( Mod1_Output(1), y1, MESH_NEWCOPY, ErrStat, ErrMsg )


   ! We fill Mod1_InputTimes with negative times, but the Mod1_Input values are identical for each of those times; this allows
   ! us to use, e.g., quadratic interpolation that effectively acts as a zeroth-order extrapolation and first-order extrapolation
   ! for the first and second time steps.  (The interpolation order in the ExtrapInput routines are determined as
   ! order = SIZE(Mod1_Input)
   DO i = 1, Mod1_interp_order + 1
      Mod1_InputTimes(i) = t_initial - (i - 1) * dt_global
      Mod1_OutputTimes(i) = t_initial - (i - 1) * dt_global
   ENDDO


   DO i = 1, Mod1_interp_order
     Call Mod1_CopyInput (Mod1_Input(i),  Mod1_Input(i+1),  MESH_NEWCOPY, Errstat, ErrMsg)
     Call Mod1_CopyOutput (Mod1_Output(i),  Mod1_Output(i+1),  MESH_NEWCOPY, Errstat, ErrMsg)
   ENDDO 

!   BD_InitInput%InputFile = 'Siemens_53_Input.inp'
   BD_InitInput%InputFile = 'GA2_Debug.inp'
   BD_InitInput%RootName  = TRIM(BD_Initinput%InputFile)
!   ALLOCATE(BD_InitInput%gravity(3)) 
   BD_InitInput%gravity(1) = 0.0D0 !-9.80665
   BD_InitInput%gravity(2) = 0.0D0 
   BD_InitInput%gravity(3) = 0.0D0
!   ALLOCATE(BD_InitInput%GlbPos(3)) 
   BD_InitInput%GlbPos(1) = 0.0D+00
   BD_InitInput%GlbPos(2) = 0.0D+01
   BD_InitInput%GlbPos(3) = 1.0D0
!   ALLOCATE(BD_InitInput%GlbRot(3,3)) 
   BD_InitInput%GlbRot(:,:) = 0.0D0
   temp_vec(1) = 0.0
   temp_vec(2) = 4.0D0*TAN((Pi_D/2.0D0)/4.0D0)
   temp_vec(3) = 0.0
   CALL BD_CrvMatrixR(temp_vec,temp_R,ErrStat,ErrMsg)
   BD_InitInput%GlbRot(1:3,1:3) = TRANSPOSE(temp_R(1:3,1:3))
!   ALLOCATE(BD_InitInput%RootDisp(3))
   BD_InitInput%RootDisp(1) = 0.0D+00
   BD_InitInput%RootDisp(2) = 0.0D+00
   BD_InitInput%RootDisp(3) = 0.1D+00
!   ALLOCATE(BD_InitInput%RootOri(3,3))
   BD_InitInput%RootOri(:,:) = 0.0D0
   temp_vec(1) = 0.0
   temp_vec(2) = 4.0D0*TAN((Pi_D/2.0D0)/4.0D0)
   temp_vec(3) = 0.0
   CALL BD_CrvMatrixR(temp_vec,temp_R,ErrStat,ErrMsg)
   BD_InitInput%RootOri(1:3,1:3) = TRANSPOSE(temp_R(1:3,1:3))
!   ALLOCATE(BD_InitInput%RootVel(6))
   BD_InitInput%RootVel(:) = 0.0D+00

   CALL BD_Init( BD_InitInput, BD_Input(1), BD_Parameter, BD_ContinuousState, BD_DiscreteState, &
                 BD_ConstraintState, BD_OtherState, BD_Output(1), dt_global, BD_InitOutput, ErrStat, ErrMsg )
   

   DO i = 1, BD_interp_order + 1
      BD_InputTimes(i) = t_initial - (i - 1) * dt_global
      BD_OutputTimes(i) = t_initial - (i - 1) * dt_global
   ENDDO

   DO i = 1, BD_interp_order
     Call BD_CopyInput (BD_Input(i),  BD_Input(i+1),  MESH_NEWCOPY, Errstat, ErrMsg)
     Call BD_CopyOutput (BD_Output(i),  BD_Output(i+1), MESH_NEWCOPY, Errstat, ErrMsg)
   ENDDO

      ! Initialize the meshes and/or allocatable arrays in inputs (u2) and outputs (y2) (required fro ExtrapInterp routines)
   CALL BD_CopyInput(BD_Input(1), u2, MESH_NEWCOPY, ErrStat, ErrMsg )
   CALL BD_CopyOutput(BD_Output(1), y2, MESH_NEWCOPY, ErrStat, ErrMsg )  
   
   ! -------------------------------------------------------------------------
   ! Initialize mesh-mapping data
   ! -------------------------------------------------------------------------
!   CALL AllocMapping( Mod1_Output(1)%PointMesh, Mod2_Input(1)%PointMesh, Map_Mod1_P_Mod2_P, ErrStat, ErrMsg )
!   CALL AllocMapping( Mod2_Output(1)%PointMesh, Mod1_Input(1)%PointMesh, Map_Mod2_P_Mod1_P, ErrStat, ErrMsg )
   
!   CALL MeshMapCreate( Mod1_Output(1)%PointMesh, Mod2_Input(1)%PointMesh, Map_Mod1_P_Mod2_P, ErrStat, ErrMsg )
!   CALL MeshMapCreate( Mod2_Output(1)%PointMesh, Mod1_Input(1)%PointMesh, Map_Mod2_P_Mod1_P, ErrStat, ErrMsg )
   ! -------------------------------------------------------------------------
   ! Time-stepping loops
   ! -------------------------------------------------------------------------

   ! write headers for output columns:
   ! write initial condition for q1
   CALL Mod1_BD_InputOutputSolve(t_global, &
                   Mod1_Input(1), Mod1_Parameter, Mod1_ContinuousState, Mod1_DiscreteState, &
                   Mod1_ConstraintState, Mod1_OtherState, Mod1_Output(1), &
                   BD_Input(1), BD_Parameter, BD_ContinuousState, BD_DiscreteState, &
                   BD_ConstraintState, BD_OtherState, BD_Output(1),  &
!                   Map_Mod1_P_Mod2_P, Map_Mod2_P_Mod1_P, &
                   ErrStat, ErrMsg)

!   CALL BD_IniAcc(BD_Input(1),BD_Output(1),BD_Parameter,BD_OtherState,ErrStat,ErrMsg)

CALL CPU_TIME(start)
   DO n_t_global = 0, n_t_final
WRITE(*,*) "Time Step: ", n_t_global
IF(n_t_global .EQ. 100) STOP
      ! Solve input-output relations; this section of code corresponds to Eq. (35) in Gasmi et al. (2013)
      ! This code will be specific to the underlying modules
IF(MOD(n_t_global,5) .EQ. 0) THEN
 CALL BD_CrvExtractCrv(BD_OutPut(1)%BldMotion%Orientation(1:3,1:3,BD_Parameter%node_elem*BD_Parameter%elem_total),temp_cc,ErrStat,ErrMsg)
      WRITE(QiDisUnit,6000) t_global,&
                            &BD_OutPut(1)%BldMotion%TranslationDisp(1:3,BD_Parameter%node_elem*BD_Parameter%elem_total),&
                            &temp_cc(1:3)
!      WRITE(QiDisUnit,6000) t_global,&
!                            &BD_ContinuousState%q(BD_Parameter%dof_total-5:BD_Parameter%dof_total)
!      WRITE(BDForce,6000) t_global,&
!                           &BD_OutPut(1)%ReactionForce%Force(1:3,1),&
!                           &BD_OutPut(1)%ReactionForce%Moment(1:3,1)
      WRITE(Mod1Disp,6000) t_global,&
                           &Mod1_OutPut(1)%PointMesh%TranslationDisp(1:3,1)
      WRITE(Mod1Vel,6000) t_global,&
                           &Mod1_OutPut(1)%PointMesh%TranslationVel(1:3,1),&
                           &Mod1_Output(1)%PointMesh%TranslationAcc(1:3,1)
ENDIF
      ! after all InputOutputSolves, we can reset the mapping flags on the meshes:
         Mod1_Input(1)%PointMesh%RemapFlag    = .FALSE. 
         Mod1_Output(1)%PointMesh%RemapFlag   = .FALSE.
         BD_Input(1)%RootMotion%RemapFlag      = .FALSE. 
         BD_Output(1)%ReactionForce%RemapFlag = .FALSE.

      ! extrapolate inputs and outputs to t + dt; will only be used by modules with an implicit dependence on input data.

      CALL Mod1_Input_ExtrapInterp(Mod1_Input, Mod1_InputTimes, u1, t_global + dt_global, ErrStat, ErrMsg)

!WRITE(*,*) 'u1%Force:',u1%PointMesh%Force(:,1)

      CALL Mod1_Output_ExtrapInterp(Mod1_Output, Mod1_OutputTimes, y1, t_global + dt_global, ErrStat, ErrMsg)

      ! Shift "window" of the Mod1_Input and Mod1_Output

      DO i = Mod1_interp_order, 1, -1
         CALL Mod1_CopyInput (Mod1_Input(i),  Mod1_Input(i+1),  MESH_UPDATECOPY, Errstat, ErrMsg)
         CALL Mod1_CopyOutput (Mod1_Output(i),  Mod1_Output(i+1),  MESH_UPDATECOPY, Errstat, ErrMsg)
         Mod1_InputTimes(i+1) = Mod1_InputTimes(i)
         Mod1_OutputTimes(i+1) = Mod1_OutputTimes(i)
      ENDDO
!WRITE(*,*) 'Mod1_InputTimes:',Mod1_InputTimes(:)

      CALL Mod1_CopyInput (u1,  Mod1_Input(1),  MESH_UPDATECOPY, Errstat, ErrMsg)
      CALL Mod1_CopyOutput (y1,  Mod1_Output(1),  MESH_UPDATECOPY, Errstat, ErrMsg)
      Mod1_InputTimes(1) = t_global + dt_global
      Mod1_OutputTimes(1) = t_global + dt_global

      ! extrapolate inputs and outputs to t + dt; will only be used by modules with an implicit dependence on input data.

      CALL BD_Input_ExtrapInterp(BD_Input, BD_InputTimes, u2, t_global + dt_global, ErrStat, ErrMsg)
      CALL BD_Output_ExtrapInterp(BD_Output, BD_OutputTimes, y2, t_global + dt_global, ErrStat, ErrMsg)

!WRITE(*,*) 'u2%Disp:',u2%RootMotion%TranslationDisp(:,1)
!WRITE(*,*) 'u2%DCM:',u2%RootMotion%Orientation(1,:,1),u2%RootMotion%Orientation(2,:,1),u2%RootMotion%Orientation(3,:,1)

      ! Shift "window" of the Mod1_Input and Mod1_Output

      DO i = BD_interp_order, 1, -1
         CALL BD_CopyInput (BD_Input(i),  BD_Input(i+1),  MESH_UPDATECOPY, Errstat, ErrMsg)
         CALL BD_CopyOutput (BD_Output(i),  BD_Output(i+1),  MESH_UPDATECOPY, Errstat, ErrMsg)
         BD_InputTimes(i+1) = BD_InputTimes(i)
         BD_OutputTimes(i+1) = BD_OutputTimes(i)
      ENDDO

      CALL BD_CopyInput (u2,  BD_Input(1),  MESH_UPDATECOPY, Errstat, ErrMsg)
      CALL BD_CopyOutput (y2, BD_Output(1),  MESH_UPDATECOPY, Errstat, ErrMsg)
      BD_InputTimes(1) = t_global + dt_global
      BD_OutputTimes(1) = t_global + dt_global

      do pc = 1, pc_max

         !----------------------------------------------------------------------------------------
         ! Module 1
         !----------------------------------------------------------------------------------------

         ! copy ContinuousState to ContinuousState_pred; don't modify ContinuousState until completion of PC iterations

         CALL Mod1_CopyContState   (Mod1_ContinuousState, Mod1_ContinuousState_pred, 0, Errstat, ErrMsg)

         CALL Mod1_CopyConstrState (Mod1_ConstraintState, Mod1_ConstraintState_pred, 0, Errstat, ErrMsg)

         CALL Mod1_CopyDiscState   (Mod1_DiscreteState,   Mod1_DiscreteState_pred,   0, Errstat, ErrMsg)

         CALL Mod1_UpdateStates( t_global, n_t_global, Mod1_Input, Mod1_InputTimes, Mod1_Parameter, Mod1_ContinuousState_pred, &
                                 Mod1_DiscreteState_pred, Mod1_ConstraintState_pred, &
                                 Mod1_OtherState, ErrStat, ErrMsg )

         !----------------------------------------------------------------------------------------
         ! Module 2
         !----------------------------------------------------------------------------------------

         ! copy ContinuousState to ContinuousState_pred; don't modify ContinuousState until completion of PC iterations

         CALL BD_CopyContState   (BD_ContinuousState, BD_ContinuousState_pred, 0, Errstat, ErrMsg)

         CALL BD_CopyConstrState (BD_ConstraintState, BD_ConstraintState_pred, 0, Errstat, ErrMsg)

         CALL BD_CopyDiscState   (BD_DiscreteState,   BD_DiscreteState_pred,   0, Errstat, ErrMsg)
         CALL BD_CopyOtherState  (BD_OtherState,   BD_OtherState_pred,   0, Errstat, ErrMsg)

!WRITE(*,*) 'BD_OtherState_pred:',BD_OtherState_pred%Acc(1:3)

         CALL BD_UpdateStates( t_global, n_t_global, BD_Input, BD_InputTimes, BD_Parameter, BD_ContinuousState_pred, &
                                 BD_DiscreteState_pred, BD_ConstraintState_pred, &
                                 BD_OtherState_pred, ErrStat, ErrMsg )
!WRITE(*,*) 'BD_ContinuousState_pred%q:'
!WRITE(*,*) BD_ContinuousState_pred%q(:)

         !-----------------------------------------------------------------------------------------
         ! If correction iteration is to be taken, solve intput-output equations; otherwise move on
         !-----------------------------------------------------------------------------------------

         if (pc .LE. pc_max) then

            call Mod1_BD_InputOutputSolve( t_global + dt_global, &
                                             Mod1_Input(1), Mod1_Parameter, Mod1_ContinuousState_pred, Mod1_DiscreteState_pred, &
                                             Mod1_ConstraintState_pred, Mod1_OtherState, Mod1_Output(1), &
                                             BD_Input(1), BD_Parameter, BD_ContinuousState_pred, BD_DiscreteState_pred, &
                                             BD_ConstraintState_pred, BD_OtherState_pred, BD_Output(1),  &
!                                             Map_Mod1_P_Mod2_P, Map_Mod2_P_Mod1_P, &
                                             ErrStat, ErrMsg)

            
            ! after all InputOutputSolves, we can reset the mapping flags on the meshes:
            Mod1_Input(1)%PointMesh%RemapFlag  = .FALSE. 
            Mod1_Output(1)%PointMesh%RemapFlag = .FALSE.
            BD_Input(1)%RootMotion%RemapFlag      = .FALSE. 
            BD_Output(1)%ReactionForce%RemapFlag = .FALSE.
            
         endif

      enddo

      ! Save all final variables

      CALL Mod1_CopyContState   (Mod1_ContinuousState_pred,  Mod1_ContinuousState, 0, Errstat, ErrMsg)
      CALL Mod1_CopyConstrState (Mod1_ConstraintState_pred,  Mod1_ConstraintState, 0, Errstat, ErrMsg)
      CALL Mod1_CopyDiscState   (Mod1_DiscreteState_pred,    Mod1_DiscreteState,   0, Errstat, ErrMsg)

      CALL BD_CopyContState   (BD_ContinuousState_pred, BD_ContinuousState, 0, Errstat, ErrMsg)
      CALL BD_CopyConstrState (BD_ConstraintState_pred, BD_ConstraintState, 0, Errstat, ErrMsg)
      CALL BD_CopyDiscState   (BD_DiscreteState_pred,   BD_DiscreteState,   0, Errstat, ErrMsg)
      CALL BD_CopyOtherState  (BD_OtherState_pred,   BD_OtherState,   0, Errstat, ErrMsg)

      ! update the global time

      t_global = REAL(n_t_global+1,DbKi) * dt_global + t_initial

      ! the following is exact solution for q_1(t) for baseline parameters in Gasmi et al. (2013)

      ! build rms_error calculation components; see Eq. (56) in Gasmi et al. (2013)

      ! print discrete q_1(t) solution to standard out


   END DO

CALL CPU_TIME(finish)

   WRITE(*,*) 'Start: ', start
   WRITE(*,*) 'Finish: ', finish
   WRITE(*,*) 'Time: ', finish-start

   ! calculate final time normalized rms error


   ! -------------------------------------------------------------------------
   ! Ending of modules
   ! -------------------------------------------------------------------------


   CALL Mod1_End(  Mod1_Input(1), Mod1_Parameter, Mod1_ContinuousState, Mod1_DiscreteState, &
                   Mod1_ConstraintState, Mod1_OtherState, Mod1_Output(1), ErrStat, ErrMsg )

   do i = 2, Mod1_interp_order+1
      CALL Mod1_DestroyInput(Mod1_Input(i), ErrStat, ErrMsg )
      CALL Mod1_DestroyOutput(Mod1_Output(i), ErrStat, ErrMsg )
   enddo

   DEALLOCATE(Mod1_InputTimes)
   DEALLOCATE(Mod1_OutputTimes)

   CALL BD_End(  BD_Input(1), BD_Parameter, BD_ContinuousState, BD_DiscreteState, &
                 BD_ConstraintState, BD_OtherState, BD_Output(1), ErrStat, ErrMsg )

   do i = 2, BD_interp_order+1
      CALL BD_DestroyInput(BD_Input(i), ErrStat, ErrMsg )
      CALL BD_DestroyOutput(BD_Output(i), ErrStat, ErrMsg )
   enddo

   DEALLOCATE(BD_InputTimes)
   DEALLOCATE(BD_Output)

   6000 FORMAT (ES12.5,6ES21.12)
   CLOSE (QiDisUnit)
!   CLOSE (BDForce)
   CLOSE (Mod1Disp)
   CLOSE (Mod1Vel)
!   CLOSE (Mod1Acc)
   ! -------------------------------------------------------------------------
   ! Deallocate arrays associated with mesh mapping
   ! -------------------------------------------------------------------------

!   CALL MeshMapDestroy(Map_Mod1_P_Mod2_P, ErrStat, ErrMsg)
!   CALL MeshMapDestroy(Map_Mod2_P_Mod1_P, ErrStat, ErrMsg)

END PROGRAM MAIN

