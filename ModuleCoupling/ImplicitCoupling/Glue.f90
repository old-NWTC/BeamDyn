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
MODULE BD1_BD_MappingModule

   USE BeamDyn
   USE BeamDyn_Types

   USE NWTC_Library

   implicit none

CONTAINS
!----------------------------------------------------------------------------------------------------------------------------------
SUBROUTINE BD1_BD_InputOutputSolve(time, &
                   BD1_Input, BD1_Parameter, BD1_ContinuousState, BD1_DiscreteState, &
                   BD1_ConstraintState, BD1_OtherState, BD1_Output, &
                   BD_Input, BD_Parameter, BD_ContinuousState, BD_DiscreteState, &
                   BD_ConstraintState, BD_OtherState, BD_Output,  &
!                   Map_Mod1_P_BD_P, Map_BD_P_Mod1_P, &
                   ErrStat, ErrMsg)
!
! Solve input-output relations for Module 1 coupled to BeamDyn; this section of code corresponds to Eq. (35) in
! Gasmi et al. (2013). This code will be specific to the underlying modules
!...................................................................................................................................


   ! Module1 Derived-types variables; see Registry_Module1.txt for details

   TYPE(BD_InputType),           INTENT(INOUT) :: BD1_Input
   TYPE(BD_ParameterType),       INTENT(IN   ) :: BD1_Parameter
   TYPE(BD_ContinuousStateType), INTENT(INOUT) :: BD1_ContinuousState
   TYPE(BD_DiscreteStateType),   INTENT(IN   ) :: BD1_DiscreteState
   TYPE(BD_ConstraintStateType), INTENT(INOUT) :: BD1_ConstraintState
   TYPE(BD_OtherStateType),      INTENT(INOUT) :: BD1_OtherState
   TYPE(BD_OutputType),          INTENT(INOUT) :: BD1_Output

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

   REAL(ReKi)                                   :: RHS(24)
   REAL(ReKi)                                   :: Coef(24,24)
   REAL(ReKi)                                   :: eps
   REAL(ReKi)                                   :: d
   REAL(ReKi)                                   :: uinc(24)
   REAL(ReKi)                                   :: temp_c0(3)
   REAL(ReKi)                                   :: temp_cc(3)
   REAL(ReKi)                                   :: temp3(3)
   REAL(ReKi)                                   :: tempBD_rr(3)
   REAL(ReKi)                                   :: tempBD1_rr(3)
   REAL(ReKi),                        PARAMETER :: TOLF = 1.0D-02
   INTEGER(IntKi)                               :: indx(24)
   INTEGER(IntKi)                               :: i
   INTEGER(IntKi)                               :: j
   INTEGER(IntKi)                               :: k
   INTEGER(IntKi),                    PARAMETER :: iter_max = 10
   TYPE(BD_OutputType)                          :: OT_tmp
   TYPE(BD_OutputType)                         :: BD1OT_tmp
   TYPE(BD_InputType)                           :: BDInput_tmp
   TYPE(BD_InputType)                          :: BD1Input_tmp
   ! Solve input-output relations; this section of code corresponds to Eq. (35) in Gasmi et al. (2013)
   ! This code will be specific to the underlying modules; could be placed in a separate routine.
   ! Note that Module2 has direct feedthrough, but Module1 does not. Thus, Module1 should be called first.

!   BD_Input%PointLoad%Force(3,BD_Parameter%node_total) = 1.0D+05*SIN(0.2*time)
   BD_Input%PointLoad%Force(3,BD_Parameter%node_total) = 0.5*(1.0D0-COS(0.2*time))*1.0D+03
   CALL BD_CopyOutput(BD_Output,OT_tmp,MESH_NEWCOPY,ErrStat,ErrMsg)
   CALL BD_CopyOutput(BD1_Output,BD1OT_tmp,MESH_NEWCOPY,ErrStat,ErrMsg)
!WRITE(*,*) 'TIME',time
   eps = 0.01D+00
   DO i=1,iter_max
!WRITE(*,*) 'i=',i

!WRITE(*,*) 'BD1_Cont%q'
!WRITE(*,*) BD1_ContinuousState%q
!WRITE(*,*) 'BD1_Cont%dqdt'
!WRITE(*,*) BD1_ContinuousState%dqdt
!WRITE(*,*) 'BD1_Other%Acc'
!WRITE(*,*) BD1_OtherState%Acc(:)
       CALL BD_CalcOutput( time, BD1_Input, BD1_Parameter, BD1_ContinuousState, BD1_DiscreteState, &
                    BD1_ConstraintState, BD1_OtherState, BD1_Output, ErrStat, ErrMsg )
!WRITE(*,*) 'BD_Cont%q'
!WRITE(*,*) BD_ContinuousState%q
!WRITE(*,*) 'BD_Cont%dqdt'
!WRITE(*,*) BD_ContinuousState%dqdt
!WRITE(*,*) 'BD_Other%Acc'
!WRITE(*,*) BD_OtherState%Acc(:)
       CALL BD_CalcOutput( time, BD_Input, BD_Parameter, BD_ContinuousState, BD_DiscreteState, &
                    BD_ConstraintState, BD_OtherState, BD_Output, ErrStat, ErrMsg )
!WRITE(*,*) 'Original BD Reaction Force:'
!WRITE(*,*) BD_Output%ReactionForce%Force(1:3,1)
!WRITE(*,*) BD_Output%ReactionForce%Moment(1:3,1)
       CALL BD_CopyInput(BD_Input,BDInput_tmp,MESH_NEWCOPY,ErrStat,ErrMsg)
       CALL BD_CopyInput(BD1_Input,BD1Input_tmp,MESH_NEWCOPY,ErrStat,ErrMsg)

       RHS(:) = 0.0D0
       RHS(1:3) = -(BD1_Input%PointLoad%Force(1:3,BD1_Parameter%node_total) - BD_Output%ReactionForce%Force(1:3,1))
       RHS(4:6) = -(BD1_Input%PointLoad%Moment(1:3,BD1_Parameter%node_total) - BD_Output%ReactionForce%Moment(1:3,1))
       RHS(7:9) = -(BD_Input%RootMotion%TranslationDisp(1:3,1) - &
                    BD1_Output%BldMotion%TranslationDisp(1:3,BD1_Parameter%node_total))
       CALL BD_CrvExtractCrv(BD_Input%RootMotion%Orientation(1:3,1:3,1),temp_cc)
       temp_c0(1:3) = MATMUL(BD_Parameter%GlbRot(1:3,1:3),BD_Parameter%uuN0(4:6,1))
       CALL BD_CrvCompose(tempBD_rr,temp_cc,temp_c0,2)
       CALL BD_CrvExtractCrv(BD1_Output%BldMotion%Orientation(1:3,1:3,BD1_Parameter%node_total),temp_cc)
       temp_c0(1:3) = MATMUL(BD1_Parameter%GlbRot(1:3,1:3),BD1_Parameter%uuN0(4:6,1))
       CALL BD_CrvCompose(tempBD1_rr,temp_cc,temp_c0,2)
       RHS(10:12) = -(tempBD_rr(1:3) - tempBD1_rr(1:3))
       RHS(13:15) = -(BD_Input%RootMotion%TranslationVel(1:3,1) - &
                      BD1_Output%BldMotion%TranslationVel(1:3,BD1_Parameter%node_total))
       RHS(16:18) = -(BD_Input%RootMotion%RotationVel(1:3,1) - &
                      BD1_Output%BldMotion%RotationVel(1:3,BD1_Parameter%node_total))
       RHS(19:21) = -(BD_Input%RootMotion%TranslationAcc(1:3,1) - &
                      BD1_Output%BldMotion%TranslationAcc(1:3,BD1_Parameter%node_total))
       RHS(22:24) = -(BD_Input%RootMotion%RotationAcc(1:3,1) - &
                      BD1_Output%BldMotion%RotationAcc(1:3,BD1_Parameter%node_total))
    
       Coef(:,:) = 0.0D0
       DO j=1,24
           Coef(j,j) = 1.0D0
       ENDDO
       DO j=1,3
           BD_Input%RootMotion%TranslationDisp(j,1) = BD_Input%RootMotion%TranslationDisp(j,1) + eps
           CALL BD_CalcOutput( time, BD_Input, BD_Parameter, BD_ContinuousState, BD_DiscreteState, &
                        BD_ConstraintState, BD_OtherState, OT_tmp, ErrStat, ErrMsg )
           Coef(1:3,6+j) = -((OT_tmp%ReactionForce%Force(1:3,1)-BD_Output%ReactionForce%Force(1:3,1))/eps)
           Coef(4:6,6+j) = -((OT_tmp%ReactionForce%Moment(1:3,1)-BD_Output%ReactionForce%Moment(1:3,1))/eps)
           CALL BD_CopyInput(BDInput_tmp,BD_Input,MESH_NEWCOPY,ErrStat,ErrMsg)
       ENDDO
       temp_c0(1:3) = MATMUL(BD_Parameter%GlbRot(1:3,1:3),BD_Parameter%uuN0(4:6,1))
       DO j=1,3
           temp_cc(:) = tempBD_rr(:)
           temp_cc(j) = temp_cc(j) + eps
           CALL BD_CrvCompose(temp3,temp_cc,temp_c0,0)
           CALL BD_CrvMatrixR(temp3,BD_Input%RootMotion%Orientation(1:3,1:3,1))
           CALL BD_CalcOutput( time, BD_Input, BD_Parameter, BD_ContinuousState, BD_DiscreteState, &
                        BD_ConstraintState, BD_OtherState, OT_tmp, ErrStat, ErrMsg )
           Coef(1:3,9+j) = -((OT_tmp%ReactionForce%Force(1:3,1)-BD_Output%ReactionForce%Force(1:3,1))/eps)
           Coef(4:6,9+j) = -((OT_tmp%ReactionForce%Moment(1:3,1)-BD_Output%ReactionForce%Moment(1:3,1))/eps)
           CALL BD_CopyInput(BDInput_tmp,BD_Input,MESH_NEWCOPY,ErrStat,ErrMsg)
       ENDDO
       DO j=1,3
           BD_Input%RootMotion%TranslationVel(j,1) = BD_Input%RootMotion%TranslationVel(j,1) + eps
           CALL BD_CalcOutput( time, BD_Input, BD_Parameter, BD_ContinuousState, BD_DiscreteState, &
                        BD_ConstraintState, BD_OtherState, OT_tmp, ErrStat, ErrMsg )
           Coef(1:3,12+j) = -((OT_tmp%ReactionForce%Force(1:3,1)-BD_Output%ReactionForce%Force(1:3,1))/eps)
           Coef(4:6,12+j) = -((OT_tmp%ReactionForce%Moment(1:3,1)-BD_Output%ReactionForce%Moment(1:3,1))/eps)
           CALL BD_CopyInput(BDInput_tmp,BD_Input,MESH_NEWCOPY,ErrStat,ErrMsg)
       ENDDO
       DO j=1,3
           BD_Input%RootMotion%RotationVel(j,1) = BD_Input%RootMotion%RotationVel(j,1) + eps
           CALL BD_CalcOutput( time, BD_Input, BD_Parameter, BD_ContinuousState, BD_DiscreteState, &
                        BD_ConstraintState, BD_OtherState, OT_tmp, ErrStat, ErrMsg )
           Coef(1:3,15+j) = -((OT_tmp%ReactionForce%Force(1:3,1)-BD_Output%ReactionForce%Force(1:3,1))/eps)
           Coef(4:6,15+j) = -((OT_tmp%ReactionForce%Moment(1:3,1)-BD_Output%ReactionForce%Moment(1:3,1))/eps)
           CALL BD_CopyInput(BDInput_tmp,BD_Input,MESH_NEWCOPY,ErrStat,ErrMsg)
       ENDDO
       DO j=1,3
           BD_Input%RootMotion%TranslationAcc(j,1) = BD_Input%RootMotion%TranslationAcc(j,1) + eps
           CALL BD_CalcOutput( time, BD_Input, BD_Parameter, BD_ContinuousState, BD_DiscreteState, &
                        BD_ConstraintState, BD_OtherState, OT_tmp, ErrStat, ErrMsg )
           Coef(1:3,18+j) = -((OT_tmp%ReactionForce%Force(1:3,1)-BD_Output%ReactionForce%Force(1:3,1))/eps)
           Coef(4:6,18+j) = -((OT_tmp%ReactionForce%Moment(1:3,1)-BD_Output%ReactionForce%Moment(1:3,1))/eps)
           CALL BD_CopyInput(BDInput_tmp,BD_Input,MESH_NEWCOPY,ErrStat,ErrMsg)
       ENDDO
       DO j=1,3
           BD_Input%RootMotion%RotationAcc(j,1) = BD_Input%RootMotion%RotationAcc(j,1) + eps
           CALL BD_CalcOutput( time, BD_Input, BD_Parameter, BD_ContinuousState, BD_DiscreteState, &
                        BD_ConstraintState, BD_OtherState, OT_tmp, ErrStat, ErrMsg )
           Coef(1:3,21+j) = -((OT_tmp%ReactionForce%Force(1:3,1)-BD_Output%ReactionForce%Force(1:3,1))/eps)
           Coef(4:6,21+j) = -((OT_tmp%ReactionForce%Moment(1:3,1)-BD_Output%ReactionForce%Moment(1:3,1))/eps)
           CALL BD_CopyInput(BDInput_tmp,BD_Input,MESH_NEWCOPY,ErrStat,ErrMsg)
       ENDDO

       DO j=1,3
           BD1_Input%PointLoad%Force(j,BD1_Parameter%node_total) = BD1_Input%PointLoad%Force(j,BD1_Parameter%node_total) + eps
           CALL BD_CalcOutput( time, BD1_Input, BD1_Parameter, BD1_ContinuousState, BD1_DiscreteState, &
                    BD1_ConstraintState, BD1_OtherState, BD1OT_tmp, ErrStat, ErrMsg )    
           DO k=1,3
               Coef(6+k,j) = -((BD1OT_tmp%BldMotion%TranslationDisp(k,BD1_Parameter%node_total) - &
                              BD1_Output%BldMotion%TranslationDisp(k,BD1_Parameter%node_total))/eps)
               Coef(12+k,j) = -((BD1OT_tmp%BldMotion%TranslationVel(k,BD1_Parameter%node_total) - &
                              BD1_Output%BldMotion%TranslationVel(k,BD1_Parameter%node_total))/eps)
               Coef(15+k,j) = -((BD1OT_tmp%BldMotion%RotationVel(k,BD1_Parameter%node_total) - &
                              BD1_Output%BldMotion%RotationVel(k,BD1_Parameter%node_total))/eps)
               Coef(18+k,j) = -((BD1OT_tmp%BldMotion%TranslationAcc(k,BD1_Parameter%node_total) - &
                              BD1_Output%BldMotion%TranslationAcc(k,BD1_Parameter%node_total))/eps)
               Coef(21+k,j) = -((BD1OT_tmp%BldMotion%RotationAcc(k,BD1_Parameter%node_total) - &
                              BD1_Output%BldMotion%RotationAcc(k,BD1_Parameter%node_total))/eps)
           ENDDO
           CALL BD_CrvExtractCrv(BD1OT_tmp%BldMotion%Orientation(1:3,1:3,BD1_Parameter%node_total),temp_cc)
           temp_c0 = MATMUL(BD1_Parameter%GlbRot(1:3,1:3),BD1_Parameter%uuN0(4:6,1))
           CALL BD_CrvCompose(temp3,temp_cc,temp_c0,2)
           DO k=1,3
               Coef(9+k,j) = -((temp3(k) - tempBD1_rr(k))/eps)
           ENDDO
           CALL BD_CopyInput(BD1Input_tmp,BD1_Input,MESH_NEWCOPY,ErrStat,ErrMsg)
       ENDDO
       DO j=1,3
           BD1_Input%PointLoad%Moment(j,BD1_Parameter%node_total) = BD1_Input%PointLoad%Moment(j,BD1_Parameter%node_total) + eps
           CALL BD_CalcOutput( time, BD1_Input, BD1_Parameter, BD1_ContinuousState, BD1_DiscreteState, &
                    BD1_ConstraintState, BD1_OtherState, BD1OT_tmp, ErrStat, ErrMsg )    
           DO k=1,3
               Coef(6+k,j+3) = -((BD1OT_tmp%BldMotion%TranslationDisp(k,BD1_Parameter%node_total) - &
                              BD1_Output%BldMotion%TranslationDisp(k,BD1_Parameter%node_total))/eps)
               Coef(12+k,j+3) = -((BD1OT_tmp%BldMotion%TranslationVel(k,BD1_Parameter%node_total) - &
                              BD1_Output%BldMotion%TranslationVel(k,BD1_Parameter%node_total))/eps)
               Coef(15+k,j+3) = -((BD1OT_tmp%BldMotion%RotationVel(k,BD1_Parameter%node_total) - &
                              BD1_Output%BldMotion%RotationVel(k,BD1_Parameter%node_total))/eps)
               Coef(18+k,j+3) = -((BD1OT_tmp%BldMotion%TranslationAcc(k,BD1_Parameter%node_total) - &
                              BD1_Output%BldMotion%TranslationAcc(k,BD1_Parameter%node_total))/eps)
               Coef(21+k,j+3) = -((BD1OT_tmp%BldMotion%RotationAcc(k,BD1_Parameter%node_total) - &
                              BD1_Output%BldMotion%RotationAcc(k,BD1_Parameter%node_total))/eps)
           ENDDO
           CALL BD_CrvExtractCrv(BD1OT_tmp%BldMotion%Orientation(1:3,1:3,BD1_Parameter%node_total),temp_cc)
           temp_c0 = MATMUL(BD1_Parameter%GlbRot(1:3,1:3),BD1_Parameter%uuN0(4:6,1))
           CALL BD_CrvCompose(temp3,temp_cc,temp_c0,2) 
           DO k=1,3
               Coef(9+k,j+3) = -((temp3(k) - tempBD1_rr(k))/eps)
           ENDDO
           CALL BD_CopyInput(BD1Input_tmp,BD1_Input,MESH_NEWCOPY,ErrStat,ErrMsg)
       ENDDO
 
       CALL ludcmp(Coef,24,indx,d)
       CALL lubksb(Coef,24,indx,RHS,uinc)

!WRITE(*,*) 'uinc:'
!WRITE(*,*) uinc(:)
!WRITE(*,*) 'BD1 Input Force'
!WRITE(*,*) BD1_Input%PointLoad%Force(1:3,BD1_Parameter%node_total),BD1_Input%PointLoad%Moment(1:3,BD1_Parameter%node_total)
!WRITE(*,*) 'BD Input Disp'
!WRITE(*,*) BD_Input%RootMotion%TranslationDisp(1:3,1),tempBD_rr
!WRITE(*,*) 'BD Input Vel'
!WRITE(*,*) BD_Input%RootMotion%TranslationVel(1:3,1),BD_Input%RootMotion%RotationVel(1:3,1)
!WRITE(*,*) 'BD Input Acc'
!WRITE(*,*) BD_Input%RootMotion%TranslationAcc(1:3,1),BD_Input%RootMotion%RotationAcc(1:3,1)
       IF(BD_Norm(uinc) .LE. TOLF) THEN
           CALL BD_DestroyInput(BDInput_tmp, ErrStat, ErrMsg )
           CALL BD_DestroyInput(BD1Input_tmp, ErrStat, ErrMsg )
           CALL BD_DestroyOutput(OT_tmp, ErrStat, ErrMsg )
           CALL BD_DestroyOutput(BD1OT_tmp, ErrStat, ErrMsg )
           RETURN
       ENDIF

       BD1_Input%PointLoad%Force(1:3,BD1_Parameter%node_total) = &
            BD1_Input%PointLoad%Force(1:3,BD1_Parameter%node_total) + uinc(1:3)
       BD1_Input%PointLoad%Moment(1:3,BD1_Parameter%node_total) = &
            BD1_Input%PointLoad%Moment(1:3,BD1_Parameter%node_total) + uinc(4:6)
       BD_Input%RootMotion%TranslationDisp(1:3,1) = BD_Input%RootMotion%TranslationDisp(1:3,1) + uinc(7:9)
       BD_Input%RootMotion%TranslationVel(1:3,1) = BD_Input%RootMotion%TranslationVel(1:3,1) + uinc(13:15)
       BD_Input%RootMotion%RotationVel(1:3,1) = BD_Input%RootMotion%RotationVel(1:3,1) + uinc(16:18)
       BD_Input%RootMotion%TranslationAcc(1:3,1) = BD_Input%RootMotion%TranslationAcc(1:3,1) + uinc(19:21)
       BD_Input%RootMotion%RotationAcc(1:3,1) = BD_Input%RootMotion%RotationAcc(1:3,1) + uinc(22:24)
       CALL BD_CrvCompose(temp3,uinc(10:12),tempBD_rr,0)
       temp_c0(:) = MATMUL(BD_Parameter%GlbRot(1:3,1:3),BD_Parameter%uuN0(4:6,1))
       CALL BD_CrvCompose(temp3,temp3,temp_c0,0)
       CALL BD_CrvMatrixR(temp3,BD_Input%RootMotion%Orientation(1:3,1:3,1))
      
       IF(i .EQ. iter_max) THEN
           WRITE(*,*) "InputOutputSolve does not converge after the maximum number of iterations"
           CALL BD_DestroyInput(BDInput_tmp, ErrStat, ErrMsg )
           CALL BD_DestroyInput(BD1Input_tmp, ErrStat, ErrMsg )
           CALL BD_DestroyOutput(OT_tmp, ErrStat, ErrMsg )
           CALL BD_DestroyOutput(BD1OT_tmp, ErrStat, ErrMsg )
           STOP 
       ENDIF

   ENDDO


END SUBROUTINE BD1_BD_InputOutputSolve

END MODULE BD1_BD_MappingModule
!----------------------------------------------------------------------------------------------------------------------------------
PROGRAM MAIN

   USE BD1_BD_MappingModule

   USE BeamDyn
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

   INTEGER(IntKi)                     :: BD1_interp_order     ! order of interpolation/extrapolation
   INTEGER(IntKi)                     :: BD_interp_order     ! order of interpolation/extrapolation

   INTEGER(IntKi)                     :: MaxPtsInMap      ! the maximum number of points in a mapping
   
   ! Module1 Derived-types variables; see Registry_Module1.txt for details

   TYPE(BD_InitInputType)           :: BD1_InitInput
   TYPE(BD_ParameterType)           :: BD1_Parameter
   TYPE(BD_ContinuousStateType)     :: BD1_ContinuousState
   TYPE(BD_ContinuousStateType)     :: BD1_ContinuousStateDeriv
   TYPE(BD_InitOutputType)          :: BD1_InitOutput
   TYPE(BD_DiscreteStateType)       :: BD1_DiscreteState
   TYPE(BD_ConstraintStateType)     :: BD1_ConstraintState
   TYPE(BD_OtherStateType)          :: BD1_OtherState

   TYPE(BD_InputType),Dimension(:),ALLOCATABLE   :: BD1_Input
   REAL(DbKi) , DIMENSION(:), ALLOCATABLE         :: BD1_InputTimes

   TYPE(BD_OutputType),Dimension(:),ALLOCATABLE  :: BD1_Output
   REAL(DbKi) , DIMENSION(:), ALLOCATABLE         :: BD1_OutputTimes

   TYPE(BD_InputType)   :: u1    ! local variable for extrapolated inputs
   TYPE(BD_OutputType)  :: y1    ! local variable for extrapolated outputs

   ! Module 1 deived data typed needed in pc-coupling; predicted states

   TYPE(BD_ContinuousStateType)     :: BD1_ContinuousState_pred
   TYPE(BD_DiscreteStateType)       :: BD1_DiscreteState_pred
   TYPE(BD_ConstraintStateType)     :: BD1_ConstraintState_pred
   TYPE(BD_OtherStateType)          :: BD1_OtherState_pred

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
   INTEGER(IntKi),PARAMETER:: QiTipDisp = 20
   INTEGER(IntKi),PARAMETER:: QiMidDisp = 30
   INTEGER(IntKi),PARAMETER:: QiMidForce = 40
   INTEGER(IntKi),PARAMETER:: QiMidVel = 50
   INTEGER(IntKi),PARAMETER:: QiMidAcc = 60

   ! -------------------------------------------------------------------------
   ! MAPPING STUFF; Likely needs to be added to ModMesh
   ! -------------------------------------------------------------------------

!   TYPE(MeshMapType) :: Map_Mod2_P_Mod1_P
!   TYPE(MeshMapType) :: Map_Mod1_P_Mod2_P

   OPEN(unit = QiTipDisp, file = 'Qi_Tip_Disp.out', status = 'REPLACE',ACTION = 'WRITE')
   OPEN(unit = QiMidDisp, file = 'Qi_Mid_Disp.out', status = 'REPLACE',ACTION = 'WRITE')
   OPEN(unit = QiMidForce, file = 'Qi_Mid_Force.out', status = 'REPLACE',ACTION = 'WRITE')
   OPEN(unit = QiMidVel, file = 'Qi_Mid_Vel.out', status = 'REPLACE',ACTION = 'WRITE')
   OPEN(unit = QiMidAcc, file = 'Qi_Mid_Acc.out', status = 'REPLACE',ACTION = 'WRITE')
   ! -------------------------------------------------------------------------
   ! Initialization of glue-code time-step variables
   ! -------------------------------------------------------------------------

   t_initial = 0.d0
   t_final   = 8.0D-00

   pc_max = 2  ! Number of predictor-corrector iterations; 1 corresponds to an explicit calculation where UpdateStates
               ! is called only once  per time step for each module; inputs and outputs are extrapolated in time and
               ! are available to modules that have an implicit dependence on other-module data

   ! specify time increment; currently, all modules will be time integrated with this increment size
   dt_global = 1.0D-03

   n_t_final = ((t_final - t_initial) / dt_global )

   t_global = t_initial

   ! define polynomial-order for ModName_Input_ExtrapInterp and ModName_Output_ExtrapInterp
   ! Must be 0, 1, or 2
   BD1_interp_order = 2
   BD_interp_order   = 2

   !Module1: allocate Input and Output arrays; used for interpolation and extrapolation
   ALLOCATE(BD1_Input(BD1_interp_order + 1))
   ALLOCATE(BD1_InputTimes(BD1_interp_order + 1)) 
   ALLOCATE(BD1_Output(BD1_interp_order + 1))
   ALLOCATE(BD1_OutputTimes(BD1_interp_order + 1))

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

   BD1_InitInput%InputFile = 'GA2_Debug_BD1.inp'
   BD1_InitInput%RootName  = TRIM(BD1_Initinput%InputFile)
   ALLOCATE(BD1_InitInput%gravity(3)) 
   BD1_InitInput%gravity(1) = 0.0D0 !-9.80665
   BD1_InitInput%gravity(2) = 0.0D0 
   BD1_InitInput%gravity(3) = 0.0D0
   ALLOCATE(BD1_InitInput%GlbPos(3)) 
   BD1_InitInput%GlbPos(1) = 0.0D+00
   BD1_InitInput%GlbPos(2) = 0.0D+00
   BD1_InitInput%GlbPos(3) = 0.0D+00
   ALLOCATE(BD1_InitInput%GlbRot(3,3)) 
   BD1_InitInput%GlbRot(:,:) = 0.0D0
   temp_vec(1) = 0.0
   temp_vec(2) = 0.0
   temp_vec(3) = 0.0
   CALL BD_CrvMatrixR(temp_vec,temp_R)
   BD1_InitInput%GlbRot(1:3,1:3) = temp_R(1:3,1:3)
   ALLOCATE(BD1_InitInput%RootDisp(3)) 
   BD1_InitInput%RootDisp(1) = 0.0D+00
   BD1_InitInput%RootDisp(2) = 0.0D+00
   BD1_InitInput%RootDisp(3) = 0.0D+00
   ALLOCATE(BD1_InitInput%RootOri(3,3)) 
   BD1_InitInput%RootOri(:,:) = 0.0D0
   temp_vec(1) = 0.0
   temp_vec(2) = 0.0
   temp_vec(3) = 0.0
   CALL BD_CrvMatrixR(temp_vec,temp_R)
   BD1_InitInput%RootOri(1:3,1:3) = temp_R(1:3,1:3)
   ALLOCATE(BD1_InitInput%RootVel(6)) 
   BD1_InitInput%RootVel(:) = 0.0D+00

   CALL BD_Init( BD1_InitInput, BD1_Input(1), BD1_Parameter, BD1_ContinuousState, BD1_DiscreteState, &
                   BD1_ConstraintState, BD1_OtherState, BD1_Output(1), dt_global, BD1_InitOutput, ErrStat, ErrMsg )

   CALL BD_CopyInput(  BD1_Input(1), u1, MESH_NEWCOPY, ErrStat, ErrMsg )
   CALL BD_CopyOutput( BD1_Output(1), y1, MESH_NEWCOPY, ErrStat, ErrMsg )


   ! We fill Mod1_InputTimes with negative times, but the Mod1_Input values are identical for each of those times; this allows
   ! us to use, e.g., quadratic interpolation that effectively acts as a zeroth-order extrapolation and first-order extrapolation
   ! for the first and second time steps.  (The interpolation order in the ExtrapInput routines are determined as
   ! order = SIZE(Mod1_Input)
   DO i = 1, BD1_interp_order + 1
      BD1_InputTimes(i) = t_initial - (i - 1) * dt_global
      BD1_OutputTimes(i) = t_initial - (i - 1) * dt_global
   ENDDO

!WRITE(*,*) 'Mod1_InputTimes:',Mod1_InputTimes(:)

   DO i = 1, BD1_interp_order
     CALL BD_CopyInput (BD1_Input(i),  BD1_Input(i+1),  MESH_NEWCOPY, Errstat, ErrMsg)
     CALL BD_CopyOutput (BD1_Output(i),  BD1_Output(i+1),  MESH_NEWCOPY, Errstat, ErrMsg)
   ENDDO 

!   BD_InitInput%InputFile = 'Siemens_53_Input.inp'
   BD_InitInput%InputFile = 'GA2_Debug_BD.inp'
   BD_InitInput%RootName  = TRIM(BD_Initinput%InputFile)
   ALLOCATE(BD_InitInput%gravity(3)) 
   BD_InitInput%gravity(1) = 0.0D0 !-9.80665
   BD_InitInput%gravity(2) = 0.0D0 
   BD_InitInput%gravity(3) = 0.0D0
   ALLOCATE(BD_InitInput%GlbPos(3)) 
   BD_InitInput%GlbPos(1) = 5.0D+00
   BD_InitInput%GlbPos(2) = 0.0D+00
   BD_InitInput%GlbPos(3) = 0.0D+00
   ALLOCATE(BD_InitInput%GlbRot(3,3)) 
   BD_InitInput%GlbRot(:,:) = 0.0D0
   temp_vec(1) = 0.0
   temp_vec(2) = 0.0
   temp_vec(3) = 0.0 !4.0D0*TAN((3.1415926D0/2.0D0)/4.0D0)
   CALL BD_CrvMatrixR(temp_vec,temp_R)
   BD_InitInput%GlbRot(1:3,1:3) = temp_R(1:3,1:3)
   ALLOCATE(BD_InitInput%RootDisp(3)) 
   BD_InitInput%RootDisp(1) = 0.0D+00
   BD_InitInput%RootDisp(2) = 0.0D+00
   BD_InitInput%RootDisp(3) = 0.0D+00
   ALLOCATE(BD_InitInput%RootOri(3,3)) 
   BD_InitInput%RootOri(:,:) = 0.0D0
   temp_vec(1) = 0.0
   temp_vec(2) = 0.0
   temp_vec(3) = 0.0
   CALL BD_CrvMatrixR(temp_vec,temp_R)
   BD_InitInput%RootOri(1:3,1:3) = temp_R(1:3,1:3)
   ALLOCATE(BD_InitInput%RootVel(6)) 
   BD_InitInput%RootVel(:) = 0.0D+00

   CALL BD_Init( BD_InitInput, BD_Input(1), BD_Parameter, BD_ContinuousState, BD_DiscreteState, &
                 BD_ConstraintState, BD_OtherState, BD_Output(1), dt_global, BD_InitOutput, ErrStat, ErrMsg )
   

   DO i = 1, BD_interp_order + 1
      BD_InputTimes(i) = t_initial - (i - 1) * dt_global
      BD_OutputTimes(i) = t_initial - (i - 1) * dt_global
   ENDDO

   DO i = 1, BD_interp_order
     CALL BD_CopyInput (BD_Input(i),  BD_Input(i+1),  MESH_NEWCOPY, Errstat, ErrMsg)
     CALL BD_CopyOutput (BD_Output(i),  BD_Output(i+1), MESH_NEWCOPY, Errstat, ErrMsg)
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
   CALL BD1_BD_InputOutputSolve(t_global, &
                   BD1_Input(1), BD1_Parameter, BD1_ContinuousState, BD1_DiscreteState, &
                   BD1_ConstraintState, BD1_OtherState, BD1_Output(1), &
                   BD_Input(1), BD_Parameter, BD_ContinuousState, BD_DiscreteState, &
                   BD_ConstraintState, BD_OtherState, BD_Output(1),  &
!                   Map_Mod1_P_Mod2_P, Map_Mod2_P_Mod1_P, &
                   ErrStat, ErrMsg)

   DO n_t_global = 0, n_t_final
WRITE(*,*) "Time Step: ", n_t_global
!IF(n_t_global .EQ. 1) STOP
      ! Solve input-output relations; this section of code corresponds to Eq. (35) in Gasmi et al. (2013)
      ! This code will be specific to the underlying modules

      CALL BD_CrvExtractCrv(BD_OutPut(1)%BldMotion%Orientation(1:3,1:3,BD_Parameter%node_total),temp_cc)
      WRITE(QiTipDisp,6000) t_global,&
                            &BD_OutPut(1)%BldMotion%TranslationDisp(1:3,BD_Parameter%node_total),&
                            &temp_cc(1:3)
      CALL BD_CrvExtractCrv(BD1_OutPut(1)%BldMotion%Orientation(1:3,1:3,BD1_Parameter%node_total),temp_cc)
      WRITE(QiMidDisp,6000) t_global,&
                            &BD1_OutPut(1)%BldMotion%TranslationDisp(1:3,BD1_Parameter%node_total),&
                            &temp_cc(1:3)
      WRITE(QiMidForce,6000) t_global,&
                             &BD_OutPut(1)%ReactionForce%Force(1:3,1),&
                             &BD_OutPut(1)%ReactionForce%Moment(1:3,1)
      WRITE(QiMidVel,6000) t_global,&
                             &BD1_OutPut(1)%BldMotion%TranslationVel(1:3,BD1_Parameter%node_total),&
                             &BD1_OutPut(1)%BldMotion%RotationVel(1:3,BD1_Parameter%node_total)
      WRITE(QiMidAcc,6000) t_global,&
                             &BD1_OutPut(1)%BldMotion%TranslationAcc(1:3,BD1_Parameter%node_total),&
                             &BD1_OutPut(1)%BldMotion%RotationAcc(1:3,BD1_Parameter%node_total)
!WRITE(*,*) 'BD_Input%Disp:',BD_Input(1)%RootMotion%TranslationDisp(:,1)
!WRITE(*,*) 'BD_Input%Disp:',BD_Input(2)%RootMotion%TranslationDisp(:,1)
!WRITE(*,*) 'BD_Input%Disp:',BD_Input(3)%RootMotion%TranslationDisp(:,1)
!WRITE(*,*) 'BD_Input%Velo:',BD_Input(1)%RootMotion%TranslationVel(:,1)
!WRITE(*,*) 'BD_Input%Velo:',BD_Input(2)%RootMotion%TranslationVel(:,1)
!WRITE(*,*) 'BD_Input%Velo:',BD_Input(3)%RootMotion%TranslationVel(:,1)
!WRITE(*,*) 'BD_Input%Acce:',BD_Input(1)%RootMotion%TranslationAcc(:,1)
!WRITE(*,*) 'BD_Input%Acce:',BD_Input(2)%RootMotion%TranslationAcc(:,1)
!WRITE(*,*) 'BD_Input%Acce:',BD_Input(3)%RootMotion%TranslationAcc(:,1)
!WRITE(*,*) 'Mod1_Input%Force:',Mod1_Input(1)%PointMesh%Force(:,1)
!WRITE(*,*) 'Mod1_Input%Force:',Mod1_Input(2)%PointMesh%Force(:,1)
!WRITE(*,*) 'Mod1_Input%Force:',Mod1_Input(3)%PointMesh%Force(:,1)

      ! after all InputOutputSolves, we can reset the mapping flags on the meshes:
         BD1_Input(1)%RootMotion%RemapFlag    = .FALSE. 
         BD1_Output(1)%ReactionForce%RemapFlag   = .FALSE.
         BD_Input(1)%RootMotion%RemapFlag      = .FALSE. 
         BD_Output(1)%ReactionForce%RemapFlag = .FALSE.

      ! extrapolate inputs and outputs to t + dt; will only be used by modules with an implicit dependence on input data.

      CALL BD_Input_ExtrapInterp(BD1_Input, BD1_InputTimes, u1, t_global + dt_global, ErrStat, ErrMsg)

!WRITE(*,*) 'u1%Force:',u1%PointMesh%Force(:,1)

      CALL BD_Output_ExtrapInterp(BD1_Output, BD1_OutputTimes, y1, t_global + dt_global, ErrStat, ErrMsg)

      ! Shift "window" of the Mod1_Input and Mod1_Output

      DO i = BD1_interp_order, 1, -1
         CALL BD_CopyInput (BD1_Input(i),  BD1_Input(i+1),  MESH_UPDATECOPY, Errstat, ErrMsg)
         CALL BD_CopyOutput (BD1_Output(i),  BD1_Output(i+1),  MESH_UPDATECOPY, Errstat, ErrMsg)
         BD1_InputTimes(i+1) = BD1_InputTimes(i)
         BD1_OutputTimes(i+1) = BD1_OutputTimes(i)
      ENDDO
!WRITE(*,*) 'Mod1_InputTimes:',Mod1_InputTimes(:)

      CALL BD_CopyInput (u1,  BD1_Input(1),  MESH_UPDATECOPY, Errstat, ErrMsg)
      CALL BD_CopyOutput (y1,  BD1_Output(1),  MESH_UPDATECOPY, Errstat, ErrMsg)
      BD1_InputTimes(1) = t_global + dt_global
      BD1_OutputTimes(1) = t_global + dt_global

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

      DO pc = 1, pc_max

         !----------------------------------------------------------------------------------------
         ! Module 1
         !----------------------------------------------------------------------------------------

         ! copy ContinuousState to ContinuousState_pred; don't modify ContinuousState until completion of PC iterations

         CALL BD_CopyContState   (BD1_ContinuousState, BD1_ContinuousState_pred, 0, Errstat, ErrMsg)

         CALL BD_CopyConstrState (BD1_ConstraintState, BD1_ConstraintState_pred, 0, Errstat, ErrMsg)

         CALL BD_CopyDiscState   (BD1_DiscreteState,   BD1_DiscreteState_pred,   0, Errstat, ErrMsg)

         CALL BD_UpdateStates( t_global, n_t_global, BD1_Input, BD1_InputTimes, BD1_Parameter, BD1_ContinuousState_pred, &
                                 BD1_DiscreteState_pred, BD1_ConstraintState_pred, &
                                 BD1_OtherState, ErrStat, ErrMsg )
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

         IF (pc .LE. pc_max) THEN

            CALL BD1_BD_InputOutputSolve( t_global + dt_global, &
                                             BD1_Input(1), BD1_Parameter, BD1_ContinuousState_pred, BD1_DiscreteState_pred, &
                                             BD1_ConstraintState_pred, BD1_OtherState, BD1_Output(1), &
                                             BD_Input(1), BD_Parameter, BD_ContinuousState_pred, BD_DiscreteState_pred, &
                                             BD_ConstraintState_pred, BD_OtherState_pred, BD_Output(1),  &
!                                             Map_Mod1_P_Mod2_P, Map_Mod2_P_Mod1_P, &
                                             ErrStat, ErrMsg)

            
            ! after all InputOutputSolves, we can reset the mapping flags on the meshes:
            BD1_Input(1)%RootMotion%RemapFlag  = .FALSE. 
            BD1_Output(1)%ReactionForce%RemapFlag = .FALSE.
            BD_Input(1)%RootMotion%RemapFlag      = .FALSE. 
            BD_Output(1)%ReactionForce%RemapFlag = .FALSE.
            
         ENDIF 

      ENDDO

      ! Save all final variables

      CALL BD_CopyContState   (BD1_ContinuousState_pred,  BD1_ContinuousState, 0, Errstat, ErrMsg)
      CALL BD_CopyConstrState (BD1_ConstraintState_pred,  BD1_ConstraintState, 0, Errstat, ErrMsg)
      CALL BD_CopyDiscState   (BD1_DiscreteState_pred,    BD1_DiscreteState,   0, Errstat, ErrMsg)
      CALL BD_CopyOtherState  (BD1_OtherState_pred,   BD1_OtherState,   0, Errstat, ErrMsg)

      CALL BD_CopyContState   (BD_ContinuousState_pred, BD_ContinuousState, 0, Errstat, ErrMsg)
      CALL BD_CopyConstrState (BD_ConstraintState_pred, BD_ConstraintState, 0, Errstat, ErrMsg)
      CALL BD_CopyDiscState   (BD_DiscreteState_pred,   BD_DiscreteState,   0, Errstat, ErrMsg)
      CALL BD_CopyOtherState  (BD_OtherState_pred,   BD_OtherState,   0, Errstat, ErrMsg)

      ! update the global time

      t_global = REAL(n_t_global+1,DbKi) * dt_global + t_initial

      ! the following is exact solution for q_1(t) for baseline parameters in Gasmi et al. (2013)

      ! build rms_error calculation components; see Eq. (56) in Gasmi et al. (2013)

      ! print discrete q_1(t) solution to standard out


   ENDDO


   ! calculate final time normalized rms error


   ! -------------------------------------------------------------------------
   ! Ending of modules
   ! -------------------------------------------------------------------------


   CALL BD_End(  BD1_Input(1), BD1_Parameter, BD1_ContinuousState, BD1_DiscreteState, &
                   BD1_ConstraintState, BD1_OtherState, BD1_Output(1), ErrStat, ErrMsg )

   DO i = 2, BD1_interp_order+1
      CALL BD_DestroyInput(BD1_Input(i), ErrStat, ErrMsg )
      CALL BD_DestroyOutput(BD1_Output(i), ErrStat, ErrMsg )
   ENDDO
   DEALLOCATE(BD1_InputTimes)
   DEALLOCATE(BD1_OutputTimes)

   CALL BD_End(  BD_Input(1), BD_Parameter, BD_ContinuousState, BD_DiscreteState, &
                 BD_ConstraintState, BD_OtherState, BD_Output(1), ErrStat, ErrMsg )

   do i = 2, BD_interp_order+1
      CALL BD_DestroyInput(BD_Input(i), ErrStat, ErrMsg )
      CALL BD_DestroyOutput(BD_Output(i), ErrStat, ErrMsg )
   enddo

   DEALLOCATE(BD_InputTimes)
   DEALLOCATE(BD_Output)

   6000 FORMAT (ES12.5,6ES21.12)
   CLOSE (QiTipDisp)
   CLOSE (QiMidDisp)
   ! -------------------------------------------------------------------------
   ! Deallocate arrays associated with mesh mapping
   ! -------------------------------------------------------------------------

!   CALL MeshMapDestroy(Map_Mod1_P_Mod2_P, ErrStat, ErrMsg)
!   CALL MeshMapDestroy(Map_Mod2_P_Mod1_P, ErrStat, ErrMsg)

END PROGRAM MAIN

