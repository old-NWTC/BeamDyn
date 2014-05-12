   SUBROUTINE BeamDyn_RK4(t,n,u,utimes,p,x,xd,z,OtherState,ErrStat,ErrMsg)
!
! This subroutine implements the fourth-order Runge-Kutta Method (RK4) for numerically integrating ordinary differential equations:
!
!   Let f(t, x) = xdot denote the time (t) derivative of the continuous states (x). 
!   Define constants k1, k2, k3, and k4 as 
!        k1 = dt * f(t        , x_t        )
!        k2 = dt * f(t + dt/2 , x_t + k1/2 )
!        k3 = dt * f(t + dt/2 , x_t + k2/2 ), and
!        k4 = dt * f(t + dt   , x_t + k3   ).
!   Then the continuous states at t = t + dt are
!        x_(t+dt) = x_t + k1/6 + k2/3 + k3/3 + k4/6 + O(dt^5)
!
! For details, see:
! Press, W. H.; Flannery, B. P.; Teukolsky, S. A.; and Vetterling, W. T. "Runge-Kutta Method" and "Adaptive Step Size Control for 
!   Runge-Kutta." ß16.1 and 16.2 in Numerical Recipes in FORTRAN: The Art of Scientific Computing, 2nd ed. Cambridge, England: 
!   Cambridge University Press, pp. 704-716, 1992.
!
!..................................................................................................................................

   REAL(DbKi),                        INTENT(IN   )  :: t           ! Current simulation time in seconds
   INTEGER(IntKi),                    INTENT(IN   )  :: n           ! time step number
   TYPE(BD_InputType),              INTENT(INOUT)  :: u(:)        ! Inputs at t
   REAL(DbKi),                        INTENT(IN   )  :: utimes(:)   ! times of input
   TYPE(BD_ParameterType),          INTENT(IN   )  :: p           ! Parameters
   TYPE(BD_ContinuousStateType),    INTENT(INOUT)  :: x           ! Continuous states at t on input at t + dt on output
   TYPE(BD_DiscreteStateType),      INTENT(IN   )  :: xd          ! Discrete states at t
   TYPE(BD_ConstraintStateType),    INTENT(IN   )  :: z           ! Constraint states at t (possibly a guess)
   TYPE(BD_OtherStateType),         INTENT(INOUT)  :: OtherState  ! Other/optimization states
   INTEGER(IntKi),                    INTENT(  OUT)  :: ErrStat     ! Error status of the operation
   CHARACTER(*),                      INTENT(  OUT)  :: ErrMsg      ! Error message if ErrStat /= ErrID_None

   ! local variables
      
   TYPE(BD_ContinuousStateType)                 :: xdot        ! time derivatives of continuous states      
   TYPE(BD_ContinuousStateType)                 :: k1          ! RK4 constant; see above
   TYPE(BD_ContinuousStateType)                 :: k2          ! RK4 constant; see above 
   TYPE(BD_ContinuousStateType)                 :: k3          ! RK4 constant; see above 
   TYPE(BD_ContinuousStateType)                 :: k4          ! RK4 constant; see above 
   TYPE(BD_ContinuousStateType)                 :: x_tmp       ! Holds temporary modification to x
   TYPE(BD_InputType)                           :: u_interp    ! interpolated value of inputs 

   ! Initialize ErrStat

   ErrStat = ErrID_None
   ErrMsg  = "" 

   CALL MeshCopy ( SrcMesh  = u(1)%PointMesh      &
                 , DestMesh = u_interp%PointMesh  &
                 , CtrlCode = MESH_NEWCOPY        &
                 , ErrStat  = ErrStat             &
                 , ErrMess  = ErrMsg               )


   ! interpolate u to find u_interp = u(t)
   CALL BD_Input_ExtrapInterp( u, utimes, u_interp, t, ErrStat, ErrMsg )

   ! find xdot at t
   CALL BeamDyn_CalcContStateDeriv( t, u_interp, p, x, xd, z, OtherState, xdot, ErrStat, ErrMsg )

   k1%q    = p%dt * xdot%q
   k1%dqdt = p%dt * xdot%dqdt
  
   x_tmp%q    = x%q    + 0.5 * k1%q
   x_tmp%dqdt = x%dqdt + 0.5 * k1%dqdt

   ! interpolate u to find u_interp = u(t + dt/2)
   CALL BD_Input_ExtrapInterp(u, utimes, u_interp, t+0.5*p%dt, ErrStat, ErrMsg)

   ! find xdot at t + dt/2
   CALL BeamDyn_CalcContStateDeriv( t + 0.5*p%dt, u_interp, p, x_tmp, xd, z, OtherState, xdot, ErrStat, ErrMsg )

   k2%q    = p%dt * xdot%q
   k2%dqdt = p%dt * xdot%dqdt

   x_tmp%q    = x%q    + 0.5 * k2%q
   x_tmp%dqdt = x%dqdt + 0.5 * k2%dqdt

   ! find xdot at t + dt/2
   CALL BeamDyn_CalcContStateDeriv( t + 0.5*p%dt, u_interp, p, x_tmp, xd, z, OtherState, xdot, ErrStat, ErrMsg )
     
   k3%q    = p%dt * xdot%q
   k3%dqdt = p%dt * xdot%dqdt

   x_tmp%q    = x%q    + k3%q
   x_tmp%dqdt = x%dqdt + k3%dqdt

   ! interpolate u to find u_interp = u(t + dt)
   CALL BD_Input_ExtrapInterp(u, utimes, u_interp, t + p%dt, ErrStat, ErrMsg)

   ! find xdot at t + dt
   CALL BeamDyn_CalcContStateDeriv( t + p%dt, u_interp, p, x_tmp, xd, z, OtherState, xdot, ErrStat, ErrMsg )

   k4%q    = p%dt * xdot%q
   k4%dqdt = p%dt * xdot%dqdt

   x%q    = x%q    +  ( k1%q    + 2. * k2%q    + 2. * k3%q    + k4%q    ) / 6.      
   x%dqdt = x%dqdt +  ( k1%dqdt + 2. * k2%dqdt + 2. * k3%dqdt + k4%dqdt ) / 6.      

   CALL MeshDestroy ( u_interp%PointMesh       &
                    , ErrStat  = ErrStat         &
                    , ErrMess  = ErrMsg           )

   END SUBROUTINE BeamDyn_RK4