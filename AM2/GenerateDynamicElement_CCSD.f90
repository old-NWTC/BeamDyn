   SUBROUTINE GenerateDynamicElement_CCSD(uuN0,uuN,vvN,Stif0,Mass0,gravity,u,&
                                        &damp_flag,beta,&
                                        &elem_total,node_elem,dof_node,ngp,RHS,MassM)
   !----------------------------------------------------------------------------------------
   ! This subroutine computes Global mass matrix and force vector for the beam.
   !----------------------------------------------------------------------------------------
   REAL(ReKi),        INTENT(IN   ):: uuN0(:,:) ! Initial position vector
   REAL(ReKi),        INTENT(IN   ):: uuN(:) ! Displacement of Mass 1: m
   REAL(ReKi),        INTENT(IN   ):: vvN(:) ! Velocity of Mass 1: m/s
   REAL(ReKi),        INTENT(IN   ):: Stif0(:,:,:) ! Element stiffness matrix
   REAL(ReKi),        INTENT(IN   ):: Mass0(:,:,:) ! Element stiffness matrix
   REAL(ReKi),        INTENT(IN   ):: gravity(:) ! Velocity of Mass 1: m/s
   TYPE(BD_InputType),INTENT(IN   ):: u           ! Inputs at t
   INTEGER(IntKi),    INTENT(IN   ):: damp_flag ! Total number of elements
   REAL(ReKi),        INTENT(IN   ):: beta(:)
   INTEGER(IntKi),    INTENT(IN   ):: elem_total ! Total number of elements
   INTEGER(IntKi),    INTENT(IN   ):: node_elem ! Node per element
   INTEGER(IntKi),    INTENT(IN   ):: dof_node ! Degrees of freedom per node
   INTEGER(IntKi),    INTENT(IN   ):: ngp ! Number of Gauss points
   REAL(ReKi),        INTENT(INOUT):: MassM(:,:) ! Mass matrix 
   REAL(ReKi),        INTENT(INOUT):: RHS(:) ! Right hand side of the equation Ax=B  

   REAL(ReKi) :: Nuu0(dof_node*node_elem) ! Nodal initial position for each element
   REAL(ReKi) :: Nuuu(dof_node*node_elem) ! Nodal displacement of Mass 1 for each element
   REAL(ReKi) :: Nrr0(3*node_elem) ! Nodal rotation parameters for initial position 
   REAL(ReKi) :: Nrrr(3*node_elem) ! Nodal rotation parameters for displacement of Mass 1
   REAL(ReKi) :: Nvvv(dof_node*node_elem) ! Nodal velocity of Mass 1: m/s for each element
   REAL(ReKi) :: EStif0_GL(6,6,node_elem-1) ! Nodal material properties for each element
   REAL(ReKi) :: EMass0_GL(6,6,node_elem-1) ! Nodal material properties for each element
   REAL(ReKi) :: DistrLoad_GL(6,node_elem-1) ! Nodal material properties for each element
   REAL(ReKi) :: elf1(dof_node*node_elem) ! Total element force (Fc, Fd, Fb)
   REAL(ReKi) :: elf2(dof_node*node_elem) ! Total element force (Fc, Fd, Fb)
   REAL(ReKi) :: elm11(dof_node*node_elem,dof_node*node_elem) ! Element mass matrix
   REAL(ReKi) :: elm12(dof_node*node_elem,dof_node*node_elem) ! Element mass matrix
   REAL(ReKi) :: elm21(dof_node*node_elem,dof_node*node_elem) ! Element mass matrix
   REAL(ReKi) :: elm22(dof_node*node_elem,dof_node*node_elem) ! Element mass matrix

   INTEGER(IntKi)                  :: dof_elem ! Degree of freedom per node
   INTEGER(IntKi)                  :: rot_elem ! Rotational degrees of freedom
   INTEGER(IntKi)                  :: nelem ! number of elements
   INTEGER(IntKi)                  :: j ! Index counter
   INTEGER(IntKi)                  :: temp_id ! Index counter
   INTEGER(IntKi)                  :: allo_stat ! Allows for an error code return

   dof_elem = dof_node * node_elem
   rot_elem = (dof_node/2) * node_elem


   DO nelem=1,elem_total
       Nuu0(:) = uuN0(:,nelem)
       CALL ElemNodalDispGL(uuN,node_elem,dof_node,nelem,Nuuu)
       temp_id = (nelem-1)*ngp
       DO j=1,ngp
           EStif0_GL(1:6,1:6,j) = Stif0(1:6,1:6,temp_id+j)
           EMass0_GL(1:6,1:6,j) = Mass0(1:6,1:6,temp_id+j)
           DistrLoad_GL(1:3,j)  = u%DistrLoad%Force(1:3,temp_id+j+1)
           DistrLoad_GL(4:6,j)  = u%DistrLoad%Moment(1:3,temp_id+j+1)
       ENDDO
       
       CALL NodalRelRotGL(Nuu0,node_elem,dof_node,Nrr0)
       CALL NodalRelRotGL(Nuuu,node_elem,dof_node,Nrrr)
       CALL ElemNodalDispGL(vvN,node_elem,dof_node,nelem,Nvvv)

       elf1(:) = 0.0D0
       elf2(:) = 0.0D0
       elm11(:,:) = 0.0D0
       elm12(:,:) = 0.0D0
       elm21(:,:) = 0.0D0
       elm22(:,:) = 0.0D0
       CALL ElementMatrix_CCSD(Nuu0,Nuuu,Nrr0,Nrrr,Nvvv,&
                              &EStif0_GL,EMass0_GL,gravity,DistrLoad_GL,&
                              &ngp,node_elem,dof_node,damp_flag,beta,&
                              &elf1,elf2,elm11,elm22)
!DO j=1,18
!WRITE(*,*) j,elm11(1,j)
!ENDDO


       CALL AssembleStiffK_AM2(nelem,node_elem,dof_elem,dof_node,&
                              &elm11,elm12,elm21,elm22,MassM)
       CALL AssembleRHS_AM2(nelem,dof_elem,node_elem,dof_node,elf1,elf2,RHS)

   ENDDO



   END SUBROUTINE GenerateDynamicElement_CCSD