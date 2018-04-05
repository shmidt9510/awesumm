MODULE user_case
  ! CASE 2/3 D channel periodic in x-direction and if dim=3 in z-derection
  USE precision
  USE wlt_vars
  USE wavelet_filters_mod
  USE elliptic_mod
  USE elliptic_vars
  USE wlt_trns_mod
  USE wlt_trns_vars
  USE wlt_trns_util_mod
  USE io_3d_vars
  USE util_mod
  USE util_vars
  USE share_consts
  USE pde
  USE variable_mapping 
  USE sizes
  USE share_kry
  USE vector_util_mod
  USE field
  USE input_file_reader
  USE debug_vars
  !
  ! case specific variables
  !
  INTEGER n_var_pressure  ! start of pressure in u array
  INTEGER n_var_temp      ! start of temperature 

  REAL (pr) :: nu     ! kinematic viscosity (1/Re)heat conduction
  REAL (pr) :: T_0    ! Temperature of the top wall + initial temperature
  REAL (pr) :: T_1    ! Temperature of the lower wall
  INTEGER, PARAMETER :: Tdim = 1
CONTAINS

  !
  ! In user_setup_pde() we setup how many variables are integrated, which are interpolated
  ! to the next times step and if any exeact solution exists to compare to during the run.
  ! We also set the variable names that are used when the result files are written.
  !
  ! The following variables must be setup in this routine:
  !
  !
  ! n_integrated     ! first n_integrated eqns will be acted on for time integration
  ! n_var_additional ! 
  ! n_var
  !
  !
  !
  !
  SUBROUTINE  user_setup_pde ( VERB ) 
    USE variable_mapping
    IMPLICIT NONE
    LOGICAL, OPTIONAL :: VERB         ! print debug info
    INTEGER :: i

    IF (verb_level.GT.0) THEN
       PRINT * ,''
       PRINT *, '**********************Setting up PDE*****************'
       PRINT * ,'CASE: incompressible Poiseulle - Couette flow channel flow '
       PRINT *, '*****************************************************'
    END IF

    !register_var( var_name, integrated, adapt, interpolate, exact, saved, req_restart, FRONT_LOAD )
    CALL register_var( 'T ',       integrated=.TRUE.,  adapt=(/.TRUE.,.TRUE./),   interpolate=(/.TRUE.,.TRUE./),   exact=(/.TRUE.,.TRUE./),   saved=.TRUE., req_restart=.TRUE., IC_MAP_NAME='T_old'  )
    CALL register_var( 'pressure ', integrated=.FALSE., adapt=(/.FALSE.,.FALSE./), interpolate=(/.FALSE.,.FALSE./), exact=(/.FALSE.,.FALSE./), saved=.TRUE., req_restart=.FALSE., IC_MAP_NAME='pressure_old' )
!    CALL register_var( 'T_old ',       integrated=.TRUE.,  adapt=(/.TRUE.,.TRUE./),   interpolate=(/.TRUE.,.TRUE./),   exact=(/.TRUE.,.TRUE./),   saved=.TRUE., req_restart=.TRUE.)
!    CALL register_var( 'pressure_old ', integrated=.FALSE., adapt=(/.FALSE.,.FALSE./), interpolate=(/.FALSE.,.FALSE./), exact=(/.FALSE.,.FALSE./), saved=.TRUE., req_restart=.FALSE.)
    CALL setup_mapping()
    CALL print_variable_registery( FULL=.TRUE.)

    n_var_pressure  = get_index('pressure ')  
    n_var_temp      = get_index('  T    ')  

    ALLOCATE ( Umn(1:n_var) )
    Umn = 0.0_pr !set up here if mean quantities are not zero and used in scales or equation

    IF (verb_level.GT.0) THEN
       PRINT *, 'n_integrated = ',n_integrated 
       PRINT *, 'n_var = ',n_var 
       PRINT *, 'n_var_exact = ',n_var_exact 
       PRINT *, '*******************Variable Names*******************'
       DO i = 1,n_var
          WRITE (*, u_variable_names_fmt) u_variable_names(i)
       END DO
       PRINT *, '****************************************************'
    END IF

  END SUBROUTINE  user_setup_pde
  !
  ! Set the exact solution for comparison to the simulated solution
  !
  ! u          - array to fill in the exact solution
  ! nlocal       - number of active wavelets   
  ! ne_local        - total number of equations
  ! t          - time of current time step 
  ! l_n_var_exact_soln_index - index into the elements of u for which we need to 
  !                            find the exact solution
  SUBROUTINE  user_exact_soln (u, nlocal,  t_local, l_n_var_exact_soln)  
    IMPLICIT NONE
    INTEGER , INTENT (IN) :: nlocal
    REAL (pr), INTENT (IN) ::  t_local
    REAL (pr), DIMENSION (nlocal,n_var_exact), INTENT (INOUT) :: u
    LOGICAL , INTENT (IN) :: l_n_var_exact_soln(n_var)
    REAL (pr) :: t_zero
    INTEGER :: i
    DOUBLE PRECISION :: DERF
    EXTERNAL :: DERF

    t_zero = 0.0_pr
    IF(t_local == 0) t_zero = 1.0e-12_pr
    DO i =1, nlocal
       u(i,n_var_temp) = T_1 - (T_1 - T_0)*DERF( (x(i,Tdim) - xyzlimits(1,Tdim))/SQRT(4.0_pr*nu*(t_local+t_zero)) ) &
            /DERF( (xyzlimits(2,Tdim) - xyzlimits(1,Tdim))/SQRT(4.0_pr*nu*t_local) )
    END DO

    !WRITE( *,'( "in user_exact_soln t_local, nwlt " , f30.20, 1x , i8.8 )' ) t_local , nlocal
   

  END SUBROUTINE  user_exact_soln

  SUBROUTINE user_initial_conditions (u, nlocal, ne_local, t_local, scl, scl_fltwt, iter)
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: nlocal, ne_local
    INTEGER  , INTENT (INOUT) :: iter ! iteration of call while adapting initial grid
    REAL (pr), DIMENSION (nlocal,ne_local), INTENT (INOUT) :: u
    REAL (pr)  :: scl(1:n_var),scl_fltwt
    REAL (pr), INTENT (IN) :: t_local
    REAL (pr) :: t_zero
    INTEGER :: i
    DOUBLE PRECISION :: DERF
    EXTERNAL :: DERF

    t_zero = 0.0_pr
    IF(t_local == 0) t_zero = 1.0e-12_pr

    IF ( IC_restart_mode .EQ. 0 ) THEN
       DO i =1, nlocal
          u(i,n_var_temp) = T_1 - (T_1 - T_0)*DERF( (x(i,Tdim) - xyzlimits(1,Tdim))/SQRT(4.0_pr*nu*(t_local+t_zero)) ) &
               /DERF( (xyzlimits(2,Tdim) - xyzlimits(1,Tdim))/SQRT(4.0_pr*nu*t_local) )
       END DO
    END IF

  END SUBROUTINE user_initial_conditions

!--********************************  
  ! Arguments
  ! u         - field on adaptive grid
  ! nlocal      - number of active points
  ! ne_local       - number of equations
  ! t         - current time
  !
  ! Global variables used
  ! ifn       - number of equations
  ! x()
  ! xO
  ! y0
  !--******************************** 
  SUBROUTINE user_algebraic_BC (Lu, u, nlocal, ne_local, jlev, meth)
    !--Defines boundary condition type
    IMPLICIT NONE
    INTEGER , INTENT (IN) :: jlev, meth, ne_local, nlocal
    REAL (pr), DIMENSION (nlocal*ne_local), INTENT (INOUT) :: Lu
    REAL (pr), DIMENSION (nlocal*ne_local), INTENT (IN)    :: u

    INTEGER :: i, ie, ii, shift
    !REAL (pr), DIMENSION (ne_nlocal,nlocal,dim) :: du, d2u !uncomment if Neuman BC are used, need to call 
    INTEGER :: face_type, nloc
    INTEGER, DIMENSION(0:dim) :: i_p_face
    INTEGER, DIMENSION(dim) :: face
    INTEGER, DIMENSION(nwlt) :: iloc

    !BC for both 2D and 3D cases

    !if Neuman BC are used, need to call 
    !CALL c_diff_fast (u, du, d2u, jlev, nlocal, meth, 10, ne_local, 1, ne_local)

    DO ie = 1, ne_local
       shift = nlocal*(ie-1)
       !--Go through all Boundary points that are specified
       i_p_face(0) = 1
       DO i=1,dim
          i_p_face(i) = i_p_face(i-1)*3
       END DO
       DO face_type = 0, 3**dim - 1
          face = INT(MOD(face_type,i_p_face(1:dim))/i_p_face(0:dim-1))-1
          IF( ANY( face(1:dim) /= 0) ) THEN ! goes only through boundary points
             CALL get_all_indices_by_face (face_type, jlev, nloc, iloc)
             IF(nloc > 0 ) THEN 
                IF( face(Tdim) == -1  ) THEN                          ! Ymin face (entire face) 
                   Lu(shift+iloc(1:nloc)) = u(shift+iloc(1:nloc))  !Dirichlet conditions
!!$                Lu(shift+iloc(1:nloc)) = du(ie,iloc(1:nloc),2)  !Neuman conditions
                ELSE IF( face(Tdim) == 1  ) THEN                      ! Ymax face (entire face) 
                   Lu(shift+iloc(1:nloc)) = u(shift+iloc(1:nloc))  !Dirichlet conditions
!!$                Lu(iloc(1:nloc),ie) = du(ie,iloc(1:nloc),2)  !Neuman conditions
                END IF
             END IF
          END IF
       END DO
    END DO

  END SUBROUTINE user_algebraic_BC

  SUBROUTINE user_algebraic_BC_diag (Lu_diag, nlocal, ne_local, jlev, meth)
    !--Defines boundary condition type
    IMPLICIT NONE
    INTEGER , INTENT (IN) :: jlev, meth, ne_local, nlocal
    REAL (pr), DIMENSION (nlocal*ne_local), INTENT (INOUT) :: Lu_diag

    INTEGER :: i, ie, ii, shift
    ! REAL (pr), DIMENSION (nlocal,dim) :: du, d2u  !uncomment if Neuman BC are used
    INTEGER :: face_type, nloc
    INTEGER, DIMENSION(0:dim) :: i_p_face
    INTEGER, DIMENSION(dim) :: face
    INTEGER, DIMENSION(nwlt) :: iloc
    
    !BC for both 2D and 3D cases

    !if Neuman BC are used, need to call 
    !CALL c_diff_diag ( du, d2u, jlev, nlocal, meth, meth, 10)

    DO ie = 1, ne_local
       shift = nlocal*(ie-1)
       !--Go through all Boundary points that are specified
       i_p_face(0) = 1
       DO i=1,dim
          i_p_face(i) = i_p_face(i-1)*3
       END DO
       DO face_type = 0, 3**dim - 1
          face = INT(MOD(face_type,i_p_face(1:dim))/i_p_face(0:dim-1))-1
          IF( ANY( face(1:dim) /= 0) ) THEN ! goes only through boundary points
             CALL get_all_indices_by_face (face_type, jlev, nloc, iloc)
             IF(nloc > 0 ) THEN 
                IF( face(Tdim) == -1 ) THEN                        ! Ymin face (entire face) 
                   Lu_diag(shift+iloc(1:nloc)) = 1.0_pr                !Dirichlet conditions
!!$                Lu_diag(shift+iloc(1:nloc)) = du(iloc(1:nloc),2)  !Neuman conditions
                ELSE IF( face(Tdim) == 1  ) THEN                    ! Ymax face (entire face) 
                   Lu_diag(shift+iloc(1:nloc)) = 1.0_pr                !Dirichlet conditions
!!$                Lu_diag(shift+iloc(1:nloc)) = du(iloc(1:nloc),2)  !Neuman conditions
                END IF
             END IF
          END IF
       END DO
    END DO
 
  END SUBROUTINE user_algebraic_BC_diag

  SUBROUTINE user_algebraic_BC_rhs (rhs, ne_local, nlocal, jlev)
    !--Sets rhs for boundary conditions
    IMPLICIT NONE
    INTEGER , INTENT (IN) :: ne_local, nlocal, jlev
    REAL (pr), DIMENSION (nlocal*ne_local), INTENT (INOUT) :: rhs

    INTEGER :: i, ie, ii, shift
    INTEGER :: face_type, nloc
    INTEGER, DIMENSION(0:dim) :: i_p_face
    INTEGER, DIMENSION(dim) :: face
    INTEGER, DIMENSION(nwlt) :: iloc

    !BC for both 2D and 3D cases

    DO ie = 1, ne_local
       shift = nlocal*(ie-1)
       !--Go through all Boundary points that are specified
       i_p_face(0) = 1
       DO i=1,dim
          i_p_face(i) = i_p_face(i-1)*3
       END DO
       DO face_type = 0, 3**dim - 1
          face = INT(MOD(face_type,i_p_face(1:dim))/i_p_face(0:dim-1))-1
          IF( ANY( face(1:dim) /= 0) ) THEN ! goes only through boundary points
             CALL get_all_indices_by_face (face_type, jlev, nloc, iloc)
             IF(nloc > 0 ) THEN 
                IF( face(Tdim) == -1  ) THEN                    ! Ymin face (entire face) 
                   rhs(shift+iloc(1:nloc)) = T_1                 
                ELSE IF( face(Tdim) == 1 ) THEN                ! Ymax face (entire face) 
                   rhs(shift+iloc(1:nloc)) = T_0
                END IF
             END IF
          END IF
       END DO
    END DO

  END SUBROUTINE user_algebraic_BC_rhs

  SUBROUTINE user_project (u, p, nlocal, meth)
    !--Makes u divergence free
    IMPLICIT NONE

    INTEGER, INTENT (IN) :: meth, nlocal
    REAL (pr), DIMENSION (nlocal),     INTENT(INOUT) :: p
    REAL (pr), DIMENSION (nlocal,n_integrated), INTENT(INOUT) :: u

    REAL (pr), DIMENSION (nlocal) :: dp
    REAL (pr), DIMENSION (nlocal) :: f

    INTEGER :: i
	INTEGER, PARAMETER :: ne = 1
    INTEGER, DIMENSION(ne) :: clip 


  END SUBROUTINE user_project

  FUNCTION user_rhs (u_integrated,p)	 
    IMPLICIT NONE
    REAL (pr), DIMENSION (ng,ne), INTENT(IN) :: u_integrated !1D flat version of integrated variables without BC 
    REAL (pr), DIMENSION (ng), INTENT(IN) :: p
    REAL (pr), DIMENSION (n) :: user_rhs


    INTEGER :: ie, shift
    INTEGER, PARAMETER :: meth=1
    REAL (pr), DIMENSION (ne,ng,dim) :: du, d2u, du_dummy
    REAL (pr), DIMENSION (ng,dim)     :: dp
    REAL (pr), DIMENSION (ng,ne) :: f

    IF (IMEXswitch .LE. 0) THEN

       DO ie = n_var_temp, n_var_temp
          shift = ng*(ie-1)
          user_rhs(shift+1:shift+ng) = 0.0_pr
       END DO

    END IF
    IF (IMEXswitch .GE. 0) THEN
       f=u_integrated**2
       CALL c_diff_fast(f, du, du_dummy, j_lev, ng, 1, 10, ne, 1, ne)
       CALL c_diff_fast(u_integrated, du, d2u, j_lev, ng, 1, 01, ne, 1, ne)
       
       !--Form right hand side of heat equations
       DO ie = n_var_temp, n_var_temp
          shift = ng*(ie-1)
          !CALL c_diff_fast(u_integrated(1:ng,ie), du, d2u, j_lev, nwlt, 1, 10, 1, 1, 1)
          user_rhs(shift+1:shift+ng) = nu*SUM(d2u(ie,:,:),2) - 0.5_pr*du(ie,:,1)
       END DO
       !---------------- Algebraic BC are set up in user_algebraic_BC and they autmatically overwrite evolution BC
    END IF

  END FUNCTION user_rhs

  ! find Jacobian of Right Hand Side of the problem
  FUNCTION user_Drhs (u, u_prev_timestep_loc, meth)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: meth
	!u_prev_timestep passed local to cast it into 2dim array to work with vector derivatives
    REAL (pr), DIMENSION (ng,ne) :: u, u_prev_timestep_loc
    REAL (pr), DIMENSION (n) :: user_Drhs
    REAL (pr), DIMENSION (ng,ne) :: f

 
    INTEGER :: ie, shift
    REAL (pr), DIMENSION (2*ne,ng,dim) :: du  ! 1st derivatives for u (in du(1:ne,:,:)) and u_prev_timestep (in du(ne+1:2*ne,:,:))
    REAL (pr), DIMENSION (ne  ,ng,dim) :: d2u ! 2nd derivatives for u 
    REAL (pr), DIMENSION (ne  ,ng,dim) :: du_dummy ! passed when only calculating 1st derivative.
    !Find batter way to do this!! du_dummy with no storage..

    IF (IMEXswitch .LE. 0) THEN

       DO ie = n_var_temp, n_var_temp
          shift = ng*(ie-1)
          user_Drhs(shift+1:shift+ng) = 0.0_pr
       END DO
       
    END IF
    IF (IMEXswitch .GE. 0) THEN
       f(:,1)=u(:,1)*u_prev_timestep_loc(:,1)
       CALL c_diff_fast(u, du(1:ne,:,:), d2u(1:ne,:,:), j_lev, ng, meth, 01, ne, 1, ne)
       CALL c_diff_fast(f, du(1:ne,:,:), du_dummy(1:ne,:,:), j_lev, ng, meth, 10, ne, 1, ne )

       DO ie = n_var_temp, n_var_temp
          shift = ng*(ie-1)
          !CALL c_diff_fast(u_integrated(1:ng,ie), du, d2u, j_lev, nwlt, 1, 11, 1, 1, 1)
          user_Drhs(shift+1:shift+ng) = nu*SUM(d2u(ie,:,:),2) - du(ie,:,1)
       END DO
       
    END IF

  END FUNCTION user_Drhs


  !
  ! Uses u_prev_timestep, which is global flat (1D ) array set in time_step_cn()
  !
  FUNCTION user_Drhs_diag (meth)
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: meth
    REAL (pr), DIMENSION (n) :: user_Drhs_diag

    INTEGER :: ie, shift
    REAL (pr), DIMENSION (ng,dim) :: du, d2u
    REAL (pr), DIMENSION (ne,ng,dim) :: du_prev_timestep
    REAL (pr), DIMENSION (ne  ,ng,dim) :: du_dummy ! passed when only calculating 1st derivative.

    IF (IMEXswitch .LE. 0) THEN
       
       DO ie = n_var_temp, n_var_temp
          shift = ng*(ie-1)
          user_Drhs_diag(shift+1:shift+ng) = 1.0_pr
       END DO

    END IF
    IF (IMEXswitch .GE. 0) THEN

       !
       ! does not rely on u so we can call it once here
       !
       CALL c_diff_diag(du, d2u, j_lev, ng, meth, meth, 11)

       DO ie = n_var_temp, n_var_temp
          shift = ng*(ie-1)
          !CALL c_diff_fast(u_integrated(1:ng,ie), du, d2u, j_lev, nwlt, 1, 11, 1, 1, 1)
          user_Drhs_diag(shift+1:shift+ng) = nu*SUM(d2u,2) - du(:,1)*u_prev_timestep(shift+1:shift+ng)
       END DO

    END IF

  END FUNCTION user_Drhs_diag

  FUNCTION user_chi (nlocal, t_local )
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: nlocal
    REAL (pr), INTENT (IN) :: t_local
    REAL (pr), DIMENSION (nlocal) :: user_chi

    user_chi = 0.0_pr
  END FUNCTION user_chi


  FUNCTION user_mapping ( xlocal, nlocal, t_local )
    USE curvilinear
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: nlocal
    REAL (pr), INTENT (IN) :: t_local
    REAL (pr), DIMENSION (nlocal,dim), INTENT(IN) :: xlocal
    REAL (pr), DIMENSION (nlocal,dim) :: user_mapping

    user_mapping(:,1:dim) = xlocal(:,1:dim)

  END FUNCTION user_mapping
  !
  ! Calculate any statitics
  !
  ! startup_flag - 0 when adapting to IC,then 1 inmain integration loop
  !
  SUBROUTINE user_stats ( u ,j_mn, startup_flag)
    IMPLICIT NONE
    INTEGER , INTENT (IN) :: startup_flag
    INTEGER , INTENT (IN) :: j_mn 
    REAL (pr), DIMENSION (nwlt,1:n_var), INTENT (IN) :: u
  END SUBROUTINE user_stats


  SUBROUTINE user_cal_force (u, n, t_local, force, drag, lift)
    !--Calculates drag and lift on obstacle using penalization formula
    IMPLICIT NONE

    INTEGER, INTENT (IN) :: n
    REAL (pr), INTENT (IN) :: t_local
    REAL (pr), DIMENSION (dim), INTENT (INOUT) :: force
    REAL (pr), INTENT (OUT) :: drag, lift
    REAL (pr), DIMENSION (n,dim) :: u
    drag = 0.0_PR
    lift = 0.0_PR

    !
    ! There is no obstacle in flow
    !

  END SUBROUTINE user_cal_force

  !
  ! Read input from "case_name"_pde.inp file
  ! case_name is in string file_name read from command line
  ! in read_command_line_input()
  ! 
  SUBROUTINE user_read_input()
    IMPLICIT NONE

  call input_real ('nu',nu,'stop', &
       ' nu: heat conduction coefficient')

  call input_real ('T_0',T_0,'stop', &
       ' T_0: initial temperature')

  call input_real ('T_1',T_1,'stop', &
       ' T_1: temperature at Xmin')

  END SUBROUTINE user_read_input



  !
  ! calculate any additional variables
  ! 
  ! arg 
  ! flag  - 0 calledwhile adapting to IC, 1 - called in main time integration loop
  !
  ! These additiona variables are calculated and left in real space.
  !
  SUBROUTINE user_additional_vars( t_local, flag )
    IMPLICIT NONE
    REAL (pr), INTENT (IN) ::  t_local
    INTEGER , INTENT(IN) :: flag ! 0- called during adaption to IC, 1 called during main integration loop
    REAL (pr) :: t_zero
    INTEGER :: i
    DOUBLE PRECISION :: DERF
    EXTERNAL :: DERF

    t_zero = 0.0_pr
    IF(t_local == 0) t_zero = 1.0e-12_pr
    DO i =1, nwlt
       u(i,n_var_pressure) = T_1 - (T_1 - T_0)*DERF( (x(i,Tdim) - xyzlimits(1,Tdim))/SQRT(4.0_pr*nu*(t_local+t_zero)) ) &
            /DERF( (xyzlimits(2,Tdim) - xyzlimits(1,Tdim))/SQRT(4.0_pr*nu*t_local) )
    END DO

  END SUBROUTINE user_additional_vars


  !
  ! calculate any additional scalar variables
  ! 
  ! arg 
  ! flag  - 0 calledwhile adapting to IC, 1 - called in main time integration loop
  !
  SUBROUTINE user_scalar_vars( flag )
    IMPLICIT NONE
    INTEGER , INTENT(IN) :: flag ! 0- called during adaption to IC, 1 called during main integration loop



  END SUBROUTINE user_scalar_vars

  !
  !************ Calculating Scales ***************************
  !
  ! Note the order of the components in the scl array
  ! correspond to u_tn, v_tn, w_tn, u_tn-1, v_tn-1, w_tn-1, u_tn-2, v_tn-2, w_tn-2
  ! 
  SUBROUTINE user_scales(flag, use_default, u, nlocal, ne_local, l_n_var_adapt , l_n_var_adapt_index, &
       scl, scl_fltwt) !add
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: flag ! 0 during initial adaptation, then 1 during time advancement
    LOGICAL , INTENT(INOUT) :: use_default
    INTEGER, INTENT (IN) :: nlocal, ne_local
    REAL (pr), DIMENSION (1:nlocal,1:ne_local), INTENT (IN) :: u
    LOGICAL , INTENT (IN) :: l_n_var_adapt(ne_local)
    INTEGER , INTENT (IN) :: l_n_var_adapt_index(1:ne_local)
    REAL (pr), DIMENSION (1:ne_local), INTENT (INOUT) :: scl
    REAL (pr) ,  INTENT(IN) ::scl_fltwt !weight for temporal filter on scl

    !
    ! Ignore the output of this routine and use default scales routine
    !
    use_default = .TRUE. 

    !
    ! NOTE: For a parallel run, synchronize scl(:) across the processors in user_scales.
    !       Use the subroutine scales of default_util.f90 as an example,
    !       subroutine parallel_global_sum of parallel.f90 may be helpful.
    !       Already synchronized: sumdA_global
  END SUBROUTINE user_scales

SUBROUTINE user_cal_cfl (use_default, u, cfl_out)
  USE precision
  USE sizes
  USE pde
  IMPLICIT NONE
  LOGICAL , INTENT(INOUT) :: use_default
  REAL (pr),                                INTENT (INOUT) :: cfl_out
  REAL (pr), DIMENSION (nwlt,n_integrated), INTENT (IN)    :: u

  INTEGER                    :: i
  REAL (pr)                  :: floor
  REAL (pr), DIMENSION (dim) :: cfl
  REAL (pr), DIMENSION(dim,nwlt) :: h_arr

  use_default = .FALSE.

  floor = 1e-12_pr
  cfl_out = floor
  
  CALL get_all_local_h (h_arr)
  
!!$  DO i = 1, nwlt
!!$     cfl(1:dim) = ABS (u(i,1:dim)+Umn(1:dim)) * dt/h_arr(1:dim,i)
!!$     cfl_out = MAX (cfl_out, MAXVAL(cfl))
!!$  END DO

  cfl_out = 1.0_pr ! no CFL condition

END SUBROUTINE user_cal_cfl

!******************************************************************************************
!************************************* SGS MODEL ROUTINES *********************************
!******************************************************************************************

!
! Intialize sgs model
! This routine is called once in the first
! iteration of the main time integration loop.
! weights and model filters have been setup for first loop when this routine is called.
!
SUBROUTINE user_init_sgs_model( )
  IMPLICIT NONE


! LDM: Giuliano

! THE sgs forcign is stored for the time step in sgs_mdl_force(1:nlocal,1:n_integrated)
! where nlocal should be nwlt.


!          print *,'initializing LDM ...'       
!          CALL sgs_mdl_force( u(:,1:n_integrated), nwlt, j_lev, .TRUE.)


END SUBROUTINE user_init_sgs_model

!
! calculate sgs model forcing term
! user_sgs_force is called int he beginning of each times step in time_adv_cn().
! THE sgs forcign is stored for the time step in sgs_mdl_force(1:nlocal,1:n_integrated)
! where nlocal should be nwlt.
! 
! Accesses u from field module, 
!          j_lev from wlt_vars module,
!
SUBROUTINE  user_sgs_force (u_loc, nlocal)
  IMPLICIT NONE

  INTEGER,                         INTENT (IN) :: nlocal
  REAL (pr), DIMENSION (nlocal,n_integrated), INTENT (INOUT) :: u_loc

END SUBROUTINE  user_sgs_force

  FUNCTION user_sound_speed (u, neq, nwlt)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: nwlt, neq
    REAL (pr), DIMENSION (nwlt,neq), INTENT(IN) :: u
    REAL (pr), DIMENSION (nwlt) :: user_sound_speed

    user_sound_speed(:) = 0.0_pr 
    !user_sound_speed(:) = SQRT( gamma*(gamma-1.0_pr)* &
    !       ( u(:,4)-0.5_pr*SUM(u(:,n_var_mom(1:dim))**2,DIM=2)/u(:,1)) ) ! pressure

  END FUNCTION user_sound_speed

  SUBROUTINE  user_pre_process
    IMPLICIT NONE
     
  END SUBROUTINE user_pre_process

  SUBROUTINE  user_post_process
    IMPLICIT NONE
    
  END SUBROUTINE user_post_process

END MODULE user_case
