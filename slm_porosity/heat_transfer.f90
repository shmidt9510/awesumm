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
  INTEGER n_var_pressure
  INTEGER n_var_temp
  INTEGER n_var_enthalpy
  INTEGER n_var_lfrac
  INTEGER n_var_porosity

  REAL (pr) :: fusion_delta         ! = liquidus - solidus
  REAL (pr) :: fusion_heat          ! the latent heat of fusion
  REAL (pr) :: convective_transfer  ! convective heat transfer coefficient
  REAL (pr) :: enthalpy_L           ! enthalpy at the liquidus temperature
  REAL (pr) :: enthalpy_S           ! enthalpy at the solidus temperature
  REAL (pr) :: power                ! laser power
  REAL (pr) :: absorb               ! absorptivity
  REAL (pr) :: scanning_speed
  REAL (pr) :: initial_porosity
  REAL (pr) :: initial_enthalpy
  REAL (pr) :: fusion_smoother                ! half-wide of step
  REAL (pr) :: conductivity_der     ! first derivative of conductivity on temperature
  REAL (pr) :: capacity_der         ! first derivative of capacity on temperature
  REAL (pr), DIMENSION(3) :: x0     ! Initial coordinates of the center of the laser beam
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
    LOGICAL, OPTIONAL :: VERB         ! debug info
    INTEGER :: i

    !register_var( var_name, integrated, adapt, interpolate, exact, saved, req_restart, FRONT_LOAD )
    CALL register_var( 'enthalpy ', &
        integrated=.TRUE., &
        adapt=(/.TRUE.,.TRUE./), &
        interpolate=(/.TRUE.,.TRUE./), &
        exact=(/.FALSE.,.FALSE./), &
        saved=.TRUE., &
        req_restart=.TRUE.)
    CALL register_var( 'temperature', &
        integrated=.FALSE., &
        adapt=(/.FALSE.,.FALSE./), &
        interpolate=(/.FALSE.,.FALSE./), &
        exact=(/.FALSE.,.FALSE./), &
        saved=.TRUE., &
        req_restart=.FALSE.)
    CALL register_var( 'liquid_fraction', &
        integrated=.FALSE., &
        adapt=(/.FALSE.,.FALSE./), &
        interpolate=(/.FALSE.,.FALSE./), &
        exact=(/.FALSE.,.FALSE./), &
        saved=.TRUE., &
        req_restart=.FALSE.)
    CALL register_var( 'porosity', &
        integrated=.FALSE., &
        adapt=(/.TRUE.,.TRUE./), &
        interpolate=(/.FALSE.,.FALSE./), &
        exact=(/.FALSE.,.FALSE./), &
        saved=.TRUE., &
        req_restart=.FALSE.)
    CALL register_var( 'pressure ', &
        integrated=.FALSE., &
        adapt=(/.FALSE.,.FALSE./), &
        interpolate=(/.FALSE.,.FALSE./), &
        exact=(/.FALSE.,.FALSE./), &
        saved=.TRUE., &
        req_restart=.FALSE.)
    CALL setup_mapping()
    CALL print_variable_registery( FULL=.TRUE.)

    n_var_pressure  = get_index('pressure ')  
    n_var_temp      = get_index('temperature')  
    n_var_enthalpy  = get_index('enthalpy')  
    n_var_lfrac     = get_index('liquid_fraction')  
    n_var_porosity  = get_index('porosity')  

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
    u = 0
  END SUBROUTINE  user_exact_soln

  SUBROUTINE user_initial_conditions (u, nlocal, ne_local, t_local, scl, scl_fltwt, iter)
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: nlocal, ne_local
    INTEGER, INTENT (INOUT) :: iter ! iteration of call while adapting initial grid
    REAL (pr), DIMENSION (nlocal,ne_local), INTENT (INOUT) :: u
    REAL (pr)  :: scl(1:n_var),scl_fltwt
    REAL (pr), INTENT (IN) :: t_local
    REAL (pr) :: t_zero
    INTEGER :: i

    IF (dim.EQ.2) x0(dim) = xyzlimits(1,dim) 
    IF ( IC_restart_mode .EQ. 0 ) THEN
       DO i = 1, nlocal
          u(i,n_var_enthalpy) = initial_enthalpy*EXP(-SUM((x(i,:)-x0)**2))*EXP(-(x(i,dim)-x0(dim))**2*power*absorb)
       END DO
       !WHERE (x(:,dim).NE.x0(dim))
       !   u(:,n_var_enthalpy) = 0
       !END WHERE
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
    IMPLICIT NONE
    INTEGER , INTENT (IN) :: jlev, meth, ne_local, nlocal
    REAL (pr), DIMENSION (nlocal*ne_local), INTENT (INOUT) :: Lu
    REAL (pr), DIMENSION (nlocal*ne_local), INTENT (IN)    :: u

    INTEGER :: i, ie, shift, face_type, nloc
    REAL (pr), DIMENSION (ne_local,nlocal,dim) :: du, d2u
    INTEGER, DIMENSION(0:dim) :: i_p_face
    INTEGER, DIMENSION(dim) :: face
    INTEGER, DIMENSION(nwlt) :: iloc

    CALL c_diff_fast (linearized_temperature(u, u_prev_timestep), du, d2u, jlev, nlocal, meth, 10, ne_local, 1, ne_local)

    DO ie = 1, ne_local
       shift = nlocal*(ie-1)
       i_p_face(0) = 1
       DO i=1,dim
          i_p_face(i) = i_p_face(i-1)*3
       END DO
       DO face_type = 0, 3**dim - 1
          face = INT(MOD(face_type,i_p_face(1:dim))/i_p_face(0:dim-1))-1
          IF( ANY( face(1:dim) /= 0) .AND. ie == n_var_enthalpy ) THEN
             CALL get_all_indices_by_face (face_type, jlev, nloc, iloc)
             iloc(1:nloc) = shift + iloc(1:nloc)
             IF( nloc > 0 ) THEN
                IF( face(dim) < 0 ) THEN    ! z=0 face
                   ! dependence on temperature should be linear; therefore, we use only u_prev_timestep here
                   CALL my_bc(Lu, linearized_temperature(u, u_prev_timestep), du(ie,:,:), temperature(u_prev_timestep), nloc, iloc)
                ELSE                        ! other faces
                   Lu(iloc(1:nloc)) = u(iloc(1:nloc))
                END IF
             END IF
          END IF
       END DO
    END DO
    !PRINT *, 'user_algebraic_BC', Lu(shift+iloc(1:nloc))

  END SUBROUTINE user_algebraic_BC

  SUBROUTINE user_algebraic_BC_diag (Lu_diag, nlocal, ne_local, jlev, meth)
    IMPLICIT NONE
    INTEGER , INTENT (IN) :: jlev, meth, ne_local, nlocal
    REAL (pr), DIMENSION (nlocal*ne_local), INTENT (INOUT) :: Lu_diag

    INTEGER :: i, ie, shift, face_type, nloc
    REAL (pr), DIMENSION (nlocal,dim) :: du, d2u
    REAL (pr), DIMENSION (nlocal) :: ones
    INTEGER, DIMENSION(0:dim) :: i_p_face
    INTEGER, DIMENSION(dim) :: face
    INTEGER, DIMENSION(nwlt) :: iloc
    
    CALL c_diff_diag (du, d2u, jlev, nlocal, meth, meth, 10)

    DO ie = 1, ne_local
       shift = nlocal*(ie-1)
       i_p_face(0) = 1
       DO i=1,dim
          i_p_face(i) = i_p_face(i-1)*3
       END DO
       DO face_type = 0, 3**dim - 1
          face = INT(MOD(face_type,i_p_face(1:dim))/i_p_face(0:dim-1))-1
          IF( ANY( face(1:dim) /= 0) .AND. ie == n_var_enthalpy ) THEN
             CALL get_all_indices_by_face (face_type, jlev, nloc, iloc)
             iloc(1:nloc) = shift + iloc(1:nloc)
             IF( nloc > 0 ) THEN 
                IF( face(dim) < 0 ) THEN    ! z=0 face
                   ones = 1.0_pr
                   CALL my_bc(Lu_diag, ones, du, temperature(u_prev_timestep), nloc, iloc)
                   Lu_diag(iloc(1:nloc)) = Lu_diag(iloc(1:nloc)) * &
                      linearized_temperature_diag(u_prev_timestep(iloc(1:nloc)))
                ELSE                        ! other faces
                   Lu_diag(iloc(1:nloc)) = 1.0_pr
                END IF
             END IF
          END IF
       END DO
    END DO
    !PRINT *, 'user_algebraic_BC_diag', Lu_diag(shift+iloc(1:nloc))
 
  END SUBROUTINE user_algebraic_BC_diag

  SUBROUTINE user_algebraic_BC_rhs (rhs, ne_local, nlocal, jlev)
    IMPLICIT NONE
    INTEGER , INTENT (IN) :: ne_local, nlocal, jlev
    REAL (pr), DIMENSION (nlocal*ne_local), INTENT (INOUT) :: rhs
    REAL (pr), DIMENSION (nlocal*ne_local) :: linearized_rhs

    INTEGER :: i, ie, shift, face_type, nloc, meth=1
    REAL (pr), DIMENSION (ne_local,nlocal,dim) :: du, d2u
    INTEGER, DIMENSION(0:dim) :: i_p_face
    INTEGER, DIMENSION(dim) :: face
    INTEGER, DIMENSION(nwlt) :: iloc

    CALL c_diff_fast (linearized_temperature_rhs(u_prev_timestep), du, d2u, jlev, nlocal, meth, 10, ne_local, 1, ne_local)

    DO ie = 1, ne_local
       shift = nlocal*(ie-1)
       i_p_face(0) = 1
       DO i=1,dim
          i_p_face(i) = i_p_face(i-1)*3
       END DO
       DO face_type = 0, 3**dim - 1
          face = INT(MOD(face_type,i_p_face(1:dim))/i_p_face(0:dim-1))-1
          IF( ANY( face(1:dim) /= 0) .AND. ie == n_var_enthalpy ) THEN
             CALL get_all_indices_by_face (face_type, jlev, nloc, iloc)
             iloc(1:nloc) = shift + iloc(1:nloc)
             IF( nloc > 0 ) THEN 
                IF( face(dim) < 0 ) THEN    ! z=0 face
                   IF( dim == 3 ) THEN
                      rhs(iloc(1:nloc)) = -exp(-(x(iloc(1:nloc), 1) - scanning_speed*t - x0(1))**2 - (x(iloc(1:nloc), 2) - x0(2))**2 ) / pi
                   ELSE
                      rhs(iloc(1:nloc)) = -exp(-(x(iloc(1:nloc), 1) - scanning_speed*t - x0(1))**2 ) / pi**.5
                   END IF
                   CALL my_bc(linearized_rhs, linearized_temperature_rhs(u_prev_timestep), du(ie,:,:), temperature(u_prev_timestep), nloc, iloc)
                   rhs(iloc(1:nloc)) = rhs(iloc(1:nloc)) * absorb * power + linearized_rhs(iloc(1:nloc))
                ELSE                        ! other faces
                   rhs(iloc(1:nloc)) = 0
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
  END SUBROUTINE user_project

  FUNCTION user_rhs (u_integrated, p)
    IMPLICIT NONE
    REAL (pr), DIMENSION (ng,ne), INTENT(IN) :: u_integrated
    REAL (pr), DIMENSION (ng), INTENT(IN) :: p
    REAL (pr), DIMENSION (n) :: user_rhs
    REAL (pr), DIMENSION (ng) :: T
    INTEGER :: ie, shift, i
    INTEGER, PARAMETER :: meth = 1
    REAL (pr), DIMENSION (ne,ng,dim) :: du, du_dummy
    REAL (pr), DIMENSION (dim,ng,dim):: d2u, d2u_dummy

    ie = n_var_enthalpy
    shift = ng*(ie-1)
    IF (IMEXswitch .LE. 0) THEN
       user_rhs(shift+1:shift+ng) = 0.0_pr
    END IF
    IF (IMEXswitch .GE. 0) THEN
       CALL c_diff_fast(enthalpy_fusion(u_integrated(:,ie)), du, du_dummy, j_lev, ng, meth, 10, ne, 1, ne)
       T = temperature(u_integrated(:,ie))
       du(ie,:,:) = du(ie,:,:) * SPREAD(conductivity(T) / capacity(T) * porosity_term(), 2, dim)
       CALL c_diff_fast(du(ie,:,:), d2u, d2u_dummy, j_lev, ng, meth, 10, dim, 1, dim)
       DO i = 1, dim
          user_rhs(shift+1:shift+ng) = user_rhs(shift+1:shift+ng) + d2u(i,:,i)
       END DO
    END IF
    !PRINT *, 'user_rhs', user_rhs

  END FUNCTION user_rhs

  FUNCTION user_Drhs (pert_u, u_prev, meth)
    IMPLICIT NONE
    INTEGER, INTENT(IN) :: meth
    REAL (pr), DIMENSION (ng,ne) :: pert_u, u_prev
    REAL (pr), DIMENSION (n) :: user_Drhs
    REAL (pr), DIMENSION (ng) :: T
    INTEGER :: ie, shift, i
    REAL (pr) :: first_der
    REAL (pr), DIMENSION (ng,2*ne) :: for_du
    REAL (pr), DIMENSION (2*ne,ng,dim) :: du, du_dummy ! der(pert_u) in du(1:ne,:,:) and der(u) in du(ne+1:2*ne,:,:)
    REAL (pr), DIMENSION (dim,ng,dim) :: d2u, d2u_dummy
    REAL (pr), DIMENSION (ng,dim) :: part1, part2

    ie = n_var_enthalpy
    shift = ng*(ie-1)
    first_der = conductivity_der - 2*capacity_der
    IF (IMEXswitch .LE. 0) THEN
       user_Drhs(shift+1:shift+ng) = 0.0_pr
    END IF
    IF (IMEXswitch .GE. 0) THEN
       for_du = RESHAPE((/ Denthalpy_fusion(pert_u(:,ie), u_prev(:,ie)), enthalpy_fusion(u_prev(:,ie)) /), SHAPE(for_du))
       CALL c_diff_fast(for_du, du, du_dummy, j_lev, ng, meth, 10, 2*ne, 1, 2*ne)
       T = temperature(u_prev(:,ie))
       part1 = du(2*ie,:,:) * SPREAD(first_der / capacity(T)**3 * Denthalpy_fusion(pert_u(:,ie), u_prev(:,ie)), 2, dim)
       part2 = du(ie,:,:) * SPREAD(conductivity(T) / capacity(T), 2, dim)
       CALL c_diff_fast((part1 + part2) * SPREAD(porosity_term(), 2, dim), d2u, d2u_dummy, j_lev, ng, meth, 10, dim, 1, dim)
       DO i = 1, dim
          user_Drhs(shift+1:shift+ng) = user_Drhs(shift+1:shift+ng) + d2u(i,:,i)
       END DO
    END IF
    !print *, 'u_prev', u_prev
    !print *, 'for_du', for_du
    !print *, 'du', du
    !print *, 'd2u', d2u
    !print *, 'porosity_term()', porosity_term()
    !print *, 'temperature', T
    !print *, 'capacity', capacity(T)**3
    !print *, 'part1', part1
    !print *, 'part2', part2
    !print *, 'pert_u(:,ie)', pert_u(1:2,ie)
    !print *, 'Denthalpy_fusion(pert_u(:,ie), u_prev(:,ie))', Denthalpy_fusion(pert_u(:,ie), u_prev(:,ie))
    !print *, 'enthalpy_fusion(u_prev(:,ie))', enthalpy_fusion(u_prev(:,ie))
    !PRINT *, 'user_Drhs', user_Drhs
    IF (user_Drhs(1).NE.user_Drhs(1)) THEN
       PRINT *, '--- NaN in user_Drhs ---'
       CALL ABORT
    END IF

  END FUNCTION user_Drhs

  FUNCTION user_Drhs_diag (meth)
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: meth
    REAL (pr), DIMENSION (n) :: user_Drhs_diag
    REAL (pr), DIMENSION (ng) :: T
    INTEGER :: ie, shift, i
    REAL (pr) :: first_der
    REAL (pr), DIMENSION (ng,dim) :: du, d2u, part1, part2
    REAL (pr), DIMENSION (ne,ng,dim) :: du_prev, du_dummy

    ie = n_var_enthalpy
    shift = ng*(ie-1)
    first_der = conductivity_der - 2*capacity_der
    IF (IMEXswitch .LE. 0) THEN
       user_Drhs_diag(shift+1:shift+ng) = 1.0_pr
    END IF
    IF (IMEXswitch .GE. 0) THEN
       CALL c_diff_fast(enthalpy_fusion(u_prev_timestep(shift+1:shift+ng)), du_prev, du_dummy, j_lev, ng, meth, 10, ne, 1, ne)
       CALL c_diff_diag(du, d2u, j_lev, ng, meth, meth, -11)
       T = temperature(u_prev_timestep(shift+1:shift+ng))
       part1 = du * du_prev(ie,:,:) * SPREAD(first_der / capacity(T)**3, 2, dim)
       part2 = d2u * SPREAD(conductivity(T) / capacity(T), 2, dim)
       user_Drhs_diag(shift+1:shift+ng) = SUM(part1 + part2, 2) &
           * Denthalpy_fusion_diag(u_prev_timestep(shift+1:shift+ng)) * porosity_term()
    END IF
    !PRINT *, 'user_Drhs_diag', user_Drhs_diag

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

  call input_real ('fusion_delta', fusion_delta, 'stop')
  call input_real ('fusion_heat', fusion_heat, 'stop')
  enthalpy_S = 1.0_pr - fusion_delta/2
  enthalpy_L = 1.0_pr + fusion_delta/2 + fusion_heat
  call input_real ('convective_transfer', convective_transfer, 'stop')
  call input_real ('power', power, 'stop')
  call input_real ('absorb', absorb, 'stop')
  call input_real ('scanning_speed', scanning_speed, 'stop')
  call input_real ('initial_porosity', initial_porosity, 'stop')
  initial_enthalpy = 2*enthalpy_S
  call input_real ('initial_enthalpy', initial_enthalpy, 'default')
  call input_real ('conductivity_der', conductivity_der, 'stop')
  call input_real ('capacity_der', capacity_der, 'stop')
  call input_real_vector ('x0', x0, 3, 'stop')
  call input_real ('fusion_smoother', fusion_smoother, 'stop')
  
  END SUBROUTINE user_read_input



  !
  ! calculate any additional variables
  ! 
  ! arg 
  ! flag  - 0 calledwhile adapting to IC, 1 - called in main time integration loop
  !
  ! These additional variables are calculated and left in real space.
  !
  SUBROUTINE user_additional_vars( t_local, flag )
    IMPLICIT NONE
    REAL (pr), INTENT (IN) ::  t_local
    INTEGER , INTENT(IN) :: flag ! 0- called during adaption to IC, 1 called during main integration loop
    IF (.NOT.flag) THEN
       u(:,n_var_porosity) = initial_porosity
    END IF
    u(:,n_var_temp) = temperature(u(:,n_var_enthalpy))
    u(:,n_var_lfrac) = liquid_fraction(u(:,n_var_enthalpy))
    u(:,n_var_porosity) = porosity(u(:,n_var_enthalpy))
    u(:,n_var_pressure) = 0.0_pr
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
    REAL (pr), DIMENSION (dim) :: cfl
    REAL (pr), DIMENSION(dim,nwlt) :: h_arr
  
    use_default = .FALSE.
    CALL get_all_local_h (h_arr)
    cfl_out = MAXVAL(dt/h_arr(1,:)*scanning_speed)
    PRINT *, "cfl_convective", cfl_out
    DO i = 1, nwlt
       cfl(1:dim) = ABS (u(i,1:dim)+Umn(1:dim)) * dt/h_arr(1:dim,i)
       cfl_out = MAX (cfl_out, MAXVAL(cfl))
    END DO
  
  !  cfl_out = 1.0_pr ! no CFL condition
    IF (u(1,1).NE.u(1,1)) THEN
      PRINT *, '--- INFINITE ---'
      CALL ABORT
    END IF
  
  END SUBROUTINE user_cal_cfl
  
  SUBROUTINE user_init_sgs_model( )
    IMPLICIT NONE
  END SUBROUTINE user_init_sgs_model
  
  SUBROUTINE  user_sgs_force (u_loc, nlocal)
    IMPLICIT NONE
    INTEGER, INTENT (IN) :: nlocal
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
  
  FUNCTION liquid_fraction (enthalpy)
    IMPLICIT NONE
    REAL (pr), INTENT(IN) :: enthalpy(:)
    REAL (pr) :: liquid_fraction(SIZE(enthalpy))
    ! liquid_fraction = 0.5+1/pi * ATAN(fusion_smoother*(enthalpy-1.2_pr))

    liquid_fraction = 1.0_pr/(1.0_pr+EXP(-2.0_pr*fusion_smoother*(enthalpy-1.2_pr)))

    ! liquid_fraction = (enthalpy - enthalpy_S) / (enthalpy_L - enthalpy_S)
    ! liquid_fraction = MAX(0.0_pr, MIN(1.0_pr, liquid_fraction))
  END FUNCTION liquid_fraction

  FUNCTION Dliquid_fraction (Denthalpy, enthalpy)
    IMPLICIT NONE
    REAL (pr), DIMENSION(:), INTENT(IN) :: enthalpy, Denthalpy
    REAL (pr) :: Dliquid_fraction(SIZE(enthalpy))
    ! Dliquid_fraction = Denthalpy *fusion_smoother/(pi*(fusion_smoother**2*(enthalpy-1.2_pr)**2+1))
    Dliquid_fraction =  Denthalpy*fusion_smoother*2.0_pr*EXP(-2.0_pr*fusion_smoother*(enthalpy-1.2_pr))&
      / ((1.0_pr+EXP(-2.0_pr*fusion_smoother*(enthalpy-1.2_pr)))**2)
    ! Dliquid_fraction = 0.0_pr
    ! WHERE (enthalpy_S < enthalpy .AND. enthalpy < enthalpy_L)
    ! Dliquid_fraction = Denthalpy / (enthalpy_L - enthalpy_S)
    ! END WHERE
  END FUNCTION Dliquid_fraction

  FUNCTION Dliquid_fraction_diag (enthalpy)
    IMPLICIT NONE
    REAL (pr), INTENT(IN) :: enthalpy(:)
    REAL (pr) :: Dliquid_fraction_diag(SIZE(enthalpy))
    ! Dliquid_fraction_diag = fusion_smoother/(pi*(fusion_smoother**2*(enthalpy-1.2_pr)**2+1))
    Dliquid_fraction_diag = fusion_smoother*2.0_pr*EXP(-2.0_pr*fusion_smoother*(enthalpy-1.2_pr))&
      / ((1.0_pr+EXP(-2.0_pr*fusion_smoother*(enthalpy-1.2_pr)))**2)
    ! Dliquid_fraction_diag = 0.0_pr
    ! WHERE (enthalpy_S < enthalpy .AND. enthalpy < enthalpy_L)
    ! Dliquid_fraction_diag = 1.0_pr / (enthalpy_L - enthalpy_S)
    ! END WHERE
  END FUNCTION Dliquid_fraction_diag
  
  FUNCTION porosity (enthalpy)
    IMPLICIT NONE
    REAL (pr), INTENT(IN) :: enthalpy(:) 
    REAL (pr) :: porosity(SIZE(enthalpy))
    porosity = MAX(0.0_pr, MIN(u(:,n_var_porosity), initial_porosity*(1.0_pr - liquid_fraction(enthalpy))))
  END FUNCTION porosity
  
  FUNCTION porosity_term ()
    IMPLICIT NONE
    REAL (pr) :: porosity_term(nwlt)
    porosity_term = (1.0_pr - u(:,n_var_porosity)) / (1.0_pr - initial_porosity)
  END FUNCTION porosity_term
  
  PURE FUNCTION conductivity (temperature)
    IMPLICIT NONE
    REAL (pr), INTENT(IN) :: temperature(:) 
    REAL (pr) :: conductivity(SIZE(temperature))
    conductivity = 1.0_pr + conductivity_der*temperature
  END FUNCTION conductivity
  
  PURE FUNCTION capacity (temperature)
    IMPLICIT NONE
    REAL (pr), INTENT(IN) :: temperature(:) 
    REAL (pr) :: capacity(SIZE(temperature))
    capacity = 1.0_pr + 2*capacity_der*temperature
  END FUNCTION capacity
  
  FUNCTION temperature (enthalpy)
    IMPLICIT NONE
    REAL (pr), INTENT(IN) :: enthalpy(:)
    REAL (pr) :: temperature(SIZE(enthalpy))
    IF (capacity_der .LE. 1e-6) THEN
        temperature = enthalpy_fusion(enthalpy)
    ELSE
        !temperature = enthalpy_fusion(enthalpy)
        temperature = (SQRT(1.0_pr + 4*capacity_der*enthalpy_fusion(enthalpy)) - 1.0_pr)/2/capacity_der
    END IF
  END FUNCTION temperature
  
  FUNCTION linearized_temperature (enthalpy, enthalpy_prev)
    IMPLICIT NONE
    REAL (pr), INTENT(IN) :: enthalpy(:), enthalpy_prev(:)
    REAL (pr) :: linearized_temperature(SIZE(enthalpy))
    linearized_temperature = enthalpy_fusion(enthalpy) / &
        (1.0_pr + 2*capacity_der*temperature(enthalpy_prev))
  END FUNCTION linearized_temperature
  
  FUNCTION linearized_temperature_diag (enthalpy_prev)
    IMPLICIT NONE
    REAL (pr), INTENT(IN) :: enthalpy_prev(:)
    REAL (pr) :: linearized_temperature_diag(SIZE(enthalpy_prev))
    linearized_temperature_diag = Denthalpy_fusion_diag(enthalpy_prev) / &
        (1.0_pr + 2*capacity_der*temperature(enthalpy_prev))
  END FUNCTION linearized_temperature_diag
  
  FUNCTION linearized_temperature_rhs (enthalpy_prev)
    IMPLICIT NONE
    REAL (pr), INTENT(IN) :: enthalpy_prev(:)
    REAL (pr) :: linearized_temperature_rhs(SIZE(enthalpy_prev))
    linearized_temperature_rhs = - temperature(enthalpy_prev) + enthalpy_fusion(enthalpy_prev) / &
        (1.0_pr + 2*capacity_der*temperature(enthalpy_prev))
  END FUNCTION linearized_temperature_rhs
  
  FUNCTION enthalpy_fusion (enthalpy)
    IMPLICIT NONE
    REAL (pr), INTENT(IN) :: enthalpy(:)
    REAL (pr) :: enthalpy_fusion(SIZE(enthalpy))
    enthalpy_fusion = enthalpy - fusion_heat*liquid_fraction(enthalpy)
  END FUNCTION enthalpy_fusion
  
  FUNCTION Denthalpy_fusion (Denthalpy, enthalpy)
    IMPLICIT NONE
    REAL (pr), DIMENSION (:), INTENT(IN) :: Denthalpy, enthalpy
    REAL (pr) :: Denthalpy_fusion(SIZE(enthalpy))
    Denthalpy_fusion = Denthalpy - fusion_heat*Dliquid_fraction(Denthalpy, enthalpy)
  END FUNCTION Denthalpy_fusion
  
  FUNCTION Denthalpy_fusion_diag (enthalpy)
    IMPLICIT NONE
    REAL (pr), INTENT(IN) :: enthalpy(:)
    REAL (pr) :: Denthalpy_fusion_diag(SIZE(enthalpy))
    Denthalpy_fusion_diag = 1.0_pr - fusion_heat*Dliquid_fraction_diag(enthalpy)
  END FUNCTION Denthalpy_fusion_diag

  SUBROUTINE my_bc (Lu, u, du, T, nloc, iloc)
    IMPLICIT NONE
    REAL (pr), INTENT(IN) :: u(:), du(:,:), T(:)
    INTEGER, INTENT(IN) :: nloc, iloc(:)
    REAL (pr), INTENT(OUT) :: Lu(:)
    REAL (pr), DIMENSION(nwlt) :: porosity_term_
    porosity_term_ = porosity_term()
    Lu(iloc(1:nloc)) = convective_transfer * u(iloc(1:nloc)) + &
        porosity_term_(iloc(1:nloc)) * conductivity(T(iloc(1:nloc))) * du(iloc(1:nloc), dim)
  END SUBROUTINE my_bc

END MODULE user_case
