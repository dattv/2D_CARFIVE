!=================================================================================================
!> CARTESIAN ADAPTIVE FIVE EQUATION MODEL
!> AUTHOR: VAN-DAT THANG
!> E-MAIL: datthangva@gmail.com
!> E-MAIL: vandatthang@gmail.com
!> SOURCE CODE LINK: https://github.com/dattv/2D_CARFIVE
!=================================================================================================  
MODULE  MODULE_INITIALCONDITION
    
    use MODULE_PRECISION
    use MODULE_NODECOORD
    use MODULE_CONSTANTS
    use MODULE_QUADTREE
    use MODULE_GENERICMETHOD
    use MODULE_CFDMAINDATA
    
    contains
!==================================================================================================
    subroutine initial_2D_simle_check_grid(nelm, tree)
    implicit none
    integer(ip),                    intent(in)              :: nelm
    type(quadtree), dimension(:),   intent(inout), target   :: tree
    
    integer(ip) :: i
    
    !call loop_on_quadtree_array(1, nelm, tree, initial_2d_dambreak_single_level)
    call loop_on_quadtree_array(1, nelm, tree, initial_2d_01_simple)
    
    return
    end subroutine  initial_2D_simle_check_grid
!================================================================================================== 
    subroutine initial_condition(nelm, tree)
    implicit none
    integer(ip),                    intent(in)              :: nelm
    type(quadtree), dimension(:),   intent(inout), target   :: tree
    
    call loop_on_quadtree_array(1, nelm, tree, initial_2d_dambreak_single_level)
    
    return
    end subroutine initial_condition
!================================================================================================== 
    subroutine initial_2d_dambreak_single_level(tree)
    implicit none
    type(quadtree), pointer, intent(inout) :: tree
    
    real(rp)    :: rho1_a1, rho2_a2, rho_u, rho_v, p, a1, a2, rho_e
    real(rp)    :: gamma_mix, pi_mix, u_vel, v_vel, rho_mix, rho
    
    ! body 
    rho = 1.e3_rp
    if (tree%pts(5)%coord(1) >= half) then
        a1 = one - tolerance; a2 = tolerance
        rho1_a1 = rho*a1; rho2_a2 = rho*a2; rho_u = 1500._rp; rho_v = zero; p = 1.e9_rp
    else
        a2 = one - tolerance; a2 = tolerance
        rho1_a1 = rho*a1; rho2_a2 = rho*a2; rho_u = 1500._rp; rho_v = zero; p = 1.e9_rp
    endif

    ! ===> NONE DIMENSIONAL <======================================================================
    rho1_a1 = rho1_a1/Rho_inf
    rho2_a2 = rho2_a2/Rho_inf
      rho_u =   rho_u/ u_inf
      rho_v =   rho_v/ u_inf
          p =       p/p_inf
          
              
    return
    end subroutine  initial_2d_dambreak_single_level
!================================================================================================== 
    subroutine initial_2d_01_simple(tree)
    implicit none
    type(quadtree), pointer, intent(inout) :: tree
    
    if (tree%pts(5)%coord(1) >= half ) then 
        tree%data%u(:) = one
        tree%data%w(:) = one
    else
        tree%data%u(:) = zero
        tree%data%w(:) = zero
    endif

    return
    end subroutine initial_2d_01_simple
!================================================================================================== 
    subroutine initial_2d_02_simple(tree)
    implicit none
    type(quadtree), pointer, intent(inout) :: tree
    
    if (tree%pts(5)%coord(1) >= 0.7_rp) then 
        tree%data%u(:) = one
        tree%data%w(:) = one
    else
        tree%data%u(:) = zero
        tree%data%w(:) = zero
    endif

    return
    end subroutine initial_2d_02_simple    
!==================================================================================================    
!==================================================================================================    
    
END MODULE  MODULE_INITIALCONDITION 