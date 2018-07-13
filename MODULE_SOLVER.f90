!=================================================================================================
!> CARTESIAN ADAPTIVE FIVE EQUATION MODEL
!> AUTHOR: VAN-DAT THANG
!> E-MAIL: datthangva@gmail.com
!> E-MAIL: vandatthang@gmail.com
!> SOURCE CODE LINK: https://github.com/dattv/2D_CARFIVE
!================================================================================================= 
MODULE MODULE_SOLVER
    
    use MODULE_PRECISION
    use MODULE_QUADTREE
    use MODULE_EXSOLVER
    
!================== DEGINS PATTERN ==================    
    
    type, abstract  :: abstract_solver
    contains
    procedure(generic_solver), deferred, pass(this) :: p_solver
    end type abstract_solver
    
    abstract interface
        subroutine generic_solver(this, first, last, tree)
        use MODULE_PRECISION
        use MODULE_QUADTREE
        import  :: abstract_solver
        class(abstract_solver), intent(in)                      :: this
        integer(ip), intent(in)                                 :: first, last
        type(quadtree), dimension(first:last), intent(inout)    :: tree
        end subroutine generic_solver
    end interface

    type :: solver_factory
    	private
        character(len = 20)             :: c_solver_type
        integer(ip)                     :: i_solver_type
        class(abstract_solver), pointer :: solver
        
    contains
    procedure   :: create_solver
    end type solver_factory
    
    type, extends(abstract_solver)  :: first_order_explicit_solver
    contains
    procedure, pass(this)   :: p_solver => first_EX_solver
    end type first_order_explicit_solver
    
    type, extends(abstract_solver)  :: second_order_explicit_solver
    contains
    procedure, pass(this)   :: p_solver => second_EX_solver
    end type second_order_explicit_solver
    
    type, extends(abstract_solver)  :: runge_kutta_explicit_solver
    contains
    procedure, pass(this)   :: p_solver => RK_EX_solver
    end type runge_kutta_explicit_solver
    
    type, extends(abstract_solver)  :: second_order_implicit_solver
    contains
    procedure, pass(this)   :: p_solver => second_IM_solver
    end type second_order_implicit_solver
    
    type, extends(abstract_solver)  :: implicit_Explicit_solver
    contains
    procedure, pass(this)   :: p_solver => IMEX_solver
    end type implicit_Explicit_solver
    
!=========================== END OF FACTORY PATTERN ==========================    
    contains
!==================================================================================================
    function create_solver(this, I_solver_type) result(ptr)
    implicit none
    class(solver_factory)           :: this
    class(abstract_solver), pointer :: ptr
    integer(ip), intent(in)         :: I_solver_type
    
    if (I_solver_type == 1) then 
        ! ===> first order explicit type <===
        allocate(first_order_explicit_solver :: this%solver)
        this%I_solver_type = 1
        this%C_solver_type = "1st explicit"
        ptr => this%solver
        
    else if (I_solver_type == 2) then
        ! ===> second order explicit type <===
        allocate(second_order_explicit_solver :: this%solver)
        this%I_solver_type = 2
        this%C_solver_type = "2st explicit"
        ptr => this%solver 
        
    else if (I_solver_type == 3) then
        ! ===> runge kutta explicit type <===
        allocate(runge_kutta_explicit_solver :: this%solver)
        this%I_solver_type = 3
        this%C_solver_type = "rk explicit"
        ptr => this%solver 
        
    else if (I_solver_type == 4) then
        ! ===> second order implicit type <===
        allocate(second_order_implicit_solver :: this%solver)
        this%I_solver_type = 4
        this%C_solver_type = "2st Implicit"
        ptr => this%solver 
        
    else if (I_solver_type == 5) then
        ! ===> implicit explicit type <===
        allocate(implicit_Explicit_solver :: this%solver)
        this%I_solver_type = 5
        this%C_solver_type = "IMEX"
        ptr => this%solver        
    else
        print*, "program have not supported this solver yet"
        print*, "Stop"
        pause
        stop
    end if
    
    return
    end function create_solver
!==================================================================================================
    subroutine first_EX_solver(this, first, last, tree) 
    implicit none
    class(first_order_explicit_solver), intent(in)          :: this
    integer(ip), intent(in)                                 :: first, last
    type(quadtree), dimension(first:last), intent(inout)    :: tree
    
    return
    end subroutine first_EX_solver
!==================================================================================================
    subroutine second_EX_solver(this, first, last, tree) 
    implicit none
    class(second_order_explicit_solver), intent(in)         :: this
    integer(ip), intent(in)                                 :: first, last
    type(quadtree), dimension(first:last), intent(inout)    :: tree
    return
    end subroutine second_EX_solver    
!==================================================================================================
    subroutine RK_EX_solver(this, first, last, tree)
    implicit none
    class(runge_kutta_explicit_solver), intent(in)          :: this
    integer(ip), intent(in)                                 :: first, last
    type(quadtree), dimension(first:last), intent(inout)    :: tree
    
    return
    end subroutine RK_EX_solver 
!================================================================================================== 
    subroutine second_IM_solver(this, first, last, tree)
    implicit none
    class(second_order_implicit_solver), intent(in)         :: this
    integer(ip), intent(in)                                 :: first, last
    type(quadtree), dimension(first:last), intent(inout)    :: tree
    
    return
    end subroutine second_IM_solver     
!==================================================================================================
    subroutine IMEX_solver(this, first, last, tree)
    implicit none
    class(implicit_Explicit_solver), intent(in)             :: this
    integer(ip), intent(in)                                 :: first, last
    type(quadtree), dimension(first:last), intent(inout)    :: tree
    
    return
    end subroutine IMEX_solver
!================================================================================================== 
    subroutine solving(first, last, tree, I_solver_type)
    implicit none
    integer(ip), intent(in)                                 :: first, last
    class(quadtree), dimension(first:last), intent(inout)   :: tree
    integer(ip), intent(in)                                 :: I_solver_type
    
    type(solver_factory)            :: factory
    class(abstract_solver), pointer :: solver => null()
    
    ! ===> create solver with I_solver_type <===
    solver => factory%create_solver(I_solver_type)
    
    ! ===> call solver <===
    call solver%p_solver(first, last, tree)
    
   
    
    return
    end subroutine solving    
!==================================================================================================    
END MODULE MODULE_SOLVER    