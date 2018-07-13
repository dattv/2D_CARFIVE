!=================================================================================================
!> CARTESIAN ADAPTIVE FIVE EQUATION MODEL
!> AUTHOR: VAN-DAT THANG
!> E-MAIL: datthangva@gmail.com
!> E-MAIL: vandatthang@gmail.com
!> SOURCE CODE LINK: https://github.com/dattv/2D_CARFIVE
!================================================================================================= 
MODULE MODULE_SOLVER
    
    use MODULE_PRECISION
    
!================== DEGINS PATTERN ==================    
    
    type, abstract  :: abstract_solver
    contains
    procedure(generic_solver), deferred, pass(this) :: p_solver
    end type abstract_solver
    
    abstract interface
        subroutine generic_solver(this)
        use MODULE_PRECISION
        import  :: abstract_solver
        class(abstract_solver), intent(in)   :: this
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
    
!    type, extends(abstract_solver)  :: second_order_explicit_solver
!    contains
!    procedure, pass(this)   :: p_solver => 2ST_EX_solver
!    end type second_order_explicit_solver
!    
!    type, extends(abstract_solver)  :: runge_kutta_explicit_solver
!    contains
!    procedure, pass(this)   :: p_solver => RK_EX_solver
!    end type runge_kutta_explicit_solver
!    
!    type, extends(abstract_solver)  :: second_order_implicit_solver
!    contains
!    procedure, pass(this)   :: p_solver => 2ST_IM_solver
!    end type second_order_implicit_solver
!    
!    type, extends(abstract_solver)  :: implicit_Explicit_solver
!    contains
!    procedure, pass(this)   :: p_solver => IMEX_solver
!    end type implicit_Explicit_solver
    
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
        
    else if (I_solver_type == 2) then 
    else if (I_solver_type == 3) then 
    else if (I_solver_type == 4) then 
    else if (I_solver_type == 5) then
    else
    end if
    
    return
    end function create_solver
!==================================================================================================
    subroutine first_EX_solver(this) 
    implicit none
    class(first_order_explicit_solver), intent(in)   :: this
    
    return
    end subroutine first_EX_solver
!==================================================================================================    
!==================================================================================================    
!==================================================================================================    
!==================================================================================================    
!==================================================================================================    
!==================================================================================================    
END MODULE MODULE_SOLVER    