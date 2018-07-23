MODULE MODULE_MUSCL
    
    use MODULE_PRECISION
    use MODULE_CONSTANTS
    use MODULE_GENERICMETHOD
    use MODULE_QUADTREE
    
!=========================== FACTORY PATTERN ===========================
    type, abstract  ::  limiter_funcs
    contains
    procedure(generic_funcs), deferred, pass(this)  :: l_funcs
    end type limiteR_funcs        
    
    abstract interface
        subroutine generic_funcs(this, ratio, omega, delta)
        use MODULE_PRECISION        
        import  :: limiter_funcs
        class(limiter_funcs), intent(in)    :: this
        real(rp), intent(in)                :: ratio, omega
        real(rp), intent(out)               :: delta
        end subroutine generic_funcs
    end interface
    
    type :: limiter_funcs_factory
    	private
        character(len = 20) :: limiter_type
        class(limiter_funcs), pointer   :: limiter
    contains
    
    procedure   :: create_limiter
    
    end type limiter_funcs_factory
    
    type, extends(limiter_funcs)    :: Godunov_first_order_upwind
    contains
    procedure, pass(this)   :: l_funcs => Godunov_first_order_upwind_limiter_funcs
    end type Godunov_first_order_upwind
    
    type, extends(limiter_funcs)    :: Upwind_second_order
    contains
    procedure, pass(this)   :: l_funcs => Upwind_second_order_limiter_funcs
    end type Upwind_second_order  
    
    type, extends(limiter_funcs)    :: Upwind_TVD_super_bee01
    contains
    procedure, pass(this)   :: l_funcs => Upwind_TVD_super_bee01_limiter_funcs
    end type Upwind_TVD_super_bee01 
    
    type, extends(limiter_funcs)    :: Upwind_TVD_VAN_LEER
    contains
    procedure, pass(this)   :: l_funcs => Upwind_TVD_VAN_LEER_limiter_funcs
    end type Upwind_TVD_VAN_LEER  
    
    type, extends(limiter_funcs)    :: Upwind_TVD_VAN_ALBADA
    contains
    procedure, pass(this)   :: l_funcs => Upwind_TVD_VAN_ALBADA_limiter_funcs
    end type Upwind_TVD_VAN_ALBADA    
    
    type, extends(limiter_funcs)    :: Upwind_TVD_super_bee
    contains
    procedure, pass(this)   :: l_funcs => Upwind_TVD_super_bee_limiter_funcs
    end type Upwind_TVD_super_bee         
!======================== END FACTORY PATTERN ==========================  
    


    contains
!==================================================================================================    
    function create_limiter(this, LIMITE) result(ptr)
    implicit none
    class(limiter_funcs_factory)    :: this
    class(limiter_funcs), pointer   :: ptr
    integer(ip), intent(in)  :: LIMITE
    
    if (LIMITE == 1) then 
        ! LIMITE = 1, Godunov's first order upwind method
        allocate(Godunov_first_order_upwind :: this%limiter)
        this%limiter_type = "first order upwind method"
        ptr => this%limiter
        
    else if (LIMITE == 2) then 
        ! LIMITE = 2, upwind second order method (non-monotone)
        allocate(Upwind_second_order :: this%limiter)
        ptr => this%limiter
        this%limiter_type = "second order upwind method"
    else if (LIMITE == 3) then 
        ! LIMITE = 3, upwind TVD, with SUPERBEE type limiter (old - reimann solver book (toro))
        allocate(Upwind_TVD_super_bee01 :: this%limiter)
        ptr => this%limiter
        this%limiter_type = "SUPERBEE type limiter old"
    else if (LIMITE == 4) then 
        ! LIMITE = 4, upwind TVD, with VAN LEER type limiter
        allocate(Upwind_TVD_VAN_LEER :: this%limiter)
        ptr => this%limiter
        this%limiter_type = "VAN LEER type limiter"
    else if (LIMITE == 5) then 
        ! LIMITE = 5, upwind TVD, with VAN ALBADA type limiter
        allocate(Upwind_TVD_VAN_ALBADA :: this%limiter)
        ptr => this%limiter
        this%limiter_type = "VAN ALBADA type limiter"
    else if (LIMITE == 6) then 
        ! LIMITE = 6, upwind TVD, with superbee type limiter (new - version, higher accuracy than old one )
        allocate(Upwind_TVD_super_bee :: this%limiter)
        ptr => this%limiter
        this%limiter_type = "superbee type limiter"
    else 
        write(*, *), "The program have not supported this limiter yet"
        write(*, *), "Stop."
        pause
        stop
    end if

    return
    end function create_limiter
!================================================================================================== 
    subroutine Godunov_first_order_upwind_limiter_funcs(this, ratio, omega, delta )
    class(Godunov_first_order_upwind), intent(in)   :: this
    real(rp), intent(in)    :: ratio, omega
    real(rp), intent(out)   :: delta
    
    delta = zero
    return
    end subroutine Godunov_first_order_upwind_limiter_funcs
!==================================================================================================  
    subroutine Upwind_second_order_limiter_funcs(this, ratio, omega, delta )
    class(Upwind_second_order), intent(in)   :: this
    real(rp), intent(in)    :: ratio, omega
    real(rp), intent(out)   :: delta
    
    delta = zero
    return
    end subroutine Upwind_second_order_limiter_funcs    
!==================================================================================================
    subroutine Upwind_TVD_super_bee01_limiter_funcs(this, ratio, omega, delta )
    class(Upwind_TVD_super_bee01), intent(in)   :: this
    real(rp), intent(in)    :: ratio, omega
    real(rp), intent(out)   :: delta
    
    real(rp)  phi, phir, denor
!
      phi             = zero
      if(ratio.ge.zero)phi = two*r
      if(ratio.ge.half)phi = one
!     
      if(ratio.ge.one)then
         denor = one - omega + (one + omega)*ratio
         phir  = two/denor
         phi   = min(phir, ratio)
         phi   = min(phi, two)
      endif
!
      delta = phi!*delta
    return
    end subroutine Upwind_TVD_super_bee01_limiter_funcs       
!==================================================================================================
    subroutine Upwind_TVD_VAN_LEER_limiter_funcs(this, ratio, omega, delta)
!   
!   purpose: to compute a van leer type slope limiter delta
!   
    implicit none
    class(Upwind_TVD_VAN_LEER), intent(in)  :: this
    real(rp), intent(in)    :: ratio, omega
    real(rp), intent(out)   :: delta
!   
!   declaration of variables
!   
    real(rp)  denor, phi, phir
!   
    phi = zero
!   
    if(ratio.ge.zero)then
       denor = one - omega + (one + omega)*ratio
       phir  = two/denor
       phi   = two*ratio/(one + ratio)
       phi   = min(phi, phir)
    endif
!   
    delta    = phi!*delta
!   
    end subroutine Upwind_TVD_VAN_LEER_limiter_funcs  
!==================================================================================================
    subroutine Upwind_TVD_VAN_ALBADA_limiter_funcs(this, ratio, omega, delta)
!
!   purpose: to compute a van albada type slope limiter delta
!
    implicit none
!
!   declaration of variables
!
    class(Upwind_TVD_VAN_ALBADA), intent(in)  :: this
    real(rp), intent(in)    :: ratio, omega
    real(rp), intent(out)   :: delta
    real(rp)    :: denor, phi, phir
!
    phi = zero
!
    if(ratio.ge.zero)then
       denor = one - omega + (one + omega)*ratio
       phir  = two/denor
       phi   = ratio*(one + ratio)/(one + ratio*ratio)
       phi   = min(phi, phir)
    endif
!
    delta    = phi!*delta
!
    end subroutine Upwind_TVD_VAN_ALBADA_limiter_funcs
!==================================================================================================
    subroutine Upwind_TVD_super_bee_limiter_funcs(this, ratio, omega, delta)
!
!   purpose: to compute a superbee type slope limiter delta
!
    implicit none
!
!   declaration of variables
!
    class(Upwind_TVD_super_bee),  intent(in)  :: this
    real(rp), intent(in)    :: ratio, omega
    real(rp), intent(out)   :: delta
    real(rp)  denor, phi, phir
!
      phi = max(zero, min(two*ratio, one), min(ratio, two))
    delta = phi
!
    end subroutine Upwind_TVD_super_bee_limiter_funcs    
!==================================================================================================    
    subroutine muscle(iVal, limite, first, last, tree)
    implicit none

    integer(ip), intent(in)                                 :: iVal, limite
    integer(ip), intent(in)                                 :: first, last
    type(quadtree), dimension(first:last), intent(inout)    :: tree
    type (limiter_funcs_factory)                            :: factory
    class (limiter_funcs), pointer                          :: limiter => null()    
    real(rp)    :: dupw, dloc
    integer(ip) :: i, j, k
    
    real(rp), parameter :: tollim = tolerance
    real(rp), parameter :: omega = third
    
    real(rp)    :: delta, ratio
    real(rp)    :: ratio_inv, delta_inv
    
!> create limiter

    limiter => factory%create_limiter(limite)
    
    call loop_on_quadtree_array(first, last,  tree, MUSCL_single)
!    ! X_component
!    do j = 1, ny
!        do i = 2, nx-1
!            dupw = w(i  ,j,iVal) - w(i-1,j,iVal)
!            dloc = w(i+1,j,iVal) - w(i  ,j,iVal)
!            
!            if(abs(dupw).le.tollim) dupw = tollim*sign(one,dupw)
!            if(abs(dloc).le.tollim) dloc = tollim*sign(one,dloc)
!            
!            !delta = half*(one + omega)*dupw + half*(one - omega)*dloc      !(old version muscle center low accuracy, dont use any more)
!
!                ratio = dupw/dloc
!            ratio_inv = one/ratio
!            
!            ! Slope limiters used are:
!            ! 
!            ! LIMITE = 1, Godunov's first order upwind method           ( both are fird order, need tobe checked, not important, we don't use it any more)
!            ! LIMITE = 2, upwind second order method (non-monotone)     ( both are fird order, need tobe checked, not important, we don't use it any more)
!            ! LIMITE = 3, upwind TVD, with SUPERBEE type limiter        (old superbee riemann solver book)
!            ! LIMITE = 4, upwind TVD, with VAN LEER type limiter        (             riemann solver book)
!            ! LIMITE = 5, upwind TVD, with VAN ALBADA type limiter      (             riemann solver book)
!            ! LIMITE = 6, upwind TVD, with MINMOD type limiter  >>>>> NOT WORK ANY MORE <<<<<
!            ! LIMITE = 7, upwind TVD, with MINMAX type limiter  >>>>> NOT WORK ANY MORE <<<<<
!            ! LIMITE = 6, new superbee (high accuracy source: https://en.wikipedia.org/wiki/Flux_limiter)
!            call limiter%l_funcs(ratio, omega, delta )    
!            delta_inv = delta/ratio
!            
!            muscl_recons(i,j,ival)%x_l = w(i,j,ival) - fourth*((one + omega)*delta_inv*dupw + (one - omega)*delta*dloc)
!            muscl_recons(i,j,ival)%x_r = w(i,j,ival) + fourth*((one - omega)*delta_inv*dupw + (one + omega)*delta*dloc)
!
!        end do
!    end do

    ! y_component
!    do j = 2, ny-1
!        do i = 1, nx
!            dupw = w(i,j  ,iVal) - w(i,j-1,iVal)
!            dloc = w(i,j+1,iVal) - w(i,j  ,iVal)
!            
!            if(abs(dupw).le.tollim)dupw=tollim*sign(one,dupw)
!            if(abs(dloc).le.tollim)dloc=tollim*sign(one,dloc)
!            
!            !delta = half*(one + omega)*dupw + half*(one - omega)*dloc
!
!                ratio = dupw/dloc
!            ratio_inv = one/ratio
!            
!            ! Slope limiters used are:
!            ! 
!            ! LIMITE = 1, Godunov's first order upwind method
!            ! LIMITE = 2, upwind second order method (non-monotone)
!            ! LIMITE = 3, upwind TVD, with SUPERBEE type limiter
!            ! LIMITE = 4, upwind TVD, with VAN LEER type limiter
!            ! LIMITE = 5, upwind TVD, with VAN ALBADA type limiter
!            ! LIMITE = 6, upwind TVD, with MINMOD type limiter
!            ! LIMITE = 7, upwind TVD, with MINMAX type limiter
!            call limiter%l_funcs(ratio, omega, delta )    ; delta_inv = delta/ratio
!
!            muscl_recons(i,j,ival)%y_l = w(i,j,ival) - fourth*((one + omega)*delta_inv*dupw + (one - omega)*delta*dloc)
!            muscl_recons(i,j,ival)%y_r = w(i,j,ival) + fourth*((one - omega)*delta_inv*dupw + (one + omega)*delta*dloc)
            
!        end do
!    end do
    
!     muscl_recons(1,:,:) = muscl_recons(2,:,:)
!    muscl_recons(nx,:,:) = muscl_recons(nx-1,:,:)
!     muscl_recons(:,1,:) = muscl_recons(:,2,:)
!    muscl_recons(:,ny,:) = muscl_recons(:,ny-1,:)
    
!> delete limiter    
    deallocate(factory%limiter)
    return
    contains

    subroutine  MUSCL_single(tree)
    implicit none
    type(quadtree), pointer, intent(inout)  :: tree
    real(rp)                                :: delta
    
    type(quadtree), pointer                 :: curr_C, adj_c
    real(rp), dimension(2)                  :: dr

    curr_c => tree
     adj_c => tree%adj_north
    
    delta = adj_c%data%w(iVal) - curr_c%data%w(iVal)
    
    dr = half*(curr_c%pts(1)%coord(:) + curr_c%pts(2)%coord(:)) - curr_c%pts(5)%coord(:)
    
    
    ! ===> compute limiter funcs 
    !call limiter%l_funcs(ratio, omega, delta)
    
    ! ===> compute muscl reconstruction
    tree%data%recons(iVal)%x_r = tree%data%w(iVal) + half*delta
    
    return
    end subroutine MUSCL_single
    
    end subroutine muscle
    
!==================================================================================================

!==================================================================================================
!
      subroutine sbslic(r, omega, delta)
!
!     purpose: to compute a superbee type slope limiter delta
!
      implicit none
!
!     declaration of variables
!
      real(rp)  delta, denor, omega, phi, phir, r
!
      phi             = zero
      if(r.ge.zero)phi = two*r
      if(r.ge.half)phi = one
!     
      if(r.ge.one)then
         denor = one - omega + (one + omega)*r
         phir  = two/denor
         phi   = min(phir, r)
         phi   = min(phi, two)
      endif

!
      delta = phi!*delta
!
      end
!
!==================================================================================================
!
      subroutine superbee(r, omega, delta)
!
!     purpose: to compute a superbee type slope limiter delta
!
      implicit none
!z
!     declaration of variables
!
      real(rp)  delta, denor, omega, phi, phir, r
!
        phi = zero
        
        phi = max(zero, min(two*r, one), min(r, two))
      delta = phi
!
      end
!
!==================================================================================================     
      
!
      subroutine vlslic(r, omega, delta)
!
!     purpose: to compute a van leer type slope limiter delta
!
      implicit none
!
!     declaration of variables
!
      real(rp)  delta, denor, omega, phi, phir, r
!
      phi = zero
!
      if(r.ge.zero)then
         denor = one - omega + (one + omega)*r
         phir  = two/denor
         phi   = two*r/(one + r)
         phi   = min(phi, phir)
      endif
!
      delta    = phi!*delta
!
      end
!
!==================================================================================================
!
      subroutine vaslic(r, omega, delta)
!
!     purpose: to compute a van albada type slope limiter delta
!
      implicit none
!
!     declaration of variables
!
      real(rp)  delta, denor, omega, phi, phir, r
!
      phi = zero
!
      if(r.ge.zero)then
         denor = one - omega + (one + omega)*r
         phir  = two/denor
         phi   = r*(one + r)/(one + r*r)
         phi   = min(phi, phir)
      endif
!
      delta    = phi!*delta
!
      end
!
!==================================================================================================
!
      subroutine mislic(r, omega, delta)
!
!     purpose: to compute a minmod type slope limiter delta
!
      implicit none
!
!     declaration of variables
!
      real(rp)  delta, denor, omega, phi, phir, r
!
      phi             = zero
      if(r.ge.zero)phi = r
!
      if(r.ge.one)then
         denor = two*(one - omega + (one + omega)*r)
         phir  = four/denor
         phi   = min(one, phir)
      endif
!
      delta    = phi!*delta
!
      end
!
!==================================================================================================
!
      subroutine minmax(dupw, dloc, delta)
!
!     purpose: to compute a minmax type slope limiter delta.
!              this is the most diffusive of all limiters
!              for centred schemes
!
      implicit none
!
!     declaration of variables
!
      real(rp)  betal, betar, delta, dloc, dupw, signo
!
      betal = one
      betar = one
      signo = half*(sign(one,dupw) + sign(one,dloc))
      delta = signo*(min(betal*abs(dupw),betar*abs(dloc)))
!
      end
!
!================================================================================================== 
END MODULE MODULE_MUSCL
