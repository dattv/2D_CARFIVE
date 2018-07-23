MODULE MODULE_CFDUTILITY
    
    use MODULE_PRECISION
    use MODULE_QUADTREE
    use MODULE_MATHUTILITY
    
    contains
!==================================================================================================    
    subroutine weighted_least_square(ival, tree, w)
    implicit none
    
    integer(ip), intent(in)             :: ival
    type(quadtree), intent(in)          :: tree
    real(rp), dimension(2), intent(out) :: w
    type(quadtree), pointer             :: cell_n, cell_s, cell_e, cell_w
    real(rp), dimension(2,2)            :: A
    real(rp), dimension(2)              :: B
    
    cell_n => tree%adj_north
    cell_e => tree%adj_east
    cell_s => tree%adj_south
    cell_w => tree%adj_west
    
    A(1,1) = compute_component(compute_alpha(cell_n%pts(5)%coord(1), tree%pts(5)%coord(1), cell_n%pts(5)%coord(2), tree%pts(5)%coord(2))**2, cell_n%pts(5)%coord(1) - tree%pts(5)%coord(1), cell_n%pts(5)%coord(1) - tree%pts(5)%coord(1)) + &
             compute_component(compute_alpha(cell_e%pts(5)%coord(1), tree%pts(5)%coord(1), cell_e%pts(5)%coord(2), tree%pts(5)%coord(2))**2, cell_e%pts(5)%coord(1) - tree%pts(5)%coord(1), cell_e%pts(5)%coord(1) - tree%pts(5)%coord(1)) + &
             compute_component(compute_alpha(cell_s%pts(5)%coord(1), tree%pts(5)%coord(1), cell_s%pts(5)%coord(2), tree%pts(5)%coord(2))**2, cell_s%pts(5)%coord(1) - tree%pts(5)%coord(1), cell_s%pts(5)%coord(1) - tree%pts(5)%coord(1)) + &
             compute_component(compute_alpha(cell_w%pts(5)%coord(1), tree%pts(5)%coord(1), cell_w%pts(5)%coord(2), tree%pts(5)%coord(2))**2, cell_w%pts(5)%coord(1) - tree%pts(5)%coord(1), cell_w%pts(5)%coord(1) - tree%pts(5)%coord(1))
    
    A(1,2) = compute_component(compute_alpha(cell_n%pts(5)%coord(1), tree%pts(5)%coord(1), cell_n%pts(5)%coord(2), tree%pts(5)%coord(2))**2, cell_n%pts(5)%coord(1) - tree%pts(5)%coord(1), cell_n%pts(5)%coord(2) - tree%pts(5)%coord(2)) + &
             compute_component(compute_alpha(cell_e%pts(5)%coord(1), tree%pts(5)%coord(1), cell_e%pts(5)%coord(2), tree%pts(5)%coord(2))**2, cell_e%pts(5)%coord(1) - tree%pts(5)%coord(1), cell_e%pts(5)%coord(2) - tree%pts(5)%coord(2)) + &
             compute_component(compute_alpha(cell_s%pts(5)%coord(1), tree%pts(5)%coord(1), cell_s%pts(5)%coord(2), tree%pts(5)%coord(2))**2, cell_s%pts(5)%coord(1) - tree%pts(5)%coord(1), cell_s%pts(5)%coord(2) - tree%pts(5)%coord(2)) + &
             compute_component(compute_alpha(cell_w%pts(5)%coord(1), tree%pts(5)%coord(1), cell_w%pts(5)%coord(2), tree%pts(5)%coord(2))**2, cell_w%pts(5)%coord(1) - tree%pts(5)%coord(1), cell_w%pts(5)%coord(2) - tree%pts(5)%coord(2))
    
    A(2,1) = A(1,2)
    
    A(2,2) = compute_component(compute_alpha(cell_n%pts(5)%coord(1), tree%pts(5)%coord(1), cell_n%pts(5)%coord(2), tree%pts(5)%coord(2))**2, cell_n%pts(5)%coord(2) - tree%pts(5)%coord(2), cell_n%pts(5)%coord(2) - tree%pts(5)%coord(2)) + &
             compute_component(compute_alpha(cell_e%pts(5)%coord(1), tree%pts(5)%coord(1), cell_e%pts(5)%coord(2), tree%pts(5)%coord(2))**2, cell_e%pts(5)%coord(2) - tree%pts(5)%coord(2), cell_e%pts(5)%coord(2) - tree%pts(5)%coord(2)) + &
             compute_component(compute_alpha(cell_s%pts(5)%coord(1), tree%pts(5)%coord(1), cell_s%pts(5)%coord(2), tree%pts(5)%coord(2))**2, cell_s%pts(5)%coord(2) - tree%pts(5)%coord(2), cell_s%pts(5)%coord(2) - tree%pts(5)%coord(2)) + &
             compute_component(compute_alpha(cell_w%pts(5)%coord(1), tree%pts(5)%coord(1), cell_w%pts(5)%coord(2), tree%pts(5)%coord(2))**2, cell_w%pts(5)%coord(2) - tree%pts(5)%coord(2), cell_w%pts(5)%coord(2) - tree%pts(5)%coord(2))
    
    B(1) = compute_component(compute_alpha(cell_n%pts(5)%coord(1), tree%pts(5)%coord(1), cell_n%pts(5)%coord(2), tree%pts(5)%coord(2))**2, cell_n%pts(5)%coord(1) - tree%pts(5)%coord(1), cell_n%data%w(ival) - tree%data%w(ival)) + &
           compute_component(compute_alpha(cell_e%pts(5)%coord(1), tree%pts(5)%coord(1), cell_e%pts(5)%coord(2), tree%pts(5)%coord(2))**2, cell_e%pts(5)%coord(1) - tree%pts(5)%coord(1), cell_e%data%w(ival) - tree%data%w(ival)) + &
           compute_component(compute_alpha(cell_s%pts(5)%coord(1), tree%pts(5)%coord(1), cell_s%pts(5)%coord(2), tree%pts(5)%coord(2))**2, cell_s%pts(5)%coord(1) - tree%pts(5)%coord(1), cell_s%data%w(ival) - tree%data%w(ival)) + &
           compute_component(compute_alpha(cell_w%pts(5)%coord(1), tree%pts(5)%coord(1), cell_w%pts(5)%coord(2), tree%pts(5)%coord(2))**2, cell_w%pts(5)%coord(1) - tree%pts(5)%coord(1), cell_w%data%w(ival) - tree%data%w(ival))
    
    B(2) = compute_component(compute_alpha(cell_n%pts(5)%coord(1), tree%pts(5)%coord(1), cell_n%pts(5)%coord(2), tree%pts(5)%coord(2))**2, cell_n%pts(5)%coord(2) - tree%pts(5)%coord(2), cell_n%data%w(ival) - tree%data%w(ival)) + &
           compute_component(compute_alpha(cell_e%pts(5)%coord(1), tree%pts(5)%coord(1), cell_e%pts(5)%coord(2), tree%pts(5)%coord(2))**2, cell_e%pts(5)%coord(2) - tree%pts(5)%coord(2), cell_e%data%w(ival) - tree%data%w(ival)) + &
           compute_component(compute_alpha(cell_s%pts(5)%coord(1), tree%pts(5)%coord(1), cell_s%pts(5)%coord(2), tree%pts(5)%coord(2))**2, cell_s%pts(5)%coord(2) - tree%pts(5)%coord(2), cell_s%data%w(ival) - tree%data%w(ival)) + &
           compute_component(compute_alpha(cell_w%pts(5)%coord(1), tree%pts(5)%coord(1), cell_w%pts(5)%coord(2), tree%pts(5)%coord(2))**2, cell_w%pts(5)%coord(2) - tree%pts(5)%coord(2), cell_w%data%w(ival) - tree%data%w(ival))
    
    w =  matmul(inv(A), B)
    return
    end subroutine weighted_least_square
!==================================================================================================    
    function compute_component(alpha, del_x, del_y) result(res)
    implicit none
    real(rp), intent(in)    :: alpha
    real(rp), intent(in)    :: del_x
    real(rp), intent(in)    :: del_y
    real(rp)                :: res
    
    res = alpha**2 * del_x * del_y
    return
    end function compute_component 
!==================================================================================================    
    function compute_alpha(xi, xc, yi, yc) result(res)
    implicit none
    real(rp), intent(in)    :: xi, xc, yi, yc
    real(rp)                :: res
    
    
    res = one/(sqrt(xi-xc)**2 + (yi-yc)**2)
    return
    end function compute_alpha
END MODULE MODULE_CFDUTILITY    