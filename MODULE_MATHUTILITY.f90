MODULE MODULE_MATHUTILITY
    
    use MODULE_PRECISION
    
    contains
    
    function inv(A) result(res)
    implicit none
    
    real(rp), dimensiON(2,2), intent(in)    :: A
    real(rp), dimension(2,2)                :: res
    
    res(1,1) =  A(2,2)/(A(1,1)*A(2,2) - A(1,2)*A(2,1))
    res(1,2) = -A(1,2)/(A(1,1)*A(2,2) - A(1,2)*A(2,1))
    res(2,1) = -A(2,1)/(A(1,1)*A(2,2) - A(1,2)*A(2,1))
    res(2,2) =  A(1,1)/(A(1,1)*A(2,2) - A(1,2)*A(2,1))
    
    return
    end function inv
    
END MODULE MODULE_MATHUTILITY    