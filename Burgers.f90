module mod_file
    !implicit none
    use, intrinsic :: iso_c_binding 
    include 'fftw3.f03'
    integer*8,parameter::N=10**3, Nh = N/2+1  
    real(kind=8),parameter::pi=4.0*atan(1.0)  , dt = 10**(-4.0)
    real(kind=8),parameter::l=2.0*pi , dx = l/N , nu =10.0**(-2.0)
    complex(kind=8),parameter:: zi = (0.0,1.0) 
    integer(kind=8) :: i , ti , MaxIter = 10**5
    real(kind=8),dimension(N)::u,x
    real(kind=8),dimension(Nh)::k
    complex(kind=8),dimension(Nh)::uk, dduk , nonlin
    integer(kind=8)::plan_forward,plan_backward 
    real(kind=4):: t1,t2
end module
program burgers
    use mod_file
    call cpu_time(t1)
    do i=1,N
        x(i) = dfloat(i-1)*dx
        u(i) = sin(x(i))
        !print*,i,u(i)
    end do
    do i=1,Nh
        k(i)= 2.0*pi*dfloat(i-1)/l
    end do

    call dfftw_plan_dft_r2c_1d(plan_forward,N,u,uk,FFTW_ESTIMATE)
    call dfftw_execute(plan_forward)
    !call dfftw_destroy_plan(plan_forward)
    uk = uk/N
    call dfftw_plan_dft_c2r_1d(plan_backward,N,uk,u,FFTW_ESTIMATE+FFTW_PRESERVE_INPUT)

    do ti = 1,MaxIter-1
        !===================================================================
        !--------------------------------------------------------------------
        call rk4
        !-------------------------------------------------------------------    
        if(mod(ti,1000)==0)then
            call dfftw_execute_dft_c2r(plan_backward,uk,u)
            do i=1,N
                write(ti,*)x(i),u(i)
            end do
            print*,'Time=',ti*dt,'Energy=',sum(u**2)*dx
        end if

        
    end do
    call dfftw_destroy_plan(plan_backward)
    call dfftw_destroy_plan(plan_forward)
    call cpu_time(t2)
    print*,'Time Elapsed=' ,t2-t1
end program burgers


subroutine rhs
    use mod_file
    !--------------- dealiasing----------------------------------------
    uk(int(Nh)-int(Nh/3.0):int(Nh))= 0
    call dfftw_execute_dft_c2r(plan_backward,uk,u)
    !-------------------------Linear Term-------------------------------------
    dduk = -k**2*uk
    !!-----------------------Nolinear Term--------------------------------------
    
    call dfftw_execute_dft_r2c(plan_forward,u*u,nonlin)
    nonlin = -zi*k*nonlin/dfloat(N) -nu*dduk
    nonlin(int(Nh)-int(Nh/3.0):int(Nh))=0

    !!--------------------------------------------------------------------------
end subroutine rhs
subroutine euler
    use mod_file
    !---------------------------------------------------------------------------
    call rhs
    uk =  exp(-nu*k**2*dt)*(uk + 0.5*nonlin*dt) 
end subroutine euler
subroutine RK4
    use mod_file
    complex(kind=8),dimension(Nh):: k1,k2,k3,k4,u_old
        u_old = uk
        call rhs
        k1 = nonlin - nu*uk*k**2 
        uk = u_old + 0.5d0*dt*k1
        call rhs
        k2 = nonlin - nu*uk*k**2 
        uk = u_old + 0.5d0*dt*k2
        call rhs
        k3 = nonlin - nu*uk*k**2 
        uk = u_old + dt*k3
        call rhs
        k4 = nonlin - nu*uk*k**2 

        uk = u_old + dt*(1d0/6.0d0)*(k1+2.0*k2+2.0*k3+k4)

end subroutine RK4
