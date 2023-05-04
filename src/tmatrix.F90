module tmatrix
    use tmatrix_extraction
    use constants
    use iso_c_binding!, only: C_FLOAT128, C_FLOAT128_COMPLEX, C_INT32_T
    implicit none

!DEC$ IF (HAVE_QUAD_PRECISION)
#define HAVE_QUAD_PRECISION .true.
#if defined(HAVE_QUAD_PRECISION)
    integer, parameter :: C_FLOAT_KND = C_FLOAT128
    integer, parameter :: C_COMPLEX_KND = C_FLOAT128_COMPLEX
!DEC$ ELSE
#else
    integer, parameter :: C_FLOAT_KND = C_DOUBLE
    integer, parameter :: C_COMPLEX_KND = C_DOUBLE_COMPLEX
!DEC$ ENDIF
#endif
contains

    subroutine get_spheroidal_tmatrix(nol, ab, rv, lambda, lnum, ri, m, res)
        real(C_FLOAT_KND) :: ab(nol), rv(nol), lambda
        integer(C_INT32_T) :: lnum, m, f, i, j, nol
        ! type(ScatteringResult) :: response
        complex(C_COMPLEX_KND) :: ri(nol), res(2*lnum, 2*lnum)

        write(*,*) 'fortran ab = ', ab, 'rv = ', rv, 'lambda = ', lambda, 'lnum = ', lnum, &
                'ri = ', ri, 'm = ', m
        ! f = 1
        ! if (ab < 1.0q0) then
        !     f = -1
        !     ab = 1q0 / ab
        ! end if

        ! response = GetSphericalTmatrix(f, 1, (/2.0q0 * PI * rv / lambda/), (/ab/), 0.5q0, &
        !         lambda, &
        ! (/ (1q0, 0q0), ri/), lnum, m, m)
        call GetSpheroidalTmatrix(nol, ab, rv, lambda, lnum, ri, m, res)
!        write( *,*) 'fortran result:'
!        do i = 1, 2*lnum
!            write(*,*) response%result_spherical_tm%uv_tmatrix(i,:,m)
!        end do
        res = transpose(res)
        ! write(*,*) 'result in fort: ', res(:,1)
        ! deallocate(response%result_spherical_tm%uv_tmatrix, response%result_spherical_te%uv_tmatrix)
    end subroutine get_spheroidal_tmatrix
end module tmatrix