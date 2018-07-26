!> @file src/util.f90
!------------------------------------------------------------------------------!
! This file is part of the PALM model system.
!
! PALM is free software: you can redistribute it and/or modify it under the
! terms of the GNU General Public License as published by the Free Software
! Foundation, either version 3 of the License, or (at your option) any later
! version.
!
! PALM is distributed in the hope that it will be useful, but WITHOUT ANY
! WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
! A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License along with
! PALM. If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 2017-2018 Leibniz Universitaet Hannover
! Copyright 2017-2018 Deutscher Wetterdienst Offenbach
!------------------------------------------------------------------------------!
!
! Current revisions:
! -----------------
! 
! 
! Former revisions:
! -----------------
! $Id: util.f90 2718 2018-01-02 08:49:38Z maronga $
! Initial revision
!
! 
!
! Authors:
! --------
! @author Eckhard Kadasch
!
! Description:
! ------------
!> The util module provides miscellaneous utility routines for INIFOR.
!------------------------------------------------------------------------------!
 MODULE util

    USE, INTRINSIC :: ISO_C_BINDING,                                           &
        ONLY :  C_CHAR, C_INT, C_PTR, C_SIZE_T
    USE defs,                                                                  &
        ONLY :  dp, PI, DATE

    IMPLICIT NONE

    TYPE, BIND(c) :: tm_struct
       INTEGER(C_INT) :: tm_sec     !< seconds after the minute [0, 61]
       INTEGER(C_INT) :: tm_min     !< minutes after the hour [0, 59]
       INTEGER(C_INT) :: tm_hour    !< hours since midnight [0, 23]
       INTEGER(C_INT) :: tm_mday    !< day of the month [1, 31]
       INTEGER(C_INT) :: tm_mon     !< month since January [0, 11]
       INTEGER(C_INT) :: tm_year    !< years since 1900
       INTEGER(C_INT) :: tm_wday    !< days since Sunday [0, 6]
       INTEGER(C_INT) :: tm_yday    !< days since January 1st [0, 356]
       INTEGER(C_INT) :: tm_isdst   !< Daylight Saving Time flag
    END TYPE

    INTERFACE

       FUNCTION strptime(string, format, timeinfo) BIND(c, NAME='strptime')
          IMPORT :: C_CHAR, C_SIZE_T, tm_struct

          IMPLICIT NONE

          CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(IN) ::  string, format
          TYPE(tm_struct), INTENT(OUT)                     ::  timeinfo

          INTEGER(C_SIZE_T)                                ::  strptime
       END FUNCTION


       FUNCTION strftime(string, string_len, format, timeinfo) BIND(c, NAME='strftime')
          IMPORT :: C_CHAR, C_SIZE_T, tm_struct

          IMPLICIT NONE

          CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(OUT) ::  string
          CHARACTER(KIND=C_CHAR), DIMENSION(*), INTENT(IN)  ::  format
          INTEGER(C_SIZE_T), INTENT(IN)                     ::  string_len
          TYPE(tm_struct), INTENT(IN)                       ::  timeinfo

          INTEGER(C_SIZE_T)                                 ::  strftime
       END FUNCTION


       FUNCTION mktime(timeinfo) BIND(c, NAME='mktime')
          IMPORT :: C_PTR, tm_struct

          IMPLICIT NONE

          TYPE(tm_struct), INTENT(IN) ::  timeinfo

          TYPE(C_PTR)                 ::  mktime
       END FUNCTION

    END INTERFACE

 CONTAINS

    CHARACTER(LEN=DATE) FUNCTION add_hours_to(date_string, hours)
       CHARACTER(LEN=DATE), INTENT(IN)          ::  date_string
       INTEGER, INTENT(IN)                      ::  hours

       CHARACTER(KIND=C_CHAR, LEN=*), PARAMETER ::  format_string = "%Y%m%d%H"
       CHARACTER(KIND=C_CHAR, LEN=DATE)         ::  c_date_string
       TYPE(C_PTR)                              ::  c_pointer
       TYPE(tm_struct)                          ::  time_info
       INTEGER                                  ::  err

       c_date_string = date_string

       ! Convert C string to C tm struct
       CALL init_tm(time_info)
       err = strptime(c_date_string, format_string, time_info)
    
       ! Manipulate and normalize C tm struct
       time_info % tm_hour = time_info % tm_hour + hours
       c_pointer = mktime(time_info)

       ! Convert back to C string
       err = strftime(c_date_string, INT(DATE, KIND=C_SIZE_T),                 &
                      format_string, time_info)

       add_hours_to = c_date_string
    END FUNCTION


    SUBROUTINE print_tm(timeinfo)
       TYPE(tm_struct), INTENT(IN) :: timeinfo

       PRINT *, "sec: ", timeinfo % tm_sec,  &  !< seconds after the minute [0, 61]
                "min: ", timeinfo % tm_min,  &  !< minutes after the hour [0, 59]
                "hr:  ", timeinfo % tm_hour, &  !< hours since midnight [0, 23]
                "day: ", timeinfo % tm_mday, &  !< day of the month [1, 31]
                "mon: ", timeinfo % tm_mon,  &  !< month since January [0, 11]
                "yr:  ", timeinfo % tm_year, &  !< years since 1900
                "wday:", timeinfo % tm_wday, &  !< days since Sunday [0, 6]
                "yday:", timeinfo % tm_yday, &  !< days since January 1st [0, 356]
                "dst: ", timeinfo % tm_isdst    !< Daylight Saving time flag
    END SUBROUTINE print_tm

    
    SUBROUTINE init_tm(timeinfo)
       TYPE(tm_struct), INTENT(INOUT) :: timeinfo

       timeinfo % tm_sec   = 0
       timeinfo % tm_min   = 0
       timeinfo % tm_hour  = 0
       timeinfo % tm_mday  = 0
       timeinfo % tm_mon   = 0
       timeinfo % tm_year  = 0
       timeinfo % tm_wday  = 0
       timeinfo % tm_yday  = 0

       ! We use UTC times, so marking Daylight Saving Time (DST) 'not available'
       ! (< 0). If this is set to 0, mktime will convert the timeinfo to DST and
       ! add one hour.
       timeinfo % tm_isdst = -1
    END SUBROUTINE init_tm


    SUBROUTINE fake_output_3d(a)

       REAL(dp), INTENT(INOUT)       ::  a(:,:,:)
       REAL(dp)                      ::  lxi, lyi
       INTEGER ::  i,j,k

       lyi = 2.0_dp * PI / (SIZE(a, 2) - 1.0_dp)
       lxi = 2.0_dp * PI / (SIZE(a, 1) - 1.0_dp)

       DO k = 1, SIZE(a, 3)
       DO j = 1, SIZE(a, 2)
       DO i = 1, SIZE(a, 1)
           a(i,j,k) = SIN(lxi * i) * COS(lyi * j) + k
       END DO
       END DO
       END DO

    END SUBROUTINE fake_output_3d


    SUBROUTINE fake_output_2d(a, offset)

       REAL(dp), INTENT(INOUT) ::  a(:,:)
       INTEGER, INTENT(IN)     ::  offset
       REAL(dp)                ::  lxi, lyi
       INTEGER                 ::  i,j

       lyi = 2.0_dp*PI / (SIZE(a, 2) - 1.0_dp)
       lxi = 2.0_dp*PI / (SIZE(a, 1) - 1.0_dp)

       a(:,:) = 1.0_dp
       DO j = 1, SIZE(a, 2)
       DO i = 1, SIZE(a, 1)
          a(i,j) = SIN(lxi * i) * COS(lyi * j) + offset
       END DO
       END DO

    END SUBROUTINE fake_output_2d


    SUBROUTINE linspace(start, stop, array)

       REAL(dp), INTENT(IN)    ::  start, stop
       REAL(dp), INTENT(INOUT) ::  array(0:)
       INTEGER                 ::  i, n

       n = UBOUND(array, 1)

       IF (n .EQ. 0)  THEN

          array(0) = start

       ELSE

          DO i = 0, n
             array(i) = start + REAL(i, dp) / n * (stop - start)
          END DO

       END IF
       
    END SUBROUTINE linspace


!------------------------------------------------------------------------------!
! Description:
! ------------
!> Reverse the order of the third (vertical) array dimension from top-down
!> (COSMO) to bottom-up (PALM)
!------------------------------------------------------------------------------!
    SUBROUTINE reverse(input_arr)

       REAL(dp), INTENT(INOUT) ::  input_arr(:,:,:)

       input_arr = input_arr(:,:,size(input_arr, 3):1:-1)

    END SUBROUTINE reverse


    SUBROUTINE deaverage(avg_1, t1, avg_2, t2, avg_3, t3)

       REAL(dp), DIMENSION(:,:,:), INTENT(IN)  ::  avg_1, avg_2
       REAL(dp), INTENT(IN)                    ::  t1, t2, t3
       REAL(dp), DIMENSION(:,:,:), INTENT(OUT) ::  avg_3

       REAL(dp)                                ::  ti
 
       ti = 1.0_dp / t3

       avg_3(:,:,:) = ti * ( t2 * avg_2(:,:,:) - t1 * avg_1(:,:,:) )

    END SUBROUTINE deaverage


    SUBROUTINE get_basic_state(z, beta, p_sl, t_sl, rd, g, p0)

       REAL(dp), INTENT(IN)  ::  z(1:)  !< height [m]
       REAL(dp), INTENT(IN)  ::  beta   !< logarithmic lapse rate, dT / d ln(p) [K]
       REAL(dp), INTENT(IN)  ::  p_sl   !< reference pressure [Pa]
       REAL(dp), INTENT(IN)  ::  t_sl   !< reference tempereature [K]
       REAL(dp), INTENT(IN)  ::  rd     !< ideal gas constant of dry air [J/kg/K]
       REAL(dp), INTENT(IN)  ::  g      !< acceleration of Earth's gravity [m/s^2]
       REAL(dp), INTENT(OUT) ::  p0(1:) !< COSMO-DE basic state pressure [Pa]
       REAL(dp) ::  root_frac, factor   !< precomputed factors

       factor = - t_sl / beta
       root_frac = (2.0_dp * beta * g) / (rd * t_sl*t_sl)

       p0(:) = p_sl * EXP(                                                     &
                  factor * ( 1.0_dp - SQRT( 1.0_dp - root_frac * z(:) ) )  &
               )

    END SUBROUTINE get_basic_state


    ! Convert a real number to a string in scientific notation
    ! showing four significant digits.
    CHARACTER(LEN=11) FUNCTION real_to_str(val, format)

        REAL(dp), INTENT(IN)                   ::  val
        CHARACTER(LEN=*), OPTIONAL, INTENT(IN) ::  format

        IF (PRESENT(format))  THEN
           WRITE(real_to_str, TRIM(format)) val
        ELSE
           WRITE(real_to_str, '(E11.4)') val
           real_to_str = ADJUSTL(real_to_str)
        END IF

    END FUNCTION real_to_str


    CHARACTER(LEN=12) FUNCTION real_to_str_f(val)

        REAL(dp), INTENT(IN) ::  val

        WRITE(real_to_str_f, '(F12.4)') val
        real_to_str_f = ADJUSTL(real_to_str_f)

    END FUNCTION real_to_str_f


    CHARACTER(LEN=10) FUNCTION str(val)

        INTEGER, INTENT(IN) ::  val

        WRITE(str, '(i10)') val
        str = ADJUSTL(str)

    END FUNCTION str


    CHARACTER(LEN=30) FUNCTION strs(vals)

        INTEGER, INTENT(IN) ::  vals(:)
        INTEGER ::  i

        strs = ''
        DO i = 1, SIZE(vals)
           strs = TRIM(strs) // TRIM(str(vals(i)))
        END DO

    END FUNCTION strs


    SUBROUTINE normalize_path(path)
        
        CHARACTER(LEN=*), INTENT(INOUT) ::  path
        INTEGER ::  n

        n = LEN_TRIM(path)

        IF (path(n:n) .NE. '/')  THEN
           path = TRIM(path) // '/'
        END IF

    END SUBROUTINE

 END MODULE

