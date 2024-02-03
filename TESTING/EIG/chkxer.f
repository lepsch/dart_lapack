      SUBROUTINE CHKXER( SRNAMT, INFOT, NOUT, LERR, OK )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      LOGICAL            LERR, OK
      CHARACTER*(*)      SRNAMT
      INTEGER            INFOT, NOUT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          LEN_TRIM
*     ..
*     .. Executable Statements ..
      IF( .NOT.LERR ) THEN
         WRITE( NOUT, FMT = 9999 )INFOT, SRNAMT( 1:LEN_TRIM( SRNAMT ) )
         OK = .FALSE.
      END IF
      LERR = .FALSE.
      RETURN
*
 9999 FORMAT( ' *** Illegal value of parameter number ', I2,
     $      ' not detected by ', A, ' ***' )
*
*     End of CHKXER
*
      END
