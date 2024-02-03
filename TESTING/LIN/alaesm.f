      SUBROUTINE ALAESM( PATH, OK, NOUT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      LOGICAL            OK
      CHARACTER*3        PATH
      int                NOUT
*     ..
*
*  =====================================================================
*
*     .. Executable Statements ..
*
      IF( OK ) THEN
         WRITE( NOUT, FMT = 9999 )PATH
      ELSE
         WRITE( NOUT, FMT = 9998 )PATH
      END IF
*
 9999 FORMAT( 1X, A3, ' routines passed the tests of the error exits'
     $       )
 9998 FORMAT( ' *** ', A3, ' routines failed the tests of the error ',
     $      'exits ***' )
      RETURN
*
*     End of ALAESM
*
      END
