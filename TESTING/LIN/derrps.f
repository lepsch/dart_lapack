      SUBROUTINE DERRPS( PATH, NUNIT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                NUNIT
      String             PATH;
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      int                NMAX
      PARAMETER          ( NMAX = 4 )
*     ..
*     .. Local Scalars ..
      int                I, INFO, J, RANK
*     ..
*     .. Local Arrays ..
      DOUBLE PRECISION   A( NMAX, NMAX ), WORK( 2*NMAX )
      int                PIV( NMAX )
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAESM, CHKXER, DPSTF2, DPSTRF
*     ..
*     .. Scalars in Common ..
      int                INFOT, NOUT
      LOGICAL            LERR, OK
      String             SRNAMT;
*     ..
*     .. Common blocks ..
      COMMON             / INFOC / INFOT, NOUT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DBLE
*     ..
*     .. Executable Statements ..
*
      NOUT = NUNIT
      WRITE( NOUT, FMT = * )
*
*     Set the variables to innocuous values.
*
      DO 110 J = 1, NMAX
         DO 100 I = 1, NMAX
            A( I, J ) = 1.D0 / DBLE( I+J )
*
  100    CONTINUE
         PIV( J ) = J
         WORK( J ) = 0.D0
         WORK( NMAX+J ) = 0.D0
*
  110 CONTINUE
      OK = .TRUE.
*
*
*        Test error exits of the routines that use the Cholesky
*        decomposition of a symmetric positive semidefinite matrix.
*
*        DPSTRF
*
      SRNAMT = 'DPSTRF'
      INFOT = 1
      CALL DPSTRF( '/', 0, A, 1, PIV, RANK, -1.D0, WORK, INFO )
      CALL CHKXER( 'DPSTRF', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL DPSTRF( 'U', -1, A, 1, PIV, RANK, -1.D0, WORK, INFO )
      CALL CHKXER( 'DPSTRF', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL DPSTRF( 'U', 2, A, 1, PIV, RANK, -1.D0, WORK, INFO )
      CALL CHKXER( 'DPSTRF', INFOT, NOUT, LERR, OK )
*
*        DPSTF2
*
      SRNAMT = 'DPSTF2'
      INFOT = 1
      CALL DPSTF2( '/', 0, A, 1, PIV, RANK, -1.D0, WORK, INFO )
      CALL CHKXER( 'DPSTF2', INFOT, NOUT, LERR, OK )
      INFOT = 2
      CALL DPSTF2( 'U', -1, A, 1, PIV, RANK, -1.D0, WORK, INFO )
      CALL CHKXER( 'DPSTF2', INFOT, NOUT, LERR, OK )
      INFOT = 4
      CALL DPSTF2( 'U', 2, A, 1, PIV, RANK, -1.D0, WORK, INFO )
      CALL CHKXER( 'DPSTF2', INFOT, NOUT, LERR, OK )
*
*
*     Print a summary line.
*
      CALL ALAESM( PATH, OK, NOUT )
*
      RETURN
*
*     End of DERRPS
*
      END
