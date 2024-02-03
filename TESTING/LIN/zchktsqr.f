      SUBROUTINE ZCHKTSQR( THRESH, TSTERR, NM, MVAL, NN, NVAL, NNB, NBVAL, NOUT )
      IMPLICIT NONE
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      LOGICAL            TSTERR
      int                NM, NN, NNB, NOUT
      DOUBLE PRECISION   THRESH
*     ..
*     .. Array Arguments ..
      int                MVAL( * ), NBVAL( * ), NVAL( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      int                NTESTS
      PARAMETER          ( NTESTS = 6 )
*     ..
*     .. Local Scalars ..
      String             PATH;
      int                I, J, K, T, M, N, NB, NFAIL, NERRS, NRUN, INB, MINMN, MB, IMB
*
*     .. Local Arrays ..
      DOUBLE PRECISION   RESULT( NTESTS )
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAERH, ALAHD, ALASUM, ZERRTSQR, ZTSQR01, XLAENV
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC  MAX, MIN
*     ..
*     .. Scalars in Common ..
      LOGICAL            LERR, OK
      String             SRNAMT;
      int                INFOT, NUNIT
*     ..
*     .. Common blocks ..
      COMMON             / INFOC / INFOT, NUNIT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
*     ..
*     .. Executable Statements ..
*
*     Initialize constants
*
      PATH( 1: 1 ) = 'Z'
      PATH( 2: 3 ) = 'TS'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
*
*     Test the error exits
*
      CALL XLAENV( 1, 0 )
      CALL XLAENV( 2, 0 )
      IF( TSTERR ) CALL ZERRTSQR( PATH, NOUT )
      INFOT = 0
*
*     Do for each value of M in MVAL.
*
      DO I = 1, NM
         M = MVAL( I )
*
*        Do for each value of N in NVAL.
*
         DO J = 1, NN
            N = NVAL( J )
              IF (MIN(M,N).NE.0) THEN
              DO INB = 1, NNB
                MB = NBVAL( INB )
                  CALL XLAENV( 1, MB )
                  DO IMB = 1, NNB
                    NB = NBVAL( IMB )
                    CALL XLAENV( 2, NB )
*
*                 Test ZGEQR and ZGEMQR
*
                    CALL ZTSQR01( 'TS', M, N, MB, NB, RESULT )
*
*                 Print information about the tests that did not
*                 pass the threshold.
*
                    DO T = 1, NTESTS
                      IF( RESULT( T ).GE.THRESH ) THEN
                        IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                         WRITE( NOUT, FMT = 9999 )M, N, MB, NB, T, RESULT( T )
                        NFAIL = NFAIL + 1
                      END IF
                    END DO
                    NRUN = NRUN + NTESTS
                  END DO
              END DO
              END IF
         END DO
      END DO
*
*     Do for each value of M in MVAL.
*
      DO I = 1, NM
         M = MVAL( I )
*
*        Do for each value of N in NVAL.
*
         DO J = 1, NN
            N = NVAL( J )
              IF (MIN(M,N).NE.0) THEN
              DO INB = 1, NNB
                MB = NBVAL( INB )
                  CALL XLAENV( 1, MB )
                  DO IMB = 1, NNB
                    NB = NBVAL( IMB )
                    CALL XLAENV( 2, NB )
*
*                 Test ZGELQ and ZGEMLQ
*
                    CALL ZTSQR01( 'SW', M, N, MB, NB, RESULT )
*
*                 Print information about the tests that did not
*                 pass the threshold.
*
                    DO T = 1, NTESTS
                      IF( RESULT( T ).GE.THRESH ) THEN
                        IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                            WRITE( NOUT, FMT = 9998 )M, N, MB, NB, T, RESULT( T )
                        NFAIL = NFAIL + 1
                      END IF
                    END DO
                    NRUN = NRUN + NTESTS
                  END DO
              END DO
              END IF
         END DO
      END DO
*
*     Print a summary of the results.
*
      CALL ALASUM( PATH, NOUT, NFAIL, NRUN, NERRS )
*
 9999 FORMAT( 'TS: M=', I5, ', N=', I5, ', MB=', I5,
     $      ', NB=', I5,' test(', I2, ')=', G12.5 )
 9998 FORMAT( 'SW: M=', I5, ', N=', I5, ', MB=', I5,
     $      ', NB=', I5,' test(', I2, ')=', G12.5 )
      RETURN
*
*     End of ZCHKTSQR
*
      END
