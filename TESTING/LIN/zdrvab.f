      SUBROUTINE ZDRVAB( DOTYPE, NM, MVAL, NNS, NSVAL, THRESH, NMAX, A, AFAC, B, X, WORK, RWORK, SWORK, IWORK, NOUT )
*
*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*
*     .. Scalar Arguments ..
      int                NM, NMAX, NNS, NOUT
      DOUBLE PRECISION   THRESH
*     ..
*     .. Array Arguments ..
      LOGICAL            DOTYPE( * )
      int                MVAL( * ), NSVAL( * ), IWORK( * )
      DOUBLE PRECISION   RWORK( * )
      COMPLEX            SWORK( * )
      COMPLEX*16         A( * ), AFAC( * ), B( * ), WORK( * ), X( * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO
      PARAMETER          ( ZERO = 0.0D+0 )
      int                NTYPES
      PARAMETER          ( NTYPES = 11 )
      int                NTESTS
      PARAMETER          ( NTESTS = 1 )
*     ..
*     .. Local Scalars ..
      LOGICAL            ZEROT
      CHARACTER          DIST, TRANS, TYPE, XTYPE
      String             PATH;
      int                I, IM, IMAT, INFO, IOFF, IRHS, IZERO, KL, KU, LDA, M, MODE, N, NERRS, NFAIL, NIMAT, NRHS, NRUN
      DOUBLE PRECISION   ANORM, CNDNUM
*     ..
*     .. Local Arrays ..
      int                ISEED( 4 ), ISEEDY( 4 )
      DOUBLE PRECISION   RESULT( NTESTS )
*     ..
*     .. Local Variables ..
      int                ITER, KASE
*     ..
*     .. External Subroutines ..
      EXTERNAL           ALAERH, ALAHD, ZGET08, ZLACPY, ZLARHS, ZLASET, ZLATB4, ZLATMS
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          DCMPLX, DBLE, MAX, MIN, SQRT
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
*     .. Data statements ..
      DATA               ISEEDY / 2006, 2007, 2008, 2009 /
*     ..
*     .. Executable Statements ..
*
*     Initialize constants and the random number seed.
*
      KASE = 0
      PATH( 1: 1 ) = 'Zomplex precision'
      PATH( 2: 3 ) = 'GE'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      DO 10 I = 1, 4
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE
*
      INFOT = 0
*
*     Do for each value of M in MVAL
*
      DO 120 IM = 1, NM
         M = MVAL( IM )
         LDA = MAX( 1, M )
*
         N = M
         NIMAT = NTYPES
         IF( M.LE.0 .OR. N.LE.0 ) NIMAT = 1
*
         DO 100 IMAT = 1, NIMAT
*
*           Do the tests only if DOTYPE( IMAT ) is true.
*
            IF( .NOT.DOTYPE( IMAT ) ) GO TO 100
*
*           Skip types 5, 6, or 7 if the matrix size is too small.
*
            ZEROT = IMAT.GE.5 .AND. IMAT.LE.7
            IF( ZEROT .AND. N.LT.IMAT-4 ) GO TO 100
*
*           Set up parameters with ZLATB4 and generate a test matrix
*           with ZLATMS.
*
            CALL ZLATB4( PATH, IMAT, M, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST )
*
            SRNAMT = 'ZLATMS'
            CALL ZLATMS( M, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, 'No packing', A, LDA, WORK, INFO )
*
*           Check error code from ZLATMS.
*
            IF( INFO.NE.0 ) THEN
               CALL ALAERH( PATH, 'ZLATMS', INFO, 0, ' ', M, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )
               GO TO 100
            END IF
*
*           For types 5-7, zero one or more columns of the matrix to
*           test that INFO is returned correctly.
*
            IF( ZEROT ) THEN
               IF( IMAT.EQ.5 ) THEN
                  IZERO = 1
               ELSE IF( IMAT.EQ.6 ) THEN
                  IZERO = MIN( M, N )
               ELSE
                  IZERO = MIN( M, N ) / 2 + 1
               END IF
               IOFF = ( IZERO-1 )*LDA
               IF( IMAT.LT.7 ) THEN
                  DO 20 I = 1, M
                     A( IOFF+I ) = ZERO
   20             CONTINUE
               ELSE
                  CALL ZLASET( 'Full', M, N-IZERO+1, DCMPLX(ZERO), DCMPLX(ZERO), A( IOFF+1 ), LDA )
               END IF
            ELSE
               IZERO = 0
            END IF
*
            DO 60 IRHS = 1, NNS
               NRHS = NSVAL( IRHS )
               XTYPE = 'N'
               TRANS = 'N'
*
               SRNAMT = 'ZLARHS'
               CALL ZLARHS( PATH, XTYPE, ' ', TRANS, N, N, KL, KU, NRHS, A, LDA, X, LDA, B, LDA, ISEED, INFO )
*
               SRNAMT = 'ZCGESV'
*
               KASE = KASE + 1
*
               CALL ZLACPY( 'Full', M, N, A, LDA, AFAC, LDA )
*
               CALL ZCGESV( N, NRHS, A, LDA, IWORK, B, LDA, X, LDA, WORK, SWORK, RWORK, ITER, INFO)
*
               IF (ITER.LT.0) THEN
                   CALL ZLACPY( 'Full', M, N, AFAC, LDA, A, LDA )
               ENDIF
*
*              Check error code from ZCGESV. This should be the same as
*              the one of DGETRF.
*
               IF( INFO.NE.IZERO ) THEN
*
                  IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )
                  NERRS = NERRS + 1
*
                  IF( INFO.NE.IZERO .AND. IZERO.NE.0 ) THEN
                     WRITE( NOUT, FMT = 9988 )'ZCGESV',INFO, IZERO,M,IMAT
                  ELSE
                     WRITE( NOUT, FMT = 9975 )'ZCGESV',INFO, M, IMAT
                  END IF
               END IF
*
*              Skip the remaining test if the matrix is singular.
*
               IF( INFO.NE.0 ) GO TO 100
*
*              Check the quality of the solution
*
               CALL ZLACPY( 'Full', N, NRHS, B, LDA, WORK, LDA )
*
               CALL ZGET08( TRANS, N, N, NRHS, A, LDA, X, LDA, WORK, LDA, RWORK, RESULT( 1 ) )
*
*              Check if the test passes the testing.
*              Print information about the tests that did not
*              pass the testing.
*
*              If iterative refinement has been used and claimed to
*              be successful (ITER>0), we want
*                NORMI(B - A*X)/(NORMI(A)*NORMI(X)*EPS*SRQT(N)) < 1
*
*              If double precision has been used (ITER<0), we want
*                NORMI(B - A*X)/(NORMI(A)*NORMI(X)*EPS) < THRES
*              (Cf. the linear solver testing routines)
*
               IF ((THRESH.LE.0.0E+00) .OR.((ITER.GE.0).AND.(N.GT.0) .AND.(RESULT(1).GE.SQRT(DBLE(N)))) .OR.((ITER.LT.0).AND.(RESULT(1).GE.THRESH))) THEN
*
                  IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) THEN
                     WRITE( NOUT, FMT = 8999 )'DGE'
                     WRITE( NOUT, FMT = '( '' Matrix types:'' )' )
                     WRITE( NOUT, FMT = 8979 )
                     WRITE( NOUT, FMT = '( '' Test ratios:'' )' )
                     WRITE( NOUT, FMT = 8960 )1
                     WRITE( NOUT, FMT = '( '' Messages:'' )' )
                  END IF
*
                  WRITE( NOUT, FMT = 9998 )TRANS, N, NRHS, IMAT, 1, RESULT( 1 )
                  NFAIL = NFAIL + 1
               END IF
               NRUN = NRUN + 1
   60       CONTINUE
  100    CONTINUE
  120 CONTINUE
*
*     Print a summary of the results.
*
      IF( NFAIL.GT.0 ) THEN
         WRITE( NOUT, FMT = 9996 )'ZCGESV', NFAIL, NRUN
      ELSE
         WRITE( NOUT, FMT = 9995 )'ZCGESV', NRUN
      END IF
      IF( NERRS.GT.0 ) THEN
         WRITE( NOUT, FMT = 9994 )NERRS
      END IF
*
 9998 FORMAT( ' TRANS=''', A1, ''', N =', I5, ', NRHS=', I3, ', type ',
     $      I2, ', test(', I2, ') =', G12.5 )
 9996 FORMAT( 1X, A6, ': ', I6, ' out of ', I6,
     $      ' tests failed to pass the threshold' )
 9995 FORMAT( /1X, 'All tests for ', A6,
     $      ' routines passed the threshold ( ', I6, ' tests run)' )
 9994 FORMAT( 6X, I6, ' error messages recorded' )
*
*     SUBNAM, INFO, INFOE, M, IMAT
*
 9988 FORMAT( ' *** ', A6, ' returned with INFO =', I5, ' instead of ',
     $      I5, / ' ==> M =', I5, ', type ',
     $      I2 )
*
*     SUBNAM, INFO, M, IMAT
*
 9975 FORMAT( ' *** Error code from ', A6, '=', I5, ' for M=', I5,
     $      ', type ', I2 )
 8999 FORMAT( / 1X, A3, ':  General dense matrices' )
 8979 FORMAT( 4X, '1. Diagonal', 24X, '7. Last n/2 columns zero', / 4X,
     $      '2. Upper triangular', 16X,
     $      '8. Random, CNDNUM = sqrt(0.1/EPS)', / 4X,
     $      '3. Lower triangular', 16X, '9. Random, CNDNUM = 0.1/EPS',
     $      / 4X, '4. Random, CNDNUM = 2', 13X,
     $      '10. Scaled near underflow', / 4X, '5. First column zero',
     $      14X, '11. Scaled near overflow', / 4X,
     $      '6. Last column zero' )
 8960 FORMAT( 3X, I2, ': norm_1( B - A * X )  / ',
     $      '( norm_1(A) * norm_1(X) * EPS * SQRT(N) ) > 1 if ITERREF',
     $      / 4x, 'or norm_1( B - A * X )  / ',
     $      '( norm_1(A) * norm_1(X) * EPS ) > THRES if DGETRF' )
      RETURN
*
*     End of ZDRVAB
*
      END
