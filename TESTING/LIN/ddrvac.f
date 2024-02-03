      SUBROUTINE DDRVAC( DOTYPE, NM, MVAL, NNS, NSVAL, THRESH, NMAX, A, AFAC, B, X, WORK, RWORK, SWORK, NOUT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                NMAX, NM, NNS, NOUT;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                MVAL( * ), NSVAL( * );
      REAL               SWORK(*)
      double             A( * ), AFAC( * ), B( * ), RWORK( * ), WORK( * ), X( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0D+0 ;
      int                NTYPES;
      const              NTYPES = 9 ;
      int                NTESTS;
      const              NTESTS = 1 ;
      // ..
      // .. Local Scalars ..
      bool               ZEROT;
      String             DIST, TYPE, UPLO, XTYPE;
      String             PATH;
      int                I, IM, IMAT, INFO, IOFF, IRHS, IUPLO, IZERO, KL, KU, LDA, MODE, N, NERRS, NFAIL, NIMAT, NRHS, NRUN;
      double             ANORM, CNDNUM;
      // ..
      // .. Local Arrays ..
      String             UPLOS( 2 );
      int                ISEED( 4 ), ISEEDY( 4 );
      double             RESULT( NTESTS );
      // ..
      // .. Local Variables ..
      int                ITER, KASE;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, DLACPY, DLARHS, DLASET, DLATB4, DLATMS, DPOT06, DSPOSV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, SQRT
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String             SRNAMT;
      int                INFOT, NUNIT;
      // ..
      // .. Common blocks ..
      COMMON             / INFOC / INFOT, NUNIT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
      // ..
      // .. Data statements ..
      DATA               ISEEDY / 1988, 1989, 1990, 1991 /
      DATA               UPLOS / 'U', 'L' /
      // ..
      // .. Executable Statements ..

      // Initialize constants and the random number seed.

      KASE = 0
      PATH( 1: 1 ) = 'double          ';
      PATH( 2: 3 ) = 'PO'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      DO 10 I = 1, 4
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE

      INFOT = 0

      // Do for each value of N in MVAL

      DO 120 IM = 1, NM
         N = MVAL( IM )
         LDA = MAX( N, 1 )
         NIMAT = NTYPES
         IF( N.LE.0 ) NIMAT = 1

         DO 110 IMAT = 1, NIMAT

            // Do the tests only if DOTYPE( IMAT ) is true.

            IF( .NOT.DOTYPE( IMAT ) ) GO TO 110

            // Skip types 3, 4, or 5 if the matrix size is too small.

            ZEROT = IMAT.GE.3 .AND. IMAT.LE.5
            IF( ZEROT .AND. N.LT.IMAT-2 ) GO TO 110

            // Do first for UPLO = 'U', then for UPLO = 'L'

            DO 100 IUPLO = 1, 2
               UPLO = UPLOS( IUPLO )

               // Set up parameters with DLATB4 and generate a test matrix
               // with DLATMS.

               CALL DLATB4( PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST )

               SRNAMT = 'DLATMS'
               CALL DLATMS( N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, UPLO, A, LDA, WORK, INFO )

               // Check error code from DLATMS.

               if ( INFO.NE.0 ) {
                  CALL ALAERH( PATH, 'DLATMS', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )
                  GO TO 100
               }

               // For types 3-5, zero one row and column of the matrix to
              t // est that INFO is returned correctly.

               if ( ZEROT ) {
                  if ( IMAT.EQ.3 ) {
                     IZERO = 1
                  } else if ( IMAT.EQ.4 ) {
                     IZERO = N
                  } else {
                     IZERO = N / 2 + 1
                  }
                  IOFF = ( IZERO-1 )*LDA

                  // Set row and column IZERO of A to 0.

                  if ( IUPLO.EQ.1 ) {
                     DO 20 I = 1, IZERO - 1
                        A( IOFF+I ) = ZERO
   20                CONTINUE
                     IOFF = IOFF + IZERO
                     DO 30 I = IZERO, N
                        A( IOFF ) = ZERO
                        IOFF = IOFF + LDA
   30                CONTINUE
                  } else {
                     IOFF = IZERO
                     DO 40 I = 1, IZERO - 1
                        A( IOFF ) = ZERO
                        IOFF = IOFF + LDA
   40                CONTINUE
                     IOFF = IOFF - IZERO
                     DO 50 I = IZERO, N
                        A( IOFF+I ) = ZERO
   50                CONTINUE
                  }
               } else {
                  IZERO = 0
               }

               DO 60 IRHS = 1, NNS
                  NRHS = NSVAL( IRHS )
                  XTYPE = 'N'

                  // Form an exact solution and set the right hand side.

                  SRNAMT = 'DLARHS'
                  CALL DLARHS( PATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A, LDA, X, LDA, B, LDA, ISEED, INFO )

                  // Compute the L*L' or U'*U factorization of the
                  // matrix and solve the system.

                  SRNAMT = 'DSPOSV '
                  KASE = KASE + 1

                  CALL DLACPY( 'All', N, N, A, LDA, AFAC, LDA)

                  CALL DSPOSV( UPLO, N, NRHS, AFAC, LDA, B, LDA, X, LDA, WORK, SWORK, ITER, INFO )

                  if (ITER.LT.0) {
                     CALL DLACPY( 'All', N, N, A, LDA, AFAC, LDA )
                  ENDIF

                  // Check error code from DSPOSV .

                  if ( INFO.NE.IZERO ) {

                     IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )
                     NERRS = NERRS + 1

                     if ( INFO.NE.IZERO .AND. IZERO.NE.0 ) {
                        WRITE( NOUT, FMT = 9988 )'DSPOSV',INFO,IZERO,N, IMAT
                     } else {
                        WRITE( NOUT, FMT = 9975 )'DSPOSV',INFO,N,IMAT
                     }
                  }

                  // Skip the remaining test if the matrix is singular.

                  IF( INFO.NE.0 ) GO TO 110

                  // Check the quality of the solution

                  CALL DLACPY( 'All', N, NRHS, B, LDA, WORK, LDA )

                  CALL DPOT06( UPLO, N, NRHS, A, LDA, X, LDA, WORK, LDA, RWORK, RESULT( 1 ) )

                  // Check if the test passes the testing.
                  // Print information about the tests that did not
                  // pass the testing.

                  // If iterative refinement has been used and claimed to
                  // be successful (ITER>0), we want
                  // NORM1(B - A*X)/(NORM1(A)*NORM1(X)*EPS*SRQT(N)) < 1

                  // If double precision has been used (ITER<0), we want
                  // NORM1(B - A*X)/(NORM1(A)*NORM1(X)*EPS) < THRES
                  // (Cf. the linear solver testing routines)

                  if ((THRESH.LE.0.0E+00) .OR.((ITER.GE.0).AND.(N.GT.0) .AND.(RESULT(1).GE.SQRT(DBLE(N)))) .OR.((ITER.LT.0).AND.(RESULT(1).GE.THRESH))) {

                     if ( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) {
                        WRITE( NOUT, FMT = 8999 )'DPO'
                        WRITE( NOUT, FMT = '( '' Matrix types:'' )' )
                        WRITE( NOUT, FMT = 8979 )
                        WRITE( NOUT, FMT = '( '' Test ratios:'' )' )
                        WRITE( NOUT, FMT = 8960 )1
                        WRITE( NOUT, FMT = '( '' Messages:'' )' )
                     }

                     WRITE( NOUT, FMT = 9998 )UPLO, N, NRHS, IMAT, 1, RESULT( 1 )

                     NFAIL = NFAIL + 1

                  }

                  NRUN = NRUN + 1

   60          CONTINUE
  100       CONTINUE
  110    CONTINUE
  120 CONTINUE

      // Print a summary of the results.

      if ( NFAIL.GT.0 ) {
         WRITE( NOUT, FMT = 9996 )'DSPOSV', NFAIL, NRUN
      } else {
         WRITE( NOUT, FMT = 9995 )'DSPOSV', NRUN
      }
      if ( NERRS.GT.0 ) {
         WRITE( NOUT, FMT = 9994 )NERRS
      }

 9998 FORMAT( ' UPLO=''', A1, ''', N =', I5, ', NRHS=', I3, ', type ',
     $      I2, ', test(', I2, ') =', G12.5 )
 9996 FORMAT( 1X, A6, ': ', I6, ' out of ', I6,
     $      ' tests failed to pass the threshold' )
 9995 FORMAT( /1X, 'All tests for ', A6,
     $      ' routines passed the threshold ( ', I6, ' tests run)' )
 9994 FORMAT( 6X, I6, ' error messages recorded' )

      // SUBNAM, INFO, INFOE, N, IMAT

 9988 FORMAT( ' *** ', A6, ' returned with INFO =', I5, ' instead of ',
     $      I5, / ' ==> N =', I5, ', type ',
     $      I2 )

      // SUBNAM, INFO, N, IMAT

 9975 FORMAT( ' *** Error code from ', A6, '=', I5, ' for M=', I5,
     $      ', type ', I2 )
 8999 FORMAT( / 1X, A3, ':  positive definite dense matrices' )
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
     $      '( norm_1(A) * norm_1(X) * EPS ) > THRES if DPOTRF' )

      RETURN

      // End of DDRVAC

      }
