      SUBROUTINE DCHKSP( DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, NMAX, A, AFAC, AINV, B, X, XACT, WORK, RWORK, IWORK, NOUT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               TSTERR;
      int                NMAX, NN, NNS, NOUT;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                IWORK( * ), NSVAL( * ), NVAL( * );
      double             A( * ), AFAC( * ), AINV( * ), B( * ), RWORK( * ), WORK( * ), X( * ), XACT( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0D+0 ;
      int                NTYPES;
      const              NTYPES = 10 ;
      int                NTESTS;
      const              NTESTS = 8 ;
      // ..
      // .. Local Scalars ..
      bool               TRFCON, ZEROT;
      String             DIST, PACKIT, TYPE, UPLO, XTYPE;
      String             PATH;
      int                I, I1, I2, IMAT, IN, INFO, IOFF, IRHS, IUPLO, IZERO, J, K, KL, KU, LDA, MODE, N, NERRS, NFAIL, NIMAT, NPP, NRHS, NRUN, NT;
      double             ANORM, CNDNUM, RCOND, RCONDC;
      // ..
      // .. Local Arrays ..
      String             UPLOS( 2 );
      int                ISEED( 4 ), ISEEDY( 4 );
      double             RESULT( NTESTS );
      // ..
      // .. External Functions ..
      bool               LSAME;
      double             DGET06, DLANSP;
      // EXTERNAL LSAME, DGET06, DLANSP
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, DCOPY, DERRSY, DGET04, DLACPY, DLARHS, DLATB4, DLATMS, DPPT02, DPPT03, DPPT05, DSPCON, DSPRFS, DSPT01, DSPTRF, DSPTRI, DSPTRS
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
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

      PATH( 1: 1 ) = 'double          ';
      PATH( 2: 3 ) = 'SP'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      DO 10 I = 1, 4
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE

      // Test the error exits

      IF( TSTERR ) CALL DERRSY( PATH, NOUT )
      INFOT = 0

      // Do for each value of N in NVAL

      DO 170 IN = 1, NN
         N = NVAL( IN )
         LDA = MAX( N, 1 )
         XTYPE = 'N'
         NIMAT = NTYPES
         IF( N.LE.0 ) NIMAT = 1

         IZERO = 0
         DO 160 IMAT = 1, NIMAT

            // Do the tests only if DOTYPE( IMAT ) is true.

            IF( .NOT.DOTYPE( IMAT ) ) GO TO 160

            // Skip types 3, 4, 5, or 6 if the matrix size is too small.

            ZEROT = IMAT.GE.3 .AND. IMAT.LE.6
            IF( ZEROT .AND. N.LT.IMAT-2 ) GO TO 160

            // Do first for UPLO = 'U', then for UPLO = 'L'

            DO 150 IUPLO = 1, 2
               UPLO = UPLOS( IUPLO )
               if ( LSAME( UPLO, 'U' ) ) {
                  PACKIT = 'C'
               } else {
                  PACKIT = 'R'
               }

               // Set up parameters with DLATB4 and generate a test matrix
               // with DLATMS.

               dlatb4(PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );

               SRNAMT = 'DLATMS'
               dlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, PACKIT, A, LDA, WORK, INFO );

               // Check error code from DLATMS.

               if ( INFO.NE.0 ) {
                  alaerh(PATH, 'DLATMS', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );
                  GO TO 150
               }

               // For types 3-6, zero one or more rows and columns of
               // the matrix to test that INFO is returned correctly.

               if ( ZEROT ) {
                  if ( IMAT.EQ.3 ) {
                     IZERO = 1
                  } else if ( IMAT.EQ.4 ) {
                     IZERO = N
                  } else {
                     IZERO = N / 2 + 1
                  }

                  if ( IMAT.LT.6 ) {

                     // Set row and column IZERO to zero.

                     if ( IUPLO.EQ.1 ) {
                        IOFF = ( IZERO-1 )*IZERO / 2
                        DO 20 I = 1, IZERO - 1
                           A( IOFF+I ) = ZERO
   20                   CONTINUE
                        IOFF = IOFF + IZERO
                        DO 30 I = IZERO, N
                           A( IOFF ) = ZERO
                           IOFF = IOFF + I
   30                   CONTINUE
                     } else {
                        IOFF = IZERO
                        DO 40 I = 1, IZERO - 1
                           A( IOFF ) = ZERO
                           IOFF = IOFF + N - I
   40                   CONTINUE
                        IOFF = IOFF - IZERO
                        DO 50 I = IZERO, N
                           A( IOFF+I ) = ZERO
   50                   CONTINUE
                     }
                  } else {
                     IOFF = 0
                     if ( IUPLO.EQ.1 ) {

                        // Set the first IZERO rows and columns to zero.

                        DO 70 J = 1, N
                           I2 = MIN( J, IZERO )
                           DO 60 I = 1, I2
                              A( IOFF+I ) = ZERO
   60                      CONTINUE
                           IOFF = IOFF + J
   70                   CONTINUE
                     } else {

                        // Set the last IZERO rows and columns to zero.

                        DO 90 J = 1, N
                           I1 = MAX( J, IZERO )
                           DO 80 I = I1, N
                              A( IOFF+I ) = ZERO
   80                      CONTINUE
                           IOFF = IOFF + N - J
   90                   CONTINUE
                     }
                  }
               } else {
                  IZERO = 0
               }

               // Compute the L*D*L' or U*D*U' factorization of the matrix.

               NPP = N*( N+1 ) / 2
               dcopy(NPP, A, 1, AFAC, 1 );
               SRNAMT = 'DSPTRF'
               dsptrf(UPLO, N, AFAC, IWORK, INFO );

               // Adjust the expected value of INFO to account for
               // pivoting.

               K = IZERO
               if ( K.GT.0 ) {
  100             CONTINUE
                  if ( IWORK( K ).LT.0 ) {
                     if ( IWORK( K ).NE.-K ) {
                        K = -IWORK( K )
                        GO TO 100
                     }
                  } else if ( IWORK( K ).NE.K ) {
                     K = IWORK( K )
                     GO TO 100
                  }
               }

               // Check error code from DSPTRF.

               IF( INFO.NE.K ) CALL ALAERH( PATH, 'DSPTRF', INFO, K, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )
               if ( INFO.NE.0 ) {
                  TRFCON = .TRUE.
               } else {
                  TRFCON = .FALSE.
               }

*+    TEST 1
               // Reconstruct matrix from factors and compute residual.

               dspt01(UPLO, N, A, AFAC, IWORK, AINV, LDA, RWORK, RESULT( 1 ) );
               NT = 1

*+    TEST 2
               // Form the inverse and compute the residual.

               if ( .NOT.TRFCON ) {
                  dcopy(NPP, AFAC, 1, AINV, 1 );
                  SRNAMT = 'DSPTRI'
                  dsptri(UPLO, N, AINV, IWORK, WORK, INFO );

               // Check error code from DSPTRI.

                  IF( INFO.NE.0 ) CALL ALAERH( PATH, 'DSPTRI', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )

                  dppt03(UPLO, N, A, AINV, WORK, LDA, RWORK, RCONDC, RESULT( 2 ) );
                  NT = 2
               }

               // Print information about the tests that did not pass
               // the threshold.

               DO 110 K = 1, NT
                  if ( RESULT( K ).GE.THRESH ) {
                     IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                      WRITE( NOUT, FMT = 9999 )UPLO, N, IMAT, K, RESULT( K )
                     NFAIL = NFAIL + 1
                  }
  110          CONTINUE
               NRUN = NRUN + NT

               // Do only the condition estimate if INFO is not 0.

               if ( TRFCON ) {
                  RCONDC = ZERO
                  GO TO 140
               }

               DO 130 IRHS = 1, NNS
                  NRHS = NSVAL( IRHS )

*+    TEST 3
               // Solve and compute residual for  A * X = B.

                  SRNAMT = 'DLARHS'
                  dlarhs(PATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO );
                  dlacpy('Full', N, NRHS, B, LDA, X, LDA );

                  SRNAMT = 'DSPTRS'
                  dsptrs(UPLO, N, NRHS, AFAC, IWORK, X, LDA, INFO );

               // Check error code from DSPTRS.

                  IF( INFO.NE.0 ) CALL ALAERH( PATH, 'DSPTRS', INFO, 0, UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT )

                  dlacpy('Full', N, NRHS, B, LDA, WORK, LDA );
                  dppt02(UPLO, N, NRHS, A, X, LDA, WORK, LDA, RWORK, RESULT( 3 ) );

*+    TEST 4
               // Check solution from generated exact solution.

                  dget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 4 ) );

*+    TESTS 5, 6, and 7
               // Use iterative refinement to improve the solution.

                  SRNAMT = 'DSPRFS'
                  dsprfs(UPLO, N, NRHS, A, AFAC, IWORK, B, LDA, X, LDA, RWORK, RWORK( NRHS+1 ), WORK, IWORK( N+1 ), INFO );

               // Check error code from DSPRFS.

                  IF( INFO.NE.0 ) CALL ALAERH( PATH, 'DSPRFS', INFO, 0, UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT )

                  dget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 5 ) )                   CALL DPPT05( UPLO, N, NRHS, A, B, LDA, X, LDA, XACT, LDA, RWORK, RWORK( NRHS+1 ), RESULT( 6 ) );

                  // Print information about the tests that did not pass
                  // the threshold.

                  DO 120 K = 3, 7
                     if ( RESULT( K ).GE.THRESH ) {
                        IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                         WRITE( NOUT, FMT = 9998 )UPLO, N, NRHS, IMAT, K, RESULT( K )
                        NFAIL = NFAIL + 1
                     }
  120             CONTINUE
                  NRUN = NRUN + 5
  130          CONTINUE

*+    TEST 8
               // Get an estimate of RCOND = 1/CNDNUM.

  140          CONTINUE
               ANORM = DLANSP( '1', UPLO, N, A, RWORK )
               SRNAMT = 'DSPCON'
               dspcon(UPLO, N, AFAC, IWORK, ANORM, RCOND, WORK, IWORK( N+1 ), INFO );

               // Check error code from DSPCON.

               IF( INFO.NE.0 ) CALL ALAERH( PATH, 'DSPCON', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )

               RESULT( 8 ) = DGET06( RCOND, RCONDC )

               // Print the test ratio if it is .GE. THRESH.

               if ( RESULT( 8 ).GE.THRESH ) {
                  IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                   WRITE( NOUT, FMT = 9999 )UPLO, N, IMAT, 8, RESULT( 8 )
                  NFAIL = NFAIL + 1
               }
               NRUN = NRUN + 1
  150       CONTINUE
  160    CONTINUE
  170 CONTINUE

      // Print a summary of the results.

      alasum(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( ' UPLO = ''', A1, ''', N =', I5, ', type ', I2, ', test ', I2, ', ratio =', G12.5 )
 9998 FORMAT( ' UPLO = ''', A1, ''', N =', I5, ', NRHS=', I3, ', type ', I2, ', test(', I2, ') =', G12.5 )
      RETURN

      // End of DCHKSP

      }
