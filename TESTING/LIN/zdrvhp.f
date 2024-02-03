      SUBROUTINE ZDRVHP( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, NMAX, A, AFAC, AINV, B, X, XACT, WORK, RWORK, IWORK, NOUT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               TSTERR;
      int                NMAX, NN, NOUT, NRHS;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                IWORK( * ), NVAL( * );
      double             RWORK( * );
      COMPLEX*16         A( * ), AFAC( * ), AINV( * ), B( * ), WORK( * ), X( * ), XACT( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      int                NTYPES, NTESTS;
      const              NTYPES = 10, NTESTS = 6 ;
      int                NFACT;
      const              NFACT = 2 ;
      // ..
      // .. Local Scalars ..
      bool               ZEROT;
      String             DIST, FACT, PACKIT, TYPE, UPLO, XTYPE;
      String             PATH;
      int                I, I1, I2, IFACT, IMAT, IN, INFO, IOFF, IUPLO, IZERO, J, K, K1, KL, KU, LDA, MODE, N, NB, NBMIN, NERRS, NFAIL, NIMAT, NPP, NRUN, NT;
      double             AINVNM, ANORM, CNDNUM, RCOND, RCONDC;
      // ..
      // .. Local Arrays ..
      String             FACTS( NFACT );
      int                ISEED( 4 ), ISEEDY( 4 );
      double             RESULT( NTESTS );
      // ..
      // .. External Functions ..
      double             DGET06, ZLANHP;
      // EXTERNAL DGET06, ZLANHP
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALADHD, ALAERH, ALASVM, XLAENV, ZCOPY, ZERRVX, ZGET04, ZHPSV, ZHPSVX, ZHPT01, ZHPTRF, ZHPTRI, ZLACPY, ZLAIPD, ZLARHS, ZLASET, ZLATB4, ZLATMS, ZPPT02, ZPPT05
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
      // .. Intrinsic Functions ..
      // INTRINSIC DCMPLX, MAX, MIN
      // ..
      // .. Data statements ..
      DATA               ISEEDY / 1988, 1989, 1990, 1991 /
      DATA               FACTS / 'F', 'N' /
      // ..
      // .. Executable Statements ..

      // Initialize constants and the random number seed.

      PATH( 1: 1 ) = 'Z'
      PATH( 2: 3 ) = 'HP'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      DO 10 I = 1, 4
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE

      // Test the error exits

      IF( TSTERR ) CALL ZERRVX( PATH, NOUT )
      INFOT = 0

      // Set the block size and minimum block size for testing.

      NB = 1
      NBMIN = 2
      CALL XLAENV( 1, NB )
      CALL XLAENV( 2, NBMIN )

      // Do for each value of N in NVAL

      DO 180 IN = 1, NN
         N = NVAL( IN )
         LDA = MAX( N, 1 )
         NPP = N*( N+1 ) / 2
         XTYPE = 'N'
         NIMAT = NTYPES
         IF( N.LE.0 ) NIMAT = 1

         DO 170 IMAT = 1, NIMAT

            // Do the tests only if DOTYPE( IMAT ) is true.

            IF( .NOT.DOTYPE( IMAT ) ) GO TO 170

            // Skip types 3, 4, 5, or 6 if the matrix size is too small.

            ZEROT = IMAT.GE.3 .AND. IMAT.LE.6
            IF( ZEROT .AND. N.LT.IMAT-2 ) GO TO 170

            // Do first for UPLO = 'U', then for UPLO = 'L'

            DO 160 IUPLO = 1, 2
               if ( IUPLO.EQ.1 ) {
                  UPLO = 'U'
                  PACKIT = 'C'
               } else {
                  UPLO = 'L'
                  PACKIT = 'R'
               }

               // Set up parameters with ZLATB4 and generate a test matrix
               // with ZLATMS.

               CALL ZLATB4( PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST )

               SRNAMT = 'ZLATMS'
               CALL ZLATMS( N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, PACKIT, A, LDA, WORK, INFO )

               // Check error code from ZLATMS.

               if ( INFO.NE.0 ) {
                  CALL ALAERH( PATH, 'ZLATMS', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )
                  GO TO 160
               }

               // For types 3-6, zero one or more rows and columns of the
               // matrix to test that INFO is returned correctly.

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

               // Set the imaginary part of the diagonals.

               if ( IUPLO.EQ.1 ) {
                  CALL ZLAIPD( N, A, 2, 1 )
               } else {
                  CALL ZLAIPD( N, A, N, -1 )
               }

               DO 150 IFACT = 1, NFACT

                  // Do first for FACT = 'F', then for other values.

                  FACT = FACTS( IFACT )

                  // Compute the condition number for comparison with
                 t // he value returned by ZHPSVX.

                  if ( ZEROT ) {
                     IF( IFACT.EQ.1 ) GO TO 150
                     RCONDC = ZERO

                  } else if ( IFACT.EQ.1 ) {

                     // Compute the 1-norm of A.

                     ANORM = ZLANHP( '1', UPLO, N, A, RWORK )

                     // Factor the matrix A.

                     CALL ZCOPY( NPP, A, 1, AFAC, 1 )
                     CALL ZHPTRF( UPLO, N, AFAC, IWORK, INFO )

                     // Compute inv(A) and take its norm.

                     CALL ZCOPY( NPP, AFAC, 1, AINV, 1 )
                     CALL ZHPTRI( UPLO, N, AINV, IWORK, WORK, INFO )
                     AINVNM = ZLANHP( '1', UPLO, N, AINV, RWORK )

                     // Compute the 1-norm condition number of A.

                     if ( ANORM.LE.ZERO .OR. AINVNM.LE.ZERO ) {
                        RCONDC = ONE
                     } else {
                        RCONDC = ( ONE / ANORM ) / AINVNM
                     }
                  }

                  // Form an exact solution and set the right hand side.

                  SRNAMT = 'ZLARHS'
                  CALL ZLARHS( PATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO )
                  XTYPE = 'C'

                  // --- Test ZHPSV  ---

                  if ( IFACT.EQ.2 ) {
                     CALL ZCOPY( NPP, A, 1, AFAC, 1 )
                     CALL ZLACPY( 'Full', N, NRHS, B, LDA, X, LDA )

                     // Factor the matrix and solve the system using ZHPSV.

                     SRNAMT = 'ZHPSV '
                     CALL ZHPSV( UPLO, N, NRHS, AFAC, IWORK, X, LDA, INFO )

                     // Adjust the expected value of INFO to account for
                     // pivoting.

                     K = IZERO
                     if ( K.GT.0 ) {
  100                   CONTINUE
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

                     // Check error code from ZHPSV .

                     if ( INFO.NE.K ) {
                        CALL ALAERH( PATH, 'ZHPSV ', INFO, K, UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT )
                        GO TO 120
                     } else if ( INFO.NE.0 ) {
                        GO TO 120
                     }

                     // Reconstruct matrix from factors and compute
                     // residual.

                     CALL ZHPT01( UPLO, N, A, AFAC, IWORK, AINV, LDA, RWORK, RESULT( 1 ) )

                     // Compute residual of the computed solution.

                     CALL ZLACPY( 'Full', N, NRHS, B, LDA, WORK, LDA )
                     CALL ZPPT02( UPLO, N, NRHS, A, X, LDA, WORK, LDA, RWORK, RESULT( 2 ) )

                     // Check solution from generated exact solution.

                     CALL ZGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) )
                     NT = 3

                     // Print information about the tests that did not pass
                    t // he threshold.

                     DO 110 K = 1, NT
                        if ( RESULT( K ).GE.THRESH ) {
                           IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, PATH )                            WRITE( NOUT, FMT = 9999 )'ZHPSV ', UPLO, N, IMAT, K, RESULT( K )
                           NFAIL = NFAIL + 1
                        }
  110                CONTINUE
                     NRUN = NRUN + NT
  120                CONTINUE
                  }

                  // --- Test ZHPSVX ---

                  IF( IFACT.EQ.2 .AND. NPP.GT.0 ) CALL ZLASET( 'Full', NPP, 1, DCMPLX( ZERO ), DCMPLX( ZERO ), AFAC, NPP )
                  CALL ZLASET( 'Full', N, NRHS, DCMPLX( ZERO ), DCMPLX( ZERO ), X, LDA )

                  // Solve the system and compute the condition number and
                  // error bounds using ZHPSVX.

                  SRNAMT = 'ZHPSVX'
                  CALL ZHPSVX( FACT, UPLO, N, NRHS, A, AFAC, IWORK, B, LDA, X, LDA, RCOND, RWORK, RWORK( NRHS+1 ), WORK, RWORK( 2*NRHS+1 ), INFO )

                  // Adjust the expected value of INFO to account for
                  // pivoting.

                  K = IZERO
                  if ( K.GT.0 ) {
  130                CONTINUE
                     if ( IWORK( K ).LT.0 ) {
                        if ( IWORK( K ).NE.-K ) {
                           K = -IWORK( K )
                           GO TO 130
                        }
                     } else if ( IWORK( K ).NE.K ) {
                        K = IWORK( K )
                        GO TO 130
                     }
                  }

                  // Check the error code from ZHPSVX.

                  if ( INFO.NE.K ) {
                     CALL ALAERH( PATH, 'ZHPSVX', INFO, K, FACT // UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT )
                     GO TO 150
                  }

                  if ( INFO.EQ.0 ) {
                     if ( IFACT.GE.2 ) {

                        // Reconstruct matrix from factors and compute
                        // residual.

                        CALL ZHPT01( UPLO, N, A, AFAC, IWORK, AINV, LDA, RWORK( 2*NRHS+1 ), RESULT( 1 ) )
                        K1 = 1
                     } else {
                        K1 = 2
                     }

                     // Compute residual of the computed solution.

                     CALL ZLACPY( 'Full', N, NRHS, B, LDA, WORK, LDA )
                     CALL ZPPT02( UPLO, N, NRHS, A, X, LDA, WORK, LDA, RWORK( 2*NRHS+1 ), RESULT( 2 ) )

                     // Check solution from generated exact solution.

                     CALL ZGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) )

                     // Check the error bounds from iterative refinement.

                     CALL ZPPT05( UPLO, N, NRHS, A, B, LDA, X, LDA, XACT, LDA, RWORK, RWORK( NRHS+1 ), RESULT( 4 ) )
                  } else {
                     K1 = 6
                  }

                  // Compare RCOND from ZHPSVX with the computed value
                  // in RCONDC.

                  RESULT( 6 ) = DGET06( RCOND, RCONDC )

                  // Print information about the tests that did not pass
                 t // he threshold.

                  DO 140 K = K1, 6
                     if ( RESULT( K ).GE.THRESH ) {
                        IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, PATH )                         WRITE( NOUT, FMT = 9998 )'ZHPSVX', FACT, UPLO, N, IMAT, K, RESULT( K )
                        NFAIL = NFAIL + 1
                     }
  140             CONTINUE
                  NRUN = NRUN + 7 - K1

  150          CONTINUE

  160       CONTINUE
  170    CONTINUE
  180 CONTINUE

      // Print a summary of the results.

      CALL ALASVM( PATH, NOUT, NFAIL, NRUN, NERRS )

 9999 FORMAT( 1X, A, ', UPLO=''', A1, ''', N =', I5, ', type ', I2,
     $      ', test ', I2, ', ratio =', G12.5 )
 9998 FORMAT( 1X, A, ', FACT=''', A1, ''', UPLO=''', A1, ''', N =', I5,
     $      ', type ', I2, ', test ', I2, ', ratio =', G12.5 )
      RETURN

      // End of ZDRVHP

      }
