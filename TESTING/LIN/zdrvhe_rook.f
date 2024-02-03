      SUBROUTINE ZDRVHE_ROOK( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, NMAX, A, AFAC, AINV, B, X, XACT, WORK, RWORK, IWORK, NOUT )

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
      const              NTYPES = 10, NTESTS = 3 ;
      int                NFACT;
      const              NFACT = 2 ;
      // ..
      // .. Local Scalars ..
      bool               ZEROT;
      String             DIST, FACT, TYPE, UPLO, XTYPE;
      String             MATPATH, PATH;
      int                I, I1, I2, IFACT, IMAT, IN, INFO, IOFF, IUPLO, IZERO, J, K, KL, KU, LDA, LWORK, MODE, N, NB, NBMIN, NERRS, NFAIL, NIMAT, NRUN, NT;
      double             AINVNM, ANORM, CNDNUM, RCONDC;
      // ..
      // .. Local Arrays ..
      String             FACTS( NFACT ), UPLOS( 2 );
      int                ISEED( 4 ), ISEEDY( 4 );
      double             RESULT( NTESTS );

      // ..
      // .. External Functions ..
      double             ZLANHE;
      // EXTERNAL ZLANHE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALADHD, ALAERH, ALASVM, XLAENV, ZERRVX, ZGET04, ZLACPY, ZLARHS, ZLATB4, ZLATMS, ZHESV_ROOK, ZHET01_ROOK, ZPOT02, ZHETRF_ROOK, ZHETRI_ROOK
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
      // INTRINSIC MAX, MIN
      // ..
      // .. Data statements ..
      DATA               ISEEDY / 1988, 1989, 1990, 1991 /
      DATA               UPLOS / 'U', 'L' / , FACTS / 'F', 'N' /
      // ..
      // .. Executable Statements ..

      // Initialize constants and the random number seed.

      // Test path

      PATH( 1: 1 ) = 'Zomplex precision'
      PATH( 2: 3 ) = 'HR'

      // Path to generate matrices

      MATPATH( 1: 1 ) = 'Zomplex precision'
      MATPATH( 2: 3 ) = 'HE'

      NRUN = 0
      NFAIL = 0
      NERRS = 0
      DO 10 I = 1, 4
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE
      LWORK = MAX( 2*NMAX, NMAX*NRHS )

      // Test the error exits

      IF( TSTERR ) CALL ZERRVX( PATH, NOUT )
      INFOT = 0

      // Set the block size and minimum block size for which the block
      // routine should be used, which will be later returned by ILAENV.

      NB = 1
      NBMIN = 2
      CALL XLAENV( 1, NB )
      CALL XLAENV( 2, NBMIN )

      // Do for each value of N in NVAL

      DO 180 IN = 1, NN
         N = NVAL( IN )
         LDA = MAX( N, 1 )
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
               UPLO = UPLOS( IUPLO )

                  // Begin generate the test matrix A.

                  // Set up parameters with ZLATB4 for the matrix generator
                  // based on the type of matrix to be generated.

                  CALL ZLATB4( MATPATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST )

                  // Generate a matrix with ZLATMS.

                  SRNAMT = 'ZLATMS'
                  CALL ZLATMS( N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, UPLO, A, LDA, WORK, INFO )

                  // Check error code from ZLATMS and handle error.

                  if ( INFO.NE.0 ) {
                     CALL ALAERH( PATH, 'ZLATMS', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT )
                     GO TO 160
                  }

                  // For types 3-6, zero one or more rows and columns of
                 t // he matrix to test that INFO is returned correctly.

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
                           IOFF = ( IZERO-1 )*LDA
                           DO 20 I = 1, IZERO - 1
                              A( IOFF+I ) = ZERO
   20                      CONTINUE
                           IOFF = IOFF + IZERO
                           DO 30 I = IZERO, N
                              A( IOFF ) = ZERO
                              IOFF = IOFF + LDA
   30                      CONTINUE
                        } else {
                           IOFF = IZERO
                           DO 40 I = 1, IZERO - 1
                              A( IOFF ) = ZERO
                              IOFF = IOFF + LDA
   40                      CONTINUE
                           IOFF = IOFF - IZERO
                           DO 50 I = IZERO, N
                              A( IOFF+I ) = ZERO
   50                      CONTINUE
                        }
                     } else {
                        if ( IUPLO.EQ.1 ) {

                        // Set the first IZERO rows and columns to zero.

                           IOFF = 0
                           DO 70 J = 1, N
                              I2 = MIN( J, IZERO )
                              DO 60 I = 1, I2
                                 A( IOFF+I ) = ZERO
   60                         CONTINUE
                              IOFF = IOFF + LDA
   70                      CONTINUE
                        } else {

                        // Set the first IZERO rows and columns to zero.

                           IOFF = 0
                           DO 90 J = 1, N
                              I1 = MAX( J, IZERO )
                              DO 80 I = I1, N
                                 A( IOFF+I ) = ZERO
   80                         CONTINUE
                              IOFF = IOFF + LDA
   90                      CONTINUE
                        }
                     }
                  } else {
                     IZERO = 0
                  }

                  // End generate the test matrix A.


               DO 150 IFACT = 1, NFACT

                  // Do first for FACT = 'F', then for other values.

                  FACT = FACTS( IFACT )

                  // Compute the condition number for comparison with
                 t // he value returned by ZHESVX_ROOK.

                  if ( ZEROT ) {
                     IF( IFACT.EQ.1 ) GO TO 150
                     RCONDC = ZERO

                  } else if ( IFACT.EQ.1 ) {

                     // Compute the 1-norm of A.

                     ANORM = ZLANHE( '1', UPLO, N, A, LDA, RWORK )

                     // Factor the matrix A.


                     CALL ZLACPY( UPLO, N, N, A, LDA, AFAC, LDA )
                     CALL ZHETRF_ROOK( UPLO, N, AFAC, LDA, IWORK, WORK, LWORK, INFO )

                     // Compute inv(A) and take its norm.

                     CALL ZLACPY( UPLO, N, N, AFAC, LDA, AINV, LDA )
                     LWORK = (N+NB+1)*(NB+3)
                     CALL ZHETRI_ROOK( UPLO, N, AINV, LDA, IWORK, WORK, INFO )
                     AINVNM = ZLANHE( '1', UPLO, N, AINV, LDA, RWORK )

                     // Compute the 1-norm condition number of A.

                     if ( ANORM.LE.ZERO .OR. AINVNM.LE.ZERO ) {
                        RCONDC = ONE
                     } else {
                        RCONDC = ( ONE / ANORM ) / AINVNM
                     }
                  }

                  // Form an exact solution and set the right hand side.

                  SRNAMT = 'ZLARHS'
                  CALL ZLARHS( MATPATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO )
                  XTYPE = 'C'

                  // --- Test ZHESV_ROOK  ---

                  if ( IFACT.EQ.2 ) {
                     CALL ZLACPY( UPLO, N, N, A, LDA, AFAC, LDA )
                     CALL ZLACPY( 'Full', N, NRHS, B, LDA, X, LDA )

                     // Factor the matrix and solve the system using
                     // ZHESV_ROOK.

                     SRNAMT = 'ZHESV_ROOK'
                     CALL ZHESV_ROOK( UPLO, N, NRHS, AFAC, LDA, IWORK, X, LDA, WORK, LWORK, INFO )

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

                     // Check error code from ZHESV_ROOK and handle error.

                     if ( INFO.NE.K ) {
                        CALL ALAERH( PATH, 'ZHESV_ROOK', INFO, K, UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT )
                        GO TO 120
                     } else if ( INFO.NE.0 ) {
                        GO TO 120
                     }

*+    TEST 1      Reconstruct matrix from factors and compute
                  // residual.

                     CALL ZHET01_ROOK( UPLO, N, A, LDA, AFAC, LDA, IWORK, AINV, LDA, RWORK, RESULT( 1 ) )

*+    TEST 2      Compute residual of the computed solution.

                     CALL ZLACPY( 'Full', N, NRHS, B, LDA, WORK, LDA )
                     CALL ZPOT02( UPLO, N, NRHS, A, LDA, X, LDA, WORK, LDA, RWORK, RESULT( 2 ) )

*+    TEST 3
                  // Check solution from generated exact solution.

                     CALL ZGET04( N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) )
                     NT = 3

                     // Print information about the tests that did not pass
                    t // he threshold.

                     DO 110 K = 1, NT
                        if ( RESULT( K ).GE.THRESH ) {
                           IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, PATH )                            WRITE( NOUT, FMT = 9999 )'ZHESV_ROOK', UPLO, N, IMAT, K, RESULT( K )
                           NFAIL = NFAIL + 1
                        }
  110                CONTINUE
                     NRUN = NRUN + NT
  120                CONTINUE
                  }

  150          CONTINUE

  160       CONTINUE
  170    CONTINUE
  180 CONTINUE

      // Print a summary of the results.

      CALL ALASVM( PATH, NOUT, NFAIL, NRUN, NERRS )

 9999 FORMAT( 1X, A, ', UPLO=''', A1, ''', N =', I5, ', type ', I2,
     $      ', test ', I2, ', ratio =', G12.5 )
      RETURN

      // End of ZDRVHE_ROOK

      }
