      SUBROUTINE DDRVSY_RK( DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, NMAX, A, AFAC, E, AINV, B, X, XACT, WORK, RWORK, IWORK, NOUT )

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
      double             A( * ), AFAC( * ), AINV( * ), B( * ), E( * ), RWORK( * ), WORK( * ), X( * ), XACT( * );
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
      String             PATH, MATPATH;
      int                I, I1, I2, IFACT, IMAT, IN, INFO, IOFF, IUPLO, IZERO, J, K, KL, KU, LDA, LWORK, MODE, N, NB, NBMIN, NERRS, NFAIL, NIMAT, NRUN, NT;
      double             AINVNM, ANORM, CNDNUM, RCONDC;
      // ..
      // .. Local Arrays ..
      String             FACTS( NFACT ), UPLOS( 2 );
      int                ISEED( 4 ), ISEEDY( 4 );
      double             RESULT( NTESTS );
      // ..
      // .. External Functions ..
      double             DLANSY;
      // EXTERNAL DLANSY
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALADHD, ALAERH, ALASVM, DERRVX, DGET04, DLACPY, DLARHS, DLATB4, DLATMS, DPOT02, DSYSV_RK, DSYT01_3, DSYTRF_RK, DSYTRI_3, XLAENV
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

      PATH( 1: 1 ) = 'double          ';
      PATH( 2: 3 ) = 'SK'

      // Path to generate matrices

      MATPATH( 1: 1 ) = 'double          ';
      MATPATH( 2: 3 ) = 'SY'

      NRUN = 0
      NFAIL = 0
      NERRS = 0
      DO 10 I = 1, 4
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE
      LWORK = MAX( 2*NMAX, NMAX*NRHS )

      // Test the error exits

      IF( TSTERR ) CALL DERRVX( PATH, NOUT )
      INFOT = 0

      // Set the block size and minimum block size for which the block
      // routine should be used, which will be later returned by ILAENV.

      NB = 1
      NBMIN = 2
      xlaenv(1, NB );
      xlaenv(2, NBMIN );

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

               // Set up parameters with DLATB4 for the matrix generator
               // based on the type of matrix to be generated.

               dlatb4(MATPATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );

               // Generate a matrix with DLATMS.

               SRNAMT = 'DLATMS'
               dlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, UPLO, A, LDA, WORK, INFO );

               // Check error code from DLATMS and handle error.

               if ( INFO.NE.0 ) {
                  alaerh(PATH, 'DLATMS', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

                  // Skip all tests for this generated matrix

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
                        IOFF = ( IZERO-1 )*LDA
                        DO 20 I = 1, IZERO - 1
                           A( IOFF+I ) = ZERO
   20                   CONTINUE
                        IOFF = IOFF + IZERO
                        DO 30 I = IZERO, N
                           A( IOFF ) = ZERO
                           IOFF = IOFF + LDA
   30                   CONTINUE
                     } else {
                        IOFF = IZERO
                        DO 40 I = 1, IZERO - 1
                           A( IOFF ) = ZERO
                           IOFF = IOFF + LDA
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
                           IOFF = IOFF + LDA
   70                   CONTINUE
                     } else {

                        // Set the last IZERO rows and columns to zero.

                        DO 90 J = 1, N
                           I1 = MAX( J, IZERO )
                           DO 80 I = I1, N
                              A( IOFF+I ) = ZERO
   80                      CONTINUE
                           IOFF = IOFF + LDA
   90                   CONTINUE
                     }
                  }
               } else {
                  IZERO = 0
               }

               // End generate the test matrix A.

               DO 150 IFACT = 1, NFACT

                  // Do first for FACT = 'F', then for other values.

                  FACT = FACTS( IFACT )

                  // Compute the condition number

                  if ( ZEROT ) {
                     IF( IFACT.EQ.1 ) GO TO 150
                     RCONDC = ZERO

                  } else if ( IFACT.EQ.1 ) {

                     // Compute the 1-norm of A.

                     ANORM = DLANSY( '1', UPLO, N, A, LDA, RWORK )

                     // Factor the matrix A.

                     dlacpy(UPLO, N, N, A, LDA, AFAC, LDA );
                     dsytrf_rk(UPLO, N, AFAC, LDA, E, IWORK, WORK, LWORK, INFO );

                     // Compute inv(A) and take its norm.

                     dlacpy(UPLO, N, N, AFAC, LDA, AINV, LDA );
                     LWORK = (N+NB+1)*(NB+3)

                     // We need to compute the inverse to compute
                     // RCONDC that is used later in TEST3.

                     dsytri_3(UPLO, N, AINV, LDA, E, IWORK, WORK, LWORK, INFO );
                     AINVNM = DLANSY( '1', UPLO, N, AINV, LDA, RWORK )

                     // Compute the 1-norm condition number of A.

                     if ( ANORM.LE.ZERO .OR. AINVNM.LE.ZERO ) {
                        RCONDC = ONE
                     } else {
                        RCONDC = ( ONE / ANORM ) / AINVNM
                     }
                  }

                  // Form an exact solution and set the right hand side.

                  SRNAMT = 'DLARHS'
                  dlarhs(MATPATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO );
                  XTYPE = 'C'

                  // --- Test DSYSV_RK  ---

                  if ( IFACT.EQ.2 ) {
                     dlacpy(UPLO, N, N, A, LDA, AFAC, LDA );
                     dlacpy('Full', N, NRHS, B, LDA, X, LDA );

                     // Factor the matrix and solve the system using
                     // DSYSV_RK.

                     SRNAMT = 'DSYSV_RK'
                     dsysv_rk(UPLO, N, NRHS, AFAC, LDA, E, IWORK, X, LDA, WORK, LWORK, INFO );

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

                     // Check error code from DSYSV_RK and handle error.

                     if ( INFO.NE.K ) {
                        alaerh(PATH, 'DSYSV_RK', INFO, K, UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );
                        GO TO 120
                     } else if ( INFO.NE.0 ) {
                        GO TO 120
                     }

*+    TEST 1      Reconstruct matrix from factors and compute
                  // residual.

                     dsyt01_3(UPLO, N, A, LDA, AFAC, LDA, E, IWORK, AINV, LDA, RWORK, RESULT( 1 ) );

*+    TEST 2      Compute residual of the computed solution.

                     dlacpy('Full', N, NRHS, B, LDA, WORK, LDA );
                     dpot02(UPLO, N, NRHS, A, LDA, X, LDA, WORK, LDA, RWORK, RESULT( 2 ) );

*+    TEST 3
                  // Check solution from generated exact solution.

                     dget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) );
                     NT = 3

                     // Print information about the tests that did not pass
                     // the threshold.

                     DO 110 K = 1, NT
                        if ( RESULT( K ).GE.THRESH ) {
                           IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALADHD( NOUT, PATH )                            WRITE( NOUT, FMT = 9999 )'DSYSV_RK', UPLO, N, IMAT, K, RESULT( K )
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

      alasvm(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( 1X, A, ', UPLO=''', A1, ''', N =', I5, ', type ', I2, ', test ', I2, ', ratio =', G12.5 )
      RETURN

      // End of DDRVSY_RK

      }
