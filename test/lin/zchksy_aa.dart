      void zchksy_aa(final Array<bool> DOTYPE_, final int NN, final Array<int> NVAL_, final int NNB, final Array<int> NBVAL_, final int NNS, final Array<int> NSVAL_, final double THRESH, final bool TSTERR, final int NMAX, final Array<double> A_, final Array<double> AFAC_, final Array<double> AINV_, final Array<double> B_, final Array<double> X_, final Array<double> XACT_, final Array<double> WORK_, final Array<double> RWORK_, final Array<int> IWORK_, final Nout NOUT,) {
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();
  final IWORK = IWORK_.having();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      bool               TSTERR;
      int                NN, NNB, NNS, NMAX, NOUT;
      double             THRESH;
      bool               DOTYPE( * );
      int                IWORK( * ), NBVAL( * ), NSVAL( * ), NVAL( * );
      double             RWORK( * );
      Complex         A( * ), AFAC( * ), AINV( * ), B( * ), WORK( * ), X( * ), XACT( * );
      // ..

      double             ZERO;
      const              ZERO = 0.0 ;
      Complex         CZERO;
      const              CZERO = ( 0.0, 0.0 ) ;
      int                NTYPES;
      const              NTYPES = 10 ;
      int                NTESTS;
      const              NTESTS = 9 ;
      bool               ZEROT;
      String             DIST, TYPE, UPLO, XTYPE;
      String             PATH, MATPATH;
      int                I, I1, I2, IMAT, IN, INB, INFO, IOFF, IRHS, IUPLO, IZERO, J, K, KL, KU, LDA, LWORK, MODE, N, NB, NERRS, NFAIL, NIMAT, NRHS, NRUN, NT;
      double             ANORM, CNDNUM;
      String             UPLOS( 2 );
      final                ISEED=Array<int>( 4 ), ISEEDY( 4 );
      final             RESULT=Array<double>( NTESTS );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, ZERRSY, ZLACPY, ZLARHS, ZLATB4, ZLATMS, ZSYT02, ZSYT01_AA, ZSYTRF_AA, ZSYTRS_AA, XLAENV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Scalars in Common ..
      bool               infoc.LERR, infoc.OK;
      String            srnamc.SRNAMT;
      int                infoc.INFOT, infoc.NUNIT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / infoc.INFOT, infoc.NUNIT, infoc.OK, infoc.LERR
      // COMMON / SRNAMC /srnamc.SRNAMT
      // ..
      // .. Data statements ..
      const ISEEDY = [ 1988, 1989, 1990, 1991 ];
      const UPLOS = [ 'U', 'L' ];

      // Initialize constants and the random number seed.

      // Test path

      PATH[1: 1] = 'Zomplex precision';
      PATH[2: 3] = 'SA';

      // Path to generate matrices

      MATPATH[1: 1] = 'Zomplex precision';
      MATPATH[2: 3] = 'SY';
      NRUN = 0;
      NFAIL = 0;
      NERRS = 0;
      for (I = 1; I <= 4; I++) { // 10
         ISEED[I] = ISEEDY( I );
      } // 10

      // Test the error exits

      if (TSTERR) zerrsy( PATH, NOUT );
      infoc.INFOT = 0;

      // Set the minimum block size for which the block routine should
      // be used, which will be later returned by ILAENV

      xlaenv(2, 2 );

      // Do for each value of N in NVAL

      for (IN = 1; IN <= NN; IN++) { // 180
         N = NVAL( IN );
         if ( N > NMAX ) {
            NFAIL = NFAIL + 1;
            WRITE(NOUT, 9995) 'M ', N, NMAX;
            GO TO 180;
         }
         LDA = max( N, 1 );
         XTYPE = 'N';
         NIMAT = NTYPES;
         if (N <= 0) NIMAT = 1;

         IZERO = 0;

         // Do for each value of matrix type IMAT

         for (IMAT = 1; IMAT <= NIMAT; IMAT++) { // 170

            // Do the tests only if DOTYPE( IMAT ) is true.

            if( !DOTYPE( IMAT ) ) GO TO 170;

            // Skip types 3, 4, 5, or 6 if the matrix size is too small.

            ZEROT = IMAT >= 3 && IMAT <= 6;
            if (ZEROT && N < IMAT-2) GO TO 170;

            // Do first for UPLO = 'U', then for UPLO = 'L'

            for (IUPLO = 1; IUPLO <= 2; IUPLO++) { // 160
               UPLO = UPLOS( IUPLO );

               // Begin generate the test matrix A.


               // Set up parameters with ZLATB4 for the matrix generator
               // based on the type of matrix to be generated.

               zlatb4(MATPATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );

               // Generate a matrix with ZLATMS.

              srnamc.SRNAMT = 'ZLATMS';
               zlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, UPLO, A, LDA, WORK, INFO );

               // Check error code from ZLATMS and handle error.

               if ( INFO != 0 ) {
                  alaerh(PATH, 'ZLATMS', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

                     // Skip all tests for this generated matrix

                  GO TO 160;
               }

               // For matrix types 3-6, zero one or more rows and
               // columns of the matrix to test that INFO is returned
               // correctly.

               if ( ZEROT ) {
                  if ( IMAT == 3 ) {
                     IZERO = 1;
                  } else if ( IMAT == 4 ) {
                     IZERO = N;
                  } else {
                     IZERO = N / 2 + 1;
                  }

                  if ( IMAT < 6 ) {

                     // Set row and column IZERO to zero.

                     if ( IUPLO == 1 ) {
                        IOFF = ( IZERO-1 )*LDA;
                        for (I = 1; I <= IZERO - 1; I++) { // 20
                           A[IOFF+I] = CZERO;
                        } // 20
                        IOFF = IOFF + IZERO;
                        for (I = IZERO; I <= N; I++) { // 30
                           A[IOFF] = CZERO;
                           IOFF = IOFF + LDA;
                        } // 30
                     } else {
                        IOFF = IZERO;
                        for (I = 1; I <= IZERO - 1; I++) { // 40
                           A[IOFF] = CZERO;
                           IOFF = IOFF + LDA;
                        } // 40
                        IOFF = IOFF - IZERO;
                        for (I = IZERO; I <= N; I++) { // 50
                           A[IOFF+I] = CZERO;
                        } // 50
                     }
                  } else {
                     if ( IUPLO == 1 ) {

                        // Set the first IZERO rows and columns to zero.

                        IOFF = 0;
                        for (J = 1; J <= N; J++) { // 70
                           I2 = min( J, IZERO );
                           for (I = 1; I <= I2; I++) { // 60
                              A[IOFF+I] = CZERO;
                           } // 60
                           IOFF = IOFF + LDA;
                        } // 70
                        IZERO = 1;
                     } else {

                        // Set the last IZERO rows and columns to zero.

                        IOFF = 0;
                        for (J = 1; J <= N; J++) { // 90
                           I1 = max( J, IZERO );
                           for (I = I1; I <= N; I++) { // 80
                              A[IOFF+I] = CZERO;
                           } // 80
                           IOFF = IOFF + LDA;
                        } // 90
                     }
                  }
               } else {
                  IZERO = 0;
               }

               // End generate the test matrix A.

               // Do for each value of NB in NBVAL

               for (INB = 1; INB <= NNB; INB++) { // 150

                  // Set the optimal blocksize, which will be later
                  // returned by ILAENV.

                  NB = NBVAL( INB );
                  xlaenv(1, NB );

                  // Copy the test matrix A into matrix AFAC which
                  // will be factorized in place. This is needed to
                  // preserve the test matrix A for subsequent tests.

                  zlacpy(UPLO, N, N, A, LDA, AFAC, LDA );

                  // Compute the L*D*L**T or U*D*U**T factorization of the
                  // matrix. IWORK stores details of the interchanges and
                  // the block structure of D. AINV is a work array for
                  // block factorization, LWORK is the length of AINV.

                 srnamc.SRNAMT = 'ZSYTRF_AA';
                  LWORK = max( 1, N*NB + N );
                  zsytrf_aa(UPLO, N, AFAC, LDA, IWORK, AINV, LWORK, INFO );

                  // Adjust the expected value of INFO to account for
                  // pivoting.

                   // IF( IZERO > 0 ) THEN
                   //    J = 1
                   //    K = IZERO
// 100                CONTINUE
                      // IF( J == K ) THEN
                      //    K = IWORK( J )
                      // ELSE IF( IWORK( J ) == K ) THEN
                      //    K = J
                      // END IF
                      // IF( J < K ) THEN
                      //    J = J + 1
                      //    GO TO 100
                      // END IF
                   // ELSE
                     K = 0;
                   // END IF

                  // Check error code from ZSYTRF and handle error.

                  if ( INFO != K ) {
                     alaerh(PATH, 'ZSYTRF_AA', INFO, K, UPLO, N, N, -1, -1, NB, IMAT, NFAIL, NERRS, NOUT );
                  }

// +    TEST 1
                  // Reconstruct matrix from factors and compute residual.

                  zsyt01_aa(UPLO, N, A, LDA, AFAC, LDA, IWORK, AINV, LDA, RWORK, RESULT( 1 ) );
                  NT = 1;


                  // Print information about the tests that did not pass
                  // the threshold.

                  for (K = 1; K <= NT; K++) { // 110
                     if ( RESULT( K ) >= THRESH ) {
                        if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                        WRITE( NOUT, FMT = 9999 )UPLO, N, NB, IMAT, K, RESULT( K );
                        NFAIL = NFAIL + 1;
                     }
                  } // 110
                  NRUN = NRUN + NT;

                  // Skip solver test if INFO is not 0.

                  if ( INFO != 0 ) {
                     GO TO 140;
                  }

                  // Do for each value of NRHS in NSVAL.

                  for (IRHS = 1; IRHS <= NNS; IRHS++) { // 130
                     NRHS = NSVAL( IRHS );

// +    TEST 2 (Using TRS)
                  // Solve and compute residual for  A * X = B.

                     // Choose a set of NRHS random solution vectors
                     // stored in XACT and set up the right hand side B

                    srnamc.SRNAMT = 'ZLARHS';
                     zlarhs(MATPATH, XTYPE, UPLO, ' ', N, N, KL, KU, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO );
                     zlacpy('Full', N, NRHS, B, LDA, X, LDA );

                    srnamc.SRNAMT = 'ZSYTRS_AA';
                     LWORK = max( 1, 3*N-2 );
                     zsytrs_aa(UPLO, N, NRHS, AFAC, LDA, IWORK, X, LDA, WORK, LWORK, INFO );

                     // Check error code from ZSYTRS and handle error.

                     if ( INFO != 0 ) {
                        if ( IZERO == 0 ) {
                           alaerh(PATH, 'ZSYTRS_AA', INFO, 0, UPLO, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );
                        }
                     } else {
                        zlacpy('Full', N, NRHS, B, LDA, WORK, LDA );

                        // Compute the residual for the solution

                        zsyt02(UPLO, N, NRHS, A, LDA, X, LDA, WORK, LDA, RWORK, RESULT( 2 ) );


                        // Print information about the tests that did not pass
                        // the threshold.

                        for (K = 2; K <= 2; K++) { // 120
                           if ( RESULT( K ) >= THRESH ) {
                              if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                              WRITE( NOUT, FMT = 9998 )UPLO, N, NRHS, IMAT, K, RESULT( K );
                              NFAIL = NFAIL + 1;
                           }
                        } // 120
                     }
                     NRUN = NRUN + 1;

                  // End do for each value of NRHS in NSVAL.

                  } // 130
                  } // 140
               } // 150
            } // 160
         } // 170
      } // 180

      // Print a summary of the results.

      alasum(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( ' UPLO = \'${.a1}\', N =${.i5}, NB =${.i4}, type ${.i2}, test ${.i2}, ratio =${.g12_5};
 9998 FORMAT( ' UPLO = \'${.a1}\', N =${.i5}, NRHS=${.i3}, type ${.i2}, test(${.i2}) =${.g12_5};
 9995 FORMAT( ' Invalid input value: ${.a4}=${.i6}; must be <=', I6 )
      }
