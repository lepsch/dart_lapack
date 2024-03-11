      void zchkpb(final Array<bool> DOTYPE_, final int NN, final Array<int> NVAL_, final int NNB, final Array<int> NBVAL_, final int NNS, final Array<int> NSVAL_, final double THRESH, final bool TSTERR, final int NMAX, final Array<double> A_, final Array<double> AFAC_, final Array<double> AINV_, final Array<double> B_, final Array<double> X_, final Array<double> XACT_, final Array<double> WORK_, final Array<double> RWORK_, final Nout NOUT,) {
  final WORK = WORK_.having();
  final RWORK = RWORK_.having();

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      bool               TSTERR;
      int                NMAX, NN, NNB, NNS, NOUT;
      double             THRESH;
      bool               DOTYPE( * );
      int                NBVAL( * ), NSVAL( * ), NVAL( * );
      double             RWORK( * );
      Complex         A( * ), AFAC( * ), AINV( * ), B( * ), WORK( * ), X( * ), XACT( * );
      // ..

      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                NTYPES, NTESTS;
      const              NTYPES = 8, NTESTS = 7 ;
      int                NBW;
      const              NBW = 4 ;
      bool               ZEROT;
      String             DIST, PACKIT, TYPE, UPLO, XTYPE;
      String             PATH;
      int                I, I1, I2, IKD, IMAT, IN, INB, INFO, IOFF, IRHS, IUPLO, IW, IZERO, K, KD, KL, KOFF, KU, LDA, LDAB, MODE, N, NB, NERRS, NFAIL, NIMAT, NKD, NRHS, NRUN;
      double             AINVNM, ANORM, CNDNUM, RCOND, RCONDC;
      final                ISEED=Array<int>( 4 ), ISEEDY( 4 ), KDVAL( NBW );
      final             RESULT=Array<double>( NTESTS );
      // ..
      // .. External Functions ..
      //- double             DGET06, ZLANGE, ZLANHB;
      // EXTERNAL DGET06, ZLANGE, ZLANHB
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, XLAENV, ZCOPY, ZERRPO, ZGET04, ZLACPY, ZLAIPD, ZLARHS, ZLASET, ZLATB4, ZLATMS, ZPBCON, ZPBRFS, ZPBT01, ZPBT02, ZPBT05, ZPBTRF, ZPBTRS, ZSWAP
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DCMPLX, MAX, MIN
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

      // Initialize constants and the random number seed.

      PATH[1: 1] = 'Zomplex precision';
      PATH[2: 3] = 'PB';
      NRUN = 0;
      NFAIL = 0;
      NERRS = 0;
      for (I = 1; I <= 4; I++) { // 10
         ISEED[I] = ISEEDY( I );
      } // 10

      // Test the error exits

      if (TSTERR) zerrpo( PATH, NOUT );
      infoc.INFOT = 0;
      KDVAL[1] = 0;

      // Do for each value of N in NVAL

      for (IN = 1; IN <= NN; IN++) { // 90
         N = NVAL( IN );
         LDA = max( N, 1 );
         XTYPE = 'N';

         // Set limits on the number of loop iterations.

         NKD = max( 1, min( N, 4 ) );
         NIMAT = NTYPES;
         if (N == 0) NIMAT = 1;

         KDVAL[2] = N + ( N+1 ) / 4;
         KDVAL[3] = ( 3*N-1 ) / 4;
         KDVAL[4] = ( N+1 ) / 4;

         for (IKD = 1; IKD <= NKD; IKD++) { // 80

            // Do for KD = 0, (5*N+1)/4, (3N-1)/4, and (N+1)/4. This order
            // makes it easier to skip redundant values for small values
            // of N.

            KD = KDVAL( IKD );
            LDAB = KD + 1;

            // Do first for UPLO = 'U', then for UPLO = 'L'

            for (IUPLO = 1; IUPLO <= 2; IUPLO++) { // 70
               KOFF = 1;
               if ( IUPLO == 1 ) {
                  UPLO = 'U';
                  KOFF = max( 1, KD+2-N );
                  PACKIT = 'Q';
               } else {
                  UPLO = 'L';
                  PACKIT = 'B';
               }

               for (IMAT = 1; IMAT <= NIMAT; IMAT++) { // 60

                  // Do the tests only if DOTYPE( IMAT ) is true.

                  if( !DOTYPE( IMAT ) ) GO TO 60;

                  // Skip types 2, 3, or 4 if the matrix size is too small.

                  ZEROT = IMAT >= 2 && IMAT <= 4;
                  if (ZEROT && N < IMAT-1) GO TO 60;

                  if ( !ZEROT || !DOTYPE( 1 ) ) {

                     // Set up parameters with ZLATB4 and generate a test
                     // matrix with ZLATMS.

                     zlatb4(PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );

                    srnamc.SRNAMT = 'ZLATMS';
                     zlatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KD, KD, PACKIT, A( KOFF ), LDAB, WORK, INFO );

                     // Check error code from ZLATMS.

                     if ( INFO != 0 ) {
                        alaerh(PATH, 'ZLATMS', INFO, 0, UPLO, N, N, KD, KD, -1, IMAT, NFAIL, NERRS, NOUT );
                        GO TO 60;
                     }
                  } else if ( IZERO > 0 ) {

                     // Use the same matrix for types 3 and 4 as for type
                     // 2 by copying back the zeroed out column,

                     IW = 2*LDA + 1;
                     if ( IUPLO == 1 ) {
                        IOFF = ( IZERO-1 )*LDAB + KD + 1;
                        zcopy(IZERO-I1, WORK( IW ), 1, A( IOFF-IZERO+I1 ), 1 );
                        IW = IW + IZERO - I1;
                        zcopy(I2-IZERO+1, WORK( IW ), 1, A( IOFF ), max( LDAB-1, 1 ) );
                     } else {
                        IOFF = ( I1-1 )*LDAB + 1;
                        zcopy(IZERO-I1, WORK( IW ), 1, A( IOFF+IZERO-I1 ), max( LDAB-1, 1 ) );
                        IOFF = ( IZERO-1 )*LDAB + 1;
                        IW = IW + IZERO - I1;
                        zcopy(I2-IZERO+1, WORK( IW ), 1, A( IOFF ), 1 );
                     }
                  }

                  // For types 2-4, zero one row and column of the matrix
                  // to test that INFO is returned correctly.

                  IZERO = 0;
                  if ( ZEROT ) {
                     if ( IMAT == 2 ) {
                        IZERO = 1;
                     } else if ( IMAT == 3 ) {
                        IZERO = N;
                     } else {
                        IZERO = N / 2 + 1;
                     }

                     // Save the zeroed out row and column in WORK(*,3)

                     IW = 2*LDA;
                     for (I = 1; I <= min( 2*KD+1, N ); I++) { // 20
                        WORK[IW+I] = ZERO;
                     } // 20
                     IW = IW + 1;
                     I1 = max( IZERO-KD, 1 );
                     I2 = min( IZERO+KD, N );

                     if ( IUPLO == 1 ) {
                        IOFF = ( IZERO-1 )*LDAB + KD + 1;
                        zswap(IZERO-I1, A( IOFF-IZERO+I1 ), 1, WORK( IW ), 1 );
                        IW = IW + IZERO - I1;
                        zswap(I2-IZERO+1, A( IOFF ), max( LDAB-1, 1 ), WORK( IW ), 1 );
                     } else {
                        IOFF = ( I1-1 )*LDAB + 1;
                        zswap(IZERO-I1, A( IOFF+IZERO-I1 ), max( LDAB-1, 1 ), WORK( IW ), 1 );
                        IOFF = ( IZERO-1 )*LDAB + 1;
                        IW = IW + IZERO - I1;
                        zswap(I2-IZERO+1, A( IOFF ), 1, WORK( IW ), 1 );
                     }
                  }

                  // Set the imaginary part of the diagonals.

                  if ( IUPLO == 1 ) {
                     zlaipd(N, A( KD+1 ), LDAB, 0 );
                  } else {
                     zlaipd(N, A( 1 ), LDAB, 0 );
                  }

                  // Do for each value of NB in NBVAL

                  for (INB = 1; INB <= NNB; INB++) { // 50
                     NB = NBVAL( INB );
                     xlaenv(1, NB );

                     // Compute the L*L' or U'*U factorization of the band
                     // matrix.

                     zlacpy('Full', KD+1, N, A, LDAB, AFAC, LDAB );
                    srnamc.SRNAMT = 'ZPBTRF';
                     zpbtrf(UPLO, N, KD, AFAC, LDAB, INFO );

                     // Check error code from ZPBTRF.

                     if ( INFO != IZERO ) {
                        alaerh(PATH, 'ZPBTRF', INFO, IZERO, UPLO, N, N, KD, KD, NB, IMAT, NFAIL, NERRS, NOUT );
                        GO TO 50;
                     }

                     // Skip the tests if INFO is not 0.

                     if (INFO != 0) GO TO 50;

// +    TEST 1
                     // Reconstruct matrix from factors and compute
                     // residual.

                     zlacpy('Full', KD+1, N, AFAC, LDAB, AINV, LDAB );
                     zpbt01(UPLO, N, KD, A, LDAB, AINV, LDAB, RWORK, RESULT( 1 ) );

                     // Print the test ratio if it is >= THRESH.

                     if ( RESULT( 1 ) >= THRESH ) {
                        if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                        WRITE( NOUT, FMT = 9999 )UPLO, N, KD, NB, IMAT, 1, RESULT( 1 );
                        NFAIL = NFAIL + 1;
                     }
                     NRUN = NRUN + 1;

                     // Only do other tests if this is the first blocksize.

                     if (INB > 1) GO TO 50;

                     // Form the inverse of A so we can get a good estimate
                     // of RCONDC = 1/(norm(A) * norm(inv(A))).

                     zlaset('Full', N, N, DCMPLX( ZERO ), DCMPLX( ONE ), AINV, LDA );
                    srnamc.SRNAMT = 'ZPBTRS';
                     zpbtrs(UPLO, N, KD, N, AFAC, LDAB, AINV, LDA, INFO );

                     // Compute RCONDC = 1/(norm(A) * norm(inv(A))).

                     ANORM = ZLANHB( '1', UPLO, N, KD, A, LDAB, RWORK );
                     AINVNM = ZLANGE( '1', N, N, AINV, LDA, RWORK );
                     if ( ANORM <= ZERO || AINVNM <= ZERO ) {
                        RCONDC = ONE;
                     } else {
                        RCONDC = ( ONE / ANORM ) / AINVNM;
                     }

                     for (IRHS = 1; IRHS <= NNS; IRHS++) { // 40
                        NRHS = NSVAL( IRHS );

// +    TEST 2
                     // Solve and compute residual for A * X = B.

                       srnamc.SRNAMT = 'ZLARHS';
                        zlarhs(PATH, XTYPE, UPLO, ' ', N, N, KD, KD, NRHS, A, LDAB, XACT, LDA, B, LDA, ISEED, INFO );
                        zlacpy('Full', N, NRHS, B, LDA, X, LDA );

                       srnamc.SRNAMT = 'ZPBTRS';
                        zpbtrs(UPLO, N, KD, NRHS, AFAC, LDAB, X, LDA, INFO );

                     // Check error code from ZPBTRS.

                        if (INFO != 0) alaerh( PATH, 'ZPBTRS', INFO, 0, UPLO, N, N, KD, KD, NRHS, IMAT, NFAIL, NERRS, NOUT );

                        zlacpy('Full', N, NRHS, B, LDA, WORK, LDA );
                        zpbt02(UPLO, N, KD, NRHS, A, LDAB, X, LDA, WORK, LDA, RWORK, RESULT( 2 ) );

// +    TEST 3
                     // Check solution from generated exact solution.

                        zget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) );

// +    TESTS 4, 5, and 6
                     // Use iterative refinement to improve the solution.

                       srnamc.SRNAMT = 'ZPBRFS';
                        zpbrfs(UPLO, N, KD, NRHS, A, LDAB, AFAC, LDAB, B, LDA, X, LDA, RWORK, RWORK( NRHS+1 ), WORK, RWORK( 2*NRHS+1 ), INFO );

                     // Check error code from ZPBRFS.

                        if (INFO != 0) alaerh( PATH, 'ZPBRFS', INFO, 0, UPLO, N, N, KD, KD, NRHS, IMAT, NFAIL, NERRS, NOUT );

                        zget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 4 ) );
                        zpbt05(UPLO, N, KD, NRHS, A, LDAB, B, LDA, X, LDA, XACT, LDA, RWORK, RWORK( NRHS+1 ), RESULT( 5 ) );

                        // Print information about the tests that did not
                        // pass the threshold.

                        for (K = 2; K <= 6; K++) { // 30
                           if ( RESULT( K ) >= THRESH ) {
                              if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                              WRITE( NOUT, FMT = 9998 )UPLO, N, KD, NRHS, IMAT, K, RESULT( K );
                              NFAIL = NFAIL + 1;
                           }
                        } // 30
                        NRUN = NRUN + 5;
                     } // 40

// +    TEST 7
                     // Get an estimate of RCOND = 1/CNDNUM.

                    srnamc.SRNAMT = 'ZPBCON';
                     zpbcon(UPLO, N, KD, AFAC, LDAB, ANORM, RCOND, WORK, RWORK, INFO );

                     // Check error code from ZPBCON.

                     if (INFO != 0) alaerh( PATH, 'ZPBCON', INFO, 0, UPLO, N, N, KD, KD, -1, IMAT, NFAIL, NERRS, NOUT );

                     RESULT[7] = DGET06( RCOND, RCONDC );

                     // Print the test ratio if it is >= THRESH.

                     if ( RESULT( 7 ) >= THRESH ) {
                        if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                        WRITE( NOUT, FMT = 9997 )UPLO, N, KD, IMAT, 7, RESULT( 7 );
                        NFAIL = NFAIL + 1;
                     }
                     NRUN = NRUN + 1;
                  } // 50
               } // 60
            } // 70
         } // 80
      } // 90

      // Print a summary of the results.

      alasum(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( ' UPLO=\'${.a1}\', N=${.i5}, KD=${.i5}, NB=${.i4}, type ${.i2}, test ${.i2}, ratio= ${.g12_5};
 9998 FORMAT( ' UPLO=\'${.a1}\', N=${.i5}, KD=${.i5}, NRHS=${.i3}, type ${.i2}, test(${.i2}) = ${.g12_5};
 9997 FORMAT( ' UPLO=\'${.a1}\', N=${.i5}, KD=${.i5},${' ' * 10} type ${.i2}, test(${.i2}) = ${.g12_5};
      }
