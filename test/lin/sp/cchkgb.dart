      void cchkgb(final int DOTYPE, final int NM, final int MVAL, final int NN, final int NVAL, final int NNB, final int NBVAL, final int NNS, final int NSVAL, final int THRESH, final int TSTERR, final int A, final int LA, final int AFAC, final int LAFAC, final int B, final int X, final int XACT, final Array<double> _WORK, final Array<double> RWORK, final Array<int> IWORK, final int NOUT,) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      bool               TSTERR;
      int                LA, LAFAC, NM, NN, NNB, NNS, NOUT;
      double               THRESH;
      bool               DOTYPE( * );
      int                IWORK( * ), MVAL( * ), NBVAL( * ), NSVAL( * ), NVAL( * );
      double               RWORK( * );
      Complex            A( * ), AFAC( * ), B( * ), WORK( * ), X( * ), XACT( * );
      // ..

      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                NTYPES, NTESTS;
      const              NTYPES = 8, NTESTS = 7 ;
      int                NBW, NTRAN;
      const              NBW = 4, NTRAN = 3 ;
      bool               TRFCON, ZEROT;
      String             DIST, NORM, TRANS, TYPE, XTYPE;
      String             PATH;
      int                I, I1, I2, IKL, IKU, IM, IMAT, IN, INB, INFO, IOFF, IRHS, ITRAN, IZERO, J, K, KL, KOFF, KU, LDA, LDAFAC, LDB, M, MODE, N, NB, NERRS, NFAIL, NIMAT, NKL, NKU, NRHS, NRUN;
      double               AINVNM, ANORM, ANORMI, ANORMO, CNDNUM, RCOND, RCONDC, RCONDI, RCONDO;
      String             TRANSS( NTRAN );
      int                ISEED( 4 ), ISEEDY( 4 ), KLVAL( NBW ), KUVAL( NBW );
      double               RESULT( NTESTS );
      // ..
      // .. External Functions ..
      //- REAL               CLANGB, CLANGE, SGET06;
      // EXTERNAL CLANGB, CLANGE, SGET06
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, CCOPY, CERRGE, CGBCON, CGBRFS, CGBT01, CGBT02, CGBT05, CGBTRF, CGBTRS, CGET04, CLACPY, CLARHS, CLASET, CLATB4, CLATMS, XLAENV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, MAX, MIN
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String            srnamc.SRNAMT;
      int                INFOT, NUNIT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / INFOT, NUNIT, OK, LERR
      // COMMON / SRNAMC /srnamc.SRNAMT
      // ..
      // .. Data statements ..
      const ISEEDY = 1988, 1989, 1990, 1991, TRANSS = 'N', 'T', 'C';

      // Initialize constants and the random number seed.

      PATH[1: 1] = 'Complex precision';
      PATH[2: 3] = 'GB';
      NRUN = 0;
      NFAIL = 0;
      NERRS = 0;
      for (I = 1; I <= 4; I++) { // 10
         ISEED[I] = ISEEDY( I );
      } // 10

      // Test the error exits

      if (TSTERR) cerrge( PATH, NOUT );
      INFOT = 0;

      // Initialize the first value for the lower and upper bandwidths.

      KLVAL[1] = 0;
      KUVAL[1] = 0;

      // Do for each value of M in MVAL

      for (IM = 1; IM <= NM; IM++) { // 160
         M = MVAL( IM );

         // Set values to use for the lower bandwidth.

         KLVAL[2] = M + ( M+1 ) / 4;

         // KLVAL( 2 ) = max( M-1, 0 )

         KLVAL[3] = ( 3*M-1 ) / 4;
         KLVAL[4] = ( M+1 ) / 4;

         // Do for each value of N in NVAL

         for (IN = 1; IN <= NN; IN++) { // 150
            N = NVAL( IN );
            XTYPE = 'N';

            // Set values to use for the upper bandwidth.

            KUVAL[2] = N + ( N+1 ) / 4;

            // KUVAL( 2 ) = max( N-1, 0 )

            KUVAL[3] = ( 3*N-1 ) / 4;
            KUVAL[4] = ( N+1 ) / 4;

            // Set limits on the number of loop iterations.

            NKL = min( M+1, 4 );
            if (N == 0) NKL = 2;
            NKU = min( N+1, 4 );
            if (M == 0) NKU = 2;
            NIMAT = NTYPES;
            if (M <= 0 || N <= 0) NIMAT = 1;

            for (IKL = 1; IKL <= NKL; IKL++) { // 140

               // Do for KL = 0, (5*M+1)/4, (3M-1)/4, and (M+1)/4. This
               // order makes it easier to skip redundant values for small
               // values of M.

               KL = KLVAL( IKL );
               for (IKU = 1; IKU <= NKU; IKU++) { // 130

                  // Do for KU = 0, (5*N+1)/4, (3N-1)/4, and (N+1)/4. This
                  // order makes it easier to skip redundant values for
                  // small values of N.

                  KU = KUVAL( IKU );

                  // Check that A and AFAC are big enough to generate this
                  // matrix.

                  LDA = KL + KU + 1;
                  LDAFAC = 2*KL + KU + 1;
                  if ( ( LDA*N ) > LA || ( LDAFAC*N ) > LAFAC ) {
                     if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                     if ( N*( KL+KU+1 ) > LA ) {
                        WRITE( NOUT, FMT = 9999 )LA, M, N, KL, KU, N*( KL+KU+1 );
                        NERRS = NERRS + 1;
                     }
                     if ( N*( 2*KL+KU+1 ) > LAFAC ) {
                        WRITE( NOUT, FMT = 9998 )LAFAC, M, N, KL, KU, N*( 2*KL+KU+1 );
                        NERRS = NERRS + 1;
                     }
                     GO TO 130;
                  }

                  for (IMAT = 1; IMAT <= NIMAT; IMAT++) { // 120

                     // Do the tests only if DOTYPE( IMAT ) is true.

                     if( !DOTYPE( IMAT ) ) GO TO 120;

                     // Skip types 2, 3, or 4 if the matrix size is too
                     // small.

                     ZEROT = IMAT >= 2 && IMAT <= 4;
                     if (ZEROT && N < IMAT-1) GO TO 120;

                     if ( !ZEROT || !DOTYPE( 1 ) ) {

                        // Set up parameters with CLATB4 and generate a
                        // test matrix with CLATMS.

                        clatb4(PATH, IMAT, M, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );

                        KOFF = max( 1, KU+2-N );
                        for (I = 1; I <= KOFF - 1; I++) { // 20
                           A[I] = ZERO;
                        } // 20
                       srnamc.SRNAMT = 'CLATMS';
                        clatms(M, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, 'Z', A( KOFF ), LDA, WORK, INFO );

                        // Check the error code from CLATMS.

                        if ( INFO != 0 ) {
                           alaerh(PATH, 'CLATMS', INFO, 0, ' ', M, N, KL, KU, -1, IMAT, NFAIL, NERRS, NOUT );
                           GO TO 120;
                        }
                     } else if ( IZERO > 0 ) {

                        // Use the same matrix for types 3 and 4 as for
                        // type 2 by copying back the zeroed out column.

                        ccopy(I2-I1+1, B, 1, A( IOFF+I1 ), 1 );
                     }

                     // For types 2, 3, and 4, zero one or more columns of
                     // the matrix to test that INFO is returned correctly.

                     IZERO = 0;
                     if ( ZEROT ) {
                        if ( IMAT == 2 ) {
                           IZERO = 1;
                        } else if ( IMAT == 3 ) {
                           IZERO = min( M, N );
                        } else {
                           IZERO = min( M, N ) / 2 + 1;
                        }
                        IOFF = ( IZERO-1 )*LDA;
                        if ( IMAT < 4 ) {

                           // Store the column to be zeroed out in B.

                           I1 = max( 1, KU+2-IZERO );
                           I2 = min( KL+KU+1, KU+1+( M-IZERO ) );
                           ccopy(I2-I1+1, A( IOFF+I1 ), 1, B, 1 );

                           for (I = I1; I <= I2; I++) { // 30
                              A[IOFF+I] = ZERO;
                           } // 30
                        } else {
                           for (J = IZERO; J <= N; J++) { // 50
                              DO 40 I = max( 1, KU+2-J ), min( KL+KU+1, KU+1+( M-J ) );
                                 A[IOFF+I] = ZERO;
                              } // 40
                              IOFF = IOFF + LDA;
                           } // 50
                        }
                     }

                     // These lines, if used in place of the calls in the
                     // loop over INB, cause the code to bomb on a Sun
                     // SPARCstation.

                      // ANORMO = CLANGB( 'O', N, KL, KU, A, LDA, RWORK )
                      // ANORMI = CLANGB( 'I', N, KL, KU, A, LDA, RWORK )

                     // Do for each blocksize in NBVAL

                     for (INB = 1; INB <= NNB; INB++) { // 110
                        NB = NBVAL( INB );
                        xlaenv(1, NB );

                        // Compute the LU factorization of the band matrix.

                        if (M > 0 && N > 0) clacpy( 'Full', KL+KU+1, N, A, LDA, AFAC( KL+1 ), LDAFAC );
                       srnamc.SRNAMT = 'CGBTRF';
                        cgbtrf(M, N, KL, KU, AFAC, LDAFAC, IWORK, INFO );

                        // Check error code from CGBTRF.

                        if (INFO != IZERO) alaerh( PATH, 'CGBTRF', INFO, IZERO, ' ', M, N, KL, KU, NB, IMAT, NFAIL, NERRS, NOUT );
                        TRFCON = false;

// +    TEST 1
                        // Reconstruct matrix from factors and compute
                        // residual.

                        cgbt01(M, N, KL, KU, A, LDA, AFAC, LDAFAC, IWORK, WORK, RESULT( 1 ) );

                        // Print information about the tests so far that
                        // did not pass the threshold.

                        if ( RESULT( 1 ) >= THRESH ) {
                           if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                           WRITE( NOUT, FMT = 9997 )M, N, KL, KU, NB, IMAT, 1, RESULT( 1 );
                           NFAIL = NFAIL + 1;
                        }
                        NRUN = NRUN + 1;

                        // Skip the remaining tests if this is not the
                        // first block size or if M != N.

                        if (INB > 1 || M != N) GO TO 110;

                        ANORMO = CLANGB( 'O', N, KL, KU, A, LDA, RWORK );
                        ANORMI = CLANGB( 'I', N, KL, KU, A, LDA, RWORK );

                        if ( INFO == 0 ) {

                           // Form the inverse of A so we can get a good
                           // estimate of CNDNUM = norm(A) * norm(inv(A)).

                           LDB = max( 1, N );
                           claset('Full', N, N, CMPLX( ZERO ), CMPLX( ONE ), WORK, LDB );
                          srnamc.SRNAMT = 'CGBTRS';
                           cgbtrs('No transpose', N, KL, KU, N, AFAC, LDAFAC, IWORK, WORK, LDB, INFO );

                           // Compute the 1-norm condition number of A.

                           AINVNM = CLANGE( 'O', N, N, WORK, LDB, RWORK );
                           if ( ANORMO <= ZERO || AINVNM <= ZERO ) {
                              RCONDO = ONE;
                           } else {
                              RCONDO = ( ONE / ANORMO ) / AINVNM;
                           }

                           // Compute the infinity-norm condition number of
                           // A.

                           AINVNM = CLANGE( 'I', N, N, WORK, LDB, RWORK );
                           if ( ANORMI <= ZERO || AINVNM <= ZERO ) {
                              RCONDI = ONE;
                           } else {
                              RCONDI = ( ONE / ANORMI ) / AINVNM;
                           }
                        } else {

                           // Do only the condition estimate if INFO != 0.

                           TRFCON = true;
                           RCONDO = ZERO;
                           RCONDI = ZERO;
                        }

                        // Skip the solve tests if the matrix is singular.

                        if (TRFCON) GO TO 90;

                        for (IRHS = 1; IRHS <= NNS; IRHS++) { // 80
                           NRHS = NSVAL( IRHS );
                           XTYPE = 'N';

                           for (ITRAN = 1; ITRAN <= NTRAN; ITRAN++) { // 70
                              TRANS = TRANSS( ITRAN );
                              if ( ITRAN == 1 ) {
                                 RCONDC = RCONDO;
                                 NORM = 'O';
                              } else {
                                 RCONDC = RCONDI;
                                 NORM = 'I';
                              }

// +    TEST 2:
                              // Solve and compute residual for op(A) * X = B.

                             srnamc.SRNAMT = 'CLARHS';
                              clarhs(PATH, XTYPE, ' ', TRANS, N, N, KL, KU, NRHS, A, LDA, XACT, LDB, B, LDB, ISEED, INFO );
                              XTYPE = 'C';
                              clacpy('Full', N, NRHS, B, LDB, X, LDB );

                             srnamc.SRNAMT = 'CGBTRS';
                              cgbtrs(TRANS, N, KL, KU, NRHS, AFAC, LDAFAC, IWORK, X, LDB, INFO );

                              // Check error code from CGBTRS.

                              if (INFO != 0) alaerh( PATH, 'CGBTRS', INFO, 0, TRANS, N, N, KL, KU, -1, IMAT, NFAIL, NERRS, NOUT );

                              clacpy('Full', N, NRHS, B, LDB, WORK, LDB );
                              cgbt02(TRANS, M, N, KL, KU, NRHS, A, LDA, X, LDB, WORK, LDB, RWORK, RESULT( 2 ) );

// +    TEST 3:
                              // Check solution from generated exact
                              // solution.

                              cget04(N, NRHS, X, LDB, XACT, LDB, RCONDC, RESULT( 3 ) );

// +    TESTS 4, 5, 6:
                              // Use iterative refinement to improve the
                              // solution.

                             srnamc.SRNAMT = 'CGBRFS';
                              cgbrfs(TRANS, N, KL, KU, NRHS, A, LDA, AFAC, LDAFAC, IWORK, B, LDB, X, LDB, RWORK, RWORK( NRHS+1 ), WORK, RWORK( 2*NRHS+1 ), INFO );

                              // Check error code from CGBRFS.

                              if (INFO != 0) alaerh( PATH, 'CGBRFS', INFO, 0, TRANS, N, N, KL, KU, NRHS, IMAT, NFAIL, NERRS, NOUT );

                              cget04(N, NRHS, X, LDB, XACT, LDB, RCONDC, RESULT( 4 ) );
                              cgbt05(TRANS, N, KL, KU, NRHS, A, LDA, B, LDB, X, LDB, XACT, LDB, RWORK, RWORK( NRHS+1 ), RESULT( 5 ) );

                              // Print information about the tests that did
                              // not pass the threshold.

                              for (K = 2; K <= 6; K++) { // 60
                                 if ( RESULT( K ) >= THRESH ) {
                                    if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                                    WRITE( NOUT, FMT = 9996 )TRANS, N, KL, KU, NRHS, IMAT, K, RESULT( K );
                                    NFAIL = NFAIL + 1;
                                 }
                              } // 60
                              NRUN = NRUN + 5;
                           } // 70
                        } // 80

// +    TEST 7:
                           // Get an estimate of RCOND = 1/CNDNUM.

                        } // 90
                        for (ITRAN = 1; ITRAN <= 2; ITRAN++) { // 100
                           if ( ITRAN == 1 ) {
                              ANORM = ANORMO;
                              RCONDC = RCONDO;
                              NORM = 'O';
                           } else {
                              ANORM = ANORMI;
                              RCONDC = RCONDI;
                              NORM = 'I';
                           }
                          srnamc.SRNAMT = 'CGBCON';
                           cgbcon(NORM, N, KL, KU, AFAC, LDAFAC, IWORK, ANORM, RCOND, WORK, RWORK, INFO );

                              // Check error code from CGBCON.

                           if (INFO != 0) alaerh( PATH, 'CGBCON', INFO, 0, NORM, N, N, KL, KU, -1, IMAT, NFAIL, NERRS, NOUT );

                           RESULT[7] = SGET06( RCOND, RCONDC );

                           // Print information about the tests that did
                           // not pass the threshold.

                           if ( RESULT( 7 ) >= THRESH ) {
                              if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                              WRITE( NOUT, FMT = 9995 )NORM, N, KL, KU, IMAT, 7, RESULT( 7 );
                              NFAIL = NFAIL + 1;
                           }
                           NRUN = NRUN + 1;
                        } // 100
                     } // 110
                  } // 120
               } // 130
            } // 140
         } // 150
      } // 160

      // Print a summary of the results.

      alasum(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( ' *** In CCHKGB, LA=${.i5} is too small for M=${.i5}, N=${.i5}, KL=${.i4}, KU=', I4, / ' ==> Increase LA to at least ${.i5}');
 9998 FORMAT( ' *** In CCHKGB, LAFAC=${.i5} is too small for M=${.i5}, N=${.i5}, KL=${.i4}, KU=', I4, / ' ==> Increase LAFAC to at least ${.i5}');
 9997 FORMAT( ' M =${.i5}, N =${.i5}, KL=${.i5}, KU=${.i5}, NB =${.i4}, type ${.i1}, test(${.i1})=${.g12_5};
 9996 FORMAT( ' TRANS=''${.a1}'', N=${.i5}, KL=${.i5}, KU=${.i5}, NRHS=${.i3}, type ${.i1}, test(${.i1})=${.g12_5};
 9995 FORMAT( ' NORM =''${.a1}'', N=${.i5}, KL=${.i5}, KU=${.i5},${' ' * 10} type ${.i1}, test(${.i1})=${.g12_5};

      }
