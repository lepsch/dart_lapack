import 'common.dart';

      void dchkge(final int DOTYPE, final int NM, final int MVAL, final int NN, final int NVAL, final int NNB, final int NBVAL, final int NNS, final int NSVAL, final int THRESH, final int TSTERR, final int NMAX, final int A, final int AFAC, final int AINV, final int B, final int X, final int XACT, final Array<double> _WORK, final Array<double> RWORK, final Array<int> IWORK, final int NOUT) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      bool               TSTERR;
      int                NM, NMAX, NN, NNB, NNS, NOUT;
      double             THRESH;
      bool               DOTYPE( * );
      int                IWORK( * ), MVAL( * ), NBVAL( * ), NSVAL( * ), NVAL( * );
      double             A( * ), AFAC( * ), AINV( * ), B( * ), RWORK( * ), WORK( * ), X( * ), XACT( * );
      // ..

      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                NTYPES;
      const              NTYPES = 11 ;
      int                NTESTS;
      const              NTESTS = 8 ;
      int                NTRAN;
      const              NTRAN = 3 ;
      bool               TRFCON, ZEROT;
      String             DIST, NORM, TRANS, TYPE, XTYPE;
      String             PATH;
      int                I, IM, IMAT, IN, INB, INFO, IOFF, IRHS, ITRAN, IZERO, K, KL, KU, LDA, LWORK, M, MODE, N, NB, NERRS, NFAIL, NIMAT, NRHS, NRUN, NT;
      double             AINVNM, ANORM, ANORMI, ANORMO, CNDNUM, DUMMY, RCOND, RCONDC, RCONDI, RCONDO;
      String             TRANSS( NTRAN );
      int                ISEED( 4 ), ISEEDY( 4 );
      double             RESULT( NTESTS );
      // ..
      // .. External Functions ..
      //- double             DGET06, DLANGE;
      // EXTERNAL DGET06, DLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, DERRGE, DGECON, DGERFS, DGET01, DGET02, DGET03, DGET04, DGET07, DGETRF, DGETRI, DGETRS, DLACPY, DLARHS, DLASET, DLATB4, DLATMS, XLAENV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Scalars in Common ..
      // bool               infoc.LERR, infoc.OK;
      // String             srnamc.SRNAMT;
      // int                infoc.INFOT, infoc.NUNIT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / infoc.INFOT, infoc.NUNIT, infoc.OK, infoc.LERR
      // COMMON / SRNAMC / srnamc.SRNAMT
      // ..
      // .. Data statements ..
      const ISEEDY = 1988, 1989, 1990, 1991, TRANSS = 'N', 'T', 'C';

      // Initialize constants and the random number seed.

      PATH = '${'Double precision'[0]}';
      PATH[2: 3] = 'GE';
      NRUN = 0;
      NFAIL = 0;
      NERRS = 0;
      for (I = 1; I <= 4; I++) { // 10
         ISEED[I] = ISEEDY( I );
      } // 10

      // Test the error exits

      xlaenv(1, 1 );
      if (TSTERR) derrge( PATH, NOUT );
      infoc.INFOT = 0;
      xlaenv(2, 2 );

      // Do for each value of M in MVAL

      for (IM = 1; IM <= NM; IM++) { // 120
         M = MVAL( IM );
         LDA = max( 1, M );

         // Do for each value of N in NVAL

         for (IN = 1; IN <= NN; IN++) { // 110
            N = NVAL( IN );
            XTYPE = 'N';
            NIMAT = NTYPES;
            if (M <= 0 || N <= 0) NIMAT = 1;

            for (IMAT = 1; IMAT <= NIMAT; IMAT++) { // 100

               // Do the tests only if DOTYPE( IMAT ) is true.

               if( !DOTYPE( IMAT ) ) GO TO 100;

               // Skip types 5, 6, or 7 if the matrix size is too small.

               ZEROT = IMAT >= 5 && IMAT <= 7;
               if (ZEROT && N < IMAT-4) GO TO 100;

               // Set up parameters with DLATB4 and generate a test matrix
               // with DLATMS.

               dlatb4(PATH, IMAT, M, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );

               srnamc.SRNAMT = 'DLATMS';
               dlatms(M, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, 'No packing', A, LDA, WORK, INFO );

               // Check error code from DLATMS.

               if ( INFO != 0 ) {
                  alaerh(PATH, 'DLATMS', INFO, 0, ' ', M, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );
                  GO TO 100;
               }

               // For types 5-7, zero one or more columns of the matrix to
               // test that INFO is returned correctly.

               if ( ZEROT ) {
                  if ( IMAT == 5 ) {
                     IZERO = 1;
                  } else if ( IMAT == 6 ) {
                     IZERO = min( M, N );
                  } else {
                     IZERO = min( M, N ) / 2 + 1;
                  }
                  IOFF = ( IZERO-1 )*LDA;
                  if ( IMAT < 7 ) {
                     for (I = 1; I <= M; I++) { // 20
                        A[IOFF+I] = ZERO;
                     } // 20
                  } else {
                     dlaset('Full', M, N-IZERO+1, ZERO, ZERO, A( IOFF+1 ), LDA );
                  }
               } else {
                  IZERO = 0;
               }

               // These lines, if used in place of the calls in the DO 60
               // loop, cause the code to bomb on a Sun SPARCstation.

                // ANORMO = dlange( 'O', M, N, A, LDA, RWORK )
                // ANORMI = dlange( 'I', M, N, A, LDA, RWORK )

               // Do for each blocksize in NBVAL

               for (INB = 1; INB <= NNB; INB++) { // 90
                  NB = NBVAL( INB );
                  xlaenv(1, NB );

                  // Compute the LU factorization of the matrix.

                  dlacpy('Full', M, N, A, LDA, AFAC, LDA );
                  srnamc.SRNAMT = 'DGETRF';
                  dgetrf(M, N, AFAC, LDA, IWORK, INFO );

                  // Check error code from DGETRF.

                  if (INFO != IZERO) alaerh( PATH, 'DGETRF', INFO, IZERO, ' ', M, N, -1, -1, NB, IMAT, NFAIL, NERRS, NOUT );
                  TRFCON = false;

// +    TEST 1
                  // Reconstruct matrix from factors and compute residual.

                  dlacpy('Full', M, N, AFAC, LDA, AINV, LDA );
                  dget01(M, N, A, LDA, AINV, LDA, IWORK, RWORK, RESULT( 1 ) );
                  NT = 1;

// +    TEST 2
                  // Form the inverse if the factorization was successful
                  // and compute the residual.

                  if ( M == N && INFO == 0 ) {
                     dlacpy('Full', N, N, AFAC, LDA, AINV, LDA );
                     srnamc.SRNAMT = 'DGETRI';
                     NRHS = NSVAL( 1 );
                     LWORK = NMAX*max( 3, NRHS );
                     dgetri(N, AINV, LDA, IWORK, WORK, LWORK, INFO );

                     // Check error code from DGETRI.

                     if (INFO != 0) alaerh( PATH, 'DGETRI', INFO, 0, ' ', N, N, -1, -1, NB, IMAT, NFAIL, NERRS, NOUT );

                     // Compute the residual for the matrix times its
                     // inverse.  Also compute the 1-norm condition number
                     // of A.

                     dget03(N, A, LDA, AINV, LDA, WORK, LDA, RWORK, RCONDO, RESULT( 2 ) );
                     ANORMO = dlange( 'O', M, N, A, LDA, RWORK );

                     // Compute the infinity-norm condition number of A.

                     ANORMI = dlange( 'I', M, N, A, LDA, RWORK );
                     AINVNM = dlange( 'I', N, N, AINV, LDA, RWORK );
                     if ( ANORMI <= ZERO || AINVNM <= ZERO ) {
                        RCONDI = ONE;
                     } else {
                        RCONDI = ( ONE / ANORMI ) / AINVNM;
                     }
                     NT = 2;
                  } else {

                     // Do only the condition estimate if INFO > 0.

                     TRFCON = true;
                     ANORMO = dlange( 'O', M, N, A, LDA, RWORK );
                     ANORMI = dlange( 'I', M, N, A, LDA, RWORK );
                     RCONDO = ZERO;
                     RCONDI = ZERO;
                  }

                  // Print information about the tests so far that did not
                  // pass the threshold.

                  for (K = 1; K <= NT; K++) { // 30
                     if ( RESULT( K ) >= THRESH ) {
                        if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                        WRITE( NOUT, FMT = 9999 )M, N, NB, IMAT, K, RESULT( K );
                        NFAIL = NFAIL + 1;
                     }
                  } // 30
                  NRUN = NRUN + NT;

                  // Skip the remaining tests if this is not the first
                  // block size or if M != N.  Skip the solve tests if
                  // the matrix is singular.

                  if (INB > 1 || M != N) GO TO 90;
                  IF( TRFCON ) GO TO 70;

                  for (IRHS = 1; IRHS <= NNS; IRHS++) { // 60
                     NRHS = NSVAL( IRHS );
                     XTYPE = 'N';

                     for (ITRAN = 1; ITRAN <= NTRAN; ITRAN++) { // 50
                        TRANS = TRANSS( ITRAN );
                        if ( ITRAN == 1 ) {
                           RCONDC = RCONDO;
                        } else {
                           RCONDC = RCONDI;
                        }

// +    TEST 3
                        // Solve and compute residual for A * X = B.

                        srnamc.SRNAMT = 'DLARHS';
                        dlarhs(PATH, XTYPE, ' ', TRANS, N, N, KL, KU, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO );
                        XTYPE = 'C';

                        dlacpy('Full', N, NRHS, B, LDA, X, LDA );
                        srnamc.SRNAMT = 'DGETRS';
                        dgetrs(TRANS, N, NRHS, AFAC, LDA, IWORK, X, LDA, INFO );

                        // Check error code from DGETRS.

                        if (INFO != 0) alaerh( PATH, 'DGETRS', INFO, 0, TRANS, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );

                        dlacpy('Full', N, NRHS, B, LDA, WORK, LDA );
                        dget02(TRANS, N, N, NRHS, A, LDA, X, LDA, WORK, LDA, RWORK, RESULT( 3 ) );

// +    TEST 4
                        // Check solution from generated exact solution.

                        dget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 4 ) );

// +    TESTS 5, 6, and 7
                        // Use iterative refinement to improve the
                        // solution.

                        srnamc.SRNAMT = 'DGERFS';
                        dgerfs(TRANS, N, NRHS, A, LDA, AFAC, LDA, IWORK, B, LDA, X, LDA, RWORK, RWORK( NRHS+1 ), WORK, IWORK( N+1 ), INFO );

                        // Check error code from DGERFS.

                        if (INFO != 0) alaerh( PATH, 'DGERFS', INFO, 0, TRANS, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );

                        dget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 5 ) );
                        dget07(TRANS, N, NRHS, A, LDA, B, LDA, X, LDA, XACT, LDA, RWORK, true , RWORK( NRHS+1 ), RESULT( 6 ) );

                        // Print information about the tests that did not
                        // pass the threshold.

                        for (K = 3; K <= 7; K++) { // 40
                           if ( RESULT( K ) >= THRESH ) {
                              if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                              WRITE( NOUT, FMT = 9998 )TRANS, N, NRHS, IMAT, K, RESULT( K );
                              NFAIL = NFAIL + 1;
                           }
                        } // 40
                        NRUN = NRUN + 5;
                     } // 50
                  } // 60

// +    TEST 8
                     // Get an estimate of RCOND = 1/CNDNUM.

                  } // 70
                  for (ITRAN = 1; ITRAN <= 2; ITRAN++) { // 80
                     if ( ITRAN == 1 ) {
                        ANORM = ANORMO;
                        RCONDC = RCONDO;
                        NORM = 'O';
                     } else {
                        ANORM = ANORMI;
                        RCONDC = RCONDI;
                        NORM = 'I';
                     }
                     srnamc.SRNAMT = 'DGECON';
                     dgecon(NORM, N, AFAC, LDA, ANORM, RCOND, WORK, IWORK( N+1 ), INFO );

                        // Check error code from DGECON.

                     if (INFO != 0) alaerh( PATH, 'DGECON', INFO, 0, NORM, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

                        // This line is needed on a Sun SPARCstation.

                     DUMMY = RCOND;

                     RESULT[8] = DGET06( RCOND, RCONDC );

                     // Print information about the tests that did not pass
                     // the threshold.

                     if ( RESULT( 8 ) >= THRESH ) {
                        if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                        WRITE( NOUT, FMT = 9997 )NORM, N, IMAT, 8, RESULT( 8 );
                        NFAIL = NFAIL + 1;
                     }
                     NRUN = NRUN + 1;
                  } // 80
               } // 90
            } // 100
         } // 110
      } // 120

      // Print a summary of the results.

      alasum(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( ' M = ${.i5}, N =${.i5}, NB =${.i4}, type ${.i2}, test(${.i2}) =${.g12_5};
 9998 FORMAT( ' TRANS=''${.a1}'', N =${.i5}, NRHS=${.i3}, type ${.i2}, test(${.i2}) =${.g12_5};
 9997 FORMAT( ' NORM =''${.a1}'', N =${.i5},${' ' * 10} type ${.i2}, test(${.i2}) =${.g12_5};
      }
