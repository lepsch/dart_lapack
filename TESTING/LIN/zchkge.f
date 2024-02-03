      SUBROUTINE ZCHKGE( DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NNS, NSVAL, THRESH, TSTERR, NMAX, A, AFAC, AINV, B, X, XACT, WORK, RWORK, IWORK, NOUT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               TSTERR;
      int                NM, NMAX, NN, NNB, NNS, NOUT;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                IWORK( * ), MVAL( * ), NBVAL( * ), NSVAL( * ), NVAL( * );
      double             RWORK( * );
      COMPLEX*16         A( * ), AFAC( * ), AINV( * ), B( * ), WORK( * ), X( * ), XACT( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      int                NTYPES;
      const              NTYPES = 11 ;
      int                NTESTS;
      const              NTESTS = 8 ;
      int                NTRAN;
      const              NTRAN = 3 ;
      // ..
      // .. Local Scalars ..
      bool               TRFCON, ZEROT;
      String             DIST, NORM, TRANS, TYPE, XTYPE;
      String             PATH;
      int                I, IM, IMAT, IN, INB, INFO, IOFF, IRHS, ITRAN, IZERO, K, KL, KU, LDA, LWORK, M, MODE, N, NB, NERRS, NFAIL, NIMAT, NRHS, NRUN, NT;
      double             AINVNM, ANORM, ANORMI, ANORMO, CNDNUM, DUMMY, RCOND, RCONDC, RCONDI, RCONDO;
      // ..
      // .. Local Arrays ..
      String             TRANSS( NTRAN );
      int                ISEED( 4 ), ISEEDY( 4 );
      double             RESULT( NTESTS );
      // ..
      // .. External Functions ..
      double             DGET06, ZLANGE;
      // EXTERNAL DGET06, ZLANGE
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, XLAENV, ZERRGE, ZGECON, ZGERFS, ZGET01, ZGET02, ZGET03, ZGET04, ZGET07, ZGETRF, ZGETRI, ZGETRS, ZLACPY, ZLARHS, ZLASET, ZLATB4, ZLATMS
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DCMPLX, MAX, MIN
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String             SRNAMT;
      int                INFOT, NUNIT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / INFOT, NUNIT, OK, LERR
      // COMMON / SRNAMC / SRNAMT
      // ..
      // .. Data statements ..
      DATA               ISEEDY / 1988, 1989, 1990, 1991 / , TRANSS / 'N', 'T', 'C' /
      // ..
      // .. Executable Statements ..

      // Initialize constants and the random number seed.

      PATH( 1: 1 ) = 'Zomplex precision'
      PATH( 2: 3 ) = 'GE'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      for (I = 1; I <= 4; I++) { // 10
         ISEED( I ) = ISEEDY( I )
      } // 10

      // Test the error exits

      xlaenv(1, 1 );
      if (TSTERR) CALL ZERRGE( PATH, NOUT );
      INFOT = 0
      xlaenv(2, 2 );

      // Do for each value of M in MVAL

      for (IM = 1; IM <= NM; IM++) { // 120
         M = MVAL( IM )
         LDA = MAX( 1, M )

         // Do for each value of N in NVAL

         for (IN = 1; IN <= NN; IN++) { // 110
            N = NVAL( IN )
            XTYPE = 'N'
            NIMAT = NTYPES
            if (M.LE.0 || N.LE.0) NIMAT = 1;

            for (IMAT = 1; IMAT <= NIMAT; IMAT++) { // 100

               // Do the tests only if DOTYPE( IMAT ) is true.

               IF( .NOT.DOTYPE( IMAT ) ) GO TO 100

               // Skip types 5, 6, or 7 if the matrix size is too small.

               ZEROT = IMAT.GE.5 && IMAT.LE.7
               if (ZEROT && N < IMAT-4) GO TO 100;

               // Set up parameters with ZLATB4 and generate a test matrix
               // with ZLATMS.

               zlatb4(PATH, IMAT, M, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );

               SRNAMT = 'ZLATMS'
               zlatms(M, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KL, KU, 'No packing', A, LDA, WORK, INFO );

               // Check error code from ZLATMS.

               if ( INFO != 0 ) {
                  alaerh(PATH, 'ZLATMS', INFO, 0, ' ', M, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );
                  GO TO 100
               }

               // For types 5-7, zero one or more columns of the matrix to
               // test that INFO is returned correctly.

               if ( ZEROT ) {
                  if ( IMAT == 5 ) {
                     IZERO = 1
                  } else if ( IMAT == 6 ) {
                     IZERO = MIN( M, N )
                  } else {
                     IZERO = MIN( M, N ) / 2 + 1
                  }
                  IOFF = ( IZERO-1 )*LDA
                  if ( IMAT < 7 ) {
                     for (I = 1; I <= M; I++) { // 20
                        A( IOFF+I ) = ZERO
                     } // 20
                  } else {
                     zlaset('Full', M, N-IZERO+1, DCMPLX( ZERO ), DCMPLX( ZERO ), A( IOFF+1 ), LDA );
                  }
               } else {
                  IZERO = 0
               }

               // These lines, if used in place of the calls in the DO 60
               // loop, cause the code to bomb on a Sun SPARCstation.

                // ANORMO = ZLANGE( 'O', M, N, A, LDA, RWORK )
                // ANORMI = ZLANGE( 'I', M, N, A, LDA, RWORK )

               // Do for each blocksize in NBVAL

               for (INB = 1; INB <= NNB; INB++) { // 90
                  NB = NBVAL( INB )
                  xlaenv(1, NB );

                  // Compute the LU factorization of the matrix.

                  zlacpy('Full', M, N, A, LDA, AFAC, LDA );
                  SRNAMT = 'ZGETRF'
                  zgetrf(M, N, AFAC, LDA, IWORK, INFO );

                  // Check error code from ZGETRF.

                  if (INFO != IZERO) CALL ALAERH( PATH, 'ZGETRF', INFO, IZERO, ' ', M, N, -1, -1, NB, IMAT, NFAIL, NERRS, NOUT );
                  TRFCON = false;

*+    TEST 1
                  // Reconstruct matrix from factors and compute residual.

                  zlacpy('Full', M, N, AFAC, LDA, AINV, LDA );
                  zget01(M, N, A, LDA, AINV, LDA, IWORK, RWORK, RESULT( 1 ) );
                  NT = 1

*+    TEST 2
                  // Form the inverse if the factorization was successful
                  // and compute the residual.

                  if ( M == N && INFO == 0 ) {
                     zlacpy('Full', N, N, AFAC, LDA, AINV, LDA );
                     SRNAMT = 'ZGETRI'
                     NRHS = NSVAL( 1 )
                     LWORK = NMAX*MAX( 3, NRHS )
                     zgetri(N, AINV, LDA, IWORK, WORK, LWORK, INFO );

                     // Check error code from ZGETRI.

                     if (INFO != 0) CALL ALAERH( PATH, 'ZGETRI', INFO, 0, ' ', N, N, -1, -1, NB, IMAT, NFAIL, NERRS, NOUT );

                     // Compute the residual for the matrix times its
                     // inverse.  Also compute the 1-norm condition number
                     // of A.

                     zget03(N, A, LDA, AINV, LDA, WORK, LDA, RWORK, RCONDO, RESULT( 2 ) );
                     ANORMO = ZLANGE( 'O', M, N, A, LDA, RWORK )

                     // Compute the infinity-norm condition number of A.

                     ANORMI = ZLANGE( 'I', M, N, A, LDA, RWORK )
                     AINVNM = ZLANGE( 'I', N, N, AINV, LDA, RWORK )
                     if ( ANORMI.LE.ZERO || AINVNM.LE.ZERO ) {
                        RCONDI = ONE
                     } else {
                        RCONDI = ( ONE / ANORMI ) / AINVNM
                     }
                     NT = 2
                  } else {

                     // Do only the condition estimate if INFO > 0.

                     TRFCON = true;
                     ANORMO = ZLANGE( 'O', M, N, A, LDA, RWORK )
                     ANORMI = ZLANGE( 'I', M, N, A, LDA, RWORK )
                     RCONDO = ZERO
                     RCONDI = ZERO
                  }

                  // Print information about the tests so far that did not
                  // pass the threshold.

                  for (K = 1; K <= NT; K++) { // 30
                     if ( RESULT( K ).GE.THRESH ) {
                        if (NFAIL == 0 && NERRS == 0) CALL ALAHD( NOUT, PATH )                         WRITE( NOUT, FMT = 9999 )M, N, NB, IMAT, K, RESULT( K );
                        NFAIL = NFAIL + 1
                     }
                  } // 30
                  NRUN = NRUN + NT

                  // Skip the remaining tests if this is not the first
                  // block size or if M != N.  Skip the solve tests if
                  // the matrix is singular.

                  if (INB.GT.1 || M != N) GO TO 90                   IF( TRFCON ) GO TO 70;

                  for (IRHS = 1; IRHS <= NNS; IRHS++) { // 60
                     NRHS = NSVAL( IRHS )
                     XTYPE = 'N'

                     for (ITRAN = 1; ITRAN <= NTRAN; ITRAN++) { // 50
                        TRANS = TRANSS( ITRAN )
                        if ( ITRAN == 1 ) {
                           RCONDC = RCONDO
                        } else {
                           RCONDC = RCONDI
                        }

*+    TEST 3
                        // Solve and compute residual for A * X = B.

                        SRNAMT = 'ZLARHS'
                        zlarhs(PATH, XTYPE, ' ', TRANS, N, N, KL, KU, NRHS, A, LDA, XACT, LDA, B, LDA, ISEED, INFO );
                        XTYPE = 'C'

                        zlacpy('Full', N, NRHS, B, LDA, X, LDA );
                        SRNAMT = 'ZGETRS'
                        zgetrs(TRANS, N, NRHS, AFAC, LDA, IWORK, X, LDA, INFO );

                        // Check error code from ZGETRS.

                        if (INFO != 0) CALL ALAERH( PATH, 'ZGETRS', INFO, 0, TRANS, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );

                        zlacpy('Full', N, NRHS, B, LDA, WORK, LDA );
                        zget02(TRANS, N, N, NRHS, A, LDA, X, LDA, WORK, LDA, RWORK, RESULT( 3 ) );

*+    TEST 4
                        // Check solution from generated exact solution.

                        zget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 4 ) );

*+    TESTS 5, 6, and 7
                        // Use iterative refinement to improve the
                        // solution.

                        SRNAMT = 'ZGERFS'
                        zgerfs(TRANS, N, NRHS, A, LDA, AFAC, LDA, IWORK, B, LDA, X, LDA, RWORK, RWORK( NRHS+1 ), WORK, RWORK( 2*NRHS+1 ), INFO );

                        // Check error code from ZGERFS.

                        if (INFO != 0) CALL ALAERH( PATH, 'ZGERFS', INFO, 0, TRANS, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );

                        zget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 5 ) );
                        zget07(TRANS, N, NRHS, A, LDA, B, LDA, X, LDA, XACT, LDA, RWORK, true , RWORK( NRHS+1 ), RESULT( 6 ) );

                        // Print information about the tests that did not
                        // pass the threshold.

                        for (K = 3; K <= 7; K++) { // 40
                           if ( RESULT( K ).GE.THRESH ) {
                              if (NFAIL == 0 && NERRS == 0) CALL ALAHD( NOUT, PATH )                               WRITE( NOUT, FMT = 9998 )TRANS, N, NRHS, IMAT, K, RESULT( K );
                              NFAIL = NFAIL + 1
                           }
                        } // 40
                        NRUN = NRUN + 5
                     } // 50
                  } // 60

*+    TEST 8
                     // Get an estimate of RCOND = 1/CNDNUM.

                  } // 70
                  for (ITRAN = 1; ITRAN <= 2; ITRAN++) { // 80
                     if ( ITRAN == 1 ) {
                        ANORM = ANORMO
                        RCONDC = RCONDO
                        NORM = 'O'
                     } else {
                        ANORM = ANORMI
                        RCONDC = RCONDI
                        NORM = 'I'
                     }
                     SRNAMT = 'ZGECON'
                     zgecon(NORM, N, AFAC, LDA, ANORM, RCOND, WORK, RWORK, INFO );

                        // Check error code from ZGECON.

                     if (INFO != 0) CALL ALAERH( PATH, 'ZGECON', INFO, 0, NORM, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

                        // This line is needed on a Sun SPARCstation.

                     DUMMY = RCOND

                     RESULT( 8 ) = DGET06( RCOND, RCONDC )

                     // Print information about the tests that did not pass
                     // the threshold.

                     if ( RESULT( 8 ).GE.THRESH ) {
                        if (NFAIL == 0 && NERRS == 0) CALL ALAHD( NOUT, PATH )                         WRITE( NOUT, FMT = 9997 )NORM, N, IMAT, 8, RESULT( 8 );
                        NFAIL = NFAIL + 1
                     }
                     NRUN = NRUN + 1
                  } // 80
               } // 90
            } // 100

         } // 110
      } // 120

      // Print a summary of the results.

      alasum(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( ' M = ', I5, ', N =', I5, ', NB =', I4, ', type ', I2, ', test(', I2, ') =', G12.5 )
 9998 FORMAT( ' TRANS=''', A1, ''', N =', I5, ', NRHS=', I3, ', type ', I2, ', test(', I2, ') =', G12.5 )
 9997 FORMAT( ' NORM =''', A1, ''', N =', I5, ',', 10X, ' type ', I2, ', test(', I2, ') =', G12.5 )
      RETURN

      // End of ZCHKGE

      }
