      void cchktp(final int DOTYPE, final int NN, final int NVAL, final int NNS, final int NSVAL, final int THRESH, final int TSTERR, final int NMAX, final int AP, final int AINVP, final int B, final int X, final int XACT, final Array<double> _WORK, final Array<double> RWORK, final int NOUT) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      bool               TSTERR;
      int                NMAX, NN, NNS, NOUT;
      double               THRESH;
      bool               DOTYPE( * );
      int                NSVAL( * ), NVAL( * );
      double               RWORK( * );
      Complex            AINVP( * ), AP( * ), B( * ), WORK( * ), X( * ), XACT( * );
      // ..

      int                NTYPE1, NTYPES;
      const              NTYPE1 = 10, NTYPES = 18 ;
      int                NTESTS;
      const              NTESTS = 9 ;
      int                NTRAN;
      const              NTRAN = 3 ;
      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      String             DIAG, NORM, TRANS, UPLO, XTYPE;
      String             PATH;
      int                I, IDIAG, IMAT, IN, INFO, IRHS, ITRAN, IUPLO, K, LAP, LDA, N, NERRS, NFAIL, NRHS, NRUN;
      double               AINVNM, ANORM, RCOND, RCONDC, RCONDI, RCONDO, SCALE;
      String             TRANSS( NTRAN ), UPLOS( 2 );
      int                ISEED( 4 ), ISEEDY( 4 );
      double               RESULT( NTESTS );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- REAL               CLANTP;
      // EXTERNAL lsame, CLANTP
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, CCOPY, CERRTR, CGET04, CLACPY, CLARHS, CLATPS, CLATTP, CTPCON, CTPRFS, CTPT01, CTPT02, CTPT03, CTPT05, CTPT06, CTPTRI, CTPTRS
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String            srnamc.SRNAMT;
      int                INFOT, IOUNIT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / INFOT, IOUNIT, OK, LERR
      // COMMON / SRNAMC /srnamc.SRNAMT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Data statements ..
      const ISEEDY = [ 1988, 1989, 1990, 1991 ];
      const UPLOS = 'U', 'L', TRANSS = 'N', 'T', 'C';

      // Initialize constants and the random number seed.

      PATH[1: 1] = 'Complex precision';
      PATH[2: 3] = 'TP';
      NRUN = 0;
      NFAIL = 0;
      NERRS = 0;
      for (I = 1; I <= 4; I++) { // 10
         ISEED[I] = ISEEDY( I );
      } // 10

      // Test the error exits

      if (TSTERR) cerrtr( PATH, NOUT );
      INFOT = 0;

      for (IN = 1; IN <= NN; IN++) { // 110

         // Do for each value of N in NVAL

         N = NVAL( IN );
         LDA = max( 1, N );
         LAP = LDA*( LDA+1 ) / 2;
         XTYPE = 'N';

         for (IMAT = 1; IMAT <= NTYPE1; IMAT++) { // 70

            // Do the tests only if DOTYPE( IMAT ) is true.

            if( !DOTYPE( IMAT ) ) GO TO 70;

            for (IUPLO = 1; IUPLO <= 2; IUPLO++) { // 60

               // Do first for UPLO = 'U', then for UPLO = 'L'

               UPLO = UPLOS( IUPLO );

               // Call CLATTP to generate a triangular test matrix.

              srnamc.SRNAMT = 'CLATTP';
               clattp(IMAT, UPLO, 'No transpose', DIAG, ISEED, N, AP, X, WORK, RWORK, INFO );

               // Set IDIAG = 1 for non-unit matrices, 2 for unit.

               if ( lsame( DIAG, 'N' ) ) {
                  IDIAG = 1;
               } else {
                  IDIAG = 2;
               }

// +    TEST 1
               // Form the inverse of A.

               if (N > 0) ccopy( LAP, AP, 1, AINVP, 1 );
              srnamc.SRNAMT = 'CTPTRI';
               ctptri(UPLO, DIAG, N, AINVP, INFO );

               // Check error code from CTPTRI.

               if (INFO != 0) alaerh( PATH, 'CTPTRI', INFO, 0, UPLO + DIAG, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

               // Compute the infinity-norm condition number of A.

               ANORM = CLANTP( 'I', UPLO, DIAG, N, AP, RWORK );
               AINVNM = CLANTP( 'I', UPLO, DIAG, N, AINVP, RWORK );
               if ( ANORM <= ZERO || AINVNM <= ZERO ) {
                  RCONDI = ONE;
               } else {
                  RCONDI = ( ONE / ANORM ) / AINVNM;
               }

               // Compute the residual for the triangular matrix times its
               // inverse.  Also compute the 1-norm condition number of A.

               ctpt01(UPLO, DIAG, N, AP, AINVP, RCONDO, RWORK, RESULT( 1 ) );

               // Print the test ratio if it is >= THRESH.

               if ( RESULT( 1 ) >= THRESH ) {
                  if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                  WRITE( NOUT, FMT = 9999 )UPLO, DIAG, N, IMAT, 1, RESULT( 1 );
                  NFAIL = NFAIL + 1;
               }
               NRUN = NRUN + 1;

               for (IRHS = 1; IRHS <= NNS; IRHS++) { // 40
                  NRHS = NSVAL( IRHS );
                  XTYPE = 'N';

                  for (ITRAN = 1; ITRAN <= NTRAN; ITRAN++) { // 30

                  // Do for op(A) = A, A**T, or A**H.

                     TRANS = TRANSS( ITRAN );
                     if ( ITRAN == 1 ) {
                        NORM = 'O';
                        RCONDC = RCONDO;
                     } else {
                        NORM = 'I';
                        RCONDC = RCONDI;
                     }

// +    TEST 2
                  // Solve and compute residual for op(A)*x = b.

                    srnamc.SRNAMT = 'CLARHS';
                     clarhs(PATH, XTYPE, UPLO, TRANS, N, N, 0, IDIAG, NRHS, AP, LAP, XACT, LDA, B, LDA, ISEED, INFO );
                     XTYPE = 'C';
                     clacpy('Full', N, NRHS, B, LDA, X, LDA );

                    srnamc.SRNAMT = 'CTPTRS';
                     ctptrs(UPLO, TRANS, DIAG, N, NRHS, AP, X, LDA, INFO );

                  // Check error code from CTPTRS.

                     if (INFO != 0) alaerh( PATH, 'CTPTRS', INFO, 0, UPLO + TRANS // DIAG, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

                     ctpt02(UPLO, TRANS, DIAG, N, NRHS, AP, X, LDA, B, LDA, WORK, RWORK, RESULT( 2 ) );

// +    TEST 3
                  // Check solution from generated exact solution.

                     cget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) );

// +    TESTS 4, 5, and 6
                  // Use iterative refinement to improve the solution and
                  // compute error bounds.

                    srnamc.SRNAMT = 'CTPRFS';
                     ctprfs(UPLO, TRANS, DIAG, N, NRHS, AP, B, LDA, X, LDA, RWORK, RWORK( NRHS+1 ), WORK, RWORK( 2*NRHS+1 ), INFO );

                  // Check error code from CTPRFS.

                     if (INFO != 0) alaerh( PATH, 'CTPRFS', INFO, 0, UPLO + TRANS // DIAG, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );

                     cget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 4 ) );
                     ctpt05(UPLO, TRANS, DIAG, N, NRHS, AP, B, LDA, X, LDA, XACT, LDA, RWORK, RWORK( NRHS+1 ), RESULT( 5 ) );

                     // Print information about the tests that did not pass
                     // the threshold.

                     for (K = 2; K <= 6; K++) { // 20
                        if ( RESULT( K ) >= THRESH ) {
                           if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                           WRITE( NOUT, FMT = 9998 )UPLO, TRANS, DIAG, N, NRHS, IMAT, K, RESULT( K );
                           NFAIL = NFAIL + 1;
                        }
                     } // 20
                     NRUN = NRUN + 5;
                  } // 30
               } // 40

// +    TEST 7
                  // Get an estimate of RCOND = 1/CNDNUM.

               for (ITRAN = 1; ITRAN <= 2; ITRAN++) { // 50
                  if ( ITRAN == 1 ) {
                     NORM = 'O';
                     RCONDC = RCONDO;
                  } else {
                     NORM = 'I';
                     RCONDC = RCONDI;
                  }
                 srnamc.SRNAMT = 'CTPCON';
                  ctpcon(NORM, UPLO, DIAG, N, AP, RCOND, WORK, RWORK, INFO );

                  // Check error code from CTPCON.

                  if (INFO != 0) alaerh( PATH, 'CTPCON', INFO, 0, NORM + UPLO // DIAG, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

                  ctpt06(RCOND, RCONDC, UPLO, DIAG, N, AP, RWORK, RESULT( 7 ) );

                  // Print the test ratio if it is >= THRESH.

                  if ( RESULT( 7 ) >= THRESH ) {
                     if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                     WRITE( NOUT, FMT = 9997 ) 'CTPCON', NORM, UPLO, DIAG, N, IMAT, 7, RESULT( 7 );
                     NFAIL = NFAIL + 1;
                  }
                  NRUN = NRUN + 1;
               } // 50
            } // 60
         } // 70

         // Use pathological test matrices to test CLATPS.

         for (IMAT = NTYPE1 + 1; IMAT <= NTYPES; IMAT++) { // 100

            // Do the tests only if DOTYPE( IMAT ) is true.

            if( !DOTYPE( IMAT ) ) GO TO 100;

            for (IUPLO = 1; IUPLO <= 2; IUPLO++) { // 90

               // Do first for UPLO = 'U', then for UPLO = 'L'

               UPLO = UPLOS( IUPLO );
               for (ITRAN = 1; ITRAN <= NTRAN; ITRAN++) { // 80

                  // Do for op(A) = A, A**T, or A**H.

                  TRANS = TRANSS( ITRAN );

                  // Call CLATTP to generate a triangular test matrix.

                 srnamc.SRNAMT = 'CLATTP';
                  clattp(IMAT, UPLO, TRANS, DIAG, ISEED, N, AP, X, WORK, RWORK, INFO );

// +    TEST 8
                  // Solve the system op(A)*x = b.

                 srnamc.SRNAMT = 'CLATPS';
                  ccopy(N, X, 1, B, 1 );
                  clatps(UPLO, TRANS, DIAG, 'N', N, AP, B, SCALE, RWORK, INFO );

                  // Check error code from CLATPS.

                  if (INFO != 0) alaerh( PATH, 'CLATPS', INFO, 0, UPLO + TRANS // DIAG // 'N', N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

                  ctpt03(UPLO, TRANS, DIAG, N, 1, AP, SCALE, RWORK, ONE, B, LDA, X, LDA, WORK, RESULT( 8 ) );

// +    TEST 9
                  // Solve op(A)*x = b again with NORMIN = 'Y'.

                  ccopy(N, X, 1, B( N+1 ), 1 );
                  clatps(UPLO, TRANS, DIAG, 'Y', N, AP, B( N+1 ), SCALE, RWORK, INFO );

                  // Check error code from CLATPS.

                  if (INFO != 0) alaerh( PATH, 'CLATPS', INFO, 0, UPLO + TRANS // DIAG // 'Y', N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

                  ctpt03(UPLO, TRANS, DIAG, N, 1, AP, SCALE, RWORK, ONE, B( N+1 ), LDA, X, LDA, WORK, RESULT( 9 ) );

                  // Print information about the tests that did not pass
                  // the threshold.

                  if ( RESULT( 8 ) >= THRESH ) {
                     if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                     WRITE( NOUT, FMT = 9996 )'CLATPS', UPLO, TRANS, DIAG, 'N', N, IMAT, 8, RESULT( 8 );
                     NFAIL = NFAIL + 1;
                  }
                  if ( RESULT( 9 ) >= THRESH ) {
                     if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                     WRITE( NOUT, FMT = 9996 )'CLATPS', UPLO, TRANS, DIAG, 'Y', N, IMAT, 9, RESULT( 9 );
                     NFAIL = NFAIL + 1;
                  }
                  NRUN = NRUN + 2;
               } // 80
            } // 90
         } // 100
      } // 110

      // Print a summary of the results.

      alasum(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( ' UPLO=''${.a1}'', DIAG=''${.a1}'', N=${.i5}, type ${.i2}, test(${.i2})= ${.g12_5};
 9998 FORMAT( ' UPLO=''${.a1}'', TRANS=''${.a1}'', DIAG=''${.a1}'', N=${.i5}'', NRHS=${.i5}, type ${.i2}, test(${.i2})= ${.g12_5};
 9997 FORMAT(' ${}( ''${.a1}''${.a1}''${.a1}'',${.i5}, ... ), type ${.i2}, test(${.i2})=${.g12_5};
 9996 FORMAT(' ${}( ''${.a1}''${.a1}''${.a1}''${.a1}'',${.i5}, ... ), type ${.i2}, test(${.i2})=${.g12_5};
      }
