import 'common.dart';

      void dchktp(DOTYPE, NN, NVAL, NNS, NSVAL, THRESH, TSTERR, NMAX, AP, AINVP, B, X, XACT, WORK, RWORK, IWORK, NOUT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      bool               TSTERR;
      int                NMAX, NN, NNS, NOUT;
      double             THRESH;
      bool               DOTYPE( * );
      int                IWORK( * ), NSVAL( * ), NVAL( * );
      double             AINVP( * ), AP( * ), B( * ), RWORK( * ), WORK( * ), X( * ), XACT( * );
      // ..

      int                NTYPE1, NTYPES;
      const              NTYPE1 = 10, NTYPES = 18 ;
      int                NTESTS;
      const              NTESTS = 9 ;
      int                NTRAN;
      const              NTRAN = 3 ;
      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      String             DIAG, NORM, TRANS, UPLO, XTYPE;
      String             PATH;
      int                I, IDIAG, IMAT, IN, INFO, IRHS, ITRAN, IUPLO, K, LAP, LDA, N, NERRS, NFAIL, NRHS, NRUN;
      double             AINVNM, ANORM, RCOND, RCONDC, RCONDI, RCONDO, SCALE;
      String             TRANSS( NTRAN ), UPLOS( 2 );
      int                ISEED( 4 ), ISEEDY( 4 );
      double             RESULT( NTESTS );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             DLANTP;
      // EXTERNAL lsame, DLANTP
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASUM, DCOPY, DERRTR, DGET04, DLACPY, DLARHS, DLATPS, DLATTP, DTPCON, DTPRFS, DTPT01, DTPT02, DTPT03, DTPT05, DTPT06, DTPTRI, DTPTRS
      // ..
      // .. Scalars in Common ..
      // bool               infoc.LERR, infoc.OK;
      // String             srnamc.SRNAMT;
      // int                infoc.INFOT, infoc.IOUNIT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / infoc.INFOT, infoc.IOUNIT, infoc.OK, infoc.LERR
      // COMMON / SRNAMC / srnamc.SRNAMT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX
      // ..
      // .. Data statements ..
      const ISEEDY = [ 1988, 1989, 1990, 1991 ];
      const UPLOS = 'U', 'L', TRANSS = 'N', 'T', 'C';

      // Initialize constants and the random number seed.

      PATH = '${'Double precision'[0]}';
      PATH[2: 3] = 'TP';
      NRUN = 0;
      NFAIL = 0;
      NERRS = 0;
      for (I = 1; I <= 4; I++) { // 10
         ISEED[I] = ISEEDY( I );
      } // 10

      // Test the error exits

      if (TSTERR) derrtr( PATH, NOUT );
      infoc.INFOT = 0;

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

               // Call DLATTP to generate a triangular test matrix.

               srnamc.SRNAMT = 'DLATTP';
               dlattp(IMAT, UPLO, 'No transpose', DIAG, ISEED, N, AP, X, WORK, INFO );

               // Set IDIAG = 1 for non-unit matrices, 2 for unit.

               if ( lsame( DIAG, 'N' ) ) {
                  IDIAG = 1;
               } else {
                  IDIAG = 2;
               }

// +    TEST 1
               // Form the inverse of A.

               if (N > 0) dcopy( LAP, AP, 1, AINVP, 1 );
               srnamc.SRNAMT = 'DTPTRI';
               dtptri(UPLO, DIAG, N, AINVP, INFO );

               // Check error code from DTPTRI.

               if (INFO != 0) alaerh( PATH, 'DTPTRI', INFO, 0, UPLO + DIAG, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

               // Compute the infinity-norm condition number of A.

               ANORM = DLANTP( 'I', UPLO, DIAG, N, AP, RWORK );
               AINVNM = DLANTP( 'I', UPLO, DIAG, N, AINVP, RWORK );
               if ( ANORM <= ZERO || AINVNM <= ZERO ) {
                  RCONDI = ONE;
               } else {
                  RCONDI = ( ONE / ANORM ) / AINVNM;
               }

               // Compute the residual for the triangular matrix times its
               // inverse.  Also compute the 1-norm condition number of A.

               dtpt01(UPLO, DIAG, N, AP, AINVP, RCONDO, RWORK, RESULT( 1 ) );

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

                     srnamc.SRNAMT = 'DLARHS';
                     dlarhs(PATH, XTYPE, UPLO, TRANS, N, N, 0, IDIAG, NRHS, AP, LAP, XACT, LDA, B, LDA, ISEED, INFO );
                     XTYPE = 'C';
                     dlacpy('Full', N, NRHS, B, LDA, X, LDA );

                     srnamc.SRNAMT = 'DTPTRS';
                     dtptrs(UPLO, TRANS, DIAG, N, NRHS, AP, X, LDA, INFO );

                  // Check error code from DTPTRS.

                     if (INFO != 0) alaerh( PATH, 'DTPTRS', INFO, 0, UPLO + TRANS // DIAG, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

                     dtpt02(UPLO, TRANS, DIAG, N, NRHS, AP, X, LDA, B, LDA, WORK, RESULT( 2 ) );

// +    TEST 3
                  // Check solution from generated exact solution.

                     dget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) );

// +    TESTS 4, 5, and 6
                  // Use iterative refinement to improve the solution and
                  // compute error bounds.

                     srnamc.SRNAMT = 'DTPRFS';
                     dtprfs(UPLO, TRANS, DIAG, N, NRHS, AP, B, LDA, X, LDA, RWORK, RWORK( NRHS+1 ), WORK, IWORK, INFO );

                  // Check error code from DTPRFS.

                     if (INFO != 0) alaerh( PATH, 'DTPRFS', INFO, 0, UPLO + TRANS // DIAG, N, N, -1, -1, NRHS, IMAT, NFAIL, NERRS, NOUT );

                     dget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 4 ) );
                     dtpt05(UPLO, TRANS, DIAG, N, NRHS, AP, B, LDA, X, LDA, XACT, LDA, RWORK, RWORK( NRHS+1 ), RESULT( 5 ) );

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

                  srnamc.SRNAMT = 'DTPCON';
                  dtpcon(NORM, UPLO, DIAG, N, AP, RCOND, WORK, IWORK, INFO );

                  // Check error code from DTPCON.

                  if (INFO != 0) alaerh( PATH, 'DTPCON', INFO, 0, NORM + UPLO // DIAG, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

                  dtpt06(RCOND, RCONDC, UPLO, DIAG, N, AP, RWORK, RESULT( 7 ) );

                  // Print the test ratio if it is >= THRESH.

                  if ( RESULT( 7 ) >= THRESH ) {
                     if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                     WRITE( NOUT, FMT = 9997 ) 'DTPCON', NORM, UPLO, DIAG, N, IMAT, 7, RESULT( 7 );
                     NFAIL = NFAIL + 1;
                  }
                  NRUN = NRUN + 1;
               } // 50
            } // 60
         } // 70

         // Use pathological test matrices to test DLATPS.

         for (IMAT = NTYPE1 + 1; IMAT <= NTYPES; IMAT++) { // 100

            // Do the tests only if DOTYPE( IMAT ) is true.

            if( !DOTYPE( IMAT ) ) GO TO 100;

            for (IUPLO = 1; IUPLO <= 2; IUPLO++) { // 90

               // Do first for UPLO = 'U', then for UPLO = 'L'

               UPLO = UPLOS( IUPLO );
               for (ITRAN = 1; ITRAN <= NTRAN; ITRAN++) { // 80

                  // Do for op(A) = A, A**T, or A**H.

                  TRANS = TRANSS( ITRAN );

                  // Call DLATTP to generate a triangular test matrix.

                  srnamc.SRNAMT = 'DLATTP';
                  dlattp(IMAT, UPLO, TRANS, DIAG, ISEED, N, AP, X, WORK, INFO );

// +    TEST 8
                  // Solve the system op(A)*x = b.

                  srnamc.SRNAMT = 'DLATPS';
                  dcopy(N, X, 1, B, 1 );
                  dlatps(UPLO, TRANS, DIAG, 'N', N, AP, B, SCALE, RWORK, INFO );

                  // Check error code from DLATPS.

                  if (INFO != 0) alaerh( PATH, 'DLATPS', INFO, 0, UPLO + TRANS // DIAG // 'N', N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

                  dtpt03(UPLO, TRANS, DIAG, N, 1, AP, SCALE, RWORK, ONE, B, LDA, X, LDA, WORK, RESULT( 8 ) );

// +    TEST 9
                  // Solve op(A)*x = b again with NORMIN = 'Y'.

                  dcopy(N, X, 1, B( N+1 ), 1 );
                  dlatps(UPLO, TRANS, DIAG, 'Y', N, AP, B( N+1 ), SCALE, RWORK, INFO );

                  // Check error code from DLATPS.

                  if (INFO != 0) alaerh( PATH, 'DLATPS', INFO, 0, UPLO + TRANS // DIAG // 'Y', N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );

                  dtpt03(UPLO, TRANS, DIAG, N, 1, AP, SCALE, RWORK, ONE, B( N+1 ), LDA, X, LDA, WORK, RESULT( 9 ) );

                  // Print information about the tests that did not pass
                  // the threshold.

                  if ( RESULT( 8 ) >= THRESH ) {
                     if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                     WRITE( NOUT, FMT = 9996 )'DLATPS', UPLO, TRANS, DIAG, 'N', N, IMAT, 8, RESULT( 8 );
                     NFAIL = NFAIL + 1;
                  }
                  if ( RESULT( 9 ) >= THRESH ) {
                     if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                     WRITE( NOUT, FMT = 9996 )'DLATPS', UPLO, TRANS, DIAG, 'Y', N, IMAT, 9, RESULT( 9 );
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
 9997 FORMAT( 1X, A, '( ''${.a1}''${.a1}''${.a1}'',${.i5}, ... ), type ${.i2}, test(${.i2})=${.g12_5};
 9996 FORMAT( 1X, A, '( ''${.a1}''${.a1}''${.a1}''${.a1}'',${.i5}, ... ), type ${.i2}, test(${.i2})=${.g12_5};
      }
