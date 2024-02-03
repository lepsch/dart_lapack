      void sdrvpb(DOTYPE, NN, NVAL, NRHS, THRESH, TSTERR, NMAX, A, AFAC, ASAV, B, BSAV, X, XACT, S, WORK, RWORK, IWORK, NOUT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               TSTERR;
      int                NMAX, NN, NOUT, NRHS;
      REAL               THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                IWORK( * ), NVAL( * );
      REAL               A( * ), AFAC( * ), ASAV( * ), B( * ), BSAV( * ), RWORK( * ), S( * ), WORK( * ), X( * ), XACT( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                NTYPES, NTESTS;
      const              NTYPES = 8, NTESTS = 6 ;
      int                NBW;
      const              NBW = 4 ;
      // ..
      // .. Local Scalars ..
      bool               EQUIL, NOFACT, PREFAC, ZEROT;
      String             DIST, EQUED, FACT, PACKIT, TYPE, UPLO, XTYPE;
      String             PATH;
      int                I, I1, I2, IEQUED, IFACT, IKD, IMAT, IN, INFO, IOFF, IUPLO, IW, IZERO, K, K1, KD, KL, KOFF, KU, LDA, LDAB, MODE, N, NB, NBMIN, NERRS, NFACT, NFAIL, NIMAT, NKD, NRUN, NT;
      REAL               AINVNM, AMAX, ANORM, CNDNUM, RCOND, RCONDC, ROLDC, SCOND;
      // ..
      // .. Local Arrays ..
      String             EQUEDS( 2 ), FACTS( 3 );
      int                ISEED( 4 ), ISEEDY( 4 ), KDVAL( NBW );
      REAL               RESULT( NTESTS );
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SGET06, SLANGE, SLANSB;
      // EXTERNAL LSAME, SGET06, SLANGE, SLANSB
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALADHD, ALAERH, ALASVM, SCOPY, SERRVX, SGET04, SLACPY, SLAQSB, SLARHS, SLASET, SLATB4, SLATMS, SPBEQU, SPBSV, SPBSVX, SPBT01, SPBT02, SPBT05, SPBTRF, SPBTRS, SSWAP, XLAENV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
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
      const ISEEDY = [ 1988, 1989, 1990, 1991 ];
      const FACTS = [ 'F', 'N', 'E' ];
      const EQUEDS = [ 'N', 'Y' ];
      // ..
      // .. Executable Statements ..

      // Initialize constants and the random number seed.

      PATH( 1: 1 ) = 'Single precision';
      PATH( 2: 3 ) = 'PB';
      NRUN = 0;
      NFAIL = 0;
      NERRS = 0;
      for (I = 1; I <= 4; I++) { // 10
         ISEED( I ) = ISEEDY( I );
      } // 10

      // Test the error exits

      if (TSTERR) CALL SERRVX( PATH, NOUT );
      INFOT = 0;
      KDVAL( 1 ) = 0;

      // Set the block size and minimum block size for testing.

      NB = 1;
      NBMIN = 2;
      xlaenv(1, NB );
      xlaenv(2, NBMIN );

      // Do for each value of N in NVAL

      for (IN = 1; IN <= NN; IN++) { // 110
         N = NVAL( IN );
         LDA = max( N, 1 );
         XTYPE = 'N';

         // Set limits on the number of loop iterations.

         NKD = max( 1, min( N, 4 ) );
         NIMAT = NTYPES;
         if (N == 0) NIMAT = 1;

         KDVAL( 2 ) = N + ( N+1 ) / 4;
         KDVAL( 3 ) = ( 3*N-1 ) / 4;
         KDVAL( 4 ) = ( N+1 ) / 4;

         for (IKD = 1; IKD <= NKD; IKD++) { // 100

            // Do for KD = 0, (5*N+1)/4, (3N-1)/4, and (N+1)/4. This order
            // makes it easier to skip redundant values for small values
            // of N.

            KD = KDVAL( IKD );
            LDAB = KD + 1;

            // Do first for UPLO = 'U', then for UPLO = 'L'

            for (IUPLO = 1; IUPLO <= 2; IUPLO++) { // 90
               KOFF = 1;
               if ( IUPLO == 1 ) {
                  UPLO = 'U';
                  PACKIT = 'Q';
                  KOFF = max( 1, KD+2-N );
               } else {
                  UPLO = 'L';
                  PACKIT = 'B';
               }

               for (IMAT = 1; IMAT <= NIMAT; IMAT++) { // 80

                  // Do the tests only if DOTYPE( IMAT ) is true.

                  if( !DOTYPE( IMAT ) ) GO TO 80;

                  // Skip types 2, 3, or 4 if the matrix size is too small.

                  ZEROT = IMAT >= 2 && IMAT <= 4;
                  if (ZEROT && N < IMAT-1) GO TO 80;

                  if ( !ZEROT || !DOTYPE( 1 ) ) {

                     // Set up parameters with SLATB4 and generate a test
                     // matrix with SLATMS.

                     slatb4(PATH, IMAT, N, N, TYPE, KL, KU, ANORM, MODE, CNDNUM, DIST );

                     SRNAMT = 'SLATMS';
                     slatms(N, N, DIST, ISEED, TYPE, RWORK, MODE, CNDNUM, ANORM, KD, KD, PACKIT, A( KOFF ), LDAB, WORK, INFO );

                     // Check error code from SLATMS.

                     if ( INFO != 0 ) {
                        alaerh(PATH, 'SLATMS', INFO, 0, UPLO, N, N, -1, -1, -1, IMAT, NFAIL, NERRS, NOUT );
                        GO TO 80;
                     }
                  } else if ( IZERO > 0 ) {

                     // Use the same matrix for types 3 and 4 as for type
                     // 2 by copying back the zeroed out column,

                     IW = 2*LDA + 1;
                     if ( IUPLO == 1 ) {
                        IOFF = ( IZERO-1 )*LDAB + KD + 1;
                        scopy(IZERO-I1, WORK( IW ), 1, A( IOFF-IZERO+I1 ), 1 );
                        IW = IW + IZERO - I1;
                        scopy(I2-IZERO+1, WORK( IW ), 1, A( IOFF ), max( LDAB-1, 1 ) );
                     } else {
                        IOFF = ( I1-1 )*LDAB + 1;
                        scopy(IZERO-I1, WORK( IW ), 1, A( IOFF+IZERO-I1 ), max( LDAB-1, 1 ) );
                        IOFF = ( IZERO-1 )*LDAB + 1;
                        IW = IW + IZERO - I1;
                        scopy(I2-IZERO+1, WORK( IW ), 1, A( IOFF ), 1 );
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
                     DO 20 I = 1, min( 2*KD+1, N );
                        WORK( IW+I ) = ZERO;
                     } // 20
                     IW = IW + 1;
                     I1 = max( IZERO-KD, 1 );
                     I2 = min( IZERO+KD, N );

                     if ( IUPLO == 1 ) {
                        IOFF = ( IZERO-1 )*LDAB + KD + 1;
                        sswap(IZERO-I1, A( IOFF-IZERO+I1 ), 1, WORK( IW ), 1 );
                        IW = IW + IZERO - I1;
                        sswap(I2-IZERO+1, A( IOFF ), max( LDAB-1, 1 ), WORK( IW ), 1 );
                     } else {
                        IOFF = ( I1-1 )*LDAB + 1;
                        sswap(IZERO-I1, A( IOFF+IZERO-I1 ), max( LDAB-1, 1 ), WORK( IW ), 1 );
                        IOFF = ( IZERO-1 )*LDAB + 1;
                        IW = IW + IZERO - I1;
                        sswap(I2-IZERO+1, A( IOFF ), 1, WORK( IW ), 1 );
                     }
                  }

                  // Save a copy of the matrix A in ASAV.

                  slacpy('Full', KD+1, N, A, LDAB, ASAV, LDAB );

                  for (IEQUED = 1; IEQUED <= 2; IEQUED++) { // 70
                     EQUED = EQUEDS( IEQUED );
                     if ( IEQUED == 1 ) {
                        NFACT = 3;
                     } else {
                        NFACT = 1;
                     }

                     for (IFACT = 1; IFACT <= NFACT; IFACT++) { // 60
                        FACT = FACTS( IFACT );
                        PREFAC = LSAME( FACT, 'F' );
                        NOFACT = LSAME( FACT, 'N' );
                        EQUIL = LSAME( FACT, 'E' );

                        if ( ZEROT ) {
                           if (PREFAC) GO TO 60;
                           RCONDC = ZERO;

                        } else if ( !LSAME( FACT, 'N' ) ) {

                           // Compute the condition number for comparison
                           // with the value returned by SPBSVX (FACT =
                           // 'N' reuses the condition number from the
                           // previous iteration with FACT = 'F').

                           slacpy('Full', KD+1, N, ASAV, LDAB, AFAC, LDAB );
                           if ( EQUIL || IEQUED > 1 ) {

                              // Compute row and column scale factors to
                              // equilibrate the matrix A.

                              spbequ(UPLO, N, KD, AFAC, LDAB, S, SCOND, AMAX, INFO );
                              if ( INFO == 0 && N > 0 ) {
                                 if (IEQUED > 1) SCOND = ZERO;

                                 // Equilibrate the matrix.

                                 slaqsb(UPLO, N, KD, AFAC, LDAB, S, SCOND, AMAX, EQUED );
                              }
                           }

                           // Save the condition number of the
                           // non-equilibrated system for use in SGET04.

                           if (EQUIL) ROLDC = RCONDC;

                           // Compute the 1-norm of A.

                           ANORM = SLANSB( '1', UPLO, N, KD, AFAC, LDAB, RWORK );

                           // Factor the matrix A.

                           spbtrf(UPLO, N, KD, AFAC, LDAB, INFO );

                           // Form the inverse of A.

                           slaset('Full', N, N, ZERO, ONE, A, LDA );
                           SRNAMT = 'SPBTRS';
                           spbtrs(UPLO, N, KD, N, AFAC, LDAB, A, LDA, INFO );

                           // Compute the 1-norm condition number of A.

                           AINVNM = SLANGE( '1', N, N, A, LDA, RWORK );
                           if ( ANORM <= ZERO || AINVNM <= ZERO ) {
                              RCONDC = ONE;
                           } else {
                              RCONDC = ( ONE / ANORM ) / AINVNM;
                           }
                        }

                        // Restore the matrix A.

                        slacpy('Full', KD+1, N, ASAV, LDAB, A, LDAB );

                        // Form an exact solution and set the right hand
                        // side.

                        SRNAMT = 'SLARHS';
                        slarhs(PATH, XTYPE, UPLO, ' ', N, N, KD, KD, NRHS, A, LDAB, XACT, LDA, B, LDA, ISEED, INFO );
                        XTYPE = 'C';
                        slacpy('Full', N, NRHS, B, LDA, BSAV, LDA );

                        if ( NOFACT ) {

                           // --- Test SPBSV  ---

                           // Compute the L*L' or U'*U factorization of the
                           // matrix and solve the system.

                           slacpy('Full', KD+1, N, A, LDAB, AFAC, LDAB );
                           slacpy('Full', N, NRHS, B, LDA, X, LDA );

                           SRNAMT = 'SPBSV ';
                           spbsv(UPLO, N, KD, NRHS, AFAC, LDAB, X, LDA, INFO );

                           // Check error code from SPBSV .

                           if ( INFO != IZERO ) {
                              alaerh(PATH, 'SPBSV ', INFO, IZERO, UPLO, N, N, KD, KD, NRHS, IMAT, NFAIL, NERRS, NOUT );
                              GO TO 40;
                           } else if ( INFO != 0 ) {
                              GO TO 40;
                           }

                           // Reconstruct matrix from factors and compute
                           // residual.

                           spbt01(UPLO, N, KD, A, LDAB, AFAC, LDAB, RWORK, RESULT( 1 ) );

                           // Compute residual of the computed solution.

                           slacpy('Full', N, NRHS, B, LDA, WORK, LDA );
                           spbt02(UPLO, N, KD, NRHS, A, LDAB, X, LDA, WORK, LDA, RWORK, RESULT( 2 ) );

                           // Check solution from generated exact solution.

                           sget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) );
                           NT = 3;

                           // Print information about the tests that did
                           // not pass the threshold.

                           for (K = 1; K <= NT; K++) { // 30
                              if ( RESULT( K ) >= THRESH ) {
                                 if (NFAIL == 0 && NERRS == 0) CALL ALADHD( NOUT, PATH );
                                 WRITE( NOUT, FMT = 9999 )'SPBSV ', UPLO, N, KD, IMAT, K, RESULT( K );
                                 NFAIL = NFAIL + 1;
                              }
                           } // 30
                           NRUN = NRUN + NT;
                           } // 40
                        }

                        // --- Test SPBSVX ---

                        if ( !PREFAC) CALL SLASET( 'Full', KD+1, N, ZERO, ZERO, AFAC, LDAB );
                        slaset('Full', N, NRHS, ZERO, ZERO, X, LDA );
                        if ( IEQUED > 1 && N > 0 ) {

                           // Equilibrate the matrix if FACT='F' and
                           // EQUED='Y'

                           slaqsb(UPLO, N, KD, A, LDAB, S, SCOND, AMAX, EQUED );
                        }

                        // Solve the system and compute the condition
                        // number and error bounds using SPBSVX.

                        SRNAMT = 'SPBSVX';
                        spbsvx(FACT, UPLO, N, KD, NRHS, A, LDAB, AFAC, LDAB, EQUED, S, B, LDA, X, LDA, RCOND, RWORK, RWORK( NRHS+1 ), WORK, IWORK, INFO );

                        // Check the error code from SPBSVX.

                        if ( INFO != IZERO ) {
                           alaerh(PATH, 'SPBSVX', INFO, IZERO, FACT // UPLO, N, N, KD, KD, NRHS, IMAT, NFAIL, NERRS, NOUT );
                           GO TO 60;
                        }

                        if ( INFO == 0 ) {
                           if ( !PREFAC ) {

                              // Reconstruct matrix from factors and
                              // compute residual.

                              spbt01(UPLO, N, KD, A, LDAB, AFAC, LDAB, RWORK( 2*NRHS+1 ), RESULT( 1 ) );
                              K1 = 1;
                           } else {
                              K1 = 2;
                           }

                           // Compute residual of the computed solution.

                           slacpy('Full', N, NRHS, BSAV, LDA, WORK, LDA );
                           spbt02(UPLO, N, KD, NRHS, ASAV, LDAB, X, LDA, WORK, LDA, RWORK( 2*NRHS+1 ), RESULT( 2 ) );

                           // Check solution from generated exact solution.

                           if( NOFACT || ( PREFAC && LSAME( EQUED, 'N' ) ) ) THEN;
                              sget04(N, NRHS, X, LDA, XACT, LDA, RCONDC, RESULT( 3 ) );
                           } else {
                              sget04(N, NRHS, X, LDA, XACT, LDA, ROLDC, RESULT( 3 ) );
                           }

                           // Check the error bounds from iterative
                           // refinement.

                           spbt05(UPLO, N, KD, NRHS, ASAV, LDAB, B, LDA, X, LDA, XACT, LDA, RWORK, RWORK( NRHS+1 ), RESULT( 4 ) );
                        } else {
                           K1 = 6;
                        }

                        // Compare RCOND from SPBSVX with the computed
                        // value in RCONDC.

                        RESULT( 6 ) = SGET06( RCOND, RCONDC );

                        // Print information about the tests that did not
                        // pass the threshold.

                        for (K = K1; K <= 6; K++) { // 50
                           if ( RESULT( K ) >= THRESH ) {
                              if (NFAIL == 0 && NERRS == 0) CALL ALADHD( NOUT, PATH );
                              if ( PREFAC ) {
                                 WRITE( NOUT, FMT = 9997 )'SPBSVX', FACT, UPLO, N, KD, EQUED, IMAT, K, RESULT( K );
                              } else {
                                 WRITE( NOUT, FMT = 9998 )'SPBSVX', FACT, UPLO, N, KD, IMAT, K, RESULT( K );
                              }
                              NFAIL = NFAIL + 1;
                           }
                        } // 50
                        NRUN = NRUN + 7 - K1;
                     } // 60
                  } // 70
               } // 80
            } // 90
         } // 100
      } // 110

      // Print a summary of the results.

      alasvm(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( 1X, A, ', UPLO=''', A1, ''', N =', I5, ', KD =', I5, ', type ', I1, ', test(', I1, ')=', G12.5 );
 9998 FORMAT( 1X, A, '( ''', A1, ''', ''', A1, ''', ', I5, ', ', I5, ', ... ), type ', I1, ', test(', I1, ')=', G12.5 );
 9997 FORMAT( 1X, A, '( ''', A1, ''', ''', A1, ''', ', I5, ', ', I5, ', ... ), EQUED=''', A1, ''', type ', I1, ', test(', I1, ')=', G12.5 );
      return;

      // End of SDRVPB

      }
