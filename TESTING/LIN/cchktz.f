      void cchktz(DOTYPE, NM, MVAL, NN, NVAL, THRESH, TSTERR, A, COPYA, S, TAU, WORK, RWORK, NOUT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               TSTERR;
      int                NM, NN, NOUT;
      REAL               THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                MVAL( * ), NVAL( * );
      REAL               S( * ), RWORK( * );
      COMPLEX            A( * ), COPYA( * ), TAU( * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      int                NTYPES;
      const              NTYPES = 3 ;
      int                NTESTS;
      const              NTESTS = 3 ;
      REAL               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      String             PATH;
      int                I, IM, IMODE, IN, INFO, K, LDA, LWORK, M, MNMIN, MODE, N, NERRS, NFAIL, NRUN;
      REAL               EPS;
      // ..
      // .. Local Arrays ..
      int                ISEED( 4 ), ISEEDY( 4 );
      REAL               RESULT( NTESTS );
      // ..
      // .. External Functions ..
      REAL               CQRT12, CRZT01, CRZT02, SLAMCH;
      // EXTERNAL CQRT12, CRZT01, CRZT02, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAHD, ALASUM, CERRTZ, CGEQR2, CLACPY, CLASET, CLATMS, CTZRZF, SLAORD
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, MAX, MIN
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String             SRNAMT;
      int                INFOT, IOUNIT;
      // ..
      // .. Common blocks ..
      // COMMON / INFOC / INFOT, IOUNIT, OK, LERR
      // COMMON / SRNAMC / SRNAMT
      // ..
      // .. Data statements ..
      const ISEEDY = [ 1988, 1989, 1990, 1991 ];
      // ..
      // .. Executable Statements ..

      // Initialize constants and the random number seed.

      PATH( 1: 1 ) = 'Complex precision';
      PATH( 2: 3 ) = 'TZ';
      NRUN = 0;
      NFAIL = 0;
      NERRS = 0;
      for (I = 1; I <= 4; I++) { // 10
         ISEED( I ) = ISEEDY( I );
      } // 10
      EPS = SLAMCH( 'Epsilon' );

      // Test the error exits

      if (TSTERR) CALL CERRTZ( PATH, NOUT );
      INFOT = 0;

      for (IM = 1; IM <= NM; IM++) { // 70

         // Do for each value of M in MVAL.

         M = MVAL( IM );
         LDA = max( 1, M );

         for (IN = 1; IN <= NN; IN++) { // 60

            // Do for each value of N in NVAL for which M <= N.

            N = NVAL( IN );
            MNMIN = min( M, N );
            LWORK = max( 1, N*N+4*M+N );

            if ( M <= N ) {
               for (IMODE = 1; IMODE <= NTYPES; IMODE++) { // 50
                  if( !DOTYPE( IMODE ) ) GO TO 50;

                  // Do for each type of singular value distribution.
                     // 0:  zero matrix
                     // 1:  one small singular value
                     // 2:  exponential distribution

                  MODE = IMODE - 1;

                  // Test CTZRZF

                  // Generate test matrix of size m by n using
                  // singular value distribution indicated by `mode'.

                  if ( MODE == 0 ) {
                     claset('Full', M, N, CMPLX( ZERO ), CMPLX( ZERO ), A, LDA );
                     for (I = 1; I <= MNMIN; I++) { // 30
                        S( I ) = ZERO;
                     } // 30
                  } else {
                     clatms(M, N, 'Uniform', ISEED, 'Nonsymmetric', S, IMODE, ONE / EPS, ONE, M, N, 'No packing', A, LDA, WORK, INFO );
                     cgeqr2(M, N, A, LDA, WORK, WORK( MNMIN+1 ), INFO );
                     claset('Lower', M-1, N, CMPLX( ZERO ), CMPLX( ZERO ), A( 2 ), LDA );
                     slaord('Decreasing', MNMIN, S, 1 );
                  }

                  // Save A and its singular values

                  clacpy('All', M, N, A, LDA, COPYA, LDA );

                  // Call CTZRZF to reduce the upper trapezoidal matrix to
                  // upper triangular form.

                  SRNAMT = 'CTZRZF';
                  ctzrzf(M, N, A, LDA, TAU, WORK, LWORK, INFO );

                  // Compute norm(svd(a) - svd(r))

                  RESULT( 1 ) = CQRT12( M, M, A, LDA, S, WORK, LWORK, RWORK );

                  // Compute norm( A - R*Q )

                  RESULT( 2 ) = CRZT01( M, N, COPYA, A, LDA, TAU, WORK, LWORK );

                  // Compute norm(Q'*Q - I).

                  RESULT( 3 ) = CRZT02( M, N, A, LDA, TAU, WORK, LWORK );

                  // Print information about the tests that did not pass
                  // the threshold.

                  for (K = 1; K <= NTESTS; K++) { // 40
                     if ( RESULT( K ) >= THRESH ) {
                        if (NFAIL == 0 && NERRS == 0) CALL ALAHD( NOUT, PATH );
                        WRITE( NOUT, FMT = 9999 )M, N, IMODE, K, RESULT( K );
                        NFAIL = NFAIL + 1;
                     }
                  } // 40
                  NRUN = NRUN + 3;
               } // 50
            }
         } // 60
      } // 70

      // Print a summary of the results.

      alasum(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( ' M =', I5, ', N =', I5, ', type ', I2, ', test ', I2, ', ratio =', G12.5 );

      // End if CCHKTZ

      }
