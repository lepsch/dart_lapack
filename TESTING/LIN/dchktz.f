      SUBROUTINE DCHKTZ( DOTYPE, NM, MVAL, NN, NVAL, THRESH, TSTERR, A, COPYA, S, TAU, WORK, NOUT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               TSTERR;
      int                NM, NN, NOUT;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                MVAL( * ), NVAL( * );
      double             A( * ), COPYA( * ), S( * ), TAU( * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      int                NTYPES;
      const              NTYPES = 3 ;
      int                NTESTS;
      const              NTESTS = 3 ;
      double             ONE, ZERO;
      const              ONE = 1.0D0, ZERO = 0.0D0 ;
      // ..
      // .. Local Scalars ..
      String             PATH;
      int                I, IM, IMODE, IN, INFO, K, LDA, LWORK, M, MNMIN, MODE, N, NERRS, NFAIL, NRUN;
      double             EPS;
      // ..
      // .. Local Arrays ..
      int                ISEED( 4 ), ISEEDY( 4 );
      double             RESULT( NTESTS );
      // ..
      // .. External Functions ..
      double             DLAMCH, DQRT12, DRZT01, DRZT02;
      // EXTERNAL DLAMCH, DQRT12, DRZT01, DRZT02
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAHD, ALASUM, DERRTZ, DGEQR2, DLACPY, DLAORD, DLASET, DLATMS, DTZRZF
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Scalars in Common ..
      bool               LERR, OK;
      String             SRNAMT;
      int                INFOT, IOUNIT;
      // ..
      // .. Common blocks ..
      COMMON             / INFOC / INFOT, IOUNIT, OK, LERR
      COMMON             / SRNAMC / SRNAMT
      // ..
      // .. Data statements ..
      DATA               ISEEDY / 1988, 1989, 1990, 1991 /
      // ..
      // .. Executable Statements ..

      // Initialize constants and the random number seed.

      PATH( 1: 1 ) = 'double          ';
      PATH( 2: 3 ) = 'TZ'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      for (I = 1; I <= 4; I++) { // 10
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE
      EPS = DLAMCH( 'Epsilon' )

      // Test the error exits

      IF( TSTERR ) CALL DERRTZ( PATH, NOUT )
      INFOT = 0

      for (IM = 1; IM <= NM; IM++) { // 70

         // Do for each value of M in MVAL.

         M = MVAL( IM )
         LDA = MAX( 1, M )

         for (IN = 1; IN <= NN; IN++) { // 60

            // Do for each value of N in NVAL for which M .LE. N.

            N = NVAL( IN )
            MNMIN = MIN( M, N )
            LWORK = MAX( 1, N*N+4*M+N, M*N+2*MNMIN+4*N )

            if ( M.LE.N ) {
               for (IMODE = 1; IMODE <= NTYPES; IMODE++) { // 50
                  IF( .NOT.DOTYPE( IMODE ) ) GO TO 50

                  // Do for each type of singular value distribution.
                     // 0:  zero matrix
                     // 1:  one small singular value
                     // 2:  exponential distribution

                  MODE = IMODE - 1

                  // Test DTZRQF

                  // Generate test matrix of size m by n using
                  // singular value distribution indicated by `mode'.

                  if ( MODE.EQ.0 ) {
                     dlaset('Full', M, N, ZERO, ZERO, A, LDA );
                     for (I = 1; I <= MNMIN; I++) { // 30
                        S( I ) = ZERO
   30                CONTINUE
                  } else {
                     dlatms(M, N, 'Uniform', ISEED, 'Nonsymmetric', S, IMODE, ONE / EPS, ONE, M, N, 'No packing', A, LDA, WORK, INFO );
                     dgeqr2(M, N, A, LDA, WORK, WORK( MNMIN+1 ), INFO )                      CALL DLASET( 'Lower', M-1, N, ZERO, ZERO, A( 2 ), LDA );
                     dlaord('Decreasing', MNMIN, S, 1 );
                  }

                  // Save A and its singular values

                  dlacpy('All', M, N, A, LDA, COPYA, LDA );

                  // Call DTZRZF to reduce the upper trapezoidal matrix to
                  // upper triangular form.

                  SRNAMT = 'DTZRZF'
                  dtzrzf(M, N, A, LDA, TAU, WORK, LWORK, INFO );

                  // Compute norm(svd(a) - svd(r))

                  RESULT( 1 ) = DQRT12( M, M, A, LDA, S, WORK, LWORK )

                  // Compute norm( A - R*Q )

                  RESULT( 2 ) = DRZT01( M, N, COPYA, A, LDA, TAU, WORK, LWORK )

                  // Compute norm(Q'*Q - I).

                  RESULT( 3 ) = DRZT02( M, N, A, LDA, TAU, WORK, LWORK )

                  // Print information about the tests that did not pass
                  // the threshold.

                  for (K = 1; K <= NTESTS; K++) { // 40
                     if ( RESULT( K ).GE.THRESH ) {
                        IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                         WRITE( NOUT, FMT = 9999 )M, N, IMODE, K, RESULT( K )
                        NFAIL = NFAIL + 1
                     }
   40             CONTINUE
                  NRUN = NRUN + 3
   50          CONTINUE
            }
   60    CONTINUE
   70 CONTINUE

      // Print a summary of the results.

      alasum(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( ' M =', I5, ', N =', I5, ', type ', I2, ', test ', I2, ', ratio =', G12.5 )

      // End if DCHKTZ

      }
