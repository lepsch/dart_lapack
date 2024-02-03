      SUBROUTINE ZCHKTZ( DOTYPE, NM, MVAL, NN, NVAL, THRESH, TSTERR, A, COPYA, S, TAU, WORK, RWORK, NOUT )

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
      double             S( * ), RWORK( * );
      COMPLEX*16         A( * ), COPYA( * ), TAU( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      int                NTYPES;
      PARAMETER          ( NTYPES = 3 )
      int                NTESTS;
      PARAMETER          ( NTESTS = 3 )
      double             ONE, ZERO;
      PARAMETER          ( ONE = 1.0D0, ZERO = 0.0D0 )
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
      double             DLAMCH, ZQRT12, ZRZT01, ZRZT02;
      // EXTERNAL DLAMCH, ZQRT12, ZRZT01, ZRZT02
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAHD, ALASUM, DLAORD, ZERRTZ, ZGEQR2, ZLACPY, ZLASET, ZLATMS, ZTZRZF
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DCMPLX, MAX, MIN
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

      PATH( 1: 1 ) = 'Zomplex precision'
      PATH( 2: 3 ) = 'TZ'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      DO 10 I = 1, 4
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE
      EPS = DLAMCH( 'Epsilon' )

      // Test the error exits

      IF( TSTERR ) CALL ZERRTZ( PATH, NOUT )
      INFOT = 0

      DO 70 IM = 1, NM

         // Do for each value of M in MVAL.

         M = MVAL( IM )
         LDA = MAX( 1, M )

         DO 60 IN = 1, NN

            // Do for each value of N in NVAL for which M .LE. N.

            N = NVAL( IN )
            MNMIN = MIN( M, N )
            LWORK = MAX( 1, N*N+4*M+N )

            IF( M.LE.N ) THEN
               DO 50 IMODE = 1, NTYPES
                  IF( .NOT.DOTYPE( IMODE ) ) GO TO 50

                  // Do for each type of singular value distribution.
                     // 0:  zero matrix
                     // 1:  one small singular value
                     // 2:  exponential distribution

                  MODE = IMODE - 1

                  // Test ZTZRQF

                  // Generate test matrix of size m by n using
                  // singular value distribution indicated by `mode'.

                  IF( MODE.EQ.0 ) THEN
                     CALL ZLASET( 'Full', M, N, DCMPLX( ZERO ), DCMPLX( ZERO ), A, LDA )
                     DO 30 I = 1, MNMIN
                        S( I ) = ZERO
   30                CONTINUE
                  ELSE
                     CALL ZLATMS( M, N, 'Uniform', ISEED, 'Nonsymmetric', S, IMODE, ONE / EPS, ONE, M, N, 'No packing', A, LDA, WORK, INFO )
                     CALL ZGEQR2( M, N, A, LDA, WORK, WORK( MNMIN+1 ), INFO )                      CALL ZLASET( 'Lower', M-1, N, DCMPLX( ZERO ), DCMPLX( ZERO ), A( 2 ), LDA )
                     CALL DLAORD( 'Decreasing', MNMIN, S, 1 )
                  END IF

                  // Save A and its singular values

                  CALL ZLACPY( 'All', M, N, A, LDA, COPYA, LDA )

                  // Call ZTZRZF to reduce the upper trapezoidal matrix to
                  // upper triangular form.

                  SRNAMT = 'ZTZRZF'
                  CALL ZTZRZF( M, N, A, LDA, TAU, WORK, LWORK, INFO )

                  // Compute norm(svd(a) - svd(r))

                  RESULT( 1 ) = ZQRT12( M, M, A, LDA, S, WORK, LWORK, RWORK )

                  // Compute norm( A - R*Q )

                  RESULT( 2 ) = ZRZT01( M, N, COPYA, A, LDA, TAU, WORK, LWORK )

                  // Compute norm(Q'*Q - I).

                  RESULT( 3 ) = ZRZT02( M, N, A, LDA, TAU, WORK, LWORK )

                  // Print information about the tests that did not pass
                 t // he threshold.

                  DO 40 K = 1, NTESTS
                     IF( RESULT( K ).GE.THRESH ) THEN
                        IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                         WRITE( NOUT, FMT = 9999 )M, N, IMODE, K, RESULT( K )
                        NFAIL = NFAIL + 1
                     END IF
   40             CONTINUE
                  NRUN = NRUN + 3
   50          CONTINUE
            END IF
   60    CONTINUE
   70 CONTINUE

      // Print a summary of the results.

      CALL ALASUM( PATH, NOUT, NFAIL, NRUN, NERRS )

 9999 FORMAT( ' M =', I5, ', N =', I5, ', type ', I2, ', test ', I2,
     $      ', ratio =', G12.5 )

      // End if ZCHKTZ

      END
