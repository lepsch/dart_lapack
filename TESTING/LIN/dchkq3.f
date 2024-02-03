      SUBROUTINE DCHKQ3( DOTYPE, NM, MVAL, NN, NVAL, NNB, NBVAL, NXVAL, THRESH, A, COPYA, S, TAU, WORK, IWORK, NOUT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                NM, NN, NNB, NOUT;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                IWORK( * ), MVAL( * ), NBVAL( * ), NVAL( * ), NXVAL( * )       double             A( * ), COPYA( * ), S( * ), TAU( * ), WORK( * );;
      // ..

*  =====================================================================

      // .. Parameters ..
      int                NTYPES;
      PARAMETER          ( NTYPES = 6 )
      int                NTESTS;
      PARAMETER          ( NTESTS = 3 )
      double             ONE, ZERO;
      PARAMETER          ( ONE = 1.0D0, ZERO = 0.0D0 )
      // ..
      // .. Local Scalars ..
      String             PATH;
      int                I, IHIGH, ILOW, IM, IMODE, IN, INB, INFO, ISTEP, K, LDA, LW, LWORK, M, MNMIN, MODE, N, NB, NERRS, NFAIL, NRUN, NX;
      double             EPS;
      // ..
      // .. Local Arrays ..
      int                ISEED( 4 ), ISEEDY( 4 );
      double             RESULT( NTESTS );
      // ..
      // .. External Functions ..
      double             DLAMCH, DQPT01, DQRT11, DQRT12;
      // EXTERNAL DLAMCH, DQPT01, DQRT11, DQRT12
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAHD, ALASUM, DGEQP3, DLACPY, DLAORD, DLASET, DLATMS, ICOPY, XLAENV
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
      PATH( 2: 3 ) = 'Q3'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      DO 10 I = 1, 4
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE
      EPS = DLAMCH( 'Epsilon' )
      INFOT = 0

      DO 90 IM = 1, NM

         // Do for each value of M in MVAL.

         M = MVAL( IM )
         LDA = MAX( 1, M )

         DO 80 IN = 1, NN

            // Do for each value of N in NVAL.

            N = NVAL( IN )
            MNMIN = MIN( M, N )
            LWORK = MAX( 1, M*MAX( M, N )+4*MNMIN+MAX( M, N ), M*N + 2*MNMIN + 4*N )

            DO 70 IMODE = 1, NTYPES
               IF( .NOT.DOTYPE( IMODE ) ) GO TO 70

               // Do for each type of matrix
                  // 1:  zero matrix
                  // 2:  one small singular value
                  // 3:  geometric distribution of singular values
                  // 4:  first n/2 columns fixed
                  // 5:  last n/2 columns fixed
                  // 6:  every second column fixed

               MODE = IMODE
               IF( IMODE.GT.3 ) MODE = 1

               // Generate test matrix of size m by n using
               // singular value distribution indicated by `mode'.

               DO 20 I = 1, N
                  IWORK( I ) = 0
   20          CONTINUE
               IF( IMODE.EQ.1 ) THEN
                  CALL DLASET( 'Full', M, N, ZERO, ZERO, COPYA, LDA )
                  DO 30 I = 1, MNMIN
                     S( I ) = ZERO
   30             CONTINUE
               ELSE
                  CALL DLATMS( M, N, 'Uniform', ISEED, 'Nonsymm', S, MODE, ONE / EPS, ONE, M, N, 'No packing', COPYA, LDA, WORK, INFO )
                  IF( IMODE.GE.4 ) THEN
                     IF( IMODE.EQ.4 ) THEN
                        ILOW = 1
                        ISTEP = 1
                        IHIGH = MAX( 1, N / 2 )
                     ELSE IF( IMODE.EQ.5 ) THEN
                        ILOW = MAX( 1, N / 2 )
                        ISTEP = 1
                        IHIGH = N
                     ELSE IF( IMODE.EQ.6 ) THEN
                        ILOW = 1
                        ISTEP = 2
                        IHIGH = N
                     END IF
                     DO 40 I = ILOW, IHIGH, ISTEP
                        IWORK( I ) = 1
   40                CONTINUE
                  END IF
                  CALL DLAORD( 'Decreasing', MNMIN, S, 1 )
               END IF

               DO 60 INB = 1, NNB

                  // Do for each pair of values (NB,NX) in NBVAL and NXVAL.

                  NB = NBVAL( INB )
                  CALL XLAENV( 1, NB )
                  NX = NXVAL( INB )
                  CALL XLAENV( 3, NX )

                  // Get a working copy of COPYA into A and a copy of
                  // vector IWORK.

                  CALL DLACPY( 'All', M, N, COPYA, LDA, A, LDA )
                  CALL ICOPY( N, IWORK( 1 ), 1, IWORK( N+1 ), 1 )

                  // Compute the QR factorization with pivoting of A

                  LW = MAX( 1, 2*N+NB*( N+1 ) )

                  // Compute the QP3 factorization of A

                  SRNAMT = 'DGEQP3'
                  CALL DGEQP3( M, N, A, LDA, IWORK( N+1 ), TAU, WORK, LW, INFO )

                  // Compute norm(svd(a) - svd(r))

                  RESULT( 1 ) = DQRT12( M, N, A, LDA, S, WORK, LWORK )

                  // Compute norm( A*P - Q*R )

                  RESULT( 2 ) = DQPT01( M, N, MNMIN, COPYA, A, LDA, TAU, IWORK( N+1 ), WORK, LWORK )

                  // Compute Q'*Q

                  RESULT( 3 ) = DQRT11( M, MNMIN, A, LDA, TAU, WORK, LWORK )

                  // Print information about the tests that did not pass
                 t // he threshold.

                  DO 50 K = 1, NTESTS
                     IF( RESULT( K ).GE.THRESH ) THEN
                        IF( NFAIL.EQ.0 .AND. NERRS.EQ.0 ) CALL ALAHD( NOUT, PATH )                         WRITE( NOUT, FMT = 9999 )'DGEQP3', M, N, NB, IMODE, K, RESULT( K )
                        NFAIL = NFAIL + 1
                     END IF
   50             CONTINUE
                  NRUN = NRUN + NTESTS

   60          CONTINUE
   70       CONTINUE
   80    CONTINUE
   90 CONTINUE

      // Print a summary of the results.

      CALL ALASUM( PATH, NOUT, NFAIL, NRUN, NERRS )

 9999 FORMAT( 1X, A, ' M =', I5, ', N =', I5, ', NB =', I4, ', type ',
     $      I2, ', test ', I2, ', ratio =', G12.5 )

      // End of DCHKQ3

      END
