      SUBROUTINE CDRVLS( DOTYPE, NM, MVAL, NN, NVAL, NNS, NSVAL, NNB, NBVAL, NXVAL, THRESH, TSTERR, A, COPYA, B, COPYB, C, S, COPYS, NOUT )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               TSTERR;
      int                NM, NN, NNB, NNS, NOUT;
      REAL               THRESH
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                MVAL( * ), NBVAL( * ), NSVAL( * ), NVAL( * ), NXVAL( * );
      REAL               COPYS( * ), S( * )
      COMPLEX            A( * ), B( * ), C( * ), COPYA( * ), COPYB( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      int                NTESTS;
      const              NTESTS = 18 ;
      int                SMLSIZ;
      const              SMLSIZ = 25 ;
      REAL               ONE, ZERO
      const              ONE = 1.0, ZERO = 0.0 ;
      COMPLEX            CONE, CZERO
      const              CONE = ( 1.0, 0.0 ), CZERO = ( 0.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      String             TRANS;
      String             PATH;
      int                CRANK, I, IM, IMB, IN, INB, INFO, INS, IRANK, ISCALE, ITRAN, ITYPE, J, K, LDA, LDB, LDWORK, LWLSY, LWORK, M, MNMIN, N, NB, NCOLS, NERRS, NFAIL, NRHS, NROWS, NRUN, RANK, MB, MMAX, NMAX, NSMAX, LIWORK, LRWORK, LWORK_CGELS, LWORK_CGELST, LWORK_CGETSLS, LWORK_CGELSS, LWORK_CGELSY,  LWORK_CGELSD, LRWORK_CGELSY, LRWORK_CGELSS, LRWORK_CGELSD;
      REAL               EPS, NORMA, NORMB, RCOND
      // ..
      // .. Local Arrays ..
      int                ISEED( 4 ), ISEEDY( 4 ), IWQ( 1 );
      REAL               RESULT( NTESTS ), RWQ( 1 )
      COMPLEX            WQ( 1 )
      // ..
      // .. Allocatable Arrays ..
      COMPLEX, ALLOCATABLE :: WORK (:)
      REAL, ALLOCATABLE :: RWORK (:), WORK2 (:)
      int    , ALLOCATABLE :: IWORK (:);
      // ..
      // .. External Functions ..
      REAL               CQRT12, CQRT14, CQRT17, SASUM, SLAMCH
      // EXTERNAL CQRT12, CQRT14, CQRT17, SASUM, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASVM, CERRLS, CGELS, CGELSD, CGELSS, CGELST, CGELSY, CGEMM, CGETSLS, CLACPY, CLARNV, CQRT13, CQRT15, CQRT16, CSSCAL, SAXPY, XLAENV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, INT, REAL, SQRT
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
      DATA               ISEEDY / 1988, 1989, 1990, 1991 /
      // ..
      // .. Executable Statements ..

      // Initialize constants and the random number seed.

      PATH( 1: 1 ) = 'Complex precision'
      PATH( 2: 3 ) = 'LS'
      NRUN = 0
      NFAIL = 0
      NERRS = 0
      for (I = 1; I <= 4; I++) { // 10
         ISEED( I ) = ISEEDY( I )
      } // 10
      EPS = SLAMCH( 'Epsilon' )

      // Threshold for rank estimation

      RCOND = SQRT( EPS ) - ( SQRT( EPS )-EPS ) / 2

      // Test the error exits

      xlaenv(9, SMLSIZ );
      if (TSTERR) CALL CERRLS( PATH, NOUT );

      // Print the header if NM = 0 or NN = 0 and THRESH = 0.

      IF( ( NM == 0 || NN == 0 ) && THRESH == ZERO ) CALL ALAHD( NOUT, PATH )
      INFOT = 0

      // Compute maximal workspace needed for all routines

      NMAX = 0
      MMAX = 0
      NSMAX = 0
      for (I = 1; I <= NM; I++) {
         if ( MVAL( I ) > MMAX ) {
            MMAX = MVAL( I )
         }
      }
      for (I = 1; I <= NN; I++) {
         if ( NVAL( I ) > NMAX ) {
            NMAX = NVAL( I )
         }
      }
      for (I = 1; I <= NNS; I++) {
         if ( NSVAL( I ) > NSMAX ) {
            NSMAX = NSVAL( I )
         }
      }
      M = MMAX
      N = NMAX
      NRHS = NSMAX
      MNMIN = MAX( MIN( M, N ), 1 )

      // Compute workspace needed for routines
      // CQRT14, CQRT17 (two side cases), CQRT15 and CQRT12

      LWORK = MAX( 1, ( M+N )*NRHS, ( N+NRHS )*( M+2 ), ( M+NRHS )*( N+2 ), MAX( M+MNMIN, NRHS*MNMIN,2*N+M ), MAX( M*N+4*MNMIN+MAX(M,N), M*N+2*MNMIN+4*N ) )
      LRWORK = 1
      LIWORK = 1

      // Iterate through all test cases and compute necessary workspace
      // sizes for ?GELS, ?GELST, ?GETSLS, ?GELSY, ?GELSS and ?GELSD
      // routines.

      for (IM = 1; IM <= NM; IM++) {
         M = MVAL( IM )
         LDA = MAX( 1, M )
         for (IN = 1; IN <= NN; IN++) {
            N = NVAL( IN )
            MNMIN = MAX(MIN( M, N ),1)
            LDB = MAX( 1, M, N )
            for (INS = 1; INS <= NNS; INS++) {
               NRHS = NSVAL( INS )
               for (IRANK = 1; IRANK <= 2; IRANK++) {
                  for (ISCALE = 1; ISCALE <= 3; ISCALE++) {
                     ITYPE = ( IRANK-1 )*3 + ISCALE
                     if ( DOTYPE( ITYPE ) ) {
                        if ( IRANK == 1 ) {
                           for (ITRAN = 1; ITRAN <= 2; ITRAN++) {
                              if ( ITRAN == 1 ) {
                                 TRANS = 'N'
                              } else {
                                 TRANS = 'C'
                              }

                              // Compute workspace needed for CGELS
                              cgels(TRANS, M, N, NRHS, A, LDA, B, LDB, WQ, -1, INFO );
                              LWORK_CGELS = INT( WQ( 1 ) )
                              // Compute workspace needed for CGELST
                              cgelst(TRANS, M, N, NRHS, A, LDA, B, LDB, WQ, -1, INFO );
                              LWORK_CGELST = INT ( WQ ( 1 ) )
                              // Compute workspace needed for CGETSLS
                              cgetsls(TRANS, M, N, NRHS, A, LDA, B, LDB, WQ, -1, INFO );
                              LWORK_CGETSLS = INT( WQ( 1 ) )
                           }
                        }
                        // Compute workspace needed for CGELSY
                        cgelsy(M, N, NRHS, A, LDA, B, LDB, IWQ, RCOND, CRANK, WQ, -1, RWQ, INFO );
                        LWORK_CGELSY = INT( WQ( 1 ) )
                        LRWORK_CGELSY = 2*N
                        // Compute workspace needed for CGELSS
                        cgelss(M, N, NRHS, A, LDA, B, LDB, S, RCOND, CRANK, WQ, -1, RWQ, INFO );
                        LWORK_CGELSS = INT( WQ( 1 ) )
                        LRWORK_CGELSS = 5*MNMIN
                        // Compute workspace needed for CGELSD
                        cgelsd(M, N, NRHS, A, LDA, B, LDB, S, RCOND, CRANK, WQ, -1, RWQ, IWQ, INFO );
                        LWORK_CGELSD = INT( WQ( 1 ) )
                        LRWORK_CGELSD = INT( RWQ ( 1 ) )
                        // Compute LIWORK workspace needed for CGELSY and CGELSD
                        LIWORK = MAX( LIWORK, N, IWQ ( 1 ) )
                        // Compute LRWORK workspace needed for CGELSY, CGELSS and CGELSD
                        LRWORK = MAX( LRWORK, LRWORK_CGELSY, LRWORK_CGELSS, LRWORK_CGELSD )
                        // Compute LWORK workspace needed for all functions
                        LWORK = MAX( LWORK, LWORK_CGELS, LWORK_CGETSLS, LWORK_CGELSY, LWORK_CGELSS, LWORK_CGELSD )
                     }
                  }
               }
            }
         }
      }

      LWLSY = LWORK

      ALLOCATE( WORK( LWORK ) )
      ALLOCATE( IWORK( LIWORK ) )
      ALLOCATE( RWORK( LRWORK ) )
      ALLOCATE( WORK2( 2 * LWORK ) )

      for (IM = 1; IM <= NM; IM++) { // 140
         M = MVAL( IM )
         LDA = MAX( 1, M )

         for (IN = 1; IN <= NN; IN++) { // 130
            N = NVAL( IN )
            MNMIN = MAX(MIN( M, N ),1)
            LDB = MAX( 1, M, N )
            MB = (MNMIN+1)

            for (INS = 1; INS <= NNS; INS++) { // 120
               NRHS = NSVAL( INS )

               for (IRANK = 1; IRANK <= 2; IRANK++) { // 110
                  for (ISCALE = 1; ISCALE <= 3; ISCALE++) { // 100
                     ITYPE = ( IRANK-1 )*3 + ISCALE
                     IF( !DOTYPE( ITYPE ) ) GO TO 100
                  // =====================================================
                        // Begin test CGELS
                  // =====================================================
                     if ( IRANK == 1 ) {

                        // Generate a matrix of scaling type ISCALE

                        cqrt13(ISCALE, M, N, COPYA, LDA, NORMA, ISEED );

                        // Loop for testing different block sizes.

                        for (INB = 1; INB <= NNB; INB++) {
                           NB = NBVAL( INB )
                           xlaenv(1, NB );
                           xlaenv(3, NXVAL( INB ) );

                           // Loop for testing non-transposed and transposed.

                           for (ITRAN = 1; ITRAN <= 2; ITRAN++) {
                              if ( ITRAN == 1 ) {
                                 TRANS = 'N'
                                 NROWS = M
                                 NCOLS = N
                              } else {
                                 TRANS = 'C'
                                 NROWS = N
                                 NCOLS = M
                              }
                              LDWORK = MAX( 1, NCOLS )

                              // Set up a consistent rhs

                              if ( NCOLS > 0 ) {
                                 clarnv(2, ISEED, NCOLS*NRHS, WORK );
                                 csscal(NCOLS*NRHS, ONE / REAL( NCOLS ), WORK, 1 );
                              }
                              cgemm(TRANS, 'No transpose', NROWS, NRHS, NCOLS, CONE, COPYA, LDA, WORK, LDWORK, CZERO, B, LDB );
                              clacpy('Full', NROWS, NRHS, B, LDB, COPYB, LDB );

                              // Solve LS or overdetermined system

                              if ( M > 0 && N > 0 ) {
                                 clacpy('Full', M, N, COPYA, LDA, A, LDA );
                                 clacpy('Full', NROWS, NRHS, COPYB, LDB, B, LDB );
                              }
                              SRNAMT = 'CGELS '
                              cgels(TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO );

                              if (INFO != 0) CALL ALAERH( PATH, 'CGELS ', INFO, 0, TRANS, M, N, NRHS, -1, NB, ITYPE, NFAIL, NERRS, NOUT );

                              // Test 1: Check correctness of results
                              // for CGELS, compute the residual:
                              // RESID = norm(B - A*X) /
                              // / ( max(m,n) * norm(A) * norm(X) * EPS )

                              if (NROWS > 0 && NRHS > 0) CALL CLACPY( 'Full', NROWS, NRHS, COPYB, LDB, C, LDB );
                              cqrt16(TRANS, M, N, NRHS, COPYA, LDA, B, LDB, C, LDB, RWORK, RESULT( 1 ) );

                              // Test 2: Check correctness of results
                              // for CGELS.

                              if ( ( ITRAN == 1 && M >= N ) || ( ITRAN == 2 && M < N ) ) {

                                 // Solving LS system

                                 RESULT( 2 ) = CQRT17( TRANS, 1, M, N, NRHS, COPYA, LDA, B, LDB, COPYB, LDB, C, WORK, LWORK )
                              } else {

                                 // Solving overdetermined system

                                 RESULT( 2 ) = CQRT14( TRANS, M, N, NRHS, COPYA, LDA, B, LDB, WORK, LWORK )
                              }

                              // Print information about the tests that
                              // did not pass the threshold.

                              for (K = 1; K <= 2; K++) {
                                 if ( RESULT( K ) >= THRESH ) {
                                    if (NFAIL == 0 && NERRS == 0) CALL ALAHD( NOUT, PATH )                                     WRITE( NOUT, FMT = 9999 )TRANS, M, N, NRHS, NB, ITYPE, K, RESULT( K );
                                    NFAIL = NFAIL + 1
                                 }
                              }
                              NRUN = NRUN + 2
                           }
                        }
                     }
                  // =====================================================
                        // End test CGELS
                  // =====================================================
                  // =====================================================
                        // Begin test CGELST
                  // =====================================================
                     if ( IRANK == 1 ) {

                        // Generate a matrix of scaling type ISCALE

                        cqrt13(ISCALE, M, N, COPYA, LDA, NORMA, ISEED );

                        // Loop for testing different block sizes.

                        for (INB = 1; INB <= NNB; INB++) {
                           NB = NBVAL( INB )
                           xlaenv(1, NB );
                           xlaenv(3, NXVAL( INB ) );

                           // Loop for testing non-transposed and transposed.

                           for (ITRAN = 1; ITRAN <= 2; ITRAN++) {
                              if ( ITRAN == 1 ) {
                                 TRANS = 'N'
                                 NROWS = M
                                 NCOLS = N
                              } else {
                                 TRANS = 'C'
                                 NROWS = N
                                 NCOLS = M
                              }
                              LDWORK = MAX( 1, NCOLS )

                              // Set up a consistent rhs

                              if ( NCOLS > 0 ) {
                                 clarnv(2, ISEED, NCOLS*NRHS, WORK );
                                 csscal(NCOLS*NRHS, ONE / REAL( NCOLS ), WORK, 1 );
                              }
                              cgemm(TRANS, 'No transpose', NROWS, NRHS, NCOLS, CONE, COPYA, LDA, WORK, LDWORK, CZERO, B, LDB );
                              clacpy('Full', NROWS, NRHS, B, LDB, COPYB, LDB );

                              // Solve LS or overdetermined system

                              if ( M > 0 && N > 0 ) {
                                 clacpy('Full', M, N, COPYA, LDA, A, LDA );
                                 clacpy('Full', NROWS, NRHS, COPYB, LDB, B, LDB );
                              }
                              SRNAMT = 'CGELST'
                              cgelst(TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO );

                              if (INFO != 0) CALL ALAERH( PATH, 'CGELST', INFO, 0, TRANS, M, N, NRHS, -1, NB, ITYPE, NFAIL, NERRS, NOUT );

                              // Test 3: Check correctness of results
                              // for CGELST, compute the residual:
                              // RESID = norm(B - A*X) /
                              // / ( max(m,n) * norm(A) * norm(X) * EPS )

                              if (NROWS > 0 && NRHS > 0) CALL CLACPY( 'Full', NROWS, NRHS, COPYB, LDB, C, LDB );
                              cqrt16(TRANS, M, N, NRHS, COPYA, LDA, B, LDB, C, LDB, RWORK, RESULT( 3 ) );

                              // Test 4: Check correctness of results
                              // for CGELST.

                              if ( ( ITRAN == 1 && M >= N ) || ( ITRAN == 2 && M < N ) ) {

                                 // Solving LS system

                                 RESULT( 4 ) = CQRT17( TRANS, 1, M, N, NRHS, COPYA, LDA, B, LDB, COPYB, LDB, C, WORK, LWORK )
                              } else {

                                 // Solving overdetermined system

                                 RESULT( 4 ) = CQRT14( TRANS, M, N, NRHS, COPYA, LDA, B, LDB, WORK, LWORK )
                              }

                              // Print information about the tests that
                              // did not pass the threshold.

                              for (K = 3; K <= 4; K++) {
                                 if ( RESULT( K ) >= THRESH ) {
                                    if (NFAIL == 0 && NERRS == 0) CALL ALAHD( NOUT, PATH )                                     WRITE( NOUT, FMT = 9999 )TRANS, M, N, NRHS, NB, ITYPE, K, RESULT( K );
                                    NFAIL = NFAIL + 1
                                 }
                              }
                              NRUN = NRUN + 2
                           }
                        }
                     }
                  // =====================================================
                        // End test CGELST
                  // =====================================================
                  // =====================================================
                        // Begin test CGELSTSLS
                  // =====================================================
                     if ( IRANK == 1 ) {

                        // Generate a matrix of scaling type ISCALE

                        cqrt13(ISCALE, M, N, COPYA, LDA, NORMA, ISEED );

                        // Loop for testing different block sizes MB.

                        for (INB = 1; INB <= NNB; INB++) {
                           MB = NBVAL( INB )
                           xlaenv(1, MB );

                           // Loop for testing different block sizes NB.

                           for (IMB = 1; IMB <= NNB; IMB++) {
                              NB = NBVAL( IMB )
                              xlaenv(2, NB );

                              // Loop for testing non-transposed
                              // and transposed.

                              for (ITRAN = 1; ITRAN <= 2; ITRAN++) {
                                 if ( ITRAN == 1 ) {
                                    TRANS = 'N'
                                    NROWS = M
                                    NCOLS = N
                                 } else {
                                    TRANS = 'C'
                                    NROWS = N
                                    NCOLS = M
                                 }
                                 LDWORK = MAX( 1, NCOLS )

                                 // Set up a consistent rhs

                                 if ( NCOLS > 0 ) {
                                    clarnv(2, ISEED, NCOLS*NRHS, WORK );
                                    cscal(NCOLS*NRHS, CONE / REAL( NCOLS ), WORK, 1 );
                                 }
                                 cgemm(TRANS, 'No transpose', NROWS, NRHS, NCOLS, CONE, COPYA, LDA, WORK, LDWORK, CZERO, B, LDB );
                                 clacpy('Full', NROWS, NRHS, B, LDB, COPYB, LDB );

                                 // Solve LS or overdetermined system

                                 if ( M > 0 && N > 0 ) {
                                    clacpy('Full', M, N, COPYA, LDA, A, LDA );
                                    clacpy('Full', NROWS, NRHS, COPYB, LDB, B, LDB );
                                 }
                                 SRNAMT = 'CGETSLS '
                                 cgetsls(TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )                                  IF( INFO != 0 ) CALL ALAERH( PATH, 'CGETSLS ', INFO, 0, TRANS, M, N, NRHS, -1, NB, ITYPE, NFAIL, NERRS, NOUT );

                              // Test 5: Check correctness of results
                              // for CGETSLS, compute the residual:
                              // RESID = norm(B - A*X) /
                              // / ( max(m,n) * norm(A) * norm(X) * EPS )

                                 if (NROWS > 0 && NRHS > 0) CALL CLACPY( 'Full', NROWS, NRHS, COPYB, LDB, C, LDB );
                                 cqrt16(TRANS, M, N, NRHS, COPYA, LDA, B, LDB, C, LDB, WORK2, RESULT( 5 ) );

                              // Test 6: Check correctness of results
                              // for CGETSLS.

                                 if ( ( ITRAN == 1 && M >= N ) || ( ITRAN == 2 && M < N ) ) {

                                    // Solving LS system, compute:
                                    // r = norm((B- A*X)**T * A) /
                                  // / (norm(A)*norm(B)*max(M,N,NRHS)*EPS)

                                    RESULT( 6 ) = CQRT17( TRANS, 1, M, N, NRHS, COPYA, LDA, B, LDB, COPYB, LDB, C, WORK, LWORK )
                                 } else {

                                    // Solving overdetermined system

                                    RESULT( 6 ) = CQRT14( TRANS, M, N, NRHS, COPYA, LDA, B, LDB, WORK, LWORK )
                                 }

                                 // Print information about the tests that
                                 // did not pass the threshold.

                                 for (K = 5; K <= 6; K++) {
                                    if ( RESULT( K ) >= THRESH ) {
                                       if (NFAIL == 0 && NERRS == 0) CALL ALAHD( NOUT, PATH )                                        WRITE( NOUT, FMT = 9997 )TRANS, M, N, NRHS, MB, NB, ITYPE, K, RESULT( K );
                                          NFAIL = NFAIL + 1
                                    }
                                 }
                                 NRUN = NRUN + 2
                              }
                           }
                        }
                     }
                  // =====================================================
                        // End test CGELSTSLS
                  // ====================================================

                     // Generate a matrix of scaling type ISCALE and rank
                     // type IRANK.

                     cqrt15(ISCALE, IRANK, M, N, NRHS, COPYA, LDA, COPYB, LDB, COPYS, RANK, NORMA, NORMB, ISEED, WORK, LWORK );

                     // workspace used: MAX(M+MIN(M,N),NRHS*MIN(M,N),2*N+M)

                     LDWORK = MAX( 1, M )

                     // Loop for testing different block sizes.

                     for (INB = 1; INB <= NNB; INB++) { // 90
                        NB = NBVAL( INB )
                        xlaenv(1, NB );
                        xlaenv(3, NXVAL( INB ) );

                        // Test CGELSY

                        // CGELSY:  Compute the minimum-norm solution
                        // X to min( norm( A * X - B ) )
                        // using the rank-revealing orthogonal
                        // factorization.

                        clacpy('Full', M, N, COPYA, LDA, A, LDA );
                        clacpy('Full', M, NRHS, COPYB, LDB, B, LDB );

                        // Initialize vector IWORK.

                        for (J = 1; J <= N; J++) { // 70
                           IWORK( J ) = 0
                        } // 70

                        SRNAMT = 'CGELSY'
                        cgelsy(M, N, NRHS, A, LDA, B, LDB, IWORK, RCOND, CRANK, WORK, LWLSY, RWORK, INFO )                         IF( INFO != 0 ) CALL ALAERH( PATH, 'CGELSY', INFO, 0, ' ', M, N, NRHS, -1, NB, ITYPE, NFAIL, NERRS, NOUT );

                        // workspace used: 2*MNMIN+NB*NB+NB*MAX(N,NRHS)

                        // Test 7:  Compute relative error in svd
                                 // workspace: M*N + 4*MIN(M,N) + MAX(M,N)

                        RESULT( 7 ) = CQRT12( CRANK, CRANK, A, LDA, COPYS, WORK, LWORK, RWORK )

                        // Test 8:  Compute error in solution
                                 // workspace:  M*NRHS + M

                        clacpy('Full', M, NRHS, COPYB, LDB, WORK, LDWORK );
                        cqrt16('No transpose', M, N, NRHS, COPYA, LDA, B, LDB, WORK, LDWORK, RWORK, RESULT( 8 ) );

                        // Test 9:  Check norm of r'*A
                                 // workspace: NRHS*(M+N)

                        RESULT( 9 ) = ZERO
                        if (M > CRANK) RESULT( 9 ) = CQRT17( 'No transpose', 1, M, N, NRHS, COPYA, LDA, B, LDB, COPYB, LDB, C, WORK, LWORK );

                        // Test 10:  Check if x is in the rowspace of A
                                 // workspace: (M+NRHS)*(N+2)

                        RESULT( 10 ) = ZERO

                        if (N > CRANK) RESULT( 10 ) = CQRT14( 'No transpose', M, N, NRHS, COPYA, LDA, B, LDB, WORK, LWORK );

                        // Test CGELSS

                        // CGELSS:  Compute the minimum-norm solution
                        // X to min( norm( A * X - B ) )
                        // using the SVD.

                        clacpy('Full', M, N, COPYA, LDA, A, LDA );
                        clacpy('Full', M, NRHS, COPYB, LDB, B, LDB );
                        SRNAMT = 'CGELSS'
                        cgelss(M, N, NRHS, A, LDA, B, LDB, S, RCOND, CRANK, WORK, LWORK, RWORK, INFO );

                        if (INFO != 0) CALL ALAERH( PATH, 'CGELSS', INFO, 0, ' ', M, N, NRHS, -1, NB, ITYPE, NFAIL, NERRS, NOUT );

                        // workspace used: 3*min(m,n) +
                                        // max(2*min(m,n),nrhs,max(m,n))

                        // Test 11:  Compute relative error in svd

                        if ( RANK > 0 ) {
                           saxpy(MNMIN, -ONE, COPYS, 1, S, 1 );
                           RESULT( 11 ) = SASUM( MNMIN, S, 1 ) / SASUM( MNMIN, COPYS, 1 ) / ( EPS*REAL( MNMIN ) )
                        } else {
                           RESULT( 11 ) = ZERO
                        }

                        // Test 12:  Compute error in solution

                        clacpy('Full', M, NRHS, COPYB, LDB, WORK, LDWORK );
                        cqrt16('No transpose', M, N, NRHS, COPYA, LDA, B, LDB, WORK, LDWORK, RWORK, RESULT( 12 ) );

                        // Test 13:  Check norm of r'*A

                        RESULT( 13 ) = ZERO
                        if (M > CRANK) RESULT( 13 ) = CQRT17( 'No transpose', 1, M, N, NRHS, COPYA, LDA, B, LDB, COPYB, LDB, C, WORK, LWORK );

                        // Test 14:  Check if x is in the rowspace of A

                        RESULT( 14 ) = ZERO
                        if (N > CRANK) RESULT( 14 ) = CQRT14( 'No transpose', M, N, NRHS, COPYA, LDA, B, LDB, WORK, LWORK );

                        // Test CGELSD

                        // CGELSD:  Compute the minimum-norm solution X
                        // to min( norm( A * X - B ) ) using a
                        // divide and conquer SVD.

                        xlaenv(9, 25 );

                        clacpy('Full', M, N, COPYA, LDA, A, LDA );
                        clacpy('Full', M, NRHS, COPYB, LDB, B, LDB );

                        SRNAMT = 'CGELSD'
                        cgelsd(M, N, NRHS, A, LDA, B, LDB, S, RCOND, CRANK, WORK, LWORK, RWORK, IWORK, INFO )                         IF( INFO != 0 ) CALL ALAERH( PATH, 'CGELSD', INFO, 0, ' ', M, N, NRHS, -1, NB, ITYPE, NFAIL, NERRS, NOUT );

                        // Test 15:  Compute relative error in svd

                        if ( RANK > 0 ) {
                           saxpy(MNMIN, -ONE, COPYS, 1, S, 1 );
                           RESULT( 15 ) = SASUM( MNMIN, S, 1 ) / SASUM( MNMIN, COPYS, 1 ) / ( EPS*REAL( MNMIN ) )
                        } else {
                           RESULT( 15 ) = ZERO
                        }

                        // Test 16:  Compute error in solution

                        clacpy('Full', M, NRHS, COPYB, LDB, WORK, LDWORK );
                        cqrt16('No transpose', M, N, NRHS, COPYA, LDA, B, LDB, WORK, LDWORK, RWORK, RESULT( 16 ) );

                        // Test 17:  Check norm of r'*A

                        RESULT( 17 ) = ZERO
                        if (M > CRANK) RESULT( 17 ) = CQRT17( 'No transpose', 1, M, N, NRHS, COPYA, LDA, B, LDB, COPYB, LDB, C, WORK, LWORK );

                        // Test 18:  Check if x is in the rowspace of A

                        RESULT( 18 ) = ZERO
                        if (N > CRANK) RESULT( 18 ) = CQRT14( 'No transpose', M, N, NRHS, COPYA, LDA, B, LDB, WORK, LWORK );

                        // Print information about the tests that did not
                        // pass the threshold.

                        for (K = 7; K <= 18; K++) { // 80
                           if ( RESULT( K ) >= THRESH ) {
                              if (NFAIL == 0 && NERRS == 0) CALL ALAHD( NOUT, PATH )                               WRITE( NOUT, FMT = 9998 )M, N, NRHS, NB, ITYPE, K, RESULT( K );
                              NFAIL = NFAIL + 1
                           }
                        } // 80
                        NRUN = NRUN + 12

                     } // 90
                  } // 100
               } // 110
            } // 120
         } // 130
      } // 140

      // Print a summary of the results.

      alasvm(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( ' TRANS=''', A1, ''', M=', I5, ', N=', I5, ', NRHS=', I4, ', NB=', I4, ', type', I2, ', test(', I2, ')=', G12.5 )
 9998 FORMAT( ' M=', I5, ', N=', I5, ', NRHS=', I4, ', NB=', I4, ', type', I2, ', test(', I2, ')=', G12.5 )
 9997 FORMAT( ' TRANS=''', A1,' M=', I5, ', N=', I5, ', NRHS=', I4, ', MB=', I4,', NB=', I4,', type', I2, ', test(', I2, ')=', G12.5 )

      DEALLOCATE( WORK )
      DEALLOCATE( RWORK )
      DEALLOCATE( IWORK )
      RETURN

      // End of CDRVLS

      }
