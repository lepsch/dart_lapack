      void zdrvls(DOTYPE, NM, MVAL, NN, NVAL, NNS, NSVAL, NNB, NBVAL, NXVAL, THRESH, TSTERR, A, COPYA, B, COPYB, C, S, COPYS, NOUT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               TSTERR;
      int                NM, NN, NNB, NNS, NOUT;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                MVAL( * ), NBVAL( * ), NSVAL( * ), NVAL( * ), NXVAL( * );
      double             COPYS( * ), S( * );
      Complex         A( * ), B( * ), C( * ), COPYA( * ), COPYB( * );
      // ..

// =====================================================================

      // .. Parameters ..
      int                NTESTS;
      const              NTESTS = 18 ;
      int                SMLSIZ;
      const              SMLSIZ = 25 ;
      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      Complex         CONE, CZERO;
      const              CONE = ( 1.0, 0.0 ), CZERO = ( 0.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      String             TRANS;
      String             PATH;
      int                CRANK, I, IM, IMB, IN, INB, INFO, INS, IRANK, ISCALE, ITRAN, ITYPE, J, K, LDA, LDB, LDWORK, LWLSY, LWORK, M, MNMIN, N, NB, NCOLS, NERRS, NFAIL, NRHS, NROWS, NRUN, RANK, MB, MMAX, NMAX, NSMAX, LIWORK, LRWORK, LWORK_ZGELS, LWORK_ZGELST, LWORK_ZGETSLS, LWORK_ZGELSS, LWORK_ZGELSY, LWORK_ZGELSD, LRWORK_ZGELSY, LRWORK_ZGELSS, LRWORK_ZGELSD;
      double             EPS, NORMA, NORMB, RCOND;
      // ..
      // .. Local Arrays ..
      int                ISEED( 4 ), ISEEDY( 4 ), IWQ( 1 );
      double             RESULT( NTESTS ), RWQ( 1 );
      Complex         WQ( 1 );
      // ..
      // .. Allocatable Arrays ..
      Complex, ALLOCATABLE :: WORK (:);
      double          , ALLOCATABLE :: RWORK (:), WORK2 (:);
      int    , ALLOCATABLE :: IWORK (:);
      // ..
      // .. External Functions ..
      //- double             DASUM, DLAMCH, ZQRT12, ZQRT14, ZQRT17;
      // EXTERNAL DASUM, DLAMCH, ZQRT12, ZQRT14, ZQRT17
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASVM, DAXPY, ZERRLS, ZGELS, ZGELSD, ZGELSS, ZGELST, ZGELSY, ZGEMM, ZGETSLS, ZLACPY, ZLARNV, ZQRT13, ZQRT15, ZQRT16, ZDSCAL, XLAENV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, MAX, MIN, INT, SQRT
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

      PATH( 1: 1 ) = 'Zomplex precision';
      PATH( 2: 3 ) = 'LS';
      NRUN = 0;
      NFAIL = 0;
      NERRS = 0;
      for (I = 1; I <= 4; I++) { // 10
         ISEED( I ) = ISEEDY( I );
      } // 10
      EPS = DLAMCH( 'Epsilon' );

      // Threshold for rank estimation

      RCOND = sqrt( EPS ) - ( sqrt( EPS )-EPS ) / 2;

      // Test the error exits

      xlaenv(9, SMLSIZ );
      if (TSTERR) zerrls( PATH, NOUT );

      // Print the header if NM = 0 or NN = 0 and THRESH = 0.

      if( ( NM == 0 || NN == 0 ) && THRESH == ZERO ) alahd( NOUT, PATH );
      INFOT = 0;

      // Compute maximal workspace needed for all routines

      NMAX = 0;
      MMAX = 0;
      NSMAX = 0;
      for (I = 1; I <= NM; I++) {
         if ( MVAL( I ) > MMAX ) {
            MMAX = MVAL( I );
         }
      }
      for (I = 1; I <= NN; I++) {
         if ( NVAL( I ) > NMAX ) {
            NMAX = NVAL( I );
         }
      }
      for (I = 1; I <= NNS; I++) {
         if ( NSVAL( I ) > NSMAX ) {
            NSMAX = NSVAL( I );
         }
      }
      M = MMAX;
      N = NMAX;
      NRHS = NSMAX;
      MNMIN = max( min( M, N ), 1 );

      // Compute workspace needed for routines
      // ZQRT14, ZQRT17 (two side cases), ZQRT15 and ZQRT12

      LWORK = max( 1, ( M+N )*NRHS, ( N+NRHS )*( M+2 ), ( M+NRHS )*( N+2 ), max( M+MNMIN, NRHS*MNMIN,2*N+M ), max( M*N+4*MNMIN+max(M,N), M*N+2*MNMIN+4*N ) );
      LRWORK = 1;
      LIWORK = 1;

      // Iterate through all test cases and compute necessary workspace
      // sizes for ?GELS, ?GELST, ?GETSLS, ?GELSY, ?GELSS and ?GELSD
      // routines.

      for (IM = 1; IM <= NM; IM++) {
         M = MVAL( IM );
         LDA = max( 1, M );
         for (IN = 1; IN <= NN; IN++) {
            N = NVAL( IN );
            MNMIN = max(min( M, N ),1);
            LDB = max( 1, M, N );
            for (INS = 1; INS <= NNS; INS++) {
               NRHS = NSVAL( INS );
               for (IRANK = 1; IRANK <= 2; IRANK++) {
                  for (ISCALE = 1; ISCALE <= 3; ISCALE++) {
                     ITYPE = ( IRANK-1 )*3 + ISCALE;
                     if ( DOTYPE( ITYPE ) ) {
                        if ( IRANK == 1 ) {
                           for (ITRAN = 1; ITRAN <= 2; ITRAN++) {
                              if ( ITRAN == 1 ) {
                                 TRANS = 'N';
                              } else {
                                 TRANS = 'C';
                              }

                              // Compute workspace needed for ZGELS
                              zgels(TRANS, M, N, NRHS, A, LDA, B, LDB, WQ, -1, INFO );
                              LWORK_ZGELS = INT ( WQ( 1 ) );
                              // Compute workspace needed for ZGELST
                              zgelst(TRANS, M, N, NRHS, A, LDA, B, LDB, WQ, -1, INFO );
                              LWORK_ZGELST = INT ( WQ ( 1 ) );
                              // Compute workspace needed for ZGETSLS
                              zgetsls(TRANS, M, N, NRHS, A, LDA, B, LDB, WQ, -1, INFO );
                              LWORK_ZGETSLS = INT( WQ( 1 ) );
                           }
                        }
                        // Compute workspace needed for ZGELSY
                        zgelsy(M, N, NRHS, A, LDA, B, LDB, IWQ, RCOND, CRANK, WQ, -1, RWQ, INFO );
                        LWORK_ZGELSY = INT( WQ( 1 ) );
                        LRWORK_ZGELSY = 2*N;
                        // Compute workspace needed for ZGELSS
                        zgelss(M, N, NRHS, A, LDA, B, LDB, S, RCOND, CRANK, WQ, -1 , RWQ, INFO );
                        LWORK_ZGELSS = INT( WQ( 1 ) );
                        LRWORK_ZGELSS = 5*MNMIN;
                        // Compute workspace needed for ZGELSD
                        zgelsd(M, N, NRHS, A, LDA, B, LDB, S, RCOND, CRANK, WQ, -1, RWQ, IWQ, INFO );
                        LWORK_ZGELSD = INT( WQ( 1 ) );
                        LRWORK_ZGELSD = INT( RWQ ( 1 ) );
                        // Compute LIWORK workspace needed for ZGELSY and ZGELSD
                        LIWORK = max( LIWORK, N, IWQ( 1 ) );
                        // Compute LRWORK workspace needed for ZGELSY, ZGELSS and ZGELSD
                        LRWORK = max( LRWORK, LRWORK_ZGELSY, LRWORK_ZGELSS, LRWORK_ZGELSD );
                        // Compute LWORK workspace needed for all functions
                        LWORK = max( LWORK, LWORK_ZGELS, LWORK_ZGELST, LWORK_ZGETSLS, LWORK_ZGELSY, LWORK_ZGELSS, LWORK_ZGELSD );
                     }
                  }
               }
            }
         }
      }

      LWLSY = LWORK;

      ALLOCATE( WORK( LWORK ) );
      ALLOCATE( WORK2( 2 * LWORK ) );
      ALLOCATE( IWORK( LIWORK ) );
      ALLOCATE( RWORK( LRWORK ) );

      for (IM = 1; IM <= NM; IM++) { // 140
         M = MVAL( IM );
         LDA = max( 1, M );

         for (IN = 1; IN <= NN; IN++) { // 130
            N = NVAL( IN );
            MNMIN = max(min( M, N ),1);
            LDB = max( 1, M, N );
            MB = (MNMIN+1);

            for (INS = 1; INS <= NNS; INS++) { // 120
               NRHS = NSVAL( INS );

               for (IRANK = 1; IRANK <= 2; IRANK++) { // 110
                  for (ISCALE = 1; ISCALE <= 3; ISCALE++) { // 100
                     ITYPE = ( IRANK-1 )*3 + ISCALE;
                     if( !DOTYPE( ITYPE ) ) GO TO 100;
                  // =====================================================
                        // Begin test ZGELS
                  // =====================================================
                     if ( IRANK == 1 ) {

                        // Generate a matrix of scaling type ISCALE

                        zqrt13(ISCALE, M, N, COPYA, LDA, NORMA, ISEED );

                        // Loop for testing different block sizes.

                        for (INB = 1; INB <= NNB; INB++) {
                           NB = NBVAL( INB );
                           xlaenv(1, NB );
                           xlaenv(3, NXVAL( INB ) );

                           // Loop for testing non-transposed and transposed.

                           for (ITRAN = 1; ITRAN <= 2; ITRAN++) {
                              if ( ITRAN == 1 ) {
                                 TRANS = 'N';
                                 NROWS = M;
                                 NCOLS = N;
                              } else {
                                 TRANS = 'C';
                                 NROWS = N;
                                 NCOLS = M;
                              }
                              LDWORK = max( 1, NCOLS );

                              // Set up a consistent rhs

                              if ( NCOLS > 0 ) {
                                 zlarnv(2, ISEED, NCOLS*NRHS, WORK );
                                 zdscal(NCOLS*NRHS, ONE / DBLE( NCOLS ), WORK, 1 );
                              }
                              zgemm(TRANS, 'No transpose', NROWS, NRHS, NCOLS, CONE, COPYA, LDA, WORK, LDWORK, CZERO, B, LDB );
                              zlacpy('Full', NROWS, NRHS, B, LDB, COPYB, LDB );

                              // Solve LS or overdetermined system

                              if ( M > 0 && N > 0 ) {
                                 zlacpy('Full', M, N, COPYA, LDA, A, LDA );
                                 zlacpy('Full', NROWS, NRHS, COPYB, LDB, B, LDB );
                              }
                              SRNAMT = 'ZGELS ';
                              zgels(TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO );

                              if (INFO != 0) alaerh( PATH, 'ZGELS ', INFO, 0, TRANS, M, N, NRHS, -1, NB, ITYPE, NFAIL, NERRS, NOUT );

                              // Test 1: Check correctness of results
                              // for ZGELS, compute the residual:
                              // RESID = norm(B - A*X) /
                              // / ( max(m,n) * norm(A) * norm(X) * EPS )

                              if (NROWS > 0 && NRHS > 0) zlacpy( 'Full', NROWS, NRHS, COPYB, LDB, C, LDB );
                              zqrt16(TRANS, M, N, NRHS, COPYA, LDA, B, LDB, C, LDB, RWORK, RESULT( 1 ) );

                              // Test 2: Check correctness of results
                              // for ZGELS.

                              if ( ( ITRAN == 1 && M >= N ) || ( ITRAN == 2 && M < N ) ) {

                                 // Solving LS system

                                 RESULT( 2 ) = ZQRT17( TRANS, 1, M, N, NRHS, COPYA, LDA, B, LDB, COPYB, LDB, C, WORK, LWORK );
                              } else {

                                 // Solving overdetermined system

                                 RESULT( 2 ) = ZQRT14( TRANS, M, N, NRHS, COPYA, LDA, B, LDB, WORK, LWORK );
                              }

                              // Print information about the tests that
                              // did not pass the threshold.

                              for (K = 1; K <= 2; K++) {
                                 if ( RESULT( K ) >= THRESH ) {
                                    if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                                    WRITE( NOUT, FMT = 9999 )TRANS, M, N, NRHS, NB, ITYPE, K, RESULT( K );
                                    NFAIL = NFAIL + 1;
                                 }
                              }
                              NRUN = NRUN + 2;
                           }
                        }
                     }
                  // =====================================================
                        // End test ZGELS
                  // =====================================================
                  // =====================================================
                        // Begin test ZGELST
                  // =====================================================
                     if ( IRANK == 1 ) {

                        // Generate a matrix of scaling type ISCALE

                        zqrt13(ISCALE, M, N, COPYA, LDA, NORMA, ISEED );

                        // Loop for testing different block sizes.

                        for (INB = 1; INB <= NNB; INB++) {
                           NB = NBVAL( INB );
                           xlaenv(1, NB );
                           xlaenv(3, NXVAL( INB ) );

                           // Loop for testing non-transposed and transposed.

                           for (ITRAN = 1; ITRAN <= 2; ITRAN++) {
                              if ( ITRAN == 1 ) {
                                 TRANS = 'N';
                                 NROWS = M;
                                 NCOLS = N;
                              } else {
                                 TRANS = 'C';
                                 NROWS = N;
                                 NCOLS = M;
                              }
                              LDWORK = max( 1, NCOLS );

                              // Set up a consistent rhs

                              if ( NCOLS > 0 ) {
                                 zlarnv(2, ISEED, NCOLS*NRHS, WORK );
                                 zdscal(NCOLS*NRHS, ONE / DBLE( NCOLS ), WORK, 1 );
                              }
                              zgemm(TRANS, 'No transpose', NROWS, NRHS, NCOLS, CONE, COPYA, LDA, WORK, LDWORK, CZERO, B, LDB );
                              zlacpy('Full', NROWS, NRHS, B, LDB, COPYB, LDB );

                              // Solve LS or overdetermined system

                              if ( M > 0 && N > 0 ) {
                                 zlacpy('Full', M, N, COPYA, LDA, A, LDA );
                                 zlacpy('Full', NROWS, NRHS, COPYB, LDB, B, LDB );
                              }
                              SRNAMT = 'ZGELST';
                              zgelst(TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO );

                              if (INFO != 0) alaerh( PATH, 'ZGELST', INFO, 0, TRANS, M, N, NRHS, -1, NB, ITYPE, NFAIL, NERRS, NOUT );

                              // Test 3: Check correctness of results
                              // for ZGELST, compute the residual:
                              // RESID = norm(B - A*X) /
                              // / ( max(m,n) * norm(A) * norm(X) * EPS )

                              if (NROWS > 0 && NRHS > 0) zlacpy( 'Full', NROWS, NRHS, COPYB, LDB, C, LDB );
                              zqrt16(TRANS, M, N, NRHS, COPYA, LDA, B, LDB, C, LDB, RWORK, RESULT( 3 ) );

                              // Test 4: Check correctness of results
                              // for ZGELST.

                              if ( ( ITRAN == 1 && M >= N ) || ( ITRAN == 2 && M < N ) ) {

                                 // Solving LS system

                                 RESULT( 4 ) = ZQRT17( TRANS, 1, M, N, NRHS, COPYA, LDA, B, LDB, COPYB, LDB, C, WORK, LWORK );
                              } else {

                                 // Solving overdetermined system

                                 RESULT( 4 ) = ZQRT14( TRANS, M, N, NRHS, COPYA, LDA, B, LDB, WORK, LWORK );
                              }

                              // Print information about the tests that
                              // did not pass the threshold.

                              for (K = 3; K <= 4; K++) {
                                 if ( RESULT( K ) >= THRESH ) {
                                    if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                                    WRITE( NOUT, FMT = 9999 )TRANS, M, N, NRHS, NB, ITYPE, K, RESULT( K );
                                    NFAIL = NFAIL + 1;
                                 }
                              }
                              NRUN = NRUN + 2;
                           }
                        }
                     }
                  // =====================================================
                        // End test ZGELST
                  // =====================================================
                  // =====================================================
                        // Begin test ZGELSTSLS
                  // =====================================================
                     if ( IRANK == 1 ) {

                        // Generate a matrix of scaling type ISCALE

                        zqrt13(ISCALE, M, N, COPYA, LDA, NORMA, ISEED );

                        // Loop for testing different block sizes MB.

                        for (INB = 1; INB <= NNB; INB++) {
                           MB = NBVAL( INB );
                           xlaenv(1, MB );

                           // Loop for testing different block sizes NB.

                           for (IMB = 1; IMB <= NNB; IMB++) {
                              NB = NBVAL( IMB );
                              xlaenv(2, NB );

                              // Loop for testing non-transposed
                              // and transposed.

                              for (ITRAN = 1; ITRAN <= 2; ITRAN++) {
                                 if ( ITRAN == 1 ) {
                                    TRANS = 'N';
                                    NROWS = M;
                                    NCOLS = N;
                                 } else {
                                    TRANS = 'C';
                                    NROWS = N;
                                    NCOLS = M;
                                 }
                                 LDWORK = max( 1, NCOLS );

                                 // Set up a consistent rhs

                                 if ( NCOLS > 0 ) {
                                    zlarnv(2, ISEED, NCOLS*NRHS, WORK );
                                    zscal(NCOLS*NRHS, CONE / DBLE( NCOLS ), WORK, 1 );
                                 }
                                 zgemm(TRANS, 'No transpose', NROWS, NRHS, NCOLS, CONE, COPYA, LDA, WORK, LDWORK, CZERO, B, LDB );
                                 zlacpy('Full', NROWS, NRHS, B, LDB, COPYB, LDB );

                                 // Solve LS or overdetermined system

                                 if ( M > 0 && N > 0 ) {
                                    zlacpy('Full', M, N, COPYA, LDA, A, LDA );
                                    zlacpy('Full', NROWS, NRHS, COPYB, LDB, B, LDB );
                                 }
                                 SRNAMT = 'ZGETSLS ';
                                 zgetsls(TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )                                  IF( INFO != 0 ) CALL ALAERH( PATH, 'ZGETSLS ', INFO, 0, TRANS, M, N, NRHS, -1, NB, ITYPE, NFAIL, NERRS, NOUT );

                              // Test 5: Check correctness of results
                              // for ZGETSLS, compute the residual:
                              // RESID = norm(B - A*X) /
                              // / ( max(m,n) * norm(A) * norm(X) * EPS )

                                 if (NROWS > 0 && NRHS > 0) zlacpy( 'Full', NROWS, NRHS, COPYB, LDB, C, LDB );
                                 zqrt16(TRANS, M, N, NRHS, COPYA, LDA, B, LDB, C, LDB, WORK2, RESULT( 5 ) );

                              // Test 6: Check correctness of results
                              // for ZGETSLS.

                                 if ( ( ITRAN == 1 && M >= N ) || ( ITRAN == 2 && M < N ) ) {

                                    // Solving LS system, compute:
                                    // r = norm((B- A*X)**T * A) /
                                  // / (norm(A)*norm(B)*max(M,N,NRHS)*EPS)

                                    RESULT( 6 ) = ZQRT17( TRANS, 1, M, N, NRHS, COPYA, LDA, B, LDB, COPYB, LDB, C, WORK, LWORK );
                                 } else {

                                    // Solving overdetermined system

                                    RESULT( 6 ) = ZQRT14( TRANS, M, N, NRHS, COPYA, LDA, B, LDB, WORK, LWORK );
                                 }

                                 // Print information about the tests that
                                 // did not pass the threshold.

                                 for (K = 5; K <= 6; K++) {
                                    if ( RESULT( K ) >= THRESH ) {
                                       if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                                       WRITE( NOUT, FMT = 9997 )TRANS, M, N, NRHS, MB, NB, ITYPE, K, RESULT( K );
                                          NFAIL = NFAIL + 1;
                                    }
                                 }
                                 NRUN = NRUN + 2;
                              }
                           }
                        }
                     }
                  // =====================================================
                        // End test ZGELSTSLS
                  // =====================================================

                     // Generate a matrix of scaling type ISCALE and rank
                     // type IRANK.

                     zqrt15(ISCALE, IRANK, M, N, NRHS, COPYA, LDA, COPYB, LDB, COPYS, RANK, NORMA, NORMB, ISEED, WORK, LWORK );

                     // workspace used: max(M+min(M,N),NRHS*min(M,N),2*N+M)

                     LDWORK = max( 1, M );

                     // Loop for testing different block sizes.

                     for (INB = 1; INB <= NNB; INB++) { // 90
                        NB = NBVAL( INB );
                        xlaenv(1, NB );
                        xlaenv(3, NXVAL( INB ) );

                        // Test ZGELSY

                        // ZGELSY:  Compute the minimum-norm solution
                        // X to min( norm( A * X - B ) )
                        // using the rank-revealing orthogonal
                        // factorization.

                        zlacpy('Full', M, N, COPYA, LDA, A, LDA );
                        zlacpy('Full', M, NRHS, COPYB, LDB, B, LDB );

                        // Initialize vector IWORK.

                        for (J = 1; J <= N; J++) { // 70
                           IWORK( J ) = 0;
                        } // 70

                        SRNAMT = 'ZGELSY';
                        zgelsy(M, N, NRHS, A, LDA, B, LDB, IWORK, RCOND, CRANK, WORK, LWLSY, RWORK, INFO )                         IF( INFO != 0 ) CALL ALAERH( PATH, 'ZGELSY', INFO, 0, ' ', M, N, NRHS, -1, NB, ITYPE, NFAIL, NERRS, NOUT );

                        // workspace used: 2*MNMIN+NB*NB+NB*max(N,NRHS)

                        // Test 7:  Compute relative error in svd
                                 // workspace: M*N + 4*min(M,N) + max(M,N)

                        RESULT( 7 ) = ZQRT12( CRANK, CRANK, A, LDA, COPYS, WORK, LWORK, RWORK );

                        // Test 8:  Compute error in solution
                                 // workspace:  M*NRHS + M

                        zlacpy('Full', M, NRHS, COPYB, LDB, WORK, LDWORK );
                        zqrt16('No transpose', M, N, NRHS, COPYA, LDA, B, LDB, WORK, LDWORK, RWORK, RESULT( 8 ) );

                        // Test 9:  Check norm of r'*A
                                 // workspace: NRHS*(M+N)

                        RESULT( 9 ) = ZERO;
                        if (M > CRANK) RESULT( 9 ) = ZQRT17( 'No transpose', 1, M, N, NRHS, COPYA, LDA, B, LDB, COPYB, LDB, C, WORK, LWORK );

                        // Test 10:  Check if x is in the rowspace of A
                                 // workspace: (M+NRHS)*(N+2)

                        RESULT( 10 ) = ZERO;

                        if (N > CRANK) RESULT( 10 ) = ZQRT14( 'No transpose', M, N, NRHS, COPYA, LDA, B, LDB, WORK, LWORK );

                        // Test ZGELSS

                        // ZGELSS:  Compute the minimum-norm solution
                        // X to min( norm( A * X - B ) )
                        // using the SVD.

                        zlacpy('Full', M, N, COPYA, LDA, A, LDA );
                        zlacpy('Full', M, NRHS, COPYB, LDB, B, LDB );
                        SRNAMT = 'ZGELSS';
                        zgelss(M, N, NRHS, A, LDA, B, LDB, S, RCOND, CRANK, WORK, LWORK, RWORK, INFO );

                        if (INFO != 0) alaerh( PATH, 'ZGELSS', INFO, 0, ' ', M, N, NRHS, -1, NB, ITYPE, NFAIL, NERRS, NOUT );

                        // workspace used: 3*min(m,n) +
                                        // max(2*min(m,n),nrhs,max(m,n))

                        // Test 11:  Compute relative error in svd

                        if ( RANK > 0 ) {
                           daxpy(MNMIN, -ONE, COPYS, 1, S, 1 );
                           RESULT( 11 ) = DASUM( MNMIN, S, 1 ) / DASUM( MNMIN, COPYS, 1 ) / ( EPS*DBLE( MNMIN ) );
                        } else {
                           RESULT( 11 ) = ZERO;
                        }

                        // Test 12:  Compute error in solution

                        zlacpy('Full', M, NRHS, COPYB, LDB, WORK, LDWORK );
                        zqrt16('No transpose', M, N, NRHS, COPYA, LDA, B, LDB, WORK, LDWORK, RWORK, RESULT( 12 ) );

                        // Test 13:  Check norm of r'*A

                        RESULT( 13 ) = ZERO;
                        if (M > CRANK) RESULT( 13 ) = ZQRT17( 'No transpose', 1, M, N, NRHS, COPYA, LDA, B, LDB, COPYB, LDB, C, WORK, LWORK );

                        // Test 14:  Check if x is in the rowspace of A

                        RESULT( 14 ) = ZERO;
                        if (N > CRANK) RESULT( 14 ) = ZQRT14( 'No transpose', M, N, NRHS, COPYA, LDA, B, LDB, WORK, LWORK );

                        // Test ZGELSD

                        // ZGELSD:  Compute the minimum-norm solution X
                        // to min( norm( A * X - B ) ) using a
                        // divide and conquer SVD.

                        xlaenv(9, 25 );

                        zlacpy('Full', M, N, COPYA, LDA, A, LDA );
                        zlacpy('Full', M, NRHS, COPYB, LDB, B, LDB );

                        SRNAMT = 'ZGELSD';
                        zgelsd(M, N, NRHS, A, LDA, B, LDB, S, RCOND, CRANK, WORK, LWORK, RWORK, IWORK, INFO )                         IF( INFO != 0 ) CALL ALAERH( PATH, 'ZGELSD', INFO, 0, ' ', M, N, NRHS, -1, NB, ITYPE, NFAIL, NERRS, NOUT );

                        // Test 15:  Compute relative error in svd

                        if ( RANK > 0 ) {
                           daxpy(MNMIN, -ONE, COPYS, 1, S, 1 );
                           RESULT( 15 ) = DASUM( MNMIN, S, 1 ) / DASUM( MNMIN, COPYS, 1 ) / ( EPS*DBLE( MNMIN ) );
                        } else {
                           RESULT( 15 ) = ZERO;
                        }

                        // Test 16:  Compute error in solution

                        zlacpy('Full', M, NRHS, COPYB, LDB, WORK, LDWORK );
                        zqrt16('No transpose', M, N, NRHS, COPYA, LDA, B, LDB, WORK, LDWORK, RWORK, RESULT( 16 ) );

                        // Test 17:  Check norm of r'*A

                        RESULT( 17 ) = ZERO;
                        if (M > CRANK) RESULT( 17 ) = ZQRT17( 'No transpose', 1, M, N, NRHS, COPYA, LDA, B, LDB, COPYB, LDB, C, WORK, LWORK );

                        // Test 18:  Check if x is in the rowspace of A

                        RESULT( 18 ) = ZERO;
                        if (N > CRANK) RESULT( 18 ) = ZQRT14( 'No transpose', M, N, NRHS, COPYA, LDA, B, LDB, WORK, LWORK );

                        // Print information about the tests that did not
                        // pass the threshold.

                        for (K = 7; K <= 18; K++) { // 80
                           if ( RESULT( K ) >= THRESH ) {
                              if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                              WRITE( NOUT, FMT = 9998 )M, N, NRHS, NB, ITYPE, K, RESULT( K );
                              NFAIL = NFAIL + 1;
                           }
                        } // 80
                        NRUN = NRUN + 12;

                     } // 90
                  } // 100
               } // 110
            } // 120
         } // 130
      } // 140

      // Print a summary of the results.

      alasvm(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( ' TRANS=''', A1, ''', M=', I5, ', N=', I5, ', NRHS=', I4, ', NB=', I4, ', type', I2, ', test(', I2, ')=', G12.5 );
 9998 FORMAT( ' M=', I5, ', N=', I5, ', NRHS=', I4, ', NB=', I4, ', type', I2, ', test(', I2, ')=', G12.5 );
 9997 FORMAT( ' TRANS=''', A1,' M=', I5, ', N=', I5, ', NRHS=', I4, ', MB=', I4,', NB=', I4,', type', I2, ', test(', I2, ')=', G12.5 );

      DEALLOCATE( WORK );
      DEALLOCATE( IWORK );
      DEALLOCATE( RWORK );
      return;
      }
