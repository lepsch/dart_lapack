      void sdrvls(DOTYPE, NM, MVAL, NN, NVAL, NNS, NSVAL, NNB, NBVAL, NXVAL, THRESH, TSTERR, A, COPYA, B, COPYB, C, S, COPYS, NOUT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               TSTERR;
      int                NM, NN, NNB, NNS, NOUT;
      double               THRESH;
      // ..
      // .. Array Arguments ..
      bool               DOTYPE( * );
      int                MVAL( * ), NBVAL( * ), NSVAL( * ), NVAL( * ), NXVAL( * );
      double               A( * ), B( * ), C( * ), COPYA( * ), COPYB( * ), COPYS( * ), S( * );
      // ..

// =====================================================================

      // .. Parameters ..
      int                NTESTS;
      const              NTESTS = 18 ;
      int                SMLSIZ;
      const              SMLSIZ = 25 ;
      double               ONE, TWO, ZERO;
      const              ONE = 1.0, TWO = 2.0, ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      String             TRANS;
      String             PATH;
      int                CRANK, I, IM, IMB, IN, INB, INFO, INS, IRANK, ISCALE, ITRAN, ITYPE, J, K, LDA, LDB, LDWORK, LWLSY, LWORK, M, MNMIN, N, NB, NCOLS, NERRS, NFAIL, NRHS, NROWS, NRUN, RANK, MB, MMAX, NMAX, NSMAX, LIWORK, LWORK_SGELS, LWORK_SGELST, LWORK_SGETSLS, LWORK_SGELSS, LWORK_SGELSY, LWORK_SGELSD;
      double               EPS, NORMA, NORMB, RCOND;
      // ..
      // .. Local Arrays ..
      int                ISEED( 4 ), ISEEDY( 4 ), IWQ( 1 );
      double               RESULT( NTESTS ), WQ( 1 );
      // ..
      // .. Allocatable Arrays ..
      double, ALLOCATABLE :: WORK (:);
      int    , ALLOCATABLE :: IWORK (:);
      // ..
      // .. External Functions ..
      //- REAL               SASUM, SLAMCH, SQRT12, SQRT14, SQRT17;
      // EXTERNAL SASUM, SLAMCH, SQRT12, SQRT14, SQRT17
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASVM, SAXPY, SERRLS, SGELS, SGELSD, SGELSS, SGELST, SGELSY, SGEMM, SGETSLS, SLACPY, SLARNV, SQRT13, SQRT15, SQRT16, SSCAL, XLAENV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC INT, MAX, MIN, REAL, SQRT
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
      // .. Data statements ..
      const ISEEDY = [ 1988, 1989, 1990, 1991 ];
      // ..
      // .. Executable Statements ..

      // Initialize constants and the random number seed.

      PATH[1: 1] = 'SINGLE PRECISION';
      PATH[2: 3] = 'LS';
      NRUN = 0;
      NFAIL = 0;
      NERRS = 0;
      for (I = 1; I <= 4; I++) { // 10
         ISEED[I] = ISEEDY( I );
      } // 10
      EPS = SLAMCH( 'Epsilon' );

      // Threshold for rank estimation

      RCOND = sqrt( EPS ) - ( sqrt( EPS )-EPS ) / 2;

      // Test the error exits

      xlaenv(2, 2 );
      xlaenv(9, SMLSIZ );
      if (TSTERR) serrls( PATH, NOUT );

      // Print the header if NM = 0 or NN = 0 and THRESH = 0.

      if( ( NM == 0 || NN == 0 ) && THRESH == ZERO ) alahd( NOUT, PATH );
      INFOT = 0;
      xlaenv(2, 2 );
      xlaenv(9, SMLSIZ );

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
      // SQRT14, SQRT17 (two side cases), SQRT15 and SQRT12

      LWORK = max( 1, ( M+N )*NRHS, ( N+NRHS )*( M+2 ), ( M+NRHS )*( N+2 ), max( M+MNMIN, NRHS*MNMIN,2*N+M ), max( M*N+4*MNMIN+max(M,N), M*N+2*MNMIN+4*N ) );
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
                                 TRANS = 'T';
                              }

                              // Compute workspace needed for SGELS
                              sgels(TRANS, M, N, NRHS, A, LDA, B, LDB, WQ( 1 ), -1, INFO );
                              LWORK_SGELS = INT ( WQ( 1 ) );
                              // Compute workspace needed for SGELST
                              sgelst(TRANS, M, N, NRHS, A, LDA, B, LDB, WQ, -1, INFO );
                              LWORK_SGELST = INT ( WQ ( 1 ) );
                              // Compute workspace needed for SGETSLS
                              sgetsls(TRANS, M, N, NRHS, A, LDA, B, LDB, WQ( 1 ), -1, INFO );
                              LWORK_SGETSLS = INT( WQ( 1 ) );
                           }
                        }
                        // Compute workspace needed for SGELSY
                        sgelsy(M, N, NRHS, A, LDA, B, LDB, IWQ, RCOND, CRANK, WQ, -1, INFO );
                        LWORK_SGELSY = INT( WQ( 1 ) );
                        // Compute workspace needed for SGELSS
                        sgelss(M, N, NRHS, A, LDA, B, LDB, S, RCOND, CRANK, WQ, -1 , INFO );
                        LWORK_SGELSS = INT( WQ( 1 ) );
                        // Compute workspace needed for SGELSD
                        sgelsd(M, N, NRHS, A, LDA, B, LDB, S, RCOND, CRANK, WQ, -1, IWQ, INFO );
                        LWORK_SGELSD = INT( WQ( 1 ) );
                        // Compute LIWORK workspace needed for SGELSY and SGELSD
                        LIWORK = max( LIWORK, N, IWQ( 1 ) );
                        // Compute LWORK workspace needed for all functions
                        LWORK = max( LWORK, LWORK_SGELS, LWORK_SGELST, LWORK_SGETSLS, LWORK_SGELSY, LWORK_SGELSS, LWORK_SGELSD );
                     }
                  }
               }
            }
         }
      }

      LWLSY = LWORK;

      ALLOCATE( WORK( LWORK ) );
      ALLOCATE( IWORK( LIWORK ) );

      for (IM = 1; IM <= NM; IM++) { // 150
         M = MVAL( IM );
         LDA = max( 1, M );

         for (IN = 1; IN <= NN; IN++) { // 140
            N = NVAL( IN );
            MNMIN = max(min( M, N ),1);
            LDB = max( 1, M, N );
            MB = (MNMIN+1);

            for (INS = 1; INS <= NNS; INS++) { // 130
               NRHS = NSVAL( INS );

               for (IRANK = 1; IRANK <= 2; IRANK++) { // 120
                  for (ISCALE = 1; ISCALE <= 3; ISCALE++) { // 110
                     ITYPE = ( IRANK-1 )*3 + ISCALE;
                     if( !DOTYPE( ITYPE ) ) GO TO 110;
                  // =====================================================
                        // Begin test SGELS
                  // =====================================================
                     if ( IRANK == 1 ) {

                        // Generate a matrix of scaling type ISCALE

                        sqrt13(ISCALE, M, N, COPYA, LDA, NORMA, ISEED );

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
                                 TRANS = 'T';
                                 NROWS = N;
                                 NCOLS = M;
                              }
                              LDWORK = max( 1, NCOLS );

                              // Set up a consistent rhs

                              if ( NCOLS > 0 ) {
                                 slarnv(2, ISEED, NCOLS*NRHS, WORK );
                                 sscal(NCOLS*NRHS, ONE / REAL( NCOLS ), WORK, 1 );
                              }
                              sgemm(TRANS, 'No transpose', NROWS, NRHS, NCOLS, ONE, COPYA, LDA, WORK, LDWORK, ZERO, B, LDB );
                              slacpy('Full', NROWS, NRHS, B, LDB, COPYB, LDB );

                              // Solve LS or overdetermined system

                              if ( M > 0 && N > 0 ) {
                                 slacpy('Full', M, N, COPYA, LDA, A, LDA );
                                 slacpy('Full', NROWS, NRHS, COPYB, LDB, B, LDB );
                              }
                             srnamc.SRNAMT = 'SGELS ';
                              sgels(TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )                               IF( INFO != 0 ) CALL ALAERH( PATH, 'SGELS ', INFO, 0, TRANS, M, N, NRHS, -1, NB, ITYPE, NFAIL, NERRS, NOUT );

                              // Test 1: Check correctness of results
                              // for SGELS, compute the residual:
                              // RESID = norm(B - A*X) /
                              // / ( max(m,n) * norm(A) * norm(X) * EPS )

                              if (NROWS > 0 && NRHS > 0) slacpy( 'Full', NROWS, NRHS, COPYB, LDB, C, LDB );
                              sqrt16(TRANS, M, N, NRHS, COPYA, LDA, B, LDB, C, LDB, WORK, RESULT( 1 ) );

                              // Test 2: Check correctness of results
                              // for SGELS.

                              if ( ( ITRAN == 1 && M >= N ) || ( ITRAN == 2 && M < N ) ) {

                                 // Solving LS system, compute:
                                 // r = norm((B- A*X)**T * A) /
                                  // / (norm(A)*norm(B)*max(M,N,NRHS)*EPS)

                                 RESULT[2] = SQRT17( TRANS, 1, M, N, NRHS, COPYA, LDA, B, LDB, COPYB, LDB, C, WORK, LWORK );
                              } else {

                                 // Solving overdetermined system

                                 RESULT[2] = SQRT14( TRANS, M, N, NRHS, COPYA, LDA, B, LDB, WORK, LWORK );
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
                        // End test SGELS
                  // =====================================================
                  // =====================================================
                        // Begin test SGELST
                  // =====================================================
                     if ( IRANK == 1 ) {

                        // Generate a matrix of scaling type ISCALE

                        sqrt13(ISCALE, M, N, COPYA, LDA, NORMA, ISEED );

                        // Loop for testing different block sizes.

                        for (INB = 1; INB <= NNB; INB++) {
                           NB = NBVAL( INB );
                           xlaenv(1, NB );

                           // Loop for testing non-transposed and transposed.

                           for (ITRAN = 1; ITRAN <= 2; ITRAN++) {
                              if ( ITRAN == 1 ) {
                                 TRANS = 'N';
                                 NROWS = M;
                                 NCOLS = N;
                              } else {
                                 TRANS = 'T';
                                 NROWS = N;
                                 NCOLS = M;
                              }
                              LDWORK = max( 1, NCOLS );

                              // Set up a consistent rhs

                              if ( NCOLS > 0 ) {
                                 slarnv(2, ISEED, NCOLS*NRHS, WORK );
                                 sscal(NCOLS*NRHS, ONE / REAL( NCOLS ), WORK, 1 );
                              }
                              sgemm(TRANS, 'No transpose', NROWS, NRHS, NCOLS, ONE, COPYA, LDA, WORK, LDWORK, ZERO, B, LDB );
                              slacpy('Full', NROWS, NRHS, B, LDB, COPYB, LDB );

                              // Solve LS or overdetermined system

                              if ( M > 0 && N > 0 ) {
                                 slacpy('Full', M, N, COPYA, LDA, A, LDA );
                                 slacpy('Full', NROWS, NRHS, COPYB, LDB, B, LDB );
                              }
                             srnamc.SRNAMT = 'SGELST';
                              sgelst(TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )                               IF( INFO != 0 ) CALL ALAERH( PATH, 'SGELST', INFO, 0, TRANS, M, N, NRHS, -1, NB, ITYPE, NFAIL, NERRS, NOUT );

                              // Test 3: Check correctness of results
                              // for SGELST, compute the residual:
                              // RESID = norm(B - A*X) /
                              // / ( max(m,n) * norm(A) * norm(X) * EPS )

                              if (NROWS > 0 && NRHS > 0) slacpy( 'Full', NROWS, NRHS, COPYB, LDB, C, LDB );
                              sqrt16(TRANS, M, N, NRHS, COPYA, LDA, B, LDB, C, LDB, WORK, RESULT( 3 ) );

                              // Test 4: Check correctness of results
                              // for SGELST.

                              if ( ( ITRAN == 1 && M >= N ) || ( ITRAN == 2 && M < N ) ) {

                                 // Solving LS system, compute:
                                 // r = norm((B- A*X)**T * A) /
                                  // / (norm(A)*norm(B)*max(M,N,NRHS)*EPS)

                                 RESULT[4] = SQRT17( TRANS, 1, M, N, NRHS, COPYA, LDA, B, LDB, COPYB, LDB, C, WORK, LWORK );
                              } else {

                                 // Solving overdetermined system

                                 RESULT[4] = SQRT14( TRANS, M, N, NRHS, COPYA, LDA, B, LDB, WORK, LWORK );
                              }

                              // Print information about the tests that
                              // did not pass the threshold.

                              for (K = 3; K <= 4; K++) {
                                 if ( RESULT( K ) >= THRESH ) {
                                    if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                                    WRITE( NOUT, FMT = 9999 ) TRANS, M, N, NRHS, NB, ITYPE, K, RESULT( K );
                                    NFAIL = NFAIL + 1;
                                 }
                              }
                              NRUN = NRUN + 2;
                           }
                        }
                     }
                  // =====================================================
                        // End test SGELST
                  // =====================================================
                  // =====================================================
                        // Begin test SGETSLS
                  // =====================================================
                     if ( IRANK == 1 ) {

                        // Generate a matrix of scaling type ISCALE

                        sqrt13(ISCALE, M, N, COPYA, LDA, NORMA, ISEED );

                        // Loop for testing different block sizes MB.

                        for (IMB = 1; IMB <= NNB; IMB++) {
                           MB = NBVAL( IMB );
                           xlaenv(1, MB );

                           // Loop for testing different block sizes NB.

                           for (INB = 1; INB <= NNB; INB++) {
                              NB = NBVAL( INB );
                              xlaenv(2, NB );

                              // Loop for testing non-transposed
                              // and transposed.

                              for (ITRAN = 1; ITRAN <= 2; ITRAN++) {
                                 if ( ITRAN == 1 ) {
                                    TRANS = 'N';
                                    NROWS = M;
                                    NCOLS = N;
                                 } else {
                                    TRANS = 'T';
                                    NROWS = N;
                                    NCOLS = M;
                                 }
                                 LDWORK = max( 1, NCOLS );

                                 // Set up a consistent rhs

                                 if ( NCOLS > 0 ) {
                                    slarnv(2, ISEED, NCOLS*NRHS, WORK );
                                    sscal(NCOLS*NRHS, ONE / REAL( NCOLS ), WORK, 1 );
                                 }
                                 sgemm(TRANS, 'No transpose', NROWS, NRHS, NCOLS, ONE, COPYA, LDA, WORK, LDWORK, ZERO, B, LDB );
                                 slacpy('Full', NROWS, NRHS, B, LDB, COPYB, LDB );

                                 // Solve LS or overdetermined system

                                 if ( M > 0 && N > 0 ) {
                                    slacpy('Full', M, N, COPYA, LDA, A, LDA );
                                    slacpy('Full', NROWS, NRHS, COPYB, LDB, B, LDB );
                                 }
                                srnamc.SRNAMT = 'SGETSLS';
                                 sgetsls(TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )                                  IF( INFO != 0 ) CALL ALAERH( PATH, 'SGETSLS', INFO, 0, TRANS, M, N, NRHS, -1, NB, ITYPE, NFAIL, NERRS, NOUT );

                              // Test 5: Check correctness of results
                              // for SGETSLS, compute the residual:
                              // RESID = norm(B - A*X) /
                              // / ( max(m,n) * norm(A) * norm(X) * EPS )

                                 if (NROWS > 0 && NRHS > 0) slacpy( 'Full', NROWS, NRHS, COPYB, LDB, C, LDB );
                                 sqrt16(TRANS, M, N, NRHS, COPYA, LDA, B, LDB, C, LDB, WORK, RESULT( 5 ) );

                              // Test 6: Check correctness of results
                              // for SGETSLS.

                                 if ( ( ITRAN == 1 && M >= N ) || ( ITRAN == 2 && M < N ) ) {

                                    // Solving LS system, compute:
                                    // r = norm((B- A*X)**T * A) /
                                  // / (norm(A)*norm(B)*max(M,N,NRHS)*EPS)

                                    RESULT[6] = SQRT17( TRANS, 1, M, N, NRHS, COPYA, LDA, B, LDB, COPYB, LDB, C, WORK, LWORK );
                                 } else {

                                    // Solving overdetermined system

                                    RESULT[6] = SQRT14( TRANS, M, N, NRHS, COPYA, LDA, B, LDB, WORK, LWORK );
                                 }

                                 // Print information about the tests that
                                 // did not pass the threshold.

                                 for (K = 5; K <= 6; K++) {
                                    if ( RESULT( K ) >= THRESH ) {
                                       if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                                       WRITE( NOUT, FMT = 9997 ) TRANS, M, N, NRHS, MB, NB, ITYPE, K, RESULT( K );
                                       NFAIL = NFAIL + 1;
                                    }
                                 }
                                 NRUN = NRUN + 2;
                              }
                           }
                        }
                     }
                  // =====================================================
                        // End test SGETSLS
                  // =====================================================

                     // Generate a matrix of scaling type ISCALE and rank
                     // type IRANK.

                     sqrt15(ISCALE, IRANK, M, N, NRHS, COPYA, LDA, COPYB, LDB, COPYS, RANK, NORMA, NORMB, ISEED, WORK, LWORK );

                     // workspace used: max(M+min(M,N),NRHS*min(M,N),2*N+M)

                     LDWORK = max( 1, M );

                     // Loop for testing different block sizes.

                     for (INB = 1; INB <= NNB; INB++) { // 100
                        NB = NBVAL( INB );
                        xlaenv(1, NB );
                        xlaenv(3, NXVAL( INB ) );

                        // Test SGELSY

                        // SGELSY:  Compute the minimum-norm solution X
                        // to min( norm( A * X - B ) )
                        // using the rank-revealing orthogonal
                        // factorization.

                        // Initialize vector IWORK.

                        for (J = 1; J <= N; J++) { // 70
                           IWORK[J] = 0;
                        } // 70

                        slacpy('Full', M, N, COPYA, LDA, A, LDA );
                        slacpy('Full', M, NRHS, COPYB, LDB, B, LDB );

                       srnamc.SRNAMT = 'SGELSY';
                        sgelsy(M, N, NRHS, A, LDA, B, LDB, IWORK, RCOND, CRANK, WORK, LWLSY, INFO )                         IF( INFO != 0 ) CALL ALAERH( PATH, 'SGELSY', INFO, 0, ' ', M, N, NRHS, -1, NB, ITYPE, NFAIL, NERRS, NOUT );

                        // Test 7:  Compute relative error in svd
                                 // workspace: M*N + 4*min(M,N) + max(M,N)

                        RESULT[7] = SQRT12( CRANK, CRANK, A, LDA, COPYS, WORK, LWORK );

                        // Test 8:  Compute error in solution
                                 // workspace:  M*NRHS + M

                        slacpy('Full', M, NRHS, COPYB, LDB, WORK, LDWORK );
                        sqrt16('No transpose', M, N, NRHS, COPYA, LDA, B, LDB, WORK, LDWORK, WORK( M*NRHS+1 ), RESULT( 8 ) );

                        // Test 9:  Check norm of r'*A
                                 // workspace: NRHS*(M+N)

                        RESULT[9] = ZERO;
                        if (M > CRANK) RESULT( 9 ) = SQRT17( 'No transpose', 1, M, N, NRHS, COPYA, LDA, B, LDB, COPYB, LDB, C, WORK, LWORK );

                        // Test 10:  Check if x is in the rowspace of A
                                 // workspace: (M+NRHS)*(N+2)

                        RESULT[10] = ZERO;

                        if (N > CRANK) RESULT( 10 ) = SQRT14( 'No transpose', M, N, NRHS, COPYA, LDA, B, LDB, WORK, LWORK );

                        // Test SGELSS

                        // SGELSS:  Compute the minimum-norm solution X
                        // to min( norm( A * X - B ) )
                        // using the SVD.

                        slacpy('Full', M, N, COPYA, LDA, A, LDA );
                        slacpy('Full', M, NRHS, COPYB, LDB, B, LDB );
                       srnamc.SRNAMT = 'SGELSS';
                        sgelss(M, N, NRHS, A, LDA, B, LDB, S, RCOND, CRANK, WORK, LWORK, INFO )                         IF( INFO != 0 ) CALL ALAERH( PATH, 'SGELSS', INFO, 0, ' ', M, N, NRHS, -1, NB, ITYPE, NFAIL, NERRS, NOUT );

                        // workspace used: 3*min(m,n) +
                                        // max(2*min(m,n),nrhs,max(m,n))

                        // Test 11:  Compute relative error in svd

                        if ( RANK > 0 ) {
                           saxpy(MNMIN, -ONE, COPYS, 1, S, 1 );
                           RESULT[11] = SASUM( MNMIN, S, 1 ) / SASUM( MNMIN, COPYS, 1 ) / ( EPS*REAL( MNMIN ) );
                        } else {
                           RESULT[11] = ZERO;
                        }

                        // Test 12:  Compute error in solution

                        slacpy('Full', M, NRHS, COPYB, LDB, WORK, LDWORK );
                        sqrt16('No transpose', M, N, NRHS, COPYA, LDA, B, LDB, WORK, LDWORK, WORK( M*NRHS+1 ), RESULT( 12 ) );

                        // Test 13:  Check norm of r'*A

                        RESULT[13] = ZERO;
                        if (M > CRANK) RESULT( 13 ) = SQRT17( 'No transpose', 1, M, N, NRHS, COPYA, LDA, B, LDB, COPYB, LDB, C, WORK, LWORK );

                        // Test 14:  Check if x is in the rowspace of A

                        RESULT[14] = ZERO;
                        if (N > CRANK) RESULT( 14 ) = SQRT14( 'No transpose', M, N, NRHS, COPYA, LDA, B, LDB, WORK, LWORK );

                        // Test SGELSD

                        // SGELSD:  Compute the minimum-norm solution X
                        // to min( norm( A * X - B ) ) using a
                        // divide and conquer SVD.

                        // Initialize vector IWORK.

                        for (J = 1; J <= N; J++) { // 80
                           IWORK[J] = 0;
                        } // 80

                        slacpy('Full', M, N, COPYA, LDA, A, LDA );
                        slacpy('Full', M, NRHS, COPYB, LDB, B, LDB );

                       srnamc.SRNAMT = 'SGELSD';
                        sgelsd(M, N, NRHS, A, LDA, B, LDB, S, RCOND, CRANK, WORK, LWORK, IWORK, INFO )                         IF( INFO != 0 ) CALL ALAERH( PATH, 'SGELSD', INFO, 0, ' ', M, N, NRHS, -1, NB, ITYPE, NFAIL, NERRS, NOUT );

                        // Test 15:  Compute relative error in svd

                        if ( RANK > 0 ) {
                           saxpy(MNMIN, -ONE, COPYS, 1, S, 1 );
                           RESULT[15] = SASUM( MNMIN, S, 1 ) / SASUM( MNMIN, COPYS, 1 ) / ( EPS*REAL( MNMIN ) );
                        } else {
                           RESULT[15] = ZERO;
                        }

                        // Test 16:  Compute error in solution

                        slacpy('Full', M, NRHS, COPYB, LDB, WORK, LDWORK );
                        sqrt16('No transpose', M, N, NRHS, COPYA, LDA, B, LDB, WORK, LDWORK, WORK( M*NRHS+1 ), RESULT( 16 ) );

                        // Test 17:  Check norm of r'*A

                        RESULT[17] = ZERO;
                        if (M > CRANK) RESULT( 17 ) = SQRT17( 'No transpose', 1, M, N, NRHS, COPYA, LDA, B, LDB, COPYB, LDB, C, WORK, LWORK );

                        // Test 18:  Check if x is in the rowspace of A

                        RESULT[18] = ZERO;
                        if (N > CRANK) RESULT( 18 ) = SQRT14( 'No transpose', M, N, NRHS, COPYA, LDA, B, LDB, WORK, LWORK );

                        // Print information about the tests that did not
                        // pass the threshold.

                        for (K = 7; K <= 18; K++) { // 90
                           if ( RESULT( K ) >= THRESH ) {
                              if (NFAIL == 0 && NERRS == 0) alahd( NOUT, PATH );
                              WRITE( NOUT, FMT = 9998 )M, N, NRHS, NB, ITYPE, K, RESULT( K );
                              NFAIL = NFAIL + 1;
                           }
                        } // 90
                        NRUN = NRUN + 12;

                     } // 100
                  } // 110
               } // 120
            } // 130
         } // 140
      } // 150

      // Print a summary of the results.

      alasvm(PATH, NOUT, NFAIL, NRUN, NERRS );

 9999 FORMAT( ' TRANS=''', A1, ''', M=', I5, ', N=', I5, ', NRHS=', I4, ', NB=', I4, ', type', I2, ', test(', I2, ')=', G12.5 );
 9998 FORMAT( ' M=', I5, ', N=', I5, ', NRHS=', I4, ', NB=', I4, ', type', I2, ', test(', I2, ')=', G12.5 );
 9997 FORMAT( ' TRANS=''', A1,' M=', I5, ', N=', I5, ', NRHS=', I4, ', MB=', I4,', NB=', I4,', type', I2, ', test(', I2, ')=', G12.5 );

      DEALLOCATE( WORK );
      DEALLOCATE( IWORK );
      return;
      }
