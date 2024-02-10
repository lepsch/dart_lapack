import 'common.dart';

      void ddrvls(DOTYPE, NM, MVAL, NN, NVAL, NNS, NSVAL, NNB, NBVAL, NXVAL, THRESH, TSTERR, A, COPYA, B, COPYB, C, S, COPYS, NOUT ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      bool               TSTERR;
      int                NM, NN, NNB, NNS, NOUT;
      double             THRESH;
      bool               DOTYPE( * );
      int                MVAL( * ), NBVAL( * ), NSVAL( * ), NVAL( * ), NXVAL( * );
      double             A( * ), B( * ), C( * ), COPYA( * ), COPYB( * ), COPYS( * ), S( * );
      // ..

      int                NTESTS;
      const              NTESTS = 18 ;
      int                SMLSIZ;
      const              SMLSIZ = 25 ;
      double             ONE, TWO, ZERO;
      const              ONE = 1.0, TWO = 2.0, ZERO = 0.0 ;
      String             TRANS;
      String             PATH;
      int                CRANK, I, IM, IMB, IN, INB, INFO, INS, IRANK, ISCALE, ITRAN, ITYPE, J, K, LDA, LDB, LDWORK, LWLSY, LWORK, M, MNMIN, N, NB, NCOLS, NERRS, NFAIL, NRHS, NROWS, NRUN, RANK, MB, MMAX, NMAX, NSMAX, LIWORK, LWORK_DGELS, LWORK_DGELST, LWORK_DGETSLS, LWORK_DGELSS, LWORK_DGELSY, LWORK_DGELSD;
      double             EPS, NORMA, NORMB, RCOND;
      int                ISEED( 4 ), ISEEDY( 4 ), IWQ( 1 );
      double             RESULT( NTESTS ), WQ( 1 );
      // ..
      // .. Allocatable Arrays ..
      double          , ALLOCATABLE :: WORK (:);
      int    , ALLOCATABLE :: IWORK (:);
      // ..
      // .. External Functions ..
      //- double             DASUM, DLAMCH, DQRT12, DQRT14, DQRT17;
      // EXTERNAL DASUM, DLAMCH, DQRT12, DQRT14, DQRT17
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAERH, ALAHD, ALASVM, DAXPY, DERRLS, DGELS, DGELSD, DGELSS, DGELST, DGELSY, DGEMM, DGETSLS, DLACPY, DLARNV, DQRT13, DQRT15, DQRT16, DSCAL, XLAENV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, INT, MAX, MIN, SQRT
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
      // .. Data statements ..
      const ISEEDY = [ 1988, 1989, 1990, 1991 ];

      // Initialize constants and the random number seed.

      PATH = '${'Double precision'[0]}';
      PATH[2: 3] = 'LS';
      NRUN = 0;
      NFAIL = 0;
      NERRS = 0;
      for (I = 1; I <= 4; I++) { // 10
         ISEED[I] = ISEEDY( I );
      } // 10
      EPS = dlamch( 'Epsilon' );

      // Threshold for rank estimation

      RCOND = sqrt( EPS ) - ( sqrt( EPS )-EPS ) / 2;

      // Test the error exits

      xlaenv(2, 2 );
      xlaenv(9, SMLSIZ );
      if (TSTERR) derrls( PATH, NOUT );

      // Print the header if NM = 0 or NN = 0 and THRESH = 0.

      if( ( NM == 0 || NN == 0 ) && THRESH == ZERO ) alahd( NOUT, PATH );
      infoc.INFOT = 0;
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
      // DQRT14, DQRT17 (two side cases), DQRT15 and DQRT12

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

                              // Compute workspace needed for DGELS
                              dgels(TRANS, M, N, NRHS, A, LDA, B, LDB, WQ, -1, INFO );
                              LWORK_DGELS = INT ( WQ ( 1 ) );
                              // Compute workspace needed for DGELST
                              dgelst(TRANS, M, N, NRHS, A, LDA, B, LDB, WQ, -1, INFO );
                              LWORK_DGELST = INT ( WQ ( 1 ) );
                              // Compute workspace needed for DGETSLS
                              dgetsls(TRANS, M, N, NRHS, A, LDA, B, LDB, WQ, -1, INFO );
                              LWORK_DGETSLS = INT( WQ ( 1 ) );
                           }
                        }
                        // Compute workspace needed for DGELSY
                        dgelsy(M, N, NRHS, A, LDA, B, LDB, IWQ, RCOND, CRANK, WQ, -1, INFO );
                        LWORK_DGELSY = INT( WQ ( 1 ) );
                        // Compute workspace needed for DGELSS
                        dgelss(M, N, NRHS, A, LDA, B, LDB, S, RCOND, CRANK, WQ, -1 , INFO );
                        LWORK_DGELSS = INT( WQ ( 1 ) );
                        // Compute workspace needed for DGELSD
                        dgelsd(M, N, NRHS, A, LDA, B, LDB, S, RCOND, CRANK, WQ, -1, IWQ, INFO );
                        LWORK_DGELSD = INT( WQ ( 1 ) );
                        // Compute LIWORK workspace needed for DGELSY and DGELSD
                        LIWORK = max( LIWORK, N, IWQ( 1 ) );
                        // Compute LWORK workspace needed for all functions
                        LWORK = max( LWORK, LWORK_DGELS, LWORK_DGELST, LWORK_DGETSLS, LWORK_DGELSY, LWORK_DGELSS, LWORK_DGELSD );
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
                  //       Begin test DGELS
                  // =====================================================
                     if ( IRANK == 1 ) {

                        // Generate a matrix of scaling type ISCALE

                        dqrt13(ISCALE, M, N, COPYA, LDA, NORMA, ISEED );

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
                                 dlarnv(2, ISEED, NCOLS*NRHS, WORK );
                                 dscal(NCOLS*NRHS, ONE / NCOLS.toDouble(), WORK, 1 );
                              }
                              dgemm(TRANS, 'No transpose', NROWS, NRHS, NCOLS, ONE, COPYA, LDA, WORK, LDWORK, ZERO, B, LDB );
                              dlacpy('Full', NROWS, NRHS, B, LDB, COPYB, LDB );

                              // Solve LS or overdetermined system

                              if ( M > 0 && N > 0 ) {
                                 dlacpy('Full', M, N, COPYA, LDA, A, LDA );
                                 dlacpy('Full', NROWS, NRHS, COPYB, LDB, B, LDB );
                              }
                              srnamc.SRNAMT = 'DGELS ';
                              dgels(TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )                               IF( INFO != 0 ) CALL ALAERH( PATH, 'DGELS ', INFO, 0, TRANS, M, N, NRHS, -1, NB, ITYPE, NFAIL, NERRS, NOUT );

                              // Test 1: Check correctness of results
                              // for DGELS, compute the residual:
                              // RESID = norm(B - A*X) /
                              // / ( max(m,n) * norm(A) * norm(X) * EPS )

                              if (NROWS > 0 && NRHS > 0) dlacpy( 'Full', NROWS, NRHS, COPYB, LDB, C, LDB );
                              dqrt16(TRANS, M, N, NRHS, COPYA, LDA, B, LDB, C, LDB, WORK, RESULT( 1 ) );

                              // Test 2: Check correctness of results
                              // for DGELS.

                              if ( ( ITRAN == 1 && M >= N ) || ( ITRAN == 2 && M < N ) ) {

                                 // Solving LS system, compute:
                                 // r = norm((B- A*X)**T * A) /
                                 //  / (norm(A)*norm(B)*max(M,N,NRHS)*EPS)

                                 RESULT[2] = DQRT17( TRANS, 1, M, N, NRHS, COPYA, LDA, B, LDB, COPYB, LDB, C, WORK, LWORK );
                              } else {

                                 // Solving overdetermined system

                                 RESULT[2] = DQRT14( TRANS, M, N, NRHS, COPYA, LDA, B, LDB, WORK, LWORK );
                              }

                              // Print information about the tests that
                              // did not pass the threshold.

                              for (K = 1; K <= 2; K++) {
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
                  //       End test DGELS
                  // =====================================================
                  // =====================================================
                  //       Begin test DGELST
                  // =====================================================
                     if ( IRANK == 1 ) {

                        // Generate a matrix of scaling type ISCALE

                        dqrt13(ISCALE, M, N, COPYA, LDA, NORMA, ISEED );

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
                                 dlarnv(2, ISEED, NCOLS*NRHS, WORK );
                                 dscal(NCOLS*NRHS, ONE / NCOLS.toDouble(), WORK, 1 );
                              }
                              dgemm(TRANS, 'No transpose', NROWS, NRHS, NCOLS, ONE, COPYA, LDA, WORK, LDWORK, ZERO, B, LDB );
                              dlacpy('Full', NROWS, NRHS, B, LDB, COPYB, LDB );

                              // Solve LS or overdetermined system

                              if ( M > 0 && N > 0 ) {
                                 dlacpy('Full', M, N, COPYA, LDA, A, LDA );
                                 dlacpy('Full', NROWS, NRHS, COPYB, LDB, B, LDB );
                              }
                              srnamc.SRNAMT = 'DGELST';
                              dgelst(TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )                               IF( INFO != 0 ) CALL ALAERH( PATH, 'DGELST', INFO, 0, TRANS, M, N, NRHS, -1, NB, ITYPE, NFAIL, NERRS, NOUT );

                              // Test 3: Check correctness of results
                              // for DGELST, compute the residual:
                              // RESID = norm(B - A*X) /
                              // / ( max(m,n) * norm(A) * norm(X) * EPS )

                              if (NROWS > 0 && NRHS > 0) dlacpy( 'Full', NROWS, NRHS, COPYB, LDB, C, LDB );
                              dqrt16(TRANS, M, N, NRHS, COPYA, LDA, B, LDB, C, LDB, WORK, RESULT( 3 ) );

                              // Test 4: Check correctness of results
                              // for DGELST.

                              if ( ( ITRAN == 1 && M >= N ) || ( ITRAN == 2 && M < N ) ) {

                                 // Solving LS system, compute:
                                 // r = norm((B- A*X)**T * A) /
                                 //  / (norm(A)*norm(B)*max(M,N,NRHS)*EPS)

                                 RESULT[4] = DQRT17( TRANS, 1, M, N, NRHS, COPYA, LDA, B, LDB, COPYB, LDB, C, WORK, LWORK );
                              } else {

                                 // Solving overdetermined system

                                 RESULT[4] = DQRT14( TRANS, M, N, NRHS, COPYA, LDA, B, LDB, WORK, LWORK );
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
                  //       End test DGELST
                  // =====================================================
                  // =====================================================
                  //       Begin test DGETSLS
                  // =====================================================
                     if ( IRANK == 1 ) {

                        // Generate a matrix of scaling type ISCALE

                        dqrt13(ISCALE, M, N, COPYA, LDA, NORMA, ISEED );

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
                                    dlarnv(2, ISEED, NCOLS*NRHS, WORK );
                                    dscal(NCOLS*NRHS, ONE / NCOLS.toDouble(), WORK, 1 );
                                 }
                                 dgemm(TRANS, 'No transpose', NROWS, NRHS, NCOLS, ONE, COPYA, LDA, WORK, LDWORK, ZERO, B, LDB );
                                 dlacpy('Full', NROWS, NRHS, B, LDB, COPYB, LDB );

                                 // Solve LS or overdetermined system

                                 if ( M > 0 && N > 0 ) {
                                    dlacpy('Full', M, N, COPYA, LDA, A, LDA );
                                    dlacpy('Full', NROWS, NRHS, COPYB, LDB, B, LDB );
                                 }
                                 srnamc.SRNAMT = 'DGETSLS';
                                 dgetsls(TRANS, M, N, NRHS, A, LDA, B, LDB, WORK, LWORK, INFO )                                  IF( INFO != 0 ) CALL ALAERH( PATH, 'DGETSLS', INFO, 0, TRANS, M, N, NRHS, -1, NB, ITYPE, NFAIL, NERRS, NOUT );

                              // Test 5: Check correctness of results
                              // for DGETSLS, compute the residual:
                              // RESID = norm(B - A*X) /
                              // / ( max(m,n) * norm(A) * norm(X) * EPS )

                                 if (NROWS > 0 && NRHS > 0) dlacpy( 'Full', NROWS, NRHS, COPYB, LDB, C, LDB );
                                 dqrt16(TRANS, M, N, NRHS, COPYA, LDA, B, LDB, C, LDB, WORK, RESULT( 5 ) );

                              // Test 6: Check correctness of results
                              // for DGETSLS.

                                 if ( ( ITRAN == 1 && M >= N ) || ( ITRAN == 2 && M < N ) ) {

                                    // Solving LS system, compute:
                                    // r = norm((B- A*X)**T * A) /
                                  // / (norm(A)*norm(B)*max(M,N,NRHS)*EPS)

                                    RESULT[6] = DQRT17( TRANS, 1, M, N, NRHS, COPYA, LDA, B, LDB, COPYB, LDB, C, WORK, LWORK );
                                 } else {

                                    // Solving overdetermined system

                                    RESULT[6] = DQRT14( TRANS, M, N, NRHS, COPYA, LDA, B, LDB, WORK, LWORK );
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
                  //       End test DGETSLS
                  // =====================================================

                     // Generate a matrix of scaling type ISCALE and rank
                     // type IRANK.

                     dqrt15(ISCALE, IRANK, M, N, NRHS, COPYA, LDA, COPYB, LDB, COPYS, RANK, NORMA, NORMB, ISEED, WORK, LWORK );

                     // workspace used: max(M+min(M,N),NRHS*min(M,N),2*N+M)

                     LDWORK = max( 1, M );

                     // Loop for testing different block sizes.

                     for (INB = 1; INB <= NNB; INB++) { // 100
                        NB = NBVAL( INB );
                        xlaenv(1, NB );
                        xlaenv(3, NXVAL( INB ) );

                        // Test DGELSY

                        // DGELSY:  Compute the minimum-norm solution X
                        // to min( norm( A * X - B ) )
                        // using the rank-revealing orthogonal
                        // factorization.

                        // Initialize vector IWORK.

                        for (J = 1; J <= N; J++) { // 70
                           IWORK[J] = 0;
                        } // 70

                        dlacpy('Full', M, N, COPYA, LDA, A, LDA );
                        dlacpy('Full', M, NRHS, COPYB, LDB, B, LDB );

                        srnamc.SRNAMT = 'DGELSY';
                        dgelsy(M, N, NRHS, A, LDA, B, LDB, IWORK, RCOND, CRANK, WORK, LWLSY, INFO )                         IF( INFO != 0 ) CALL ALAERH( PATH, 'DGELSY', INFO, 0, ' ', M, N, NRHS, -1, NB, ITYPE, NFAIL, NERRS, NOUT );

                        // Test 7:  Compute relative error in svd
                        //          workspace: M*N + 4*min(M,N) + max(M,N)

                        RESULT[7] = DQRT12( CRANK, CRANK, A, LDA, COPYS, WORK, LWORK );

                        // Test 8:  Compute error in solution
                        //          workspace:  M*NRHS + M

                        dlacpy('Full', M, NRHS, COPYB, LDB, WORK, LDWORK );
                        dqrt16('No transpose', M, N, NRHS, COPYA, LDA, B, LDB, WORK, LDWORK, WORK( M*NRHS+1 ), RESULT( 8 ) );

                        // Test 9:  Check norm of r'*A
                        //          workspace: NRHS*(M+N)

                        RESULT[9] = ZERO;
                        if (M > CRANK) RESULT( 9 ) = DQRT17( 'No transpose', 1, M, N, NRHS, COPYA, LDA, B, LDB, COPYB, LDB, C, WORK, LWORK );

                        // Test 10:  Check if x is in the rowspace of A
                        //          workspace: (M+NRHS)*(N+2)

                        RESULT[10] = ZERO;

                        if (N > CRANK) RESULT( 10 ) = DQRT14( 'No transpose', M, N, NRHS, COPYA, LDA, B, LDB, WORK, LWORK );

                        // Test DGELSS

                        // DGELSS:  Compute the minimum-norm solution X
                        // to min( norm( A * X - B ) )
                        // using the SVD.

                        dlacpy('Full', M, N, COPYA, LDA, A, LDA );
                        dlacpy('Full', M, NRHS, COPYB, LDB, B, LDB );
                        srnamc.SRNAMT = 'DGELSS';
                        dgelss(M, N, NRHS, A, LDA, B, LDB, S, RCOND, CRANK, WORK, LWORK, INFO )                         IF( INFO != 0 ) CALL ALAERH( PATH, 'DGELSS', INFO, 0, ' ', M, N, NRHS, -1, NB, ITYPE, NFAIL, NERRS, NOUT );

                        // workspace used: 3*min(m,n) +
                        //                 max(2*min(m,n),nrhs,max(m,n))

                        // Test 11:  Compute relative error in svd

                        if ( RANK > 0 ) {
                           daxpy(MNMIN, -ONE, COPYS, 1, S, 1 );
                           RESULT[11] = dasum( MNMIN, S, 1 ) / dasum( MNMIN, COPYS, 1 ) / ( EPS*MNMIN.toDouble() );
                        } else {
                           RESULT[11] = ZERO;
                        }

                        // Test 12:  Compute error in solution

                        dlacpy('Full', M, NRHS, COPYB, LDB, WORK, LDWORK );
                        dqrt16('No transpose', M, N, NRHS, COPYA, LDA, B, LDB, WORK, LDWORK, WORK( M*NRHS+1 ), RESULT( 12 ) );

                        // Test 13:  Check norm of r'*A

                        RESULT[13] = ZERO;
                        if (M > CRANK) RESULT( 13 ) = DQRT17( 'No transpose', 1, M, N, NRHS, COPYA, LDA, B, LDB, COPYB, LDB, C, WORK, LWORK );

                        // Test 14:  Check if x is in the rowspace of A

                        RESULT[14] = ZERO;
                        if (N > CRANK) RESULT( 14 ) = DQRT14( 'No transpose', M, N, NRHS, COPYA, LDA, B, LDB, WORK, LWORK );

                        // Test DGELSD

                        // DGELSD:  Compute the minimum-norm solution X
                        // to min( norm( A * X - B ) ) using a
                        // divide and conquer SVD.

                        // Initialize vector IWORK.

                        for (J = 1; J <= N; J++) { // 80
                           IWORK[J] = 0;
                        } // 80

                        dlacpy('Full', M, N, COPYA, LDA, A, LDA );
                        dlacpy('Full', M, NRHS, COPYB, LDB, B, LDB );

                        srnamc.SRNAMT = 'DGELSD';
                        dgelsd(M, N, NRHS, A, LDA, B, LDB, S, RCOND, CRANK, WORK, LWORK, IWORK, INFO )                         IF( INFO != 0 ) CALL ALAERH( PATH, 'DGELSD', INFO, 0, ' ', M, N, NRHS, -1, NB, ITYPE, NFAIL, NERRS, NOUT );

                        // Test 15:  Compute relative error in svd

                        if ( RANK > 0 ) {
                           daxpy(MNMIN, -ONE, COPYS, 1, S, 1 );
                           RESULT[15] = dasum( MNMIN, S, 1 ) / dasum( MNMIN, COPYS, 1 ) / ( EPS*MNMIN.toDouble() );
                        } else {
                           RESULT[15] = ZERO;
                        }

                        // Test 16:  Compute error in solution

                        dlacpy('Full', M, NRHS, COPYB, LDB, WORK, LDWORK );
                        dqrt16('No transpose', M, N, NRHS, COPYA, LDA, B, LDB, WORK, LDWORK, WORK( M*NRHS+1 ), RESULT( 16 ) );

                        // Test 17:  Check norm of r'*A

                        RESULT[17] = ZERO;
                        if (M > CRANK) RESULT( 17 ) = DQRT17( 'No transpose', 1, M, N, NRHS, COPYA, LDA, B, LDB, COPYB, LDB, C, WORK, LWORK );

                        // Test 18:  Check if x is in the rowspace of A

                        RESULT[18] = ZERO;
                        if (N > CRANK) RESULT( 18 ) = DQRT14( 'No transpose', M, N, NRHS, COPYA, LDA, B, LDB, WORK, LWORK );

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

 9999 FORMAT( ' TRANS=''${.a1}'', M=${.i5}, N=${.i5}, NRHS=${.i4}, NB=${.i4}, type${.i2}, test(${.i2})=${.g12_5};
 9998 FORMAT( ' M=${.i5}, N=${.i5}, NRHS=${.i4}, NB=${.i4}, type${.i2}, test(${.i2})=${.g12_5};
 9997 FORMAT( ' TRANS=''${.a1} M=${.i5}, N=${.i5}, NRHS=${.i4}, MB=', I4,', NB=', I4,', type${.i2}, test(${.i2})=${.g12_5};

      DEALLOCATE( WORK );
      DEALLOCATE( IWORK );
      }
