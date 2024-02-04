      void ddrvrf3(NOUT, NN, NVAL, THRESH, A, LDA, ARF, B1, B2, D_WORK_DLANGE, D_WORK_DGEQRF, TAU ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, NN, NOUT;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      int                NVAL( NN );
      double             A( LDA, * ), ARF( * ), B1( LDA, * ), B2( LDA, * ), D_WORK_DGEQRF( * ), D_WORK_DLANGE( * ), TAU( * );
      // ..

// =====================================================================
      // ..
      // .. Parameters ..
      double             ZERO, ONE;
      const              ZERO = ( 0.0, 0.0 ) , ONE  = ( 1.0, 0.0 ) ;
      int                NTESTS;
      const              NTESTS = 1 ;
      // ..
      // .. Local Scalars ..
      String             UPLO, CFORM, DIAG, TRANS, SIDE;
      int                I, IFORM, IIM, IIN, INFO, IUPLO, J, M, N, NA, NFAIL, NRUN, ISIDE, IDIAG, IALPHA, ITRANS;
      double             EPS, ALPHA;
      // ..
      // .. Local Arrays ..
      String             UPLOS( 2 ), FORMS( 2 ), TRANSS( 2 ), DIAGS( 2 ), SIDES( 2 );
      int                ISEED( 4 ), ISEEDY( 4 );
      double             RESULT( NTESTS );
      // ..
      // .. External Functions ..
      //- bool               lsame;
      //- double             DLAMCH, DLANGE, DLARND;
      // EXTERNAL DLAMCH, DLANGE, DLARND, lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL DTRTTF, DGEQRF, DGEQLF, DTFSM, DTRSM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SQRT
      // ..
      // .. Scalars in Common ..
      String             SRNAMT;
      // ..
      // .. Common blocks ..
      // COMMON / SRNAMC / SRNAMT
      // ..
      // .. Data statements ..
      const ISEEDY = [ 1988, 1989, 1990, 1991 ];
      const UPLOS = [ 'U', 'L' ];
      const FORMS = [ 'N', 'T' ];
      const SIDES = [ 'L', 'R' ];
      const TRANSS = [ 'N', 'T' ];
      const DIAGS = [ 'N', 'U' ];
      // ..
      // .. Executable Statements ..

      // Initialize constants and the random number seed.

      NRUN = 0;
      NFAIL = 0;
      INFO = 0;
      for (I = 1; I <= 4; I++) { // 10
         ISEED[I] = ISEEDY( I );
      } // 10
      EPS = DLAMCH( 'Precision' );

      for (IIM = 1; IIM <= NN; IIM++) { // 170

         M = NVAL( IIM );

         for (IIN = 1; IIN <= NN; IIN++) { // 160

            N = NVAL( IIN );

            for (IFORM = 1; IFORM <= 2; IFORM++) { // 150

               CFORM = FORMS( IFORM );

               for (IUPLO = 1; IUPLO <= 2; IUPLO++) { // 140

                  UPLO = UPLOS( IUPLO );

                  for (ISIDE = 1; ISIDE <= 2; ISIDE++) { // 130

                     SIDE = SIDES( ISIDE );

                     for (ITRANS = 1; ITRANS <= 2; ITRANS++) { // 120

                        TRANS = TRANSS( ITRANS );

                        for (IDIAG = 1; IDIAG <= 2; IDIAG++) { // 110

                           DIAG = DIAGS( IDIAG );

                           for (IALPHA = 1; IALPHA <= 3; IALPHA++) { // 100

                              if ( IALPHA == 1 ) {
                                 ALPHA = ZERO;
                              } else if ( IALPHA == 2 ) {
                                 ALPHA = ONE;
                              } else {
                                 ALPHA = DLARND( 2, ISEED );
                              }

                              // All the parameters are set:
                                 // CFORM, SIDE, UPLO, TRANS, DIAG, M, N,
                                 // and ALPHA
                              // READY TO TEST!

                              NRUN = NRUN + 1;

                              if ( ISIDE == 1 ) {

                                 // The case ISIDE == 1 is when SIDE == 'L'
                                 // -> A is M-by-M ( B is M-by-N )

                                 NA = M;

                              } else {

                                 // The case ISIDE == 2 is when SIDE == 'R'
                                 // -> A is N-by-N ( B is M-by-N )

                                 NA = N;

                              }

                              // Generate A our NA--by--NA triangular
                              // matrix.
                              // Our test is based on forward error so we
                              // do want A to be well conditioned! To get
                              // a well-conditioned triangular matrix, we
                              // take the R factor of the QR/LQ factorization
                              // of a random matrix.

                              for (J = 1; J <= NA; J++) {
                                 for (I = 1; I <= NA; I++) {
                                    A[I, J] = DLARND( 2, ISEED );
                                 }
                              }

                              if ( IUPLO == 1 ) {

                                 // The case IUPLO == 1 is when SIDE == 'U'
                                 // -> QR factorization.

                                 SRNAMT = 'DGEQRF';
                                 dgeqrf(NA, NA, A, LDA, TAU, D_WORK_DGEQRF, LDA, INFO );

                                 // Forcing main diagonal of test matrix to
                                 // be unit makes it ill-conditioned for
                                 // some test cases

                                 if ( lsame( DIAG, 'U' ) ) {
                                    for (J = 1; J <= NA; J++) {
                                       for (I = 1; I <= J; I++) {
                                          A[I, J] = A( I, J ) / ( 2.0 * A( J, J ) );
                                       }
                                    }
                                 }

                              } else {

                                 // The case IUPLO == 2 is when SIDE == 'L'
                                 // -> QL factorization.

                                 SRNAMT = 'DGELQF';
                                 dgelqf(NA, NA, A, LDA, TAU, D_WORK_DGEQRF, LDA, INFO );

                                 // Forcing main diagonal of test matrix to
                                 // be unit makes it ill-conditioned for
                                 // some test cases

                                 if ( lsame( DIAG, 'U' ) ) {
                                    for (I = 1; I <= NA; I++) {
                                       for (J = 1; J <= I; J++) {
                                          A[I, J] = A( I, J ) / ( 2.0 * A( I, I ) );
                                       }
                                    }
                                 }

                              }

                              // Store a copy of A in RFP format (in ARF).

                              SRNAMT = 'DTRTTF';
                              dtrttf(CFORM, UPLO, NA, A, LDA, ARF, INFO );

                              // Generate B1 our M--by--N right-hand side
                              // and store a copy in B2.

                              for (J = 1; J <= N; J++) {
                                 for (I = 1; I <= M; I++) {
                                    B1[I, J] = DLARND( 2, ISEED );
                                    B2[I, J] = B1( I, J );
                                 }
                              }

                              // Solve op( A ) X = B or X op( A ) = B
                              // with DTRSM

                              SRNAMT = 'DTRSM';
                              dtrsm(SIDE, UPLO, TRANS, DIAG, M, N, ALPHA, A, LDA, B1, LDA );

                              // Solve op( A ) X = B or X op( A ) = B
                              // with DTFSM

                              SRNAMT = 'DTFSM';
                              dtfsm(CFORM, SIDE, UPLO, TRANS, DIAG, M, N, ALPHA, ARF, B2, LDA );

                              // Check that the result agrees.

                              for (J = 1; J <= N; J++) {
                                 for (I = 1; I <= M; I++) {
                                    B1[I, J] = B2( I, J ) - B1( I, J );
                                 }
                              }

                              RESULT[1] = DLANGE( 'I', M, N, B1, LDA, D_WORK_DLANGE );

                              RESULT[1] = RESULT( 1 ) / sqrt( EPS ) / max( max( M, N ), 1 );

                              if ( RESULT( 1 ) >= THRESH ) {
                                 if ( NFAIL == 0 ) {
                                    WRITE( NOUT, * );
                                    WRITE( NOUT, FMT = 9999 );
                                 }
                                 WRITE( NOUT, FMT = 9997 ) 'DTFSM', CFORM, SIDE, UPLO, TRANS, DIAG, M, N, RESULT( 1 );
                                 NFAIL = NFAIL + 1;
                              }

                           } // 100
                        } // 110
                     } // 120
                  } // 130
               } // 140
            } // 150
         } // 160
      } // 170

      // Print a summary of the results.

      if ( NFAIL == 0 ) {
         WRITE( NOUT, FMT = 9996 ) 'DTFSM', NRUN;
      } else {
         WRITE( NOUT, FMT = 9995 ) 'DTFSM', NFAIL, NRUN;
      }

 9999 FORMAT( 1X, ' *** Error(s) or Failure(s) while testing DTFSM ***');
 9997 FORMAT( 1X, '     Failure in ',A5,', CFORM=''',A1,''',', ' SIDE=''',A1,''',',' UPLO=''',A1,''',',' TRANS=''',A1,''',', ' DIAG=''',A1,''',',' M=',I3,', N =', I3,', test=',G12.5);
 9996 FORMAT( 1X, 'All tests for ',A5,' auxiliary routine passed the ', 'threshold ( ',I5,' tests run)');
 9995 FORMAT( 1X, A6, ' auxiliary routine: ',I5,' out of ',I5, ' tests failed to pass the threshold');

      return;
      }