      SUBROUTINE CDRVRF3( NOUT, NN, NVAL, THRESH, A, LDA, ARF, B1, B2, S_WORK_CLANGE, C_WORK_CGEQRF, TAU )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, NN, NOUT;
      REAL               THRESH
      // ..
      // .. Array Arguments ..
      int                NVAL( NN );
      REAL               S_WORK_CLANGE( * )
      COMPLEX            A( LDA, * ), ARF( * ), B1( LDA, * ), B2( LDA, * )
      COMPLEX            C_WORK_CGEQRF( * ), TAU( * )
      // ..

*  =====================================================================
      // ..
      // .. Parameters ..
      COMPLEX            ZERO, ONE
      const              ZERO = ( 0.0E+0, 0.0E+0 ) , ONE  = ( 1.0E+0, 0.0E+0 ) ;
      int                NTESTS;
      const              NTESTS = 1 ;
      // ..
      // .. Local Scalars ..
      String             UPLO, CFORM, DIAG, TRANS, SIDE;
      int                I, IFORM, IIM, IIN, INFO, IUPLO, J, M, N, NA, NFAIL, NRUN, ISIDE, IDIAG, IALPHA, ITRANS;
      COMPLEX            ALPHA
      REAL               EPS
      // ..
      // .. Local Arrays ..
      String             UPLOS( 2 ), FORMS( 2 ), TRANSS( 2 ), DIAGS( 2 ), SIDES( 2 );
      int                ISEED( 4 ), ISEEDY( 4 );
      REAL               RESULT( NTESTS )
      // ..
      // .. External Functions ..
      bool               LSAME;
      REAL               SLAMCH, CLANGE
      COMPLEX            CLARND
      // EXTERNAL SLAMCH, CLARND, CLANGE, LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL CTRTTF, CGEQRF, CGEQLF, CTFSM, CTRSM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SQRT
      // ..
      // .. Scalars in Common ..
      String             SRNAMT;
      // ..
      // .. Common blocks ..
      COMMON             / SRNAMC / SRNAMT
      // ..
      // .. Data statements ..
      DATA               ISEEDY / 1988, 1989, 1990, 1991 /
      DATA               UPLOS  / 'U', 'L' /
      DATA               FORMS  / 'N', 'C' /
      DATA               SIDES  / 'L', 'R' /
      DATA               TRANSS / 'N', 'C' /
      DATA               DIAGS  / 'N', 'U' /
      // ..
      // .. Executable Statements ..

      // Initialize constants and the random number seed.

      NRUN = 0
      NFAIL = 0
      INFO = 0
      for (I = 1; I <= 4; I++) { // 10
         ISEED( I ) = ISEEDY( I )
   10 CONTINUE
      EPS = SLAMCH( 'Precision' )

      for (IIM = 1; IIM <= NN; IIM++) { // 170

         M = NVAL( IIM )

         for (IIN = 1; IIN <= NN; IIN++) { // 160

            N = NVAL( IIN )

            for (IFORM = 1; IFORM <= 2; IFORM++) { // 150

               CFORM = FORMS( IFORM )

               for (IUPLO = 1; IUPLO <= 2; IUPLO++) { // 140

                  UPLO = UPLOS( IUPLO )

                  for (ISIDE = 1; ISIDE <= 2; ISIDE++) { // 130

                     SIDE = SIDES( ISIDE )

                     for (ITRANS = 1; ITRANS <= 2; ITRANS++) { // 120

                        TRANS = TRANSS( ITRANS )

                        for (IDIAG = 1; IDIAG <= 2; IDIAG++) { // 110

                           DIAG = DIAGS( IDIAG )

                           for (IALPHA = 1; IALPHA <= 3; IALPHA++) { // 100

                              if ( IALPHA.EQ.1 ) {
                                 ALPHA = ZERO
                              } else if ( IALPHA.EQ.2 ) {
                                 ALPHA = ONE
                              } else {
                                 ALPHA = CLARND( 4, ISEED )
                              }

                              // All the parameters are set:
                                 // CFORM, SIDE, UPLO, TRANS, DIAG, M, N,
                                 // and ALPHA
                              // READY TO TEST!

                              NRUN = NRUN + 1

                              if ( ISIDE.EQ.1 ) {

                                 // The case ISIDE.EQ.1 is when SIDE.EQ.'L'
                                 // -> A is M-by-M ( B is M-by-N )

                                 NA = M

                              } else {

                                 // The case ISIDE.EQ.2 is when SIDE.EQ.'R'
                                 // -> A is N-by-N ( B is M-by-N )

                                 NA = N

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
                                    A( I, J ) = CLARND( 4, ISEED )
                                 END DO
                              END DO

                              if ( IUPLO.EQ.1 ) {

                                 // The case IUPLO.EQ.1 is when SIDE.EQ.'U'
                                 // -> QR factorization.

                                 SRNAMT = 'CGEQRF'
                                 cgeqrf(NA, NA, A, LDA, TAU, C_WORK_CGEQRF, LDA, INFO );

                                 // Forcing main diagonal of test matrix to
                                 // be unit makes it ill-conditioned for
                                 // some test cases

                                 if ( LSAME( DIAG, 'U' ) ) {
                                    for (J = 1; J <= NA; J++) {
                                       for (I = 1; I <= J; I++) {
                                          A( I, J ) = A( I, J ) / ( 2.0 * A( J, J ) )
                                       END DO
                                    END DO
                                 }

                              } else {

                                 // The case IUPLO.EQ.2 is when SIDE.EQ.'L'
                                 // -> QL factorization.

                                 SRNAMT = 'CGELQF'
                                 cgelqf(NA, NA, A, LDA, TAU, C_WORK_CGEQRF, LDA, INFO );

                                 // Forcing main diagonal of test matrix to
                                 // be unit makes it ill-conditioned for
                                 // some test cases

                                 if ( LSAME( DIAG, 'U' ) ) {
                                    for (I = 1; I <= NA; I++) {
                                       for (J = 1; J <= I; J++) {
                                          A( I, J ) = A( I, J ) / ( 2.0 * A( I, I ) )
                                       END DO
                                    END DO
                                 }

                              }

                              // After the QR factorization, the diagonal
                              // of A is made of real numbers, we multiply
                              // by a random complex number of absolute
                              // value 1.0E+00.

                              for (J = 1; J <= NA; J++) {
                                 A( J, J ) = A( J, J ) * CLARND( 5, ISEED )
                              END DO

                              // Store a copy of A in RFP format (in ARF).

                              SRNAMT = 'CTRTTF'
                              ctrttf(CFORM, UPLO, NA, A, LDA, ARF, INFO );

                              // Generate B1 our M--by--N right-hand side
                              // and store a copy in B2.

                              for (J = 1; J <= N; J++) {
                                 for (I = 1; I <= M; I++) {
                                    B1( I, J ) = CLARND( 4, ISEED )
                                    B2( I, J ) = B1( I, J )
                                 END DO
                              END DO

                              // Solve op( A ) X = B or X op( A ) = B
                              // with CTRSM

                              SRNAMT = 'CTRSM'
                              ctrsm(SIDE, UPLO, TRANS, DIAG, M, N, ALPHA, A, LDA, B1, LDA );

                              // Solve op( A ) X = B or X op( A ) = B
                              // with CTFSM

                              SRNAMT = 'CTFSM'
                              ctfsm(CFORM, SIDE, UPLO, TRANS, DIAG, M, N, ALPHA, ARF, B2, LDA );

                              // Check that the result agrees.

                              for (J = 1; J <= N; J++) {
                                 for (I = 1; I <= M; I++) {
                                    B1( I, J ) = B2( I, J ) - B1( I, J )
                                 END DO
                              END DO

                              RESULT( 1 ) = CLANGE( 'I', M, N, B1, LDA, S_WORK_CLANGE )

                              RESULT( 1 ) = RESULT( 1 ) / SQRT( EPS ) / MAX ( MAX( M, N ), 1 )

                              if ( RESULT( 1 ).GE.THRESH ) {
                                 if ( NFAIL.EQ.0 ) {
                                    WRITE( NOUT, * )
                                    WRITE( NOUT, FMT = 9999 )
                                 }
                                 WRITE( NOUT, FMT = 9997 ) 'CTFSM', CFORM, SIDE, UPLO, TRANS, DIAG, M, N, RESULT( 1 )
                                 NFAIL = NFAIL + 1
                              }

  100                      CONTINUE
  110                   CONTINUE
  120                CONTINUE
  130             CONTINUE
  140          CONTINUE
  150       CONTINUE
  160    CONTINUE
  170 CONTINUE

      // Print a summary of the results.

      if ( NFAIL.EQ.0 ) {
         WRITE( NOUT, FMT = 9996 ) 'CTFSM', NRUN
      } else {
         WRITE( NOUT, FMT = 9995 ) 'CTFSM', NFAIL, NRUN
      }

 9999 FORMAT( 1X, ' *** Error(s) or Failure(s) while testing CTFSM ***')
 9997 FORMAT( 1X, '     Failure in ',A5,', CFORM=''',A1,''',', ' SIDE=''',A1,''',',' UPLO=''',A1,''',',' TRANS=''',A1,''',', ' DIAG=''',A1,''',',' M=',I3,', N =', I3,', test=',G12.5)
 9996 FORMAT( 1X, 'All tests for ',A5,' auxiliary routine passed the ', 'threshold ( ',I5,' tests run)')
 9995 FORMAT( 1X, A6, ' auxiliary routine: ',I5,' out of ',I5, ' tests failed to pass the threshold')

      RETURN

      // End of CDRVRF3

      }
