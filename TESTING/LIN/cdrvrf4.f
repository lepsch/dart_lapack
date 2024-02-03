      SUBROUTINE CDRVRF4( NOUT, NN, NVAL, THRESH, C1, C2, LDC, CRF, A, LDA, S_WORK_CLANGE );

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LDC, NN, NOUT;
      REAL               THRESH;
      // ..
      // .. Array Arguments ..
      int                NVAL( NN );
      REAL               S_WORK_CLANGE( * );
      COMPLEX            A( LDA, * ), C1( LDC, * ), C2( LDC, *), CRF( * );
      // ..

*  =====================================================================
      // ..
      // .. Parameters ..
      REAL               ZERO, ONE;
      const              ZERO = 0.0, ONE  = 1.0 ;
      int                NTESTS;
      const              NTESTS = 1 ;
      // ..
      // .. Local Scalars ..
      String             UPLO, CFORM, TRANS;
      int                I, IFORM, IIK, IIN, INFO, IUPLO, J, K, N, NFAIL, NRUN, IALPHA, ITRANS;
      REAL               ALPHA, BETA, EPS, NORMA, NORMC;
      // ..
      // .. Local Arrays ..
      String             UPLOS( 2 ), FORMS( 2 ), TRANSS( 2 );
      int                ISEED( 4 ), ISEEDY( 4 );
      REAL               RESULT( NTESTS );
      // ..
      // .. External Functions ..
      REAL               SLAMCH, SLARND, CLANGE;
      COMPLEX            CLARND;
      // EXTERNAL SLAMCH, SLARND, CLANGE, CLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL CHERK, CHFRK, CTFTTR, CTRTTF
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. Scalars in Common ..
      String             SRNAMT;
      // ..
      // .. Common blocks ..
      // COMMON / SRNAMC / SRNAMT
      // ..
      // .. Data statements ..
      DATA               ISEEDY / 1988, 1989, 1990, 1991 /;
      DATA               UPLOS  / 'U', 'L' /;
      DATA               FORMS  / 'N', 'C' /;
      DATA               TRANSS / 'N', 'C' /;
      // ..
      // .. Executable Statements ..

      // Initialize constants and the random number seed.

      NRUN = 0;
      NFAIL = 0;
      INFO = 0;
      for (I = 1; I <= 4; I++) { // 10
         ISEED( I ) = ISEEDY( I );
      } // 10
      EPS = SLAMCH( 'Precision' );

      for (IIN = 1; IIN <= NN; IIN++) { // 150

         N = NVAL( IIN );

         for (IIK = 1; IIK <= NN; IIK++) { // 140

            K = NVAL( IIN );

            for (IFORM = 1; IFORM <= 2; IFORM++) { // 130

               CFORM = FORMS( IFORM );

               for (IUPLO = 1; IUPLO <= 2; IUPLO++) { // 120

                  UPLO = UPLOS( IUPLO );

                  for (ITRANS = 1; ITRANS <= 2; ITRANS++) { // 110

                     TRANS = TRANSS( ITRANS );

                     for (IALPHA = 1; IALPHA <= 4; IALPHA++) { // 100

                        if ( IALPHA == 1) {
                           ALPHA = ZERO;
                           BETA = ZERO;
                        } else if ( IALPHA == 2) {
                           ALPHA = ONE;
                           BETA = ZERO;
                        } else if ( IALPHA == 3) {
                           ALPHA = ZERO;
                           BETA = ONE;
                        } else {
                           ALPHA = SLARND( 2, ISEED );
                           BETA = SLARND( 2, ISEED );
                        }

                        // All the parameters are set:
                           // CFORM, UPLO, TRANS, M, N,
                           // ALPHA, and BETA
                        // READY TO TEST!

                        NRUN = NRUN + 1;

                        if ( ITRANS == 1 ) {

                           // In this case we are NOTRANS, so A is N-by-K

                           for (J = 1; J <= K; J++) {
                              for (I = 1; I <= N; I++) {
                                 A( I, J) = CLARND( 4, ISEED );
                              }
                           }

                           NORMA = CLANGE( 'I', N, K, A, LDA, S_WORK_CLANGE );

                        } else {

                           // In this case we are TRANS, so A is K-by-N

                           for (J = 1; J <= N; J++) {
                              for (I = 1; I <= K; I++) {
                                 A( I, J) = CLARND( 4, ISEED );
                              }
                           }

                           NORMA = CLANGE( 'I', K, N, A, LDA, S_WORK_CLANGE );

                        }


                        // Generate C1 our N--by--N Hermitian matrix.
                        // Make sure C2 has the same upper/lower part,
                        // (the one that we do not touch), so
                        // copy the initial C1 in C2 in it.

                        for (J = 1; J <= N; J++) {
                           for (I = 1; I <= N; I++) {
                              C1( I, J) = CLARND( 4, ISEED );
                              C2(I,J) = C1(I,J);
                           }
                        }

                        // (See comment later on for why we use CLANGE and
                        // not CLANHE for C1.)

                        NORMC = CLANGE( 'I', N, N, C1, LDC, S_WORK_CLANGE );

                        SRNAMT = 'CTRTTF';
                        ctrttf(CFORM, UPLO, N, C1, LDC, CRF, INFO );

                        // call zherk the BLAS routine -> gives C1

                        SRNAMT = 'CHERK ';
                        cherk(UPLO, TRANS, N, K, ALPHA, A, LDA, BETA, C1, LDC );

                        // call zhfrk the RFP routine -> gives CRF

                        SRNAMT = 'CHFRK ';
                        chfrk(CFORM, UPLO, TRANS, N, K, ALPHA, A, LDA, BETA, CRF );

                        // convert CRF in full format -> gives C2

                        SRNAMT = 'CTFTTR';
                        ctfttr(CFORM, UPLO, N, CRF, C2, LDC, INFO );

                        // compare C1 and C2

                        for (J = 1; J <= N; J++) {
                           for (I = 1; I <= N; I++) {
                              C1(I,J) = C1(I,J)-C2(I,J);
                           }
                        }

                        // Yes, C1 is Hermitian so we could call CLANHE,
                        // but we want to check the upper part that is
                        // supposed to be unchanged and the diagonal that
                        // is supposed to be real -> CLANGE

                        RESULT(1) = CLANGE( 'I', N, N, C1, LDC, S_WORK_CLANGE )                         RESULT(1) = RESULT(1) / MAX( ABS( ALPHA ) * NORMA * NORMA + ABS( BETA ) * NORMC, ONE ) / MAX( N , 1 ) / EPS;

                        if ( RESULT(1) >= THRESH ) {
                           if ( NFAIL == 0 ) {
                              WRITE( NOUT, * );
                              WRITE( NOUT, FMT = 9999 );
                           }
                           WRITE( NOUT, FMT = 9997 ) 'CHFRK', CFORM, UPLO, TRANS, N, K, RESULT(1);
                           NFAIL = NFAIL + 1;
                        }

                     } // 100
                  } // 110
               } // 120
            } // 130
         } // 140
      } // 150

      // Print a summary of the results.

      if ( NFAIL == 0 ) {
         WRITE( NOUT, FMT = 9996 ) 'CHFRK', NRUN;
      } else {
         WRITE( NOUT, FMT = 9995 ) 'CHFRK', NFAIL, NRUN;
      }

 9999 FORMAT( 1X, ' *** Error(s) or Failure(s) while testing CHFRK ***');
 9997 FORMAT( 1X, '     Failure in ',A5,', CFORM=''',A1,''',', ' UPLO=''',A1,''',',' TRANS=''',A1,''',', ' N=',I3,', K =', I3, ', test=',G12.5);
 9996 FORMAT( 1X, 'All tests for ',A5,' auxiliary routine passed the ', 'threshold ( ',I5,' tests run)');
 9995 FORMAT( 1X, A6, ' auxiliary routine: ',I5,' out of ',I5, ' tests failed to pass the threshold');

      RETURN;

      // End of CDRVRF4

      }
