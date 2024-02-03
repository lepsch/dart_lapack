      SUBROUTINE ZCKGQR( NM, MVAL, NP, PVAL, NN, NVAL, NMATS, ISEED, THRESH, NMAX, A, AF, AQ, AR, TAUA, B, BF, BZ, BT, BWK, TAUB, WORK, RWORK, NIN, NOUT, INFO )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, NIN, NM, NMATS, NMAX, NN, NOUT, NP;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      int                ISEED( 4 ), MVAL( * ), NVAL( * ), PVAL( * );
      double             RWORK( * );
      COMPLEX*16         A( * ), AF( * ), AQ( * ), AR( * ), B( * ), BF( * ), BT( * ), BWK( * ), BZ( * ), TAUA( * ), TAUB( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      int                NTESTS;
      const              NTESTS = 7 ;
      int                NTYPES;
      const              NTYPES = 8 ;
      // ..
      // .. Local Scalars ..
      bool               FIRSTT;
      String             DISTA, DISTB, TYPE;
      String             PATH;
      int                I, IINFO, IM, IMAT, IN, IP, KLA, KLB, KUA, KUB, LDA, LDB, LWORK, M, MODEA, MODEB, N, NFAIL, NRUN, NT, P;
      double             ANORM, BNORM, CNDNMA, CNDNMB;
      // ..
      // .. Local Arrays ..
      bool               DOTYPE( NTYPES );
      double             RESULT( NTESTS );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAHDG, ALAREQ, ALASUM, DLATB9, ZGQRTS, ZGRQTS, ZLATMS
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS
      // ..
      // .. Executable Statements ..

      // Initialize constants.

      PATH( 1: 3 ) = 'GQR'
      INFO = 0
      NRUN = 0
      NFAIL = 0
      FIRSTT = .TRUE.
      alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );
      LDA = NMAX
      LDB = NMAX
      LWORK = NMAX*NMAX

      // Do for each value of M in MVAL.

      DO 60 IM = 1, NM
         M = MVAL( IM )

         // Do for each value of P in PVAL.

         DO 50 IP = 1, NP
            P = PVAL( IP )

            // Do for each value of N in NVAL.

            DO 40 IN = 1, NN
               N = NVAL( IN )

               DO 30 IMAT = 1, NTYPES

                  // Do the tests only if DOTYPE( IMAT ) is true.

                  IF( .NOT.DOTYPE( IMAT ) ) GO TO 30

                  // Test ZGGRQF

                  // Set up parameters with DLATB9 and generate test
                  // matrices A and B with ZLATMS.

                  dlatb9('GRQ', IMAT, M, P, N, TYPE, KLA, KUA, KLB, KUB, ANORM, BNORM, MODEA, MODEB, CNDNMA, CNDNMB, DISTA, DISTB );

                  zlatms(M, N, DISTA, ISEED, TYPE, RWORK, MODEA, CNDNMA, ANORM, KLA, KUA, 'No packing', A, LDA, WORK, IINFO );
                  if ( IINFO.NE.0 ) {
                     WRITE( NOUT, FMT = 9999 )IINFO
                     INFO = ABS( IINFO )
                     GO TO 30
                  }

                  zlatms(P, N, DISTB, ISEED, TYPE, RWORK, MODEB, CNDNMB, BNORM, KLB, KUB, 'No packing', B, LDB, WORK, IINFO );
                  if ( IINFO.NE.0 ) {
                     WRITE( NOUT, FMT = 9999 )IINFO
                     INFO = ABS( IINFO )
                     GO TO 30
                  }

                  NT = 4

                  zgrqts(M, P, N, A, AF, AQ, AR, LDA, TAUA, B, BF, BZ, BT, BWK, LDB, TAUB, WORK, LWORK, RWORK, RESULT );

                  // Print information about the tests that did not
                  // pass the threshold.

                  DO 10 I = 1, NT
                     if ( RESULT( I ).GE.THRESH ) {
                        if ( NFAIL.EQ.0 .AND. FIRSTT ) {
                           FIRSTT = .FALSE.
                           alahdg(NOUT, 'GRQ' );
                        }
                        WRITE( NOUT, FMT = 9998 )M, P, N, IMAT, I, RESULT( I )
                        NFAIL = NFAIL + 1
                     }
   10             CONTINUE
                  NRUN = NRUN + NT

                  // Test ZGGQRF

                  // Set up parameters with DLATB9 and generate test
                  // matrices A and B with ZLATMS.

                  dlatb9('GQR', IMAT, M, P, N, TYPE, KLA, KUA, KLB, KUB, ANORM, BNORM, MODEA, MODEB, CNDNMA, CNDNMB, DISTA, DISTB );

                  zlatms(N, M, DISTA, ISEED, TYPE, RWORK, MODEA, CNDNMA, ANORM, KLA, KUA, 'No packing', A, LDA, WORK, IINFO );
                  if ( IINFO.NE.0 ) {
                     WRITE( NOUT, FMT = 9999 )IINFO
                     INFO = ABS( IINFO )
                     GO TO 30
                  }

                  zlatms(N, P, DISTB, ISEED, TYPE, RWORK, MODEA, CNDNMA, BNORM, KLB, KUB, 'No packing', B, LDB, WORK, IINFO );
                  if ( IINFO.NE.0 ) {
                     WRITE( NOUT, FMT = 9999 )IINFO
                     INFO = ABS( IINFO )
                     GO TO 30
                  }

                  NT = 4

                  zgqrts(N, M, P, A, AF, AQ, AR, LDA, TAUA, B, BF, BZ, BT, BWK, LDB, TAUB, WORK, LWORK, RWORK, RESULT );

                  // Print information about the tests that did not
                  // pass the threshold.

                  DO 20 I = 1, NT
                     if ( RESULT( I ).GE.THRESH ) {
                        if ( NFAIL.EQ.0 .AND. FIRSTT ) {
                           FIRSTT = .FALSE.
                           alahdg(NOUT, PATH );
                        }
                        WRITE( NOUT, FMT = 9997 )N, M, P, IMAT, I, RESULT( I )
                        NFAIL = NFAIL + 1
                     }
   20             CONTINUE
                  NRUN = NRUN + NT

   30          CONTINUE
   40       CONTINUE
   50    CONTINUE
   60 CONTINUE

      // Print a summary of the results.

      alasum(PATH, NOUT, NFAIL, NRUN, 0 );

 9999 FORMAT( ' ZLATMS in ZCKGQR:    INFO = ', I5 )
 9998 FORMAT( ' M=', I4, ' P=', I4, ', N=', I4, ', type ', I2, ', test ', I2, ', ratio=', G13.6 )
 9997 FORMAT( ' N=', I4, ' M=', I4, ', P=', I4, ', type ', I2, ', test ', I2, ', ratio=', G13.6 )
      RETURN

      // End of ZCKGQR

      }
