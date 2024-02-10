      void cckgqr(NM, MVAL, NP, PVAL, NN, NVAL, NMATS, final Array<int> ISEED, THRESH, NMAX, A, AF, AQ, AR, TAUA, B, BF, BZ, BT, BWK, TAUB, final Array<double> _WORK, final Array<double> RWORK, NIN, NOUT, final Box<int> INFO ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, NIN, NM, NMATS, NMAX, NN, NOUT, NP;
      double               THRESH;
      int                ISEED( 4 ), MVAL( * ), NVAL( * ), PVAL( * );
      double               RWORK( * );
      Complex            A( * ), AF( * ), AQ( * ), AR( * ), B( * ), BF( * ), BT( * ), BWK( * ), BZ( * ), TAUA( * ), TAUB( * ), WORK( * );
      // ..

      int                NTESTS;
      const              NTESTS = 7 ;
      int                NTYPES;
      const              NTYPES = 8 ;
      bool               FIRSTT;
      String             DISTA, DISTB, TYPE;
      String             PATH;
      int                I, IINFO, IM, IMAT, IN, IP, KLA, KLB, KUA, KUB, LDA, LDB, LWORK, M, MODEA, MODEB, N, NFAIL, NRUN, NT, P;
      double               ANORM, BNORM, CNDNMA, CNDNMB;
      bool               DOTYPE( NTYPES );
      double               RESULT( NTESTS );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAHDG, ALAREQ, ALASUM, CGQRTS, CGRQTS, CLATMS, SLATB9
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS

      // Initialize constants.

      PATH[1: 3] = 'GQR';
      INFO = 0;
      NRUN = 0;
      NFAIL = 0;
      FIRSTT = true;
      alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );
      LDA = NMAX;
      LDB = NMAX;
      LWORK = NMAX*NMAX;

      // Do for each value of M in MVAL.

      for (IM = 1; IM <= NM; IM++) { // 60
         M = MVAL( IM );

         // Do for each value of P in PVAL.

         for (IP = 1; IP <= NP; IP++) { // 50
            P = PVAL( IP );

            // Do for each value of N in NVAL.

            for (IN = 1; IN <= NN; IN++) { // 40
               N = NVAL( IN );

               for (IMAT = 1; IMAT <= NTYPES; IMAT++) { // 30

                  // Do the tests only if DOTYPE( IMAT ) is true.

                  if( !DOTYPE( IMAT ) ) GO TO 30;

                  // Test CGGRQF

                  // Set up parameters with SLATB9 and generate test
                  // matrices A and B with CLATMS.

                  slatb9('GRQ', IMAT, M, P, N, TYPE, KLA, KUA, KLB, KUB, ANORM, BNORM, MODEA, MODEB, CNDNMA, CNDNMB, DISTA, DISTB );

                  clatms(M, N, DISTA, ISEED, TYPE, RWORK, MODEA, CNDNMA, ANORM, KLA, KUA, 'No packing', A, LDA, WORK, IINFO );
                  if ( IINFO != 0 ) {
                     WRITE( NOUT, FMT = 9999 )IINFO;
                     INFO = ( IINFO ).abs();
                     GO TO 30;
                  }

                  clatms(P, N, DISTB, ISEED, TYPE, RWORK, MODEB, CNDNMB, BNORM, KLB, KUB, 'No packing', B, LDB, WORK, IINFO );
                  if ( IINFO != 0 ) {
                     WRITE( NOUT, FMT = 9999 )IINFO;
                     INFO = ( IINFO ).abs();
                     GO TO 30;
                  }

                  NT = 4;

                  cgrqts(M, P, N, A, AF, AQ, AR, LDA, TAUA, B, BF, BZ, BT, BWK, LDB, TAUB, WORK, LWORK, RWORK, RESULT );

                  // Print information about the tests that did not
                  // pass the threshold.

                  for (I = 1; I <= NT; I++) { // 10
                     if ( RESULT( I ) >= THRESH ) {
                        if ( NFAIL == 0 && FIRSTT ) {
                           FIRSTT = false;
                           alahdg(NOUT, 'GRQ' );
                        }
                        WRITE( NOUT, FMT = 9998 )M, P, N, IMAT, I, RESULT( I );
                        NFAIL = NFAIL + 1;
                     }
                  } // 10
                  NRUN = NRUN + NT;

                  // Test CGGQRF

                  // Set up parameters with SLATB9 and generate test
                  // matrices A and B with CLATMS.

                  slatb9('GQR', IMAT, M, P, N, TYPE, KLA, KUA, KLB, KUB, ANORM, BNORM, MODEA, MODEB, CNDNMA, CNDNMB, DISTA, DISTB );

                  clatms(N, M, DISTA, ISEED, TYPE, RWORK, MODEA, CNDNMA, ANORM, KLA, KUA, 'No packing', A, LDA, WORK, IINFO );
                  if ( IINFO != 0 ) {
                     WRITE( NOUT, FMT = 9999 )IINFO;
                     INFO = ( IINFO ).abs();
                     GO TO 30;
                  }

                  clatms(N, P, DISTB, ISEED, TYPE, RWORK, MODEA, CNDNMA, BNORM, KLB, KUB, 'No packing', B, LDB, WORK, IINFO );
                  if ( IINFO != 0 ) {
                     WRITE( NOUT, FMT = 9999 )IINFO;
                     INFO = ( IINFO ).abs();
                     GO TO 30;
                  }

                  NT = 4;

                  cgqrts(N, M, P, A, AF, AQ, AR, LDA, TAUA, B, BF, BZ, BT, BWK, LDB, TAUB, WORK, LWORK, RWORK, RESULT );

                  // Print information about the tests that did not
                  // pass the threshold.

                  for (I = 1; I <= NT; I++) { // 20
                     if ( RESULT( I ) >= THRESH ) {
                        if ( NFAIL == 0 && FIRSTT ) {
                           FIRSTT = false;
                           alahdg(NOUT, PATH );
                        }
                        WRITE( NOUT, FMT = 9997 )N, M, P, IMAT, I, RESULT( I );
                        NFAIL = NFAIL + 1;
                     }
                  } // 20
                  NRUN = NRUN + NT;

               } // 30
            } // 40
         } // 50
      } // 60

      // Print a summary of the results.

      alasum(PATH, NOUT, NFAIL, NRUN, 0 );

 9999 FORMAT( ' CLATMS in CCKGQR:    INFO = ${.i5}');
 9998 FORMAT( ' M=${.i4} P=${.i4}, N=${.i4}, type ${.i2}, test ${.i2}, ratio=${.g13_6};
 9997 FORMAT( ' N=${.i4} M=${.i4}, P=${.i4}, type ${.i2}, test ${.i2}, ratio=${.g13_6};
      }
