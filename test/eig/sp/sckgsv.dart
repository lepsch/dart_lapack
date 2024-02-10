      void sckgsv(NM, MVAL, PVAL, NVAL, NMATS, final Array<int> ISEED, THRESH, NMAX, A, AF, B, BF, U, V, Q, ALPHA, BETA, R, final Array<int> IWORK, WORK, RWORK, NIN, NOUT, final Box<int> INFO ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, NIN, NM, NMATS, NMAX, NOUT;
      double               THRESH;
      int                ISEED( 4 ), IWORK( * ), MVAL( * ), NVAL( * ), PVAL( * );
      double               A( * ), AF( * ), ALPHA( * ), B( * ), BETA( * ), BF( * ), Q( * ), R( * ), RWORK( * ), U( * ), V( * ), WORK( * );
      // ..

      int                NTESTS;
      const              NTESTS = 12 ;
      int                NTYPES;
      const              NTYPES = 8 ;
      bool               FIRSTT;
      String             DISTA, DISTB, TYPE;
      String             PATH;
      int                I, IINFO, IM, IMAT, KLA, KLB, KUA, KUB, LDA, LDB, LDQ, LDR, LDU, LDV, LWORK, M, MODEA, MODEB, N, NFAIL, NRUN, NT, P;
      double               ANORM, BNORM, CNDNMA, CNDNMB;
      bool               DOTYPE( NTYPES );
      double               RESULT( NTESTS );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAHDG, ALAREQ, ALASUM, SGSVTS3, SLATB9, SLATMS
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS

      // Initialize constants and the random number seed.

      PATH[1: 3] = 'GSV';
      INFO = 0;
      NRUN = 0;
      NFAIL = 0;
      FIRSTT = true;
      alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );
      LDA = NMAX;
      LDB = NMAX;
      LDU = NMAX;
      LDV = NMAX;
      LDQ = NMAX;
      LDR = NMAX;
      LWORK = NMAX*NMAX;

      // Do for each value of M in MVAL.

      for (IM = 1; IM <= NM; IM++) { // 30
         M = MVAL( IM );
         P = PVAL( IM );
         N = NVAL( IM );

         for (IMAT = 1; IMAT <= NTYPES; IMAT++) { // 20

            // Do the tests only if DOTYPE( IMAT ) is true.

            if( !DOTYPE( IMAT ) ) GO TO 20;

            // Set up parameters with SLATB9 and generate test
            // matrices A and B with SLATMS.

            slatb9(PATH, IMAT, M, P, N, TYPE, KLA, KUA, KLB, KUB, ANORM, BNORM, MODEA, MODEB, CNDNMA, CNDNMB, DISTA, DISTB );

            // Generate M by N matrix A

            slatms(M, N, DISTA, ISEED, TYPE, RWORK, MODEA, CNDNMA, ANORM, KLA, KUA, 'No packing', A, LDA, WORK, IINFO );
            if ( IINFO != 0 ) {
               WRITE( NOUT, FMT = 9999 )IINFO;
               INFO = ( IINFO ).abs();
               GO TO 20;
            }

            slatms(P, N, DISTB, ISEED, TYPE, RWORK, MODEB, CNDNMB, BNORM, KLB, KUB, 'No packing', B, LDB, WORK, IINFO );
            if ( IINFO != 0 ) {
               WRITE( NOUT, FMT = 9999 )IINFO;
               INFO = ( IINFO ).abs();
               GO TO 20;
            }

            NT = 6;

            sgsvts3(M, P, N, A, AF, LDA, B, BF, LDB, U, LDU, V, LDV, Q, LDQ, ALPHA, BETA, R, LDR, IWORK, WORK, LWORK, RWORK, RESULT );

            // Print information about the tests that did not
            // pass the threshold.

            for (I = 1; I <= NT; I++) { // 10
               if ( RESULT( I ) >= THRESH ) {
                  if ( NFAIL == 0 && FIRSTT ) {
                     FIRSTT = false;
                     alahdg(NOUT, PATH );
                  }
                  WRITE( NOUT, FMT = 9998 )M, P, N, IMAT, I, RESULT( I );
                  NFAIL = NFAIL + 1;
               }
            } // 10
            NRUN = NRUN + NT;
         } // 20
      } // 30

      // Print a summary of the results.

      alasum(PATH, NOUT, NFAIL, NRUN, 0 );

 9999 FORMAT( ' SLATMS in SCKGSV   INFO = ${.i5}');
 9998 FORMAT( ' M=${.i4} P=${.i4}, N=${.i4}, type ${.i2}, test ${.i2}, ratio=${.g13_6};
      }
