      void cckgsv(NM, MVAL, PVAL, NVAL, NMATS, ISEED, THRESH, NMAX, A, AF, B, BF, U, V, Q, ALPHA, BETA, R, IWORK, WORK, RWORK, NIN, NOUT, INFO ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, NIN, NM, NMATS, NMAX, NOUT;
      double               THRESH;
      int                ISEED( 4 ), IWORK( * ), MVAL( * ), NVAL( * ), PVAL( * );
      double               ALPHA( * ), BETA( * ), RWORK( * );
      Complex            A( * ), AF( * ), B( * ), BF( * ), Q( * ), R( * ), U( * ), V( * ), WORK( * );
      // ..

      int                NTESTS;
      const              NTESTS = 12 ;
      int                NTYPES;
      const              NTYPES = 8 ;
      bool               FIRSTT;
      String             DISTA, DISTB, TYPE;
      String             PATH;
      int                I, IINFO, IM, IMAT, KLA, KLB, KUA, KUB, LDA, LDB, LDQ, LDR, LDU, LDV, LWORK, M, MODEA, MODEB, N, NFAIL, NRUN, NT, P, K, L;
      double               ANORM, BNORM, CNDNMA, CNDNMB;
      bool               DOTYPE( NTYPES );
      double               RESULT( NTESTS );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAHDG, ALAREQ, ALASUM, CLATMS, SLATB9, CGSVTS3
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

      // Specific cases

      // Test: https://github.com/Reference-LAPACK/lapack/issues/411#issue-608776973

      M = 6;
      P = 6;
      N = 6;
      A[1:M*N] = CMPLX(1.0, 0.0);
      B[1:M*N] = CMPLX(0.0, 0.0);
      B[1+0*M] = CMPLX(9e19, 0.0);
      B[2+1*M] = CMPLX(9e18, 0.0);
      B[3+2*M] = CMPLX(9e17, 0.0);
      B[4+3*M] = CMPLX(9e16, 0.0);
      B[5+4*M] = CMPLX(9e15, 0.0);
      B[6+5*M] = CMPLX(9e14, 0.0);
      cggsvd3('N','N','N', M, P, N, K, L, A, M, B, M, ALPHA, BETA, U, 1, V, 1, Q, 1, WORK, M*N, RWORK, IWORK, INFO);

      // Print information there is a NAN in BETA
      for (I = 1; I <= L; I++) { // 40
         if ( BETA(I) != BETA(I) ) {
            INFO = -I;
            break;
         }
      } // 40
      if ( INFO < 0 ) {
         if ( NFAIL == 0 && FIRSTT ) {
            FIRSTT = false;
            alahdg(NOUT, PATH );
         }
         WRITE( NOUT, FMT = 9997 ) -INFO;
         NFAIL = NFAIL + 1;
      }
      NRUN = NRUN + 1;
      INFO = 0;

      // Do for each value of M in MVAL.

      for (IM = 1; IM <= NM; IM++) { // 30
         M = MVAL( IM );
         P = PVAL( IM );
         N = NVAL( IM );

         for (IMAT = 1; IMAT <= NTYPES; IMAT++) { // 20

            // Do the tests only if DOTYPE( IMAT ) is true.

            if( !DOTYPE( IMAT ) ) GO TO 20;

            // Set up parameters with SLATB9 and generate test
            // matrices A and B with CLATMS.

            slatb9(PATH, IMAT, M, P, N, TYPE, KLA, KUA, KLB, KUB, ANORM, BNORM, MODEA, MODEB, CNDNMA, CNDNMB, DISTA, DISTB );

            // Generate M by N matrix A

            clatms(M, N, DISTA, ISEED, TYPE, RWORK, MODEA, CNDNMA, ANORM, KLA, KUA, 'No packing', A, LDA, WORK, IINFO );
            if ( IINFO != 0 ) {
               WRITE( NOUT, FMT = 9999 )IINFO;
               INFO = ( IINFO ).abs();
               GO TO 20;
            }

            // Generate P by N matrix B

            clatms(P, N, DISTB, ISEED, TYPE, RWORK, MODEB, CNDNMB, BNORM, KLB, KUB, 'No packing', B, LDB, WORK, IINFO );
            if ( IINFO != 0 ) {
               WRITE( NOUT, FMT = 9999 )IINFO;
               INFO = ( IINFO ).abs();
               GO TO 20;
            }

            NT = 6;

            cgsvts3(M, P, N, A, AF, LDA, B, BF, LDB, U, LDU, V, LDV, Q, LDQ, ALPHA, BETA, R, LDR, IWORK, WORK, LWORK, RWORK, RESULT );

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

 9999 FORMAT( ' CLATMS in CCKGSV   INFO = ${.i5}');
 9998 FORMAT( ' M=${.i4} P=${.i4}, N=${.i4}, type ${.i2}, test ${.i2}, ratio=${.g13_6};
 9997 FORMAT( ' FOUND NaN in BETA(', I4,')' );
      return;
      }
