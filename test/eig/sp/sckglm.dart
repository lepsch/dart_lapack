      void sckglm(NN, MVAL, PVAL, NVAL, NMATS, final Array<int> ISEED, THRESH, NMAX, A, AF, B, BF, X, final Array<double> _WORK, final Array<double> RWORK, NIN, NOUT, final Box<int> INFO ) {

// -- LAPACK test routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      int                INFO, NIN, NMATS, NMAX, NN, NOUT;
      double               THRESH;
      int                ISEED( 4 ), MVAL( * ), NVAL( * ), PVAL( * );
      double               A( * ), AF( * ), B( * ), BF( * ), RWORK( * ), WORK( * ), X( * );
      // ..

      int                NTYPES;
      const              NTYPES = 8 ;
      bool               FIRSTT;
      String             DISTA, DISTB, TYPE;
      String             PATH;
      int                I, IINFO, IK, IMAT, KLA, KLB, KUA, KUB, LDA, LDB, LWORK, M, MODEA, MODEB, N, NFAIL, NRUN, P;
      double               ANORM, BNORM, CNDNMA, CNDNMB, RESID;
      bool               DOTYPE( NTYPES );
      // ..
      // .. External Functions ..
      //- REAL               SLARND;
      // EXTERNAL SLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAHDG, ALAREQ, ALASUM, SGLMTS, SLATB9, SLATMS
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS

      // Initialize constants.

      PATH[1: 3] = 'GLM';
      INFO = 0;
      NRUN = 0;
      NFAIL = 0;
      FIRSTT = true;
      alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );
      LDA = NMAX;
      LDB = NMAX;
      LWORK = NMAX*NMAX;

      // Check for valid input values.

      for (IK = 1; IK <= NN; IK++) { // 10
         M = MVAL( IK );
         P = PVAL( IK );
         N = NVAL( IK );
         if ( M > N || N > M+P ) {
            if ( FIRSTT ) {
               WRITE( NOUT, FMT = * );
               FIRSTT = false;
            }
            WRITE( NOUT, FMT = 9997 )M, P, N;
         }
      } // 10
      FIRSTT = true;

      // Do for each value of M in MVAL.

      for (IK = 1; IK <= NN; IK++) { // 40
         M = MVAL( IK );
         P = PVAL( IK );
         N = NVAL( IK );
         if (M > N || N > M+P) GO TO 40;

         for (IMAT = 1; IMAT <= NTYPES; IMAT++) { // 30

            // Do the tests only if DOTYPE( IMAT ) is true.

            if( !DOTYPE( IMAT ) ) GO TO 30;

            // Set up parameters with SLATB9 and generate test
            // matrices A and B with SLATMS.

            slatb9(PATH, IMAT, M, P, N, TYPE, KLA, KUA, KLB, KUB, ANORM, BNORM, MODEA, MODEB, CNDNMA, CNDNMB, DISTA, DISTB );

            slatms(N, M, DISTA, ISEED, TYPE, RWORK, MODEA, CNDNMA, ANORM, KLA, KUA, 'No packing', A, LDA, WORK, IINFO );
            if ( IINFO != 0 ) {
               WRITE( NOUT, FMT = 9999 )IINFO;
               INFO = ( IINFO ).abs();
               GO TO 30;
            }

            slatms(N, P, DISTB, ISEED, TYPE, RWORK, MODEB, CNDNMB, BNORM, KLB, KUB, 'No packing', B, LDB, WORK, IINFO );
            if ( IINFO != 0 ) {
               WRITE( NOUT, FMT = 9999 )IINFO;
               INFO = ( IINFO ).abs();
               GO TO 30;
            }

            // Generate random left hand side vector of GLM

            for (I = 1; I <= N; I++) { // 20
               X[I] = SLARND( 2, ISEED );
            } // 20

            sglmts(N, M, P, A, AF, LDA, B, BF, LDB, X, X( NMAX+1 ), X( 2*NMAX+1 ), X( 3*NMAX+1 ), WORK, LWORK, RWORK, RESID );

            // Print information about the tests that did not
            // pass the threshold.

            if ( RESID >= THRESH ) {
               if ( NFAIL == 0 && FIRSTT ) {
                  FIRSTT = false;
                  alahdg(NOUT, PATH );
               }
               WRITE( NOUT, FMT = 9998 )N, M, P, IMAT, 1, RESID;
               NFAIL = NFAIL + 1;
            }
            NRUN = NRUN + 1;

         } // 30
      } // 40

      // Print a summary of the results.

      alasum(PATH, NOUT, NFAIL, NRUN, 0 );

 9999 FORMAT( ' SLATMS in SCKGLM INFO = ${.i5}');
 9998 FORMAT( ' N=${.i4} M=${.i4}, P=${.i4}, type ${.i2}, test ${.i2}, ratio=${.g13_6};
 9997 FORMAT( ' *** Invalid input  for GLM:  M = ${.i6}, P = ${.i6}, N = ${.i6};\n     must satisfy M <= N <= M+P  (this set of values will be skipped)' )
      }
