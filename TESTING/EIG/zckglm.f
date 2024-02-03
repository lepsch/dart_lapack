      SUBROUTINE ZCKGLM( NN, NVAL, MVAL, PVAL, NMATS, ISEED, THRESH, NMAX, A, AF, B, BF, X, WORK, RWORK, NIN, NOUT, INFO )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, NIN, NMATS, NMAX, NN, NOUT;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      int                ISEED( 4 ), MVAL( * ), NVAL( * ), PVAL( * );
      double             RWORK( * );
      COMPLEX*16         A( * ), AF( * ), B( * ), BF( * ), WORK( * ), X( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      int                NTYPES;
      const              NTYPES = 8 ;
      // ..
      // .. Local Scalars ..
      bool               FIRSTT;
      String             DISTA, DISTB, TYPE;
      String             PATH;
      int                I, IINFO, IK, IMAT, KLA, KLB, KUA, KUB, LDA, LDB, LWORK, M, MODEA, MODEB, N, NFAIL, NRUN, P;
      double             ANORM, BNORM, CNDNMA, CNDNMB, RESID;
      // ..
      // .. Local Arrays ..
      bool               DOTYPE( NTYPES );
      // ..
      // .. External Functions ..
      COMPLEX*16         ZLARND
      // EXTERNAL ZLARND
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAHDG, ALAREQ, ALASUM, DLATB9, ZGLMTS, ZLATMS
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS
      // ..
      // .. Executable Statements ..

      // Initialize constants.

      PATH( 1: 3 ) = 'GLM'
      INFO = 0
      NRUN = 0
      NFAIL = 0
      FIRSTT = true;
      alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );
      LDA = NMAX
      LDB = NMAX
      LWORK = NMAX*NMAX

      // Check for valid input values.

      for (IK = 1; IK <= NN; IK++) { // 10
         M = MVAL( IK )
         P = PVAL( IK )
         N = NVAL( IK )
         if ( M.GT.N .OR. N.GT.M+P ) {
            if ( FIRSTT ) {
               WRITE( NOUT, FMT = * )
               FIRSTT = false;
            }
            WRITE( NOUT, FMT = 9997 )M, P, N
         }
      } // 10
      FIRSTT = true;

      // Do for each value of M in MVAL.

      for (IK = 1; IK <= NN; IK++) { // 40
         M = MVAL( IK )
         P = PVAL( IK )
         N = NVAL( IK )
         if (M.GT.N .OR. N.GT.M+P) GO TO 40;

         for (IMAT = 1; IMAT <= NTYPES; IMAT++) { // 30

            // Do the tests only if DOTYPE( IMAT ) is true.

            IF( .NOT.DOTYPE( IMAT ) ) GO TO 30

            // Set up parameters with DLATB9 and generate test
            // matrices A and B with ZLATMS.

            dlatb9(PATH, IMAT, M, P, N, TYPE, KLA, KUA, KLB, KUB, ANORM, BNORM, MODEA, MODEB, CNDNMA, CNDNMB, DISTA, DISTB );

            zlatms(N, M, DISTA, ISEED, TYPE, RWORK, MODEA, CNDNMA, ANORM, KLA, KUA, 'No packing', A, LDA, WORK, IINFO );
            if ( IINFO.NE.0 ) {
               WRITE( NOUT, FMT = 9999 )IINFO
               INFO = ABS( IINFO )
               GO TO 30
            }

            zlatms(N, P, DISTB, ISEED, TYPE, RWORK, MODEB, CNDNMB, BNORM, KLB, KUB, 'No packing', B, LDB, WORK, IINFO );
            if ( IINFO.NE.0 ) {
               WRITE( NOUT, FMT = 9999 )IINFO
               INFO = ABS( IINFO )
               GO TO 30
            }

            // Generate random left hand side vector of GLM

            for (I = 1; I <= N; I++) { // 20
               X( I ) = ZLARND( 2, ISEED )
            } // 20

            zglmts(N, M, P, A, AF, LDA, B, BF, LDB, X, X( NMAX+1 ), X( 2*NMAX+1 ), X( 3*NMAX+1 ), WORK, LWORK, RWORK, RESID );

            // Print information about the tests that did not
            // pass the threshold.

            if ( RESID.GE.THRESH ) {
               if ( NFAIL == 0 .AND. FIRSTT ) {
                  FIRSTT = false;
                  alahdg(NOUT, PATH );
               }
               WRITE( NOUT, FMT = 9998 )N, M, P, IMAT, 1, RESID
               NFAIL = NFAIL + 1
            }
            NRUN = NRUN + 1

         } // 30
      } // 40

      // Print a summary of the results.

      alasum(PATH, NOUT, NFAIL, NRUN, 0 );

 9999 FORMAT( ' ZLATMS in ZCKGLM INFO = ', I5 )
 9998 FORMAT( ' N=', I4, ' M=', I4, ', P=', I4, ', type ', I2, ', test ', I2, ', ratio=', G13.6 )
 9997 FORMAT( ' *** Invalid input  for GLM:  M = ', I6, ', P = ', I6, ', N = ', I6, ';', / '     must satisfy M <= N <= M+P  ', '(this set of values will be skipped)' )
      RETURN

      // End of ZCKGLM

      }
