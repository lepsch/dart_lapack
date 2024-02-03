      SUBROUTINE DCKLSE( NN, MVAL, PVAL, NVAL, NMATS, ISEED, THRESH, NMAX, A, AF, B, BF, X, WORK, RWORK, NIN, NOUT, INFO )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, NIN, NMATS, NMAX, NN, NOUT;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      int                ISEED( 4 ), MVAL( * ), NVAL( * ), PVAL( * );
      double             A( * ), AF( * ), B( * ), BF( * ), RWORK( * ), WORK( * ), X( * );
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
      int                I, IINFO, IK, IMAT, KLA, KLB, KUA, KUB, LDA, LDB, LWORK, M, MODEA, MODEB, N, NFAIL, NRUN, NT, P;
      double             ANORM, BNORM, CNDNMA, CNDNMB;
      // ..
      // .. Local Arrays ..
      bool               DOTYPE( NTYPES );
      double             RESULT( NTESTS );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAHDG, ALAREQ, ALASUM, DLARHS, DLATB9, DLATMS, DLSETS
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX
      // ..
      // .. Executable Statements ..

      // Initialize constants and the random number seed.

      PATH( 1: 3 ) = 'LSE'
      INFO = 0
      NRUN = 0
      NFAIL = 0
      FIRSTT = .TRUE.
      alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );
      LDA = NMAX
      LDB = NMAX
      LWORK = NMAX*NMAX

      // Check for valid input values.

      for (IK = 1; IK <= NN; IK++) { // 10
         M = MVAL( IK )
         P = PVAL( IK )
         N = NVAL( IK )
         if ( P.GT.N .OR. N.GT.M+P ) {
            if ( FIRSTT ) {
               WRITE( NOUT, FMT = * )
               FIRSTT = .FALSE.
            }
            WRITE( NOUT, FMT = 9997 )M, P, N
         }
   10 CONTINUE
      FIRSTT = .TRUE.

      // Do for each value of M in MVAL.

      for (IK = 1; IK <= NN; IK++) { // 40
         M = MVAL( IK )
         P = PVAL( IK )
         N = NVAL( IK )
         IF( P.GT.N .OR. N.GT.M+P ) GO TO 40

         for (IMAT = 1; IMAT <= NTYPES; IMAT++) { // 30

            // Do the tests only if DOTYPE( IMAT ) is true.

            IF( .NOT.DOTYPE( IMAT ) ) GO TO 30

            // Set up parameters with DLATB9 and generate test
            // matrices A and B with DLATMS.

            dlatb9(PATH, IMAT, M, P, N, TYPE, KLA, KUA, KLB, KUB, ANORM, BNORM, MODEA, MODEB, CNDNMA, CNDNMB, DISTA, DISTB );

            dlatms(M, N, DISTA, ISEED, TYPE, RWORK, MODEA, CNDNMA, ANORM, KLA, KUA, 'No packing', A, LDA, WORK, IINFO );
            if ( IINFO.NE.0 ) {
               WRITE( NOUT, FMT = 9999 )IINFO
               INFO = ABS( IINFO )
               GO TO 30
            }

            dlatms(P, N, DISTB, ISEED, TYPE, RWORK, MODEB, CNDNMB, BNORM, KLB, KUB, 'No packing', B, LDB, WORK, IINFO );
            if ( IINFO.NE.0 ) {
               WRITE( NOUT, FMT = 9999 )IINFO
               INFO = ABS( IINFO )
               GO TO 30
            }

            // Generate the right-hand sides C and D for the LSE.

            dlarhs('DGE', 'New solution', 'Upper', 'N', M, N, MAX( M-1, 0 ), MAX( N-1, 0 ), 1, A, LDA, X( 4*NMAX+1 ), MAX( N, 1 ), X, MAX( M, 1 ), ISEED, IINFO );

            dlarhs('DGE', 'Computed', 'Upper', 'N', P, N, MAX( P-1, 0 ), MAX( N-1, 0 ), 1, B, LDB, X( 4*NMAX+1 ), MAX( N, 1 ), X( 2*NMAX+1 ), MAX( P, 1 ), ISEED, IINFO );

            NT = 2

            dlsets(M, P, N, A, AF, LDA, B, BF, LDB, X, X( NMAX+1 ), X( 2*NMAX+1 ), X( 3*NMAX+1 ), X( 4*NMAX+1 ), WORK, LWORK, RWORK, RESULT( 1 ) );

            // Print information about the tests that did not
            // pass the threshold.

            for (I = 1; I <= NT; I++) { // 20
               if ( RESULT( I ).GE.THRESH ) {
                  if ( NFAIL.EQ.0 .AND. FIRSTT ) {
                     FIRSTT = .FALSE.
                     alahdg(NOUT, PATH );
                  }
                  WRITE( NOUT, FMT = 9998 )M, P, N, IMAT, I, RESULT( I )
                  NFAIL = NFAIL + 1
               }
   20       CONTINUE
            NRUN = NRUN + NT

   30    CONTINUE
   40 CONTINUE

      // Print a summary of the results.

      alasum(PATH, NOUT, NFAIL, NRUN, 0 );

 9999 FORMAT( ' DLATMS in DCKLSE   INFO = ', I5 )
 9998 FORMAT( ' M=', I4, ' P=', I4, ', N=', I4, ', type ', I2, ', test ', I2, ', ratio=', G13.6 )
 9997 FORMAT( ' *** Invalid input  for LSE:  M = ', I6, ', P = ', I6, ', N = ', I6, ';', / '     must satisfy P <= N <= P+M  ', '(this set of values will be skipped)' )
      RETURN

      // End of DCKLSE

      }
