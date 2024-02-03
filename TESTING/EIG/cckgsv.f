      SUBROUTINE CCKGSV( NM, MVAL, PVAL, NVAL, NMATS, ISEED, THRESH, NMAX, A, AF, B, BF, U, V, Q, ALPHA, BETA, R, IWORK, WORK, RWORK, NIN, NOUT, INFO )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, NIN, NM, NMATS, NMAX, NOUT;
      REAL               THRESH
      // ..
      // .. Array Arguments ..
      int                ISEED( 4 ), IWORK( * ), MVAL( * ), NVAL( * ), PVAL( * );
      REAL               ALPHA( * ), BETA( * ), RWORK( * )
      COMPLEX            A( * ), AF( * ), B( * ), BF( * ), Q( * ), R( * ), U( * ), V( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      int                NTESTS;
      const              NTESTS = 12 ;
      int                NTYPES;
      const              NTYPES = 8 ;
      // ..
      // .. Local Scalars ..
      bool               FIRSTT;
      String             DISTA, DISTB, TYPE;
      String             PATH;
      int                I, IINFO, IM, IMAT, KLA, KLB, KUA, KUB, LDA, LDB, LDQ, LDR, LDU, LDV, LWORK, M, MODEA, MODEB, N, NFAIL, NRUN, NT, P, K, L;
      REAL               ANORM, BNORM, CNDNMA, CNDNMB
      // ..
      // .. Local Arrays ..
      bool               DOTYPE( NTYPES );
      REAL               RESULT( NTESTS )
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAHDG, ALAREQ, ALASUM, CLATMS, SLATB9, CGSVTS3
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS
      // ..
      // .. Executable Statements ..

      // Initialize constants and the random number seed.

      PATH( 1: 3 ) = 'GSV'
      INFO = 0
      NRUN = 0
      NFAIL = 0
      FIRSTT = .TRUE.
      CALL ALAREQ( PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT )
      LDA = NMAX
      LDB = NMAX
      LDU = NMAX
      LDV = NMAX
      LDQ = NMAX
      LDR = NMAX
      LWORK = NMAX*NMAX

      // Specific cases

      // Test: https://github.com/Reference-LAPACK/lapack/issues/411#issue-608776973

      M = 6
      P = 6
      N = 6
      A(1:M*N) = CMPLX(1.E0, 0.E0)
      B(1:M*N) = CMPLX(0.E0, 0.E0)
      B(1+0*M) = CMPLX(9.E19, 0.E0)
      B(2+1*M) = CMPLX(9.E18, 0.E0)
      B(3+2*M) = CMPLX(9.E17, 0.E0)
      B(4+3*M) = CMPLX(9.E16, 0.E0)
      B(5+4*M) = CMPLX(9.E15, 0.E0)
      B(6+5*M) = CMPLX(9.E14, 0.E0)
      CALL CGGSVD3('N','N','N', M, P, N, K, L, A, M, B, M, ALPHA, BETA, U, 1, V, 1, Q, 1, WORK, M*N, RWORK, IWORK, INFO)

      // Print information there is a NAN in BETA
      DO 40 I = 1, L
         IF( BETA(I).NE.BETA(I) ) THEN
            INFO = -I
            EXIT
         END IF
   40 CONTINUE
      IF( INFO.LT.0 ) THEN
         IF( NFAIL.EQ.0 .AND. FIRSTT ) THEN
            FIRSTT = .FALSE.
            CALL ALAHDG( NOUT, PATH )
         END IF
         WRITE( NOUT, FMT = 9997 ) -INFO
         NFAIL = NFAIL + 1
      END IF
      NRUN = NRUN + 1
      INFO = 0

      // Do for each value of M in MVAL.

      DO 30 IM = 1, NM
         M = MVAL( IM )
         P = PVAL( IM )
         N = NVAL( IM )

         DO 20 IMAT = 1, NTYPES

            // Do the tests only if DOTYPE( IMAT ) is true.

            IF( .NOT.DOTYPE( IMAT ) ) GO TO 20

            // Set up parameters with SLATB9 and generate test
            // matrices A and B with CLATMS.

            CALL SLATB9( PATH, IMAT, M, P, N, TYPE, KLA, KUA, KLB, KUB, ANORM, BNORM, MODEA, MODEB, CNDNMA, CNDNMB, DISTA, DISTB )

            // Generate M by N matrix A

            CALL CLATMS( M, N, DISTA, ISEED, TYPE, RWORK, MODEA, CNDNMA, ANORM, KLA, KUA, 'No packing', A, LDA, WORK, IINFO )
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUT, FMT = 9999 )IINFO
               INFO = ABS( IINFO )
               GO TO 20
            END IF

            // Generate P by N matrix B

            CALL CLATMS( P, N, DISTB, ISEED, TYPE, RWORK, MODEB, CNDNMB, BNORM, KLB, KUB, 'No packing', B, LDB, WORK, IINFO )
            IF( IINFO.NE.0 ) THEN
               WRITE( NOUT, FMT = 9999 )IINFO
               INFO = ABS( IINFO )
               GO TO 20
            END IF

            NT = 6

            CALL CGSVTS3( M, P, N, A, AF, LDA, B, BF, LDB, U, LDU, V, LDV, Q, LDQ, ALPHA, BETA, R, LDR, IWORK, WORK, LWORK, RWORK, RESULT )

            // Print information about the tests that did not
            // pass the threshold.

            DO 10 I = 1, NT
               IF( RESULT( I ).GE.THRESH ) THEN
                  IF( NFAIL.EQ.0 .AND. FIRSTT ) THEN
                     FIRSTT = .FALSE.
                     CALL ALAHDG( NOUT, PATH )
                  END IF
                  WRITE( NOUT, FMT = 9998 )M, P, N, IMAT, I, RESULT( I )
                  NFAIL = NFAIL + 1
               END IF
   10       CONTINUE
            NRUN = NRUN + NT

   20    CONTINUE
   30 CONTINUE

      // Print a summary of the results.

      CALL ALASUM( PATH, NOUT, NFAIL, NRUN, 0 )

 9999 FORMAT( ' CLATMS in CCKGSV   INFO = ', I5 )
 9998 FORMAT( ' M=', I4, ' P=', I4, ', N=', I4, ', type ', I2,
     $      ', test ', I2, ', ratio=', G13.6 )
 9997 FORMAT( ' FOUND NaN in BETA(', I4,')' )
      RETURN

      // End of CCKGSV

      }
