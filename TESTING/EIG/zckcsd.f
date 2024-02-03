      SUBROUTINE ZCKCSD( NM, MVAL, PVAL, QVAL, NMATS, ISEED, THRESH, MMAX, X, XF, U1, U2, V1T, V2T, THETA, IWORK, WORK, RWORK, NIN, NOUT, INFO )

*  -- LAPACK test routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, NIN, NM, NMATS, MMAX, NOUT;
      double             THRESH;
      // ..
      // .. Array Arguments ..
      int                ISEED( 4 ), IWORK( * ), MVAL( * ), PVAL( * ), QVAL( * );
      double             RWORK( * ), THETA( * );
      COMPLEX*16         U1( * ), U2( * ), V1T( * ), V2T( * ), WORK( * ), X( * ), XF( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      int                NTESTS;
      const              NTESTS = 15 ;
      int                NTYPES;
      const              NTYPES = 4 ;
      double             GAPDIGIT, ORTH, REALONE, REALZERO, TEN;
      const              GAPDIGIT = 18.0D0, ORTH = 1.0D-12, REALONE = 1.0D0, REALZERO = 0.0D0, TEN = 10.0D0 ;
      COMPLEX*16         ONE, ZERO
      const              ONE = (1.0D0,0.0D0), ZERO = (0.0D0,0.0D0) ;
      double             PIOVER2;
      const     PIOVER2 = 1.57079632679489661923132169163975144210D0 ;
      // ..
      // .. Local Scalars ..
      bool               FIRSTT;
      String             PATH;
      int                I, IINFO, IM, IMAT, J, LDU1, LDU2, LDV1T, LDV2T, LDX, LWORK, M, NFAIL, NRUN, NT, P, Q, R;
      // ..
      // .. Local Arrays ..
      bool               DOTYPE( NTYPES );
      double             RESULT( NTESTS );
      // ..
      // .. External Subroutines ..
      // EXTERNAL ALAHDG, ALAREQ, ALASUM, ZCSDTS, ZLACSG, ZLAROR, ZLASET, ZDROT
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MIN
      // ..
      // .. External Functions ..
      double             DLARAN, DLARND;
      // EXTERNAL DLARAN, DLARND
      // ..
      // .. Executable Statements ..

      // Initialize constants and the random number seed.

      PATH( 1: 3 ) = 'CSD'
      INFO = 0
      NRUN = 0
      NFAIL = 0
      FIRSTT = .TRUE.
      alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );
      LDX = MMAX
      LDU1 = MMAX
      LDU2 = MMAX
      LDV1T = MMAX
      LDV2T = MMAX
      LWORK = MMAX*MMAX

      // Do for each value of M in MVAL.

      for (IM = 1; IM <= NM; IM++) { // 30
         M = MVAL( IM )
         P = PVAL( IM )
         Q = QVAL( IM )

         for (IMAT = 1; IMAT <= NTYPES; IMAT++) { // 20

            // Do the tests only if DOTYPE( IMAT ) is true.

            IF( .NOT.DOTYPE( IMAT ) ) GO TO 20

            // Generate X

            if ( IMAT.EQ.1 ) {
               zlaror('L', 'I', M, M, X, LDX, ISEED, WORK, IINFO );
               if ( M .NE. 0 .AND. IINFO .NE. 0 ) {
                  WRITE( NOUT, FMT = 9999 ) M, IINFO
                  INFO = ABS( IINFO )
                  GO TO 20
               }
            } else if ( IMAT.EQ.2 ) {
               R = MIN( P, M-P, Q, M-Q )
               for (I = 1; I <= R; I++) {
                  THETA(I) = PIOVER2 * DLARND( 1, ISEED )
               END DO
               zlacsg(M, P, Q, THETA, ISEED, X, LDX, WORK );
               for (I = 1; I <= M; I++) {
                  for (J = 1; J <= M; J++) {
                     X(I+(J-1)*LDX) = X(I+(J-1)*LDX) + ORTH*DLARND(2,ISEED)
                  END DO
               END DO
            } else if ( IMAT.EQ.3 ) {
               R = MIN( P, M-P, Q, M-Q )
               DO I = 1, R+1
                  THETA(I) = TEN**(-DLARND(1,ISEED)*GAPDIGIT)
               END DO
               DO I = 2, R+1
                  THETA(I) = THETA(I-1) + THETA(I)
               END DO
               for (I = 1; I <= R; I++) {
                  THETA(I) = PIOVER2 * THETA(I) / THETA(R+1)
               END DO
               zlacsg(M, P, Q, THETA, ISEED, X, LDX, WORK );
            } else {
               zlaset('F', M, M, ZERO, ONE, X, LDX );
               for (I = 1; I <= M; I++) {
                  J = INT( DLARAN( ISEED ) * M ) + 1
                  if ( J .NE. I ) {
                     zdrot(M, X(1+(I-1)*LDX), 1, X(1+(J-1)*LDX), 1, REALZERO, REALONE );
                  }
               END DO
            }

            NT = 15

            zcsdts(M, P, Q, X, XF, LDX, U1, LDU1, U2, LDU2, V1T, LDV1T, V2T, LDV2T, THETA, IWORK, WORK, LWORK, RWORK, RESULT );

            // Print information about the tests that did not
            // pass the threshold.

            for (I = 1; I <= NT; I++) { // 10
               if ( RESULT( I ).GE.THRESH ) {
                  if ( NFAIL.EQ.0 .AND. FIRSTT ) {
                     FIRSTT = .FALSE.
                     alahdg(NOUT, PATH );
                  }
                  WRITE( NOUT, FMT = 9998 )M, P, Q, IMAT, I, RESULT( I )
                  NFAIL = NFAIL + 1
               }
            } // 10
            NRUN = NRUN + NT
         } // 20
      } // 30

      // Print a summary of the results.

      alasum(PATH, NOUT, NFAIL, NRUN, 0 );

 9999 FORMAT( ' ZLAROR in ZCKCSD: M = ', I5, ', INFO = ', I15 )
 9998 FORMAT( ' M=', I4, ' P=', I4, ', Q=', I4, ', type ', I2, ', test ', I2, ', ratio=', G13.6 )
      RETURN

      // End of ZCKCSD

      }



      SUBROUTINE ZLACSG( M, P, Q, THETA, ISEED, X, LDX, WORK )
      IMPLICIT NONE

      int                LDX, M, P, Q;
      int                ISEED( 4 );
      double             THETA( * );
      COMPLEX*16         WORK( * ), X( LDX, * )

      COMPLEX*16         ONE, ZERO
      const              ONE = (1.0D0,0.0D0), ZERO = (0.0D0,0.0D0) ;

      int                I, INFO, R;

      R = MIN( P, M-P, Q, M-Q )

      zlaset('Full', M, M, ZERO, ZERO, X, LDX );

      DO I = 1, MIN(P,Q)-R
         X(I,I) = ONE
      END DO
      for (I = 1; I <= R; I++) {
         X(MIN(P,Q)-R+I,MIN(P,Q)-R+I) = DCMPLX( COS(THETA(I)), 0.0D0 )
      END DO
      DO I = 1, MIN(P,M-Q)-R
         X(P-I+1,M-I+1) = -ONE
      END DO
      for (I = 1; I <= R; I++) {
         X(P-(MIN(P,M-Q)-R)+1-I,M-(MIN(P,M-Q)-R)+1-I) = DCMPLX( -SIN(THETA(R-I+1)), 0.0D0 )
      END DO
      DO I = 1, MIN(M-P,Q)-R
         X(M-I+1,Q-I+1) = ONE
      END DO
      for (I = 1; I <= R; I++) {
         X(M-(MIN(M-P,Q)-R)+1-I,Q-(MIN(M-P,Q)-R)+1-I) = DCMPLX( SIN(THETA(R-I+1)), 0.0D0 )
      END DO
      DO I = 1, MIN(M-P,M-Q)-R
         X(P+I,Q+I) = ONE
      END DO
      for (I = 1; I <= R; I++) {
         X(P+(MIN(M-P,M-Q)-R)+I,Q+(MIN(M-P,M-Q)-R)+I) = DCMPLX( COS(THETA(I)), 0.0D0 )
      END DO
      zlaror('Left', 'No init', P, M, X, LDX, ISEED, WORK, INFO );
      zlaror('Left', 'No init', M-P, M, X(P+1,1), LDX, ISEED, WORK, INFO )       CALL ZLAROR( 'Right', 'No init', M, Q, X, LDX, ISEED, WORK, INFO )       CALL ZLAROR( 'Right', 'No init', M, M-Q, X(1,Q+1), LDX, ISEED, WORK, INFO );

      }
