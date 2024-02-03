      SUBROUTINE DCKCSD( NM, MVAL, PVAL, QVAL, NMATS, ISEED, THRESH, MMAX, X, XF, U1, U2, V1T, V2T, THETA, IWORK, WORK, RWORK, NIN, NOUT, INFO );

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
      double             U1( * ), U2( * ), V1T( * ), V2T( * ), WORK( * ), X( * ), XF( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      int                NTESTS;
      const              NTESTS = 15 ;
      int                NTYPES;
      const              NTYPES = 4 ;
      double             GAPDIGIT, ONE, ORTH, TEN, ZERO;
      const              GAPDIGIT = 18.0, ONE = 1.0, ORTH = 1.0e-12, TEN = 10.0, ZERO = 0.0 ;
      double             PIOVER2;
      const     PIOVER2 = 1.57079632679489661923132169163975144210 ;
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
      // EXTERNAL ALAHDG, ALAREQ, ALASUM, DCSDTS, DLACSG, DLAROR, DLASET, DROT
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

      PATH( 1: 3 ) = 'CSD';
      INFO = 0;
      NRUN = 0;
      NFAIL = 0;
      FIRSTT = true;
      alareq(PATH, NMATS, DOTYPE, NTYPES, NIN, NOUT );
      LDX = MMAX;
      LDU1 = MMAX;
      LDU2 = MMAX;
      LDV1T = MMAX;
      LDV2T = MMAX;
      LWORK = MMAX*MMAX;

      // Do for each value of M in MVAL.

      for (IM = 1; IM <= NM; IM++) { // 30
         M = MVAL( IM );
         P = PVAL( IM );
         Q = QVAL( IM );

         for (IMAT = 1; IMAT <= NTYPES; IMAT++) { // 20

            // Do the tests only if DOTYPE( IMAT ) is true.

            IF( !DOTYPE( IMAT ) ) GO TO 20;

            // Generate X

            if ( IMAT == 1 ) {
               dlaror('L', 'I', M, M, X, LDX, ISEED, WORK, IINFO );
               if ( M != 0 && IINFO != 0 ) {
                  WRITE( NOUT, FMT = 9999 ) M, IINFO;
                  INFO = ABS( IINFO );
                  GO TO 20;
               }
            } else if ( IMAT == 2 ) {
               R = MIN( P, M-P, Q, M-Q );
               for (I = 1; I <= R; I++) {
                  THETA(I) = PIOVER2 * DLARND( 1, ISEED );
               }
               dlacsg(M, P, Q, THETA, ISEED, X, LDX, WORK );
               for (I = 1; I <= M; I++) {
                  for (J = 1; J <= M; J++) {
                     X(I+(J-1)*LDX) = X(I+(J-1)*LDX) + ORTH*DLARND(2,ISEED);
                  }
               }
            } else if ( IMAT == 3 ) {
               R = MIN( P, M-P, Q, M-Q );
               for (I = 1; I <= R+1; I++) {
                  THETA(I) = TEN**(-DLARND(1,ISEED)*GAPDIGIT);
               }
               for (I = 2; I <= R+1; I++) {
                  THETA(I) = THETA(I-1) + THETA(I);
               }
               for (I = 1; I <= R; I++) {
                  THETA(I) = PIOVER2 * THETA(I) / THETA(R+1);
               }
               dlacsg(M, P, Q, THETA, ISEED, X, LDX, WORK );
            } else {
               dlaset('F', M, M, ZERO, ONE, X, LDX );
               for (I = 1; I <= M; I++) {
                  J = INT( DLARAN( ISEED ) * M ) + 1;
                  if ( J != I ) {
                     drot(M, X(1+(I-1)*LDX), 1, X(1+(J-1)*LDX), 1, ZERO, ONE );
                  }
               }
            }

            NT = 15;

            dcsdts(M, P, Q, X, XF, LDX, U1, LDU1, U2, LDU2, V1T, LDV1T, V2T, LDV2T, THETA, IWORK, WORK, LWORK, RWORK, RESULT );

            // Print information about the tests that did not
            // pass the threshold.

            for (I = 1; I <= NT; I++) { // 10
               if ( RESULT( I ) >= THRESH ) {
                  if ( NFAIL == 0 && FIRSTT ) {
                     FIRSTT = false;
                     alahdg(NOUT, PATH );
                  }
                  WRITE( NOUT, FMT = 9998 )M, P, Q, IMAT, I, RESULT( I );
                  NFAIL = NFAIL + 1;
               }
            } // 10
            NRUN = NRUN + NT;
         } // 20
      } // 30

      // Print a summary of the results.

      alasum(PATH, NOUT, NFAIL, NRUN, 0 );

 9999 FORMAT( ' DLAROR in DCKCSD: M = ', I5, ', INFO = ', I15 );
 9998 FORMAT( ' M=', I4, ' P=', I4, ', Q=', I4, ', type ', I2, ', test ', I2, ', ratio=', G13.6 );
      RETURN;

      // End of DCKCSD

      }



      SUBROUTINE DLACSG( M, P, Q, THETA, ISEED, X, LDX, WORK );
      IMPLICIT NONE;

      int                LDX, M, P, Q;
      int                ISEED( 4 );
      double             THETA( * );
      double             WORK( * ), X( LDX, * );

      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;

      int                I, INFO, R;

      R = MIN( P, M-P, Q, M-Q );

      dlaset('Full', M, M, ZERO, ZERO, X, LDX );

      DO I = 1, MIN(P,Q)-R;
         X(I,I) = ONE;
      }
      for (I = 1; I <= R; I++) {
         X(MIN(P,Q)-R+I,MIN(P,Q)-R+I) = COS(THETA(I));
      }
      DO I = 1, MIN(P,M-Q)-R;
         X(P-I+1,M-I+1) = -ONE;
      }
      for (I = 1; I <= R; I++) {
         X(P-(MIN(P,M-Q)-R)+1-I,M-(MIN(P,M-Q)-R)+1-I) = -SIN(THETA(R-I+1));
      }
      DO I = 1, MIN(M-P,Q)-R;
         X(M-I+1,Q-I+1) = ONE;
      }
      for (I = 1; I <= R; I++) {
         X(M-(MIN(M-P,Q)-R)+1-I,Q-(MIN(M-P,Q)-R)+1-I) = SIN(THETA(R-I+1));
      }
      DO I = 1, MIN(M-P,M-Q)-R;
         X(P+I,Q+I) = ONE;
      }
      for (I = 1; I <= R; I++) {
         X(P+(MIN(M-P,M-Q)-R)+I,Q+(MIN(M-P,M-Q)-R)+I) = COS(THETA(I));
      }
      dlaror('Left', 'No init', P, M, X, LDX, ISEED, WORK, INFO );
      dlaror('Left', 'No init', M-P, M, X(P+1,1), LDX, ISEED, WORK, INFO )       CALL DLAROR( 'Right', 'No init', M, Q, X, LDX, ISEED, WORK, INFO )       CALL DLAROR( 'Right', 'No init', M, M-Q, X(1,Q+1), LDX, ISEED, WORK, INFO );

      }
