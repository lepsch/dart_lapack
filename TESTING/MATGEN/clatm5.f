      SUBROUTINE CLATM5( PRTYPE, M, N, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE, F, LDF, R, LDR, L, LDL, ALPHA, QBLCKA, QBLCKB )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LDB, LDC, LDD, LDE, LDF, LDL, LDR, M, N, PRTYPE, QBLCKA, QBLCKB;
      REAL               ALPHA
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), B( LDB, * ), C( LDC, * ), D( LDD, * ), E( LDE, * ), F( LDF, * ), L( LDL, * ), R( LDR, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX            ONE, TWO, ZERO, HALF, TWENTY
      const              ONE = ( 1.0E+0, 0.0E+0 ), TWO = ( 2.0E+0, 0.0E+0 ), ZERO = ( 0.0E+0, 0.0E+0 ), HALF = ( 0.5E+0, 0.0E+0 ), TWENTY = ( 2.0E+1, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, J, K;
      COMPLEX            IMEPS, REEPS
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CMPLX, MOD, SIN
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMM
      // ..
      // .. Executable Statements ..

      if ( PRTYPE.EQ.1 ) {
         for (I = 1; I <= M; I++) { // 20
            for (J = 1; J <= M; J++) { // 10
               if ( I.EQ.J ) {
                  A( I, J ) = ONE
                  D( I, J ) = ONE
               } else if ( I.EQ.J-1 ) {
                  A( I, J ) = -ONE
                  D( I, J ) = ZERO
               } else {
                  A( I, J ) = ZERO
                  D( I, J ) = ZERO
               }
   10       CONTINUE
   20    CONTINUE

         for (I = 1; I <= N; I++) { // 40
            for (J = 1; J <= N; J++) { // 30
               if ( I.EQ.J ) {
                  B( I, J ) = ONE - ALPHA
                  E( I, J ) = ONE
               } else if ( I.EQ.J-1 ) {
                  B( I, J ) = ONE
                  E( I, J ) = ZERO
               } else {
                  B( I, J ) = ZERO
                  E( I, J ) = ZERO
               }
   30       CONTINUE
   40    CONTINUE

         for (I = 1; I <= M; I++) { // 60
            for (J = 1; J <= N; J++) { // 50
               R( I, J ) = ( HALF-SIN( CMPLX( I / J ) ) )*TWENTY
               L( I, J ) = R( I, J )
   50       CONTINUE
   60    CONTINUE

      } else if ( PRTYPE.EQ.2 .OR. PRTYPE.EQ.3 ) {
         for (I = 1; I <= M; I++) { // 80
            for (J = 1; J <= M; J++) { // 70
               if ( I.LE.J ) {
                  A( I, J ) = ( HALF-SIN( CMPLX( I ) ) )*TWO
                  D( I, J ) = ( HALF-SIN( CMPLX( I*J ) ) )*TWO
               } else {
                  A( I, J ) = ZERO
                  D( I, J ) = ZERO
               }
   70       CONTINUE
   80    CONTINUE

         for (I = 1; I <= N; I++) { // 100
            for (J = 1; J <= N; J++) { // 90
               if ( I.LE.J ) {
                  B( I, J ) = ( HALF-SIN( CMPLX( I+J ) ) )*TWO
                  E( I, J ) = ( HALF-SIN( CMPLX( J ) ) )*TWO
               } else {
                  B( I, J ) = ZERO
                  E( I, J ) = ZERO
               }
   90       CONTINUE
  100    CONTINUE

         for (I = 1; I <= M; I++) { // 120
            for (J = 1; J <= N; J++) { // 110
               R( I, J ) = ( HALF-SIN( CMPLX( I*J ) ) )*TWENTY
               L( I, J ) = ( HALF-SIN( CMPLX( I+J ) ) )*TWENTY
  110       CONTINUE
  120    CONTINUE

         if ( PRTYPE.EQ.3 ) {
            IF( QBLCKA.LE.1 ) QBLCKA = 2
            DO 130 K = 1, M - 1, QBLCKA
               A( K+1, K+1 ) = A( K, K )
               A( K+1, K ) = -SIN( A( K, K+1 ) )
  130       CONTINUE

            IF( QBLCKB.LE.1 ) QBLCKB = 2
            DO 140 K = 1, N - 1, QBLCKB
               B( K+1, K+1 ) = B( K, K )
               B( K+1, K ) = -SIN( B( K, K+1 ) )
  140       CONTINUE
         }

      } else if ( PRTYPE.EQ.4 ) {
         for (I = 1; I <= M; I++) { // 160
            for (J = 1; J <= M; J++) { // 150
               A( I, J ) = ( HALF-SIN( CMPLX( I*J ) ) )*TWENTY
               D( I, J ) = ( HALF-SIN( CMPLX( I+J ) ) )*TWO
  150       CONTINUE
  160    CONTINUE

         for (I = 1; I <= N; I++) { // 180
            for (J = 1; J <= N; J++) { // 170
               B( I, J ) = ( HALF-SIN( CMPLX( I+J ) ) )*TWENTY
               E( I, J ) = ( HALF-SIN( CMPLX( I*J ) ) )*TWO
  170       CONTINUE
  180    CONTINUE

         for (I = 1; I <= M; I++) { // 200
            for (J = 1; J <= N; J++) { // 190
               R( I, J ) = ( HALF-SIN( CMPLX( J / I ) ) )*TWENTY
               L( I, J ) = ( HALF-SIN( CMPLX( I*J ) ) )*TWO
  190       CONTINUE
  200    CONTINUE

      } else if ( PRTYPE.GE.5 ) {
         REEPS = HALF*TWO*TWENTY / ALPHA
         IMEPS = ( HALF-TWO ) / ALPHA
         for (I = 1; I <= M; I++) { // 220
            for (J = 1; J <= N; J++) { // 210
               R( I, J ) = ( HALF-SIN( CMPLX( I*J ) ) )*ALPHA / TWENTY
               L( I, J ) = ( HALF-SIN( CMPLX( I+J ) ) )*ALPHA / TWENTY
  210       CONTINUE
  220    CONTINUE

         for (I = 1; I <= M; I++) { // 230
            D( I, I ) = ONE
  230    CONTINUE

         for (I = 1; I <= M; I++) { // 240
            if ( I.LE.4 ) {
               A( I, I ) = ONE
               IF( I.GT.2 ) A( I, I ) = ONE + REEPS
               if ( MOD( I, 2 ).NE.0 .AND. I.LT.M ) {
                  A( I, I+1 ) = IMEPS
               } else if ( I.GT.1 ) {
                  A( I, I-1 ) = -IMEPS
               }
            } else if ( I.LE.8 ) {
               if ( I.LE.6 ) {
                  A( I, I ) = REEPS
               } else {
                  A( I, I ) = -REEPS
               }
               if ( MOD( I, 2 ).NE.0 .AND. I.LT.M ) {
                  A( I, I+1 ) = ONE
               } else if ( I.GT.1 ) {
                  A( I, I-1 ) = -ONE
               }
            } else {
               A( I, I ) = ONE
               if ( MOD( I, 2 ).NE.0 .AND. I.LT.M ) {
                  A( I, I+1 ) = IMEPS*2
               } else if ( I.GT.1 ) {
                  A( I, I-1 ) = -IMEPS*2
               }
            }
  240    CONTINUE

         for (I = 1; I <= N; I++) { // 250
            E( I, I ) = ONE
            if ( I.LE.4 ) {
               B( I, I ) = -ONE
               IF( I.GT.2 ) B( I, I ) = ONE - REEPS
               if ( MOD( I, 2 ).NE.0 .AND. I.LT.N ) {
                  B( I, I+1 ) = IMEPS
               } else if ( I.GT.1 ) {
                  B( I, I-1 ) = -IMEPS
               }
            } else if ( I.LE.8 ) {
               if ( I.LE.6 ) {
                  B( I, I ) = REEPS
               } else {
                  B( I, I ) = -REEPS
               }
               if ( MOD( I, 2 ).NE.0 .AND. I.LT.N ) {
                  B( I, I+1 ) = ONE + IMEPS
               } else if ( I.GT.1 ) {
                  B( I, I-1 ) = -ONE - IMEPS
               }
            } else {
               B( I, I ) = ONE - REEPS
               if ( MOD( I, 2 ).NE.0 .AND. I.LT.N ) {
                  B( I, I+1 ) = IMEPS*2
               } else if ( I.GT.1 ) {
                  B( I, I-1 ) = -IMEPS*2
               }
            }
  250    CONTINUE
      }

      // Compute rhs (C, F)

      cgemm('N', 'N', M, N, M, ONE, A, LDA, R, LDR, ZERO, C, LDC );
      cgemm('N', 'N', M, N, N, -ONE, L, LDL, B, LDB, ONE, C, LDC );
      cgemm('N', 'N', M, N, M, ONE, D, LDD, R, LDR, ZERO, F, LDF );
      cgemm('N', 'N', M, N, N, -ONE, L, LDL, E, LDE, ONE, F, LDF );

      // End of CLATM5

      }
