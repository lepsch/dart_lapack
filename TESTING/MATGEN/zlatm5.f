      SUBROUTINE ZLATM5( PRTYPE, M, N, A, LDA, B, LDB, C, LDC, D, LDD, E, LDE, F, LDF, R, LDR, L, LDL, ALPHA, QBLCKA, QBLCKB )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                LDA, LDB, LDC, LDD, LDE, LDF, LDL, LDR, M, N, PRTYPE, QBLCKA, QBLCKB;
      double             ALPHA;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * ), C( LDC, * ), D( LDD, * ), E( LDE, * ), F( LDF, * ), L( LDL, * ), R( LDR, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         ONE, TWO, ZERO, HALF, TWENTY
      const              ONE = ( 1.0D+0, 0.0D+0 ), TWO = ( 2.0D+0, 0.0D+0 ), ZERO = ( 0.0D+0, 0.0D+0 ), HALF = ( 0.5D+0, 0.0D+0 ), TWENTY = ( 2.0D+1, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      int                I, J, K;
      COMPLEX*16         IMEPS, REEPS
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DCMPLX, MOD, SIN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZGEMM
      // ..
      // .. Executable Statements ..

      if ( PRTYPE == 1 ) {
         for (I = 1; I <= M; I++) { // 20
            for (J = 1; J <= M; J++) { // 10
               if ( I == J ) {
                  A( I, J ) = ONE
                  D( I, J ) = ONE
               } else if ( I == J-1 ) {
                  A( I, J ) = -ONE
                  D( I, J ) = ZERO
               } else {
                  A( I, J ) = ZERO
                  D( I, J ) = ZERO
               }
            } // 10
         } // 20

         for (I = 1; I <= N; I++) { // 40
            for (J = 1; J <= N; J++) { // 30
               if ( I == J ) {
                  B( I, J ) = ONE - ALPHA
                  E( I, J ) = ONE
               } else if ( I == J-1 ) {
                  B( I, J ) = ONE
                  E( I, J ) = ZERO
               } else {
                  B( I, J ) = ZERO
                  E( I, J ) = ZERO
               }
            } // 30
         } // 40

         for (I = 1; I <= M; I++) { // 60
            for (J = 1; J <= N; J++) { // 50
               R( I, J ) = ( HALF-SIN( DCMPLX( I / J ) ) )*TWENTY
               L( I, J ) = R( I, J )
            } // 50
         } // 60

      } else if ( PRTYPE == 2 || PRTYPE == 3 ) {
         for (I = 1; I <= M; I++) { // 80
            for (J = 1; J <= M; J++) { // 70
               if ( I.LE.J ) {
                  A( I, J ) = ( HALF-SIN( DCMPLX( I ) ) )*TWO
                  D( I, J ) = ( HALF-SIN( DCMPLX( I*J ) ) )*TWO
               } else {
                  A( I, J ) = ZERO
                  D( I, J ) = ZERO
               }
            } // 70
         } // 80

         for (I = 1; I <= N; I++) { // 100
            for (J = 1; J <= N; J++) { // 90
               if ( I.LE.J ) {
                  B( I, J ) = ( HALF-SIN( DCMPLX( I+J ) ) )*TWO
                  E( I, J ) = ( HALF-SIN( DCMPLX( J ) ) )*TWO
               } else {
                  B( I, J ) = ZERO
                  E( I, J ) = ZERO
               }
            } // 90
         } // 100

         for (I = 1; I <= M; I++) { // 120
            for (J = 1; J <= N; J++) { // 110
               R( I, J ) = ( HALF-SIN( DCMPLX( I*J ) ) )*TWENTY
               L( I, J ) = ( HALF-SIN( DCMPLX( I+J ) ) )*TWENTY
            } // 110
         } // 120

         if ( PRTYPE == 3 ) {
            if (QBLCKA.LE.1) QBLCKA = 2;
            DO 130 K = 1, M - 1, QBLCKA
               A( K+1, K+1 ) = A( K, K )
               A( K+1, K ) = -SIN( A( K, K+1 ) )
            } // 130

            if (QBLCKB.LE.1) QBLCKB = 2;
            DO 140 K = 1, N - 1, QBLCKB
               B( K+1, K+1 ) = B( K, K )
               B( K+1, K ) = -SIN( B( K, K+1 ) )
            } // 140
         }

      } else if ( PRTYPE == 4 ) {
         for (I = 1; I <= M; I++) { // 160
            for (J = 1; J <= M; J++) { // 150
               A( I, J ) = ( HALF-SIN( DCMPLX( I*J ) ) )*TWENTY
               D( I, J ) = ( HALF-SIN( DCMPLX( I+J ) ) )*TWO
            } // 150
         } // 160

         for (I = 1; I <= N; I++) { // 180
            for (J = 1; J <= N; J++) { // 170
               B( I, J ) = ( HALF-SIN( DCMPLX( I+J ) ) )*TWENTY
               E( I, J ) = ( HALF-SIN( DCMPLX( I*J ) ) )*TWO
            } // 170
         } // 180

         for (I = 1; I <= M; I++) { // 200
            for (J = 1; J <= N; J++) { // 190
               R( I, J ) = ( HALF-SIN( DCMPLX( J / I ) ) )*TWENTY
               L( I, J ) = ( HALF-SIN( DCMPLX( I*J ) ) )*TWO
            } // 190
         } // 200

      } else if ( PRTYPE >= 5 ) {
         REEPS = HALF*TWO*TWENTY / ALPHA
         IMEPS = ( HALF-TWO ) / ALPHA
         for (I = 1; I <= M; I++) { // 220
            for (J = 1; J <= N; J++) { // 210
               R( I, J ) = ( HALF-SIN( DCMPLX( I*J ) ) )*ALPHA / TWENTY
               L( I, J ) = ( HALF-SIN( DCMPLX( I+J ) ) )*ALPHA / TWENTY
            } // 210
         } // 220

         for (I = 1; I <= M; I++) { // 230
            D( I, I ) = ONE
         } // 230

         for (I = 1; I <= M; I++) { // 240
            if ( I.LE.4 ) {
               A( I, I ) = ONE
               if (I > 2) A( I, I ) = ONE + REEPS;
               if ( MOD( I, 2 ) != 0 && I < M ) {
                  A( I, I+1 ) = IMEPS
               } else if ( I > 1 ) {
                  A( I, I-1 ) = -IMEPS
               }
            } else if ( I.LE.8 ) {
               if ( I.LE.6 ) {
                  A( I, I ) = REEPS
               } else {
                  A( I, I ) = -REEPS
               }
               if ( MOD( I, 2 ) != 0 && I < M ) {
                  A( I, I+1 ) = ONE
               } else if ( I > 1 ) {
                  A( I, I-1 ) = -ONE
               }
            } else {
               A( I, I ) = ONE
               if ( MOD( I, 2 ) != 0 && I < M ) {
                  A( I, I+1 ) = IMEPS*2
               } else if ( I > 1 ) {
                  A( I, I-1 ) = -IMEPS*2
               }
            }
         } // 240

         for (I = 1; I <= N; I++) { // 250
            E( I, I ) = ONE
            if ( I.LE.4 ) {
               B( I, I ) = -ONE
               if (I > 2) B( I, I ) = ONE - REEPS;
               if ( MOD( I, 2 ) != 0 && I < N ) {
                  B( I, I+1 ) = IMEPS
               } else if ( I > 1 ) {
                  B( I, I-1 ) = -IMEPS
               }
            } else if ( I.LE.8 ) {
               if ( I.LE.6 ) {
                  B( I, I ) = REEPS
               } else {
                  B( I, I ) = -REEPS
               }
               if ( MOD( I, 2 ) != 0 && I < N ) {
                  B( I, I+1 ) = ONE + IMEPS
               } else if ( I > 1 ) {
                  B( I, I-1 ) = -ONE - IMEPS
               }
            } else {
               B( I, I ) = ONE - REEPS
               if ( MOD( I, 2 ) != 0 && I < N ) {
                  B( I, I+1 ) = IMEPS*2
               } else if ( I > 1 ) {
                  B( I, I-1 ) = -IMEPS*2
               }
            }
         } // 250
      }

      // Compute rhs (C, F)

      zgemm('N', 'N', M, N, M, ONE, A, LDA, R, LDR, ZERO, C, LDC );
      zgemm('N', 'N', M, N, N, -ONE, L, LDL, B, LDB, ONE, C, LDC );
      zgemm('N', 'N', M, N, M, ONE, D, LDD, R, LDR, ZERO, F, LDF );
      zgemm('N', 'N', M, N, N, -ONE, L, LDL, E, LDE, ONE, F, LDF );

      // End of ZLATM5

      }
