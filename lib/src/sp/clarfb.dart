      void clarfb(SIDE, TRANS, DIRECT, STOREV, M, N, K, V, LDV, T, LDT, C, LDC, WORK, LDWORK ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             DIRECT, SIDE, STOREV, TRANS;
      int                K, LDC, LDT, LDV, LDWORK, M, N;
      // ..
      // .. Array Arguments ..
      Complex            C( LDC, * ), T( LDT, * ), V( LDV, * ), WORK( LDWORK, * );
      // ..

// =====================================================================

      // .. Parameters ..
      Complex            ONE;
      const              ONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      String             TRANST;
      int                I, J;
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. External Subroutines ..
      // EXTERNAL CCOPY, CGEMM, CLACGV, CTRMM
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if (M <= 0 || N <= 0) return;

      if ( lsame( TRANS, 'N' ) ) {
         TRANST = 'C';
      } else {
         TRANST = 'N';
      }

      if ( lsame( STOREV, 'C' ) ) {

         if ( lsame( DIRECT, 'F' ) ) {

            // Let  V =  ( V1 )    (first K rows)
                      // ( V2 )
            // where  V1  is unit lower triangular.

            if ( lsame( SIDE, 'L' ) ) {

               // Form  H * C  or  H**H * C  where  C = ( C1 )
                                                     // ( C2 )

               // W := C**H * V  =  (C1**H * V1 + C2**H * V2)  (stored in WORK)

               // W := C1**H

               for (J = 1; J <= K; J++) { // 10
                  ccopy(N, C( J, 1 ), LDC, WORK( 1, J ), 1 );
                  clacgv(N, WORK( 1, J ), 1 );
               } // 10

               // W := W * V1

               ctrmm('Right', 'Lower', 'No transpose', 'Unit', N, K, ONE, V, LDV, WORK, LDWORK );
               if ( M > K ) {

                  // W := W + C2**H *V2

                  cgemm('Conjugate transpose', 'No transpose', N, K, M-K, ONE, C( K+1, 1 ), LDC, V( K+1, 1 ), LDV, ONE, WORK, LDWORK );
               }

               // W := W * T**H  or  W * T

               ctrmm('Right', 'Upper', TRANST, 'Non-unit', N, K, ONE, T, LDT, WORK, LDWORK );

               // C := C - V * W**H

               if ( M > K ) {

                  // C2 := C2 - V2 * W**H

                  cgemm('No transpose', 'Conjugate transpose', M-K, N, K, -ONE, V( K+1, 1 ), LDV, WORK, LDWORK, ONE, C( K+1, 1 ), LDC );
               }

               // W := W * V1**H

               ctrmm('Right', 'Lower', 'Conjugate transpose', 'Unit', N, K, ONE, V, LDV, WORK, LDWORK );

               // C1 := C1 - W**H

               for (J = 1; J <= K; J++) { // 30
                  for (I = 1; I <= N; I++) { // 20
                     C[J, I] = C( J, I ) - CONJG( WORK( I, J ) );
                  } // 20
               } // 30

            } else if ( lsame( SIDE, 'R' ) ) {

               // Form  C * H  or  C * H**H  where  C = ( C1  C2 )

               // W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)

               // W := C1

               for (J = 1; J <= K; J++) { // 40
                  ccopy(M, C( 1, J ), 1, WORK( 1, J ), 1 );
               } // 40

               // W := W * V1

               ctrmm('Right', 'Lower', 'No transpose', 'Unit', M, K, ONE, V, LDV, WORK, LDWORK );
               if ( N > K ) {

                  // W := W + C2 * V2

                  cgemm('No transpose', 'No transpose', M, K, N-K, ONE, C( 1, K+1 ), LDC, V( K+1, 1 ), LDV, ONE, WORK, LDWORK );
               }

               // W := W * T  or  W * T**H

               ctrmm('Right', 'Upper', TRANS, 'Non-unit', M, K, ONE, T, LDT, WORK, LDWORK );

               // C := C - W * V**H

               if ( N > K ) {

                  // C2 := C2 - W * V2**H

                  cgemm('No transpose', 'Conjugate transpose', M, N-K, K, -ONE, WORK, LDWORK, V( K+1, 1 ), LDV, ONE, C( 1, K+1 ), LDC );
               }

               // W := W * V1**H

               ctrmm('Right', 'Lower', 'Conjugate transpose', 'Unit', M, K, ONE, V, LDV, WORK, LDWORK );

               // C1 := C1 - W

               for (J = 1; J <= K; J++) { // 60
                  for (I = 1; I <= M; I++) { // 50
                     C[I, J] = C( I, J ) - WORK( I, J );
                  } // 50
               } // 60
            }

         } else {

            // Let  V =  ( V1 )
                      // ( V2 )    (last K rows)
            // where  V2  is unit upper triangular.

            if ( lsame( SIDE, 'L' ) ) {

               // Form  H * C  or  H**H * C  where  C = ( C1 )
                                                   // ( C2 )

               // W := C**H * V  =  (C1**H * V1 + C2**H * V2)  (stored in WORK)

               // W := C2**H

               for (J = 1; J <= K; J++) { // 70
                  ccopy(N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 );
                  clacgv(N, WORK( 1, J ), 1 );
               } // 70

               // W := W * V2

               ctrmm('Right', 'Upper', 'No transpose', 'Unit', N, K, ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK );
               if ( M > K ) {

                  // W := W + C1**H * V1

                  cgemm('Conjugate transpose', 'No transpose', N, K, M-K, ONE, C, LDC, V, LDV, ONE, WORK, LDWORK );
               }

               // W := W * T**H  or  W * T

               ctrmm('Right', 'Lower', TRANST, 'Non-unit', N, K, ONE, T, LDT, WORK, LDWORK );

               // C := C - V * W**H

               if ( M > K ) {

                  // C1 := C1 - V1 * W**H

                  cgemm('No transpose', 'Conjugate transpose', M-K, N, K, -ONE, V, LDV, WORK, LDWORK, ONE, C, LDC );
               }

               // W := W * V2**H

               ctrmm('Right', 'Upper', 'Conjugate transpose', 'Unit', N, K, ONE, V( M-K+1, 1 ), LDV, WORK, LDWORK );

               // C2 := C2 - W**H

               for (J = 1; J <= K; J++) { // 90
                  for (I = 1; I <= N; I++) { // 80
                     C[M-K+J, I] = C( M-K+J, I ) - CONJG( WORK( I, J ) );
                  } // 80
               } // 90

            } else if ( lsame( SIDE, 'R' ) ) {

               // Form  C * H  or  C * H**H  where  C = ( C1  C2 )

               // W := C * V  =  (C1*V1 + C2*V2)  (stored in WORK)

               // W := C2

               for (J = 1; J <= K; J++) { // 100
                  ccopy(M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 );
               } // 100

               // W := W * V2

               ctrmm('Right', 'Upper', 'No transpose', 'Unit', M, K, ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK );
               if ( N > K ) {

                  // W := W + C1 * V1

                  cgemm('No transpose', 'No transpose', M, K, N-K, ONE, C, LDC, V, LDV, ONE, WORK, LDWORK );
               }

               // W := W * T  or  W * T**H

               ctrmm('Right', 'Lower', TRANS, 'Non-unit', M, K, ONE, T, LDT, WORK, LDWORK );

               // C := C - W * V**H

               if ( N > K ) {

                  // C1 := C1 - W * V1**H

                  cgemm('No transpose', 'Conjugate transpose', M, N-K, K, -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC );
               }

               // W := W * V2**H

               ctrmm('Right', 'Upper', 'Conjugate transpose', 'Unit', M, K, ONE, V( N-K+1, 1 ), LDV, WORK, LDWORK );

               // C2 := C2 - W

               for (J = 1; J <= K; J++) { // 120
                  for (I = 1; I <= M; I++) { // 110
                     C[I, N-K+J] = C( I, N-K+J ) - WORK( I, J );
                  } // 110
               } // 120
            }
         }

      } else if ( lsame( STOREV, 'R' ) ) {

         if ( lsame( DIRECT, 'F' ) ) {

            // Let  V =  ( V1  V2 )    (V1: first K columns)
            // where  V1  is unit upper triangular.

            if ( lsame( SIDE, 'L' ) ) {

               // Form  H * C  or  H**H * C  where  C = ( C1 )
                                                     // ( C2 )

               // W := C**H * V**H  =  (C1**H * V1**H + C2**H * V2**H) (stored in WORK)

               // W := C1**H

               for (J = 1; J <= K; J++) { // 130
                  ccopy(N, C( J, 1 ), LDC, WORK( 1, J ), 1 );
                  clacgv(N, WORK( 1, J ), 1 );
               } // 130

               // W := W * V1**H

               ctrmm('Right', 'Upper', 'Conjugate transpose', 'Unit', N, K, ONE, V, LDV, WORK, LDWORK );
               if ( M > K ) {

                  // W := W + C2**H * V2**H

                  cgemm('Conjugate transpose', 'Conjugate transpose', N, K, M-K, ONE, C( K+1, 1 ), LDC, V( 1, K+1 ), LDV, ONE, WORK, LDWORK );
               }

               // W := W * T**H  or  W * T

               ctrmm('Right', 'Upper', TRANST, 'Non-unit', N, K, ONE, T, LDT, WORK, LDWORK );

               // C := C - V**H * W**H

               if ( M > K ) {

                  // C2 := C2 - V2**H * W**H

                  cgemm('Conjugate transpose', 'Conjugate transpose', M-K, N, K, -ONE, V( 1, K+1 ), LDV, WORK, LDWORK, ONE, C( K+1, 1 ), LDC );
               }

               // W := W * V1

               ctrmm('Right', 'Upper', 'No transpose', 'Unit', N, K, ONE, V, LDV, WORK, LDWORK );

               // C1 := C1 - W**H

               for (J = 1; J <= K; J++) { // 150
                  for (I = 1; I <= N; I++) { // 140
                     C[J, I] = C( J, I ) - CONJG( WORK( I, J ) );
                  } // 140
               } // 150

            } else if ( lsame( SIDE, 'R' ) ) {

               // Form  C * H  or  C * H**H  where  C = ( C1  C2 )

               // W := C * V**H  =  (C1*V1**H + C2*V2**H)  (stored in WORK)

               // W := C1

               for (J = 1; J <= K; J++) { // 160
                  ccopy(M, C( 1, J ), 1, WORK( 1, J ), 1 );
               } // 160

               // W := W * V1**H

               ctrmm('Right', 'Upper', 'Conjugate transpose', 'Unit', M, K, ONE, V, LDV, WORK, LDWORK );
               if ( N > K ) {

                  // W := W + C2 * V2**H

                  cgemm('No transpose', 'Conjugate transpose', M, K, N-K, ONE, C( 1, K+1 ), LDC, V( 1, K+1 ), LDV, ONE, WORK, LDWORK );
               }

               // W := W * T  or  W * T**H

               ctrmm('Right', 'Upper', TRANS, 'Non-unit', M, K, ONE, T, LDT, WORK, LDWORK );

               // C := C - W * V

               if ( N > K ) {

                  // C2 := C2 - W * V2

                  cgemm('No transpose', 'No transpose', M, N-K, K, -ONE, WORK, LDWORK, V( 1, K+1 ), LDV, ONE, C( 1, K+1 ), LDC );
               }

               // W := W * V1

               ctrmm('Right', 'Upper', 'No transpose', 'Unit', M, K, ONE, V, LDV, WORK, LDWORK );

               // C1 := C1 - W

               for (J = 1; J <= K; J++) { // 180
                  for (I = 1; I <= M; I++) { // 170
                     C[I, J] = C( I, J ) - WORK( I, J );
                  } // 170
               } // 180

            }

         } else {

            // Let  V =  ( V1  V2 )    (V2: last K columns)
            // where  V2  is unit lower triangular.

            if ( lsame( SIDE, 'L' ) ) {

               // Form  H * C  or  H**H * C  where  C = ( C1 )
                                                     // ( C2 )

               // W := C**H * V**H  =  (C1**H * V1**H + C2**H * V2**H) (stored in WORK)

               // W := C2**H

               for (J = 1; J <= K; J++) { // 190
                  ccopy(N, C( M-K+J, 1 ), LDC, WORK( 1, J ), 1 );
                  clacgv(N, WORK( 1, J ), 1 );
               } // 190

               // W := W * V2**H

               ctrmm('Right', 'Lower', 'Conjugate transpose', 'Unit', N, K, ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK );
               if ( M > K ) {

                  // W := W + C1**H * V1**H

                  cgemm('Conjugate transpose', 'Conjugate transpose', N, K, M-K, ONE, C, LDC, V, LDV, ONE, WORK, LDWORK );
               }

               // W := W * T**H  or  W * T

               ctrmm('Right', 'Lower', TRANST, 'Non-unit', N, K, ONE, T, LDT, WORK, LDWORK );

               // C := C - V**H * W**H

               if ( M > K ) {

                  // C1 := C1 - V1**H * W**H

                  cgemm('Conjugate transpose', 'Conjugate transpose', M-K, N, K, -ONE, V, LDV, WORK, LDWORK, ONE, C, LDC );
               }

               // W := W * V2

               ctrmm('Right', 'Lower', 'No transpose', 'Unit', N, K, ONE, V( 1, M-K+1 ), LDV, WORK, LDWORK );

               // C2 := C2 - W**H

               for (J = 1; J <= K; J++) { // 210
                  for (I = 1; I <= N; I++) { // 200
                     C[M-K+J, I] = C( M-K+J, I ) - CONJG( WORK( I, J ) );
                  } // 200
               } // 210

            } else if ( lsame( SIDE, 'R' ) ) {

               // Form  C * H  or  C * H**H  where  C = ( C1  C2 )

               // W := C * V**H  =  (C1*V1**H + C2*V2**H)  (stored in WORK)

               // W := C2

               for (J = 1; J <= K; J++) { // 220
                  ccopy(M, C( 1, N-K+J ), 1, WORK( 1, J ), 1 );
               } // 220

               // W := W * V2**H

               ctrmm('Right', 'Lower', 'Conjugate transpose', 'Unit', M, K, ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK );
               if ( N > K ) {

                  // W := W + C1 * V1**H

                  cgemm('No transpose', 'Conjugate transpose', M, K, N-K, ONE, C, LDC, V, LDV, ONE, WORK, LDWORK );
               }

               // W := W * T  or  W * T**H

               ctrmm('Right', 'Lower', TRANS, 'Non-unit', M, K, ONE, T, LDT, WORK, LDWORK );

               // C := C - W * V

               if ( N > K ) {

                  // C1 := C1 - W * V1

                  cgemm('No transpose', 'No transpose', M, N-K, K, -ONE, WORK, LDWORK, V, LDV, ONE, C, LDC );
               }

               // W := W * V2

               ctrmm('Right', 'Lower', 'No transpose', 'Unit', M, K, ONE, V( 1, N-K+1 ), LDV, WORK, LDWORK );

               // C1 := C1 - W

               for (J = 1; J <= K; J++) { // 240
                  for (I = 1; I <= M; I++) { // 230
                     C[I, N-K+J] = C( I, N-K+J ) - WORK( I, J );
                  } // 230
               } // 240

            }

         }
      }

      return;
      }