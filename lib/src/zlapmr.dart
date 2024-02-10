      void zlapmr(FORWRD, M, N, final Matrix<double> X, final int LDX, K ) {

// -- LAPACK auxiliary routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      bool               FORWRD;
      int                LDX, M, N;
      int                K( * );
      Complex         X( LDX, * );
      // ..

// =====================================================================

      // .. Local Scalars ..
      int                I, IN, J, JJ;
      Complex         TEMP;

      if (M <= 1) return;

      for (I = 1; I <= M; I++) { // 10
         K[I] = -K( I );
      } // 10

      if ( FORWRD ) {

         // Forward permutation

         for (I = 1; I <= M; I++) { // 50

            if( K( I ) > 0 ) GO TO 40;

            J = I;
            K[J] = -K( J );
            IN = K( J );

            } // 20
            if( K( IN ) > 0 ) GO TO 40;

            for (JJ = 1; JJ <= N; JJ++) { // 30
               TEMP = X( J, JJ );
               X[J][JJ] = X( IN, JJ );
               X[IN][JJ] = TEMP;
            } // 30

            K[IN] = -K( IN );
            J = IN;
            IN = K( IN );
            GO TO 20;

            } // 40

         } // 50

      } else {

         // Backward permutation

         for (I = 1; I <= M; I++) { // 90

            if( K( I ) > 0 ) GO TO 80;

            K[I] = -K( I );
            J = K( I );
            } // 60
            if (J == I) GO TO 80;

            for (JJ = 1; JJ <= N; JJ++) { // 70
               TEMP = X( I, JJ );
               X[I][JJ] = X( J, JJ );
               X[J][JJ] = TEMP;
            } // 70

            K[J] = -K( J );
            J = K( J );
            GO TO 60;

            } // 80

         } // 90

      }

      }
