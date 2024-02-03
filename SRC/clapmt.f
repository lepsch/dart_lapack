      SUBROUTINE CLAPMT( FORWRD, M, N, X, LDX, K )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               FORWRD;
      int                LDX, M, N;
      // ..
      // .. Array Arguments ..
      int                K( * );
      COMPLEX            X( LDX, * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      int                I, II, J, IN;
      COMPLEX            TEMP
      // ..
      // .. Executable Statements ..

      if (N <= 1) RETURN;

      for (I = 1; I <= N; I++) { // 10
         K( I ) = -K( I )
      } // 10

      if ( FORWRD ) {

         // Forward permutation

         for (I = 1; I <= N; I++) { // 60

            IF( K( I ) > 0 ) GO TO 40

            J = I
            K( J ) = -K( J )
            IN = K( J )

            } // 20
            IF( K( IN ) > 0 ) GO TO 40

            for (II = 1; II <= M; II++) { // 30
               TEMP = X( II, J )
               X( II, J ) = X( II, IN )
               X( II, IN ) = TEMP
            } // 30

            K( IN ) = -K( IN )
            J = IN
            IN = K( IN )
            GO TO 20

            } // 40

         } // 60

      } else {

         // Backward permutation

         for (I = 1; I <= N; I++) { // 110

            IF( K( I ) > 0 ) GO TO 100

            K( I ) = -K( I )
            J = K( I )
            } // 80
            if (J == I) GO TO 100;

            for (II = 1; II <= M; II++) { // 90
               TEMP = X( II, I )
               X( II, I ) = X( II, J )
               X( II, J ) = TEMP
            } // 90

            K( J ) = -K( J )
            J = K( J )
            GO TO 80

            } // 100

         } // 110

      }

      RETURN

      // End of CLAPMT

      }
