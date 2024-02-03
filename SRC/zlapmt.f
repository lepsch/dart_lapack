      SUBROUTINE ZLAPMT( FORWRD, M, N, X, LDX, K )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      bool               FORWRD;
      int                LDX, M, N;
      // ..
      // .. Array Arguments ..
      int                K( * );
      COMPLEX*16         X( LDX, * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      int                I, II, IN, J;
      COMPLEX*16         TEMP
      // ..
      // .. Executable Statements ..

      IF( N.LE.1 ) RETURN

      for (I = 1; I <= N; I++) { // 10
         K( I ) = -K( I )
   10 CONTINUE

      if ( FORWRD ) {

         // Forward permutation

         for (I = 1; I <= N; I++) { // 50

            IF( K( I ).GT.0 ) GO TO 40

            J = I
            K( J ) = -K( J )
            IN = K( J )

   20       CONTINUE
            IF( K( IN ).GT.0 ) GO TO 40

            for (II = 1; II <= M; II++) { // 30
               TEMP = X( II, J )
               X( II, J ) = X( II, IN )
               X( II, IN ) = TEMP
   30       CONTINUE

            K( IN ) = -K( IN )
            J = IN
            IN = K( IN )
            GO TO 20

   40       CONTINUE

   50    CONTINUE

      } else {

         // Backward permutation

         for (I = 1; I <= N; I++) { // 90

            IF( K( I ).GT.0 ) GO TO 80

            K( I ) = -K( I )
            J = K( I )
   60       CONTINUE
            IF( J.EQ.I ) GO TO 80

            for (II = 1; II <= M; II++) { // 70
               TEMP = X( II, I )
               X( II, I ) = X( II, J )
               X( II, J ) = TEMP
   70       CONTINUE

            K( J ) = -K( J )
            J = K( J )
            GO TO 60

   80       CONTINUE

   90    CONTINUE

      }

      RETURN

      // End of ZLAPMT

      }
