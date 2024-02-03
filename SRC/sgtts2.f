      SUBROUTINE SGTTS2( ITRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                ITRANS, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      REAL               B( LDB, * ), D( * ), DL( * ), DU( * ), DU2( * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      int                I, IP, J;
      REAL               TEMP
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      IF( N.EQ.0 .OR. NRHS.EQ.0 ) RETURN

      if ( ITRANS.EQ.0 ) {

         // Solve A*X = B using the LU factorization of A,
         // overwriting each right hand side vector with its solution.

         if ( NRHS.LE.1 ) {
            J = 1
   10       CONTINUE

            // Solve L*x = b.

            DO 20 I = 1, N - 1
               IP = IPIV( I )
               TEMP = B( I+1-IP+I, J ) - DL( I )*B( IP, J )
               B( I, J ) = B( IP, J )
               B( I+1, J ) = TEMP
   20       CONTINUE

            // Solve U*x = b.

            B( N, J ) = B( N, J ) / D( N )
            IF( N.GT.1 ) B( N-1, J ) = ( B( N-1, J )-DU( N-1 )*B( N, J ) ) / D( N-1 )
            DO 30 I = N - 2, 1, -1
               B( I, J ) = ( B( I, J )-DU( I )*B( I+1, J )-DU2( I )* B( I+2, J ) ) / D( I )
   30       CONTINUE
            if ( J.LT.NRHS ) {
               J = J + 1
               GO TO 10
            }
         } else {
            for (J = 1; J <= NRHS; J++) { // 60

               // Solve L*x = b.

               DO 40 I = 1, N - 1
                  if ( IPIV( I ).EQ.I ) {
                     B( I+1, J ) = B( I+1, J ) - DL( I )*B( I, J )
                  } else {
                     TEMP = B( I, J )
                     B( I, J ) = B( I+1, J )
                     B( I+1, J ) = TEMP - DL( I )*B( I, J )
                  }
   40          CONTINUE

               // Solve U*x = b.

               B( N, J ) = B( N, J ) / D( N )
               IF( N.GT.1 ) B( N-1, J ) = ( B( N-1, J )-DU( N-1 )*B( N, J ) ) / D( N-1 )
               DO 50 I = N - 2, 1, -1
                  B( I, J ) = ( B( I, J )-DU( I )*B( I+1, J )-DU2( I )* B( I+2, J ) ) / D( I )
   50          CONTINUE
   60       CONTINUE
         }
      } else {

         // Solve A**T * X = B.

         if ( NRHS.LE.1 ) {

            // Solve U**T*x = b.

            J = 1
   70       CONTINUE
            B( 1, J ) = B( 1, J ) / D( 1 )
            IF( N.GT.1 ) B( 2, J ) = ( B( 2, J )-DU( 1 )*B( 1, J ) ) / D( 2 )
            for (I = 3; I <= N; I++) { // 80
               B( I, J ) = ( B( I, J )-DU( I-1 )*B( I-1, J )-DU2( I-2 )* B( I-2, J ) ) / D( I )
   80       CONTINUE

            // Solve L**T*x = b.

            DO 90 I = N - 1, 1, -1
               IP = IPIV( I )
               TEMP = B( I, J ) - DL( I )*B( I+1, J )
               B( I, J ) = B( IP, J )
               B( IP, J ) = TEMP
   90       CONTINUE
            if ( J.LT.NRHS ) {
               J = J + 1
               GO TO 70
            }

         } else {
            for (J = 1; J <= NRHS; J++) { // 120

               // Solve U**T*x = b.

               B( 1, J ) = B( 1, J ) / D( 1 )
               IF( N.GT.1 ) B( 2, J ) = ( B( 2, J )-DU( 1 )*B( 1, J ) ) / D( 2 )
               for (I = 3; I <= N; I++) { // 100
                  B( I, J ) = ( B( I, J )-DU( I-1 )*B( I-1, J )- DU2( I-2 )*B( I-2, J ) ) / D( I )
  100          CONTINUE
               DO 110 I = N - 1, 1, -1
                  if ( IPIV( I ).EQ.I ) {
                     B( I, J ) = B( I, J ) - DL( I )*B( I+1, J )
                  } else {
                     TEMP = B( I+1, J )
                     B( I+1, J ) = B( I, J ) - DL( I )*TEMP
                     B( I, J ) = TEMP
                  }
  110          CONTINUE
  120       CONTINUE
         }
      }

      // End of SGTTS2

      }
