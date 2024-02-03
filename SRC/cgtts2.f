      SUBROUTINE CGTTS2( ITRANS, N, NRHS, DL, D, DU, DU2, IPIV, B, LDB )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                ITRANS, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX            B( LDB, * ), D( * ), DL( * ), DU( * ), DU2( * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      int                I, J;
      COMPLEX            TEMP
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if (N == 0 || NRHS == 0) RETURN;

      if ( ITRANS == 0 ) {

         // Solve A*X = B using the LU factorization of A,
         // overwriting each right hand side vector with its solution.

         if ( NRHS <= 1 ) {
            J = 1
            } // 10

            // Solve L*x = b.

            for (I = 1; I <= N - 1; I++) { // 20
               if ( IPIV( I ) == I ) {
                  B( I+1, J ) = B( I+1, J ) - DL( I )*B( I, J )
               } else {
                  TEMP = B( I, J )
                  B( I, J ) = B( I+1, J )
                  B( I+1, J ) = TEMP - DL( I )*B( I, J )
               }
            } // 20

            // Solve U*x = b.

            B( N, J ) = B( N, J ) / D( N )
            if (N > 1) B( N-1, J ) = ( B( N-1, J )-DU( N-1 )*B( N, J ) ) / D( N-1 );
            DO 30 I = N - 2, 1, -1
               B( I, J ) = ( B( I, J )-DU( I )*B( I+1, J )-DU2( I )* B( I+2, J ) ) / D( I )
            } // 30
            if ( J < NRHS ) {
               J = J + 1
               GO TO 10
            }
         } else {
            for (J = 1; J <= NRHS; J++) { // 60

            // Solve L*x = b.

               for (I = 1; I <= N - 1; I++) { // 40
                  if ( IPIV( I ) == I ) {
                     B( I+1, J ) = B( I+1, J ) - DL( I )*B( I, J )
                  } else {
                     TEMP = B( I, J )
                     B( I, J ) = B( I+1, J )
                     B( I+1, J ) = TEMP - DL( I )*B( I, J )
                  }
               } // 40

            // Solve U*x = b.

               B( N, J ) = B( N, J ) / D( N )
               if (N > 1) B( N-1, J ) = ( B( N-1, J )-DU( N-1 )*B( N, J ) ) / D( N-1 );
               DO 50 I = N - 2, 1, -1
                  B( I, J ) = ( B( I, J )-DU( I )*B( I+1, J )-DU2( I )* B( I+2, J ) ) / D( I )
               } // 50
            } // 60
         }
      } else if ( ITRANS == 1 ) {

         // Solve A**T * X = B.

         if ( NRHS <= 1 ) {
            J = 1
            } // 70

            // Solve U**T * x = b.

            B( 1, J ) = B( 1, J ) / D( 1 )
            if (N > 1) B( 2, J ) = ( B( 2, J )-DU( 1 )*B( 1, J ) ) / D( 2 );
            for (I = 3; I <= N; I++) { // 80
               B( I, J ) = ( B( I, J )-DU( I-1 )*B( I-1, J )-DU2( I-2 )* B( I-2, J ) ) / D( I )
            } // 80

            // Solve L**T * x = b.

            DO 90 I = N - 1, 1, -1
               if ( IPIV( I ) == I ) {
                  B( I, J ) = B( I, J ) - DL( I )*B( I+1, J )
               } else {
                  TEMP = B( I+1, J )
                  B( I+1, J ) = B( I, J ) - DL( I )*TEMP
                  B( I, J ) = TEMP
               }
            } // 90
            if ( J < NRHS ) {
               J = J + 1
               GO TO 70
            }
         } else {
            for (J = 1; J <= NRHS; J++) { // 120

            // Solve U**T * x = b.

               B( 1, J ) = B( 1, J ) / D( 1 )
               if (N > 1) B( 2, J ) = ( B( 2, J )-DU( 1 )*B( 1, J ) ) / D( 2 );
               for (I = 3; I <= N; I++) { // 100
                  B( I, J ) = ( B( I, J )-DU( I-1 )*B( I-1, J )- DU2( I-2 )*B( I-2, J ) ) / D( I )
               } // 100

            // Solve L**T * x = b.

               DO 110 I = N - 1, 1, -1
                  if ( IPIV( I ) == I ) {
                     B( I, J ) = B( I, J ) - DL( I )*B( I+1, J )
                  } else {
                     TEMP = B( I+1, J )
                     B( I+1, J ) = B( I, J ) - DL( I )*TEMP
                     B( I, J ) = TEMP
                  }
               } // 110
            } // 120
         }
      } else {

         // Solve A**H * X = B.

         if ( NRHS <= 1 ) {
            J = 1
            } // 130

            // Solve U**H * x = b.

            B( 1, J ) = B( 1, J ) / CONJG( D( 1 ) )
            if (N > 1) B( 2, J ) = ( B( 2, J )-CONJG( DU( 1 ) )*B( 1, J ) ) / CONJG( D( 2 ) );
            for (I = 3; I <= N; I++) { // 140
               B( I, J ) = ( B( I, J )-CONJG( DU( I-1 ) )*B( I-1, J )- CONJG( DU2( I-2 ) )*B( I-2, J ) ) / CONJG( D( I ) )
            } // 140

            // Solve L**H * x = b.

            DO 150 I = N - 1, 1, -1
               if ( IPIV( I ) == I ) {
                  B( I, J ) = B( I, J ) - CONJG( DL( I ) )*B( I+1, J )
               } else {
                  TEMP = B( I+1, J )
                  B( I+1, J ) = B( I, J ) - CONJG( DL( I ) )*TEMP
                  B( I, J ) = TEMP
               }
            } // 150
            if ( J < NRHS ) {
               J = J + 1
               GO TO 130
            }
         } else {
            for (J = 1; J <= NRHS; J++) { // 180

            // Solve U**H * x = b.

               B( 1, J ) = B( 1, J ) / CONJG( D( 1 ) )
               if (N > 1) B( 2, J ) = ( B( 2, J )-CONJG( DU( 1 ) )*B( 1, J ) ) / CONJG( D( 2 ) );
               for (I = 3; I <= N; I++) { // 160
                  B( I, J ) = ( B( I, J )-CONJG( DU( I-1 ) )* B( I-1, J )-CONJG( DU2( I-2 ) )* B( I-2, J ) ) / CONJG( D( I ) )
               } // 160

            // Solve L**H * x = b.

               DO 170 I = N - 1, 1, -1
                  if ( IPIV( I ) == I ) {
                     B( I, J ) = B( I, J ) - CONJG( DL( I ) )* B( I+1, J )
                  } else {
                     TEMP = B( I+1, J )
                     B( I+1, J ) = B( I, J ) - CONJG( DL( I ) )*TEMP
                     B( I, J ) = TEMP
                  }
               } // 170
            } // 180
         }
      }

      // End of CGTTS2

      }
