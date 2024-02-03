      SUBROUTINE CPTTS2( IUPLO, N, NRHS, D, E, B, LDB )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                IUPLO, LDB, N, NRHS;
      // ..
      // .. Array Arguments ..
      REAL               D( * )
      COMPLEX            B( LDB, * ), E( * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      int                I, J;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CSSCAL
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG
      // ..
      // .. Executable Statements ..

      // Quick return if possible

      if ( N.LE.1 ) {
         if (N.EQ.1) CALL CSSCAL( NRHS, 1. / D( 1 ), B, LDB );
         RETURN
      }

      if ( IUPLO.EQ.1 ) {

         // Solve A * X = B using the factorization A = U**H *D*U,
         // overwriting each right hand side vector with its solution.

         if ( NRHS.LE.2 ) {
            J = 1
            } // 5

            // Solve U**H * x = b.

            for (I = 2; I <= N; I++) { // 10
               B( I, J ) = B( I, J ) - B( I-1, J )*CONJG( E( I-1 ) )
            } // 10

            // Solve D * U * x = b.

            for (I = 1; I <= N; I++) { // 20
               B( I, J ) = B( I, J ) / D( I )
            } // 20
            DO 30 I = N - 1, 1, -1
               B( I, J ) = B( I, J ) - B( I+1, J )*E( I )
            } // 30
            if ( J.LT.NRHS ) {
               J = J + 1
               GO TO 5
            }
         } else {
            for (J = 1; J <= NRHS; J++) { // 60

               // Solve U**H * x = b.

               for (I = 2; I <= N; I++) { // 40
                  B( I, J ) = B( I, J ) - B( I-1, J )*CONJG( E( I-1 ) )
               } // 40

               // Solve D * U * x = b.

               B( N, J ) = B( N, J ) / D( N )
               DO 50 I = N - 1, 1, -1
                  B( I, J ) = B( I, J ) / D( I ) - B( I+1, J )*E( I )
               } // 50
            } // 60
         }
      } else {

         // Solve A * X = B using the factorization A = L*D*L**H,
         // overwriting each right hand side vector with its solution.

         if ( NRHS.LE.2 ) {
            J = 1
            } // 65

            // Solve L * x = b.

            for (I = 2; I <= N; I++) { // 70
               B( I, J ) = B( I, J ) - B( I-1, J )*E( I-1 )
            } // 70

            // Solve D * L**H * x = b.

            for (I = 1; I <= N; I++) { // 80
               B( I, J ) = B( I, J ) / D( I )
            } // 80
            DO 90 I = N - 1, 1, -1
               B( I, J ) = B( I, J ) - B( I+1, J )*CONJG( E( I ) )
            } // 90
            if ( J.LT.NRHS ) {
               J = J + 1
               GO TO 65
            }
         } else {
            for (J = 1; J <= NRHS; J++) { // 120

               // Solve L * x = b.

               for (I = 2; I <= N; I++) { // 100
                  B( I, J ) = B( I, J ) - B( I-1, J )*E( I-1 )
               } // 100

               // Solve D * L**H * x = b.

               B( N, J ) = B( N, J ) / D( N )
               DO 110 I = N - 1, 1, -1
                  B( I, J ) = B( I, J ) / D( I ) - B( I+1, J )*CONJG( E( I ) )
               } // 110
            } // 120
         }
      }

      RETURN

      // End of CPTTS2

      }
