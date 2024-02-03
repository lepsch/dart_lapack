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
         IF( N.EQ.1 ) CALL CSSCAL( NRHS, 1. / D( 1 ), B, LDB )
         RETURN
      }

      if ( IUPLO.EQ.1 ) {

         // Solve A * X = B using the factorization A = U**H *D*U,
         // overwriting each right hand side vector with its solution.

         if ( NRHS.LE.2 ) {
            J = 1
    5       CONTINUE

            // Solve U**H * x = b.

            DO 10 I = 2, N
               B( I, J ) = B( I, J ) - B( I-1, J )*CONJG( E( I-1 ) )
   10       CONTINUE

            // Solve D * U * x = b.

            DO 20 I = 1, N
               B( I, J ) = B( I, J ) / D( I )
   20       CONTINUE
            DO 30 I = N - 1, 1, -1
               B( I, J ) = B( I, J ) - B( I+1, J )*E( I )
   30       CONTINUE
            if ( J.LT.NRHS ) {
               J = J + 1
               GO TO 5
            }
         } else {
            DO 60 J = 1, NRHS

               // Solve U**H * x = b.

               DO 40 I = 2, N
                  B( I, J ) = B( I, J ) - B( I-1, J )*CONJG( E( I-1 ) )
   40          CONTINUE

               // Solve D * U * x = b.

               B( N, J ) = B( N, J ) / D( N )
               DO 50 I = N - 1, 1, -1
                  B( I, J ) = B( I, J ) / D( I ) - B( I+1, J )*E( I )
   50          CONTINUE
   60       CONTINUE
         }
      } else {

         // Solve A * X = B using the factorization A = L*D*L**H,
         // overwriting each right hand side vector with its solution.

         if ( NRHS.LE.2 ) {
            J = 1
   65       CONTINUE

            // Solve L * x = b.

            DO 70 I = 2, N
               B( I, J ) = B( I, J ) - B( I-1, J )*E( I-1 )
   70       CONTINUE

            // Solve D * L**H * x = b.

            DO 80 I = 1, N
               B( I, J ) = B( I, J ) / D( I )
   80       CONTINUE
            DO 90 I = N - 1, 1, -1
               B( I, J ) = B( I, J ) - B( I+1, J )*CONJG( E( I ) )
   90       CONTINUE
            if ( J.LT.NRHS ) {
               J = J + 1
               GO TO 65
            }
         } else {
            DO 120 J = 1, NRHS

               // Solve L * x = b.

               DO 100 I = 2, N
                  B( I, J ) = B( I, J ) - B( I-1, J )*E( I-1 )
  100          CONTINUE

               // Solve D * L**H * x = b.

               B( N, J ) = B( N, J ) / D( N )
               DO 110 I = N - 1, 1, -1
                  B( I, J ) = B( I, J ) / D( I ) - B( I+1, J )*CONJG( E( I ) )
  110          CONTINUE
  120       CONTINUE
         }
      }

      RETURN

      // End of CPTTS2

      }
