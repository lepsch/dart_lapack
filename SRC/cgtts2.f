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

      IF( N.EQ.0 .OR. NRHS.EQ.0 ) RETURN

      if ( ITRANS.EQ.0 ) {

         // Solve A*X = B using the LU factorization of A,
         // overwriting each right hand side vector with its solution.

         if ( NRHS.LE.1 ) {
            J = 1
   10       CONTINUE

            // Solve L*x = b.

            DO 20 I = 1, N - 1
               if ( IPIV( I ).EQ.I ) {
                  B( I+1, J ) = B( I+1, J ) - DL( I )*B( I, J )
               } else {
                  TEMP = B( I, J )
                  B( I, J ) = B( I+1, J )
                  B( I+1, J ) = TEMP - DL( I )*B( I, J )
               }
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
            DO 60 J = 1, NRHS

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
      } else if ( ITRANS.EQ.1 ) {

         // Solve A**T * X = B.

         if ( NRHS.LE.1 ) {
            J = 1
   70       CONTINUE

            // Solve U**T * x = b.

            B( 1, J ) = B( 1, J ) / D( 1 )
            IF( N.GT.1 ) B( 2, J ) = ( B( 2, J )-DU( 1 )*B( 1, J ) ) / D( 2 )
            DO 80 I = 3, N
               B( I, J ) = ( B( I, J )-DU( I-1 )*B( I-1, J )-DU2( I-2 )* B( I-2, J ) ) / D( I )
   80       CONTINUE

            // Solve L**T * x = b.

            DO 90 I = N - 1, 1, -1
               if ( IPIV( I ).EQ.I ) {
                  B( I, J ) = B( I, J ) - DL( I )*B( I+1, J )
               } else {
                  TEMP = B( I+1, J )
                  B( I+1, J ) = B( I, J ) - DL( I )*TEMP
                  B( I, J ) = TEMP
               }
   90       CONTINUE
            if ( J.LT.NRHS ) {
               J = J + 1
               GO TO 70
            }
         } else {
            DO 120 J = 1, NRHS

            // Solve U**T * x = b.

               B( 1, J ) = B( 1, J ) / D( 1 )
               IF( N.GT.1 ) B( 2, J ) = ( B( 2, J )-DU( 1 )*B( 1, J ) ) / D( 2 )
               DO 100 I = 3, N
                  B( I, J ) = ( B( I, J )-DU( I-1 )*B( I-1, J )- DU2( I-2 )*B( I-2, J ) ) / D( I )
  100          CONTINUE

            // Solve L**T * x = b.

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
      } else {

         // Solve A**H * X = B.

         if ( NRHS.LE.1 ) {
            J = 1
  130       CONTINUE

            // Solve U**H * x = b.

            B( 1, J ) = B( 1, J ) / CONJG( D( 1 ) )
            IF( N.GT.1 ) B( 2, J ) = ( B( 2, J )-CONJG( DU( 1 ) )*B( 1, J ) ) / CONJG( D( 2 ) )
            DO 140 I = 3, N
               B( I, J ) = ( B( I, J )-CONJG( DU( I-1 ) )*B( I-1, J )- CONJG( DU2( I-2 ) )*B( I-2, J ) ) / CONJG( D( I ) )
  140       CONTINUE

            // Solve L**H * x = b.

            DO 150 I = N - 1, 1, -1
               if ( IPIV( I ).EQ.I ) {
                  B( I, J ) = B( I, J ) - CONJG( DL( I ) )*B( I+1, J )
               } else {
                  TEMP = B( I+1, J )
                  B( I+1, J ) = B( I, J ) - CONJG( DL( I ) )*TEMP
                  B( I, J ) = TEMP
               }
  150       CONTINUE
            if ( J.LT.NRHS ) {
               J = J + 1
               GO TO 130
            }
         } else {
            DO 180 J = 1, NRHS

            // Solve U**H * x = b.

               B( 1, J ) = B( 1, J ) / CONJG( D( 1 ) )
               IF( N.GT.1 ) B( 2, J ) = ( B( 2, J )-CONJG( DU( 1 ) )*B( 1, J ) ) / CONJG( D( 2 ) )
               DO 160 I = 3, N
                  B( I, J ) = ( B( I, J )-CONJG( DU( I-1 ) )* B( I-1, J )-CONJG( DU2( I-2 ) )* B( I-2, J ) ) / CONJG( D( I ) )
  160          CONTINUE

            // Solve L**H * x = b.

               DO 170 I = N - 1, 1, -1
                  if ( IPIV( I ).EQ.I ) {
                     B( I, J ) = B( I, J ) - CONJG( DL( I ) )* B( I+1, J )
                  } else {
                     TEMP = B( I+1, J )
                     B( I+1, J ) = B( I, J ) - CONJG( DL( I ) )*TEMP
                     B( I, J ) = TEMP
                  }
  170          CONTINUE
  180       CONTINUE
         }
      }

      // End of CGTTS2

      }
