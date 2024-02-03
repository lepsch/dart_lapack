      SUBROUTINE CLAGTM( TRANS, N, NRHS, ALPHA, DL, D, DU, X, LDX, BETA, B, LDB )

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             TRANS;
      int                LDB, LDX, N, NRHS;
      REAL               ALPHA, BETA
      // ..
      // .. Array Arguments ..
      COMPLEX            B( LDB, * ), D( * ), DL( * ), DU( * ), X( LDX, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      int                I, J;
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG
      // ..
      // .. Executable Statements ..

      IF( N.EQ.0 ) RETURN

      // Multiply B by BETA if BETA.NE.1.

      if ( BETA.EQ.ZERO ) {
         DO 20 J = 1, NRHS
            DO 10 I = 1, N
               B( I, J ) = ZERO
   10       CONTINUE
   20    CONTINUE
      } else if ( BETA.EQ.-ONE ) {
         DO 40 J = 1, NRHS
            DO 30 I = 1, N
               B( I, J ) = -B( I, J )
   30       CONTINUE
   40    CONTINUE
      }

      if ( ALPHA.EQ.ONE ) {
         if ( LSAME( TRANS, 'N' ) ) {

            // Compute B := B + A*X

            DO 60 J = 1, NRHS
               if ( N.EQ.1 ) {
                  B( 1, J ) = B( 1, J ) + D( 1 )*X( 1, J )
               } else {
                  B( 1, J ) = B( 1, J ) + D( 1 )*X( 1, J ) + DU( 1 )*X( 2, J )                   B( N, J ) = B( N, J ) + DL( N-1 )*X( N-1, J ) + D( N )*X( N, J )
                  DO 50 I = 2, N - 1
                     B( I, J ) = B( I, J ) + DL( I-1 )*X( I-1, J ) + D( I )*X( I, J ) + DU( I )*X( I+1, J )
   50             CONTINUE
               }
   60       CONTINUE
         } else if ( LSAME( TRANS, 'T' ) ) {

            // Compute B := B + A**T * X

            DO 80 J = 1, NRHS
               if ( N.EQ.1 ) {
                  B( 1, J ) = B( 1, J ) + D( 1 )*X( 1, J )
               } else {
                  B( 1, J ) = B( 1, J ) + D( 1 )*X( 1, J ) + DL( 1 )*X( 2, J )                   B( N, J ) = B( N, J ) + DU( N-1 )*X( N-1, J ) + D( N )*X( N, J )
                  DO 70 I = 2, N - 1
                     B( I, J ) = B( I, J ) + DU( I-1 )*X( I-1, J ) + D( I )*X( I, J ) + DL( I )*X( I+1, J )
   70             CONTINUE
               }
   80       CONTINUE
         } else if ( LSAME( TRANS, 'C' ) ) {

            // Compute B := B + A**H * X

            DO 100 J = 1, NRHS
               if ( N.EQ.1 ) {
                  B( 1, J ) = B( 1, J ) + CONJG( D( 1 ) )*X( 1, J )
               } else {
                  B( 1, J ) = B( 1, J ) + CONJG( D( 1 ) )*X( 1, J ) + CONJG( DL( 1 ) )*X( 2, J )                   B( N, J ) = B( N, J ) + CONJG( DU( N-1 ) )* X( N-1, J ) + CONJG( D( N ) )*X( N, J )
                  DO 90 I = 2, N - 1
                     B( I, J ) = B( I, J ) + CONJG( DU( I-1 ) )* X( I-1, J ) + CONJG( D( I ) )* X( I, J ) + CONJG( DL( I ) )* X( I+1, J )
   90             CONTINUE
               }
  100       CONTINUE
         }
      } else if ( ALPHA.EQ.-ONE ) {
         if ( LSAME( TRANS, 'N' ) ) {

            // Compute B := B - A*X

            DO 120 J = 1, NRHS
               if ( N.EQ.1 ) {
                  B( 1, J ) = B( 1, J ) - D( 1 )*X( 1, J )
               } else {
                  B( 1, J ) = B( 1, J ) - D( 1 )*X( 1, J ) - DU( 1 )*X( 2, J )                   B( N, J ) = B( N, J ) - DL( N-1 )*X( N-1, J ) - D( N )*X( N, J )
                  DO 110 I = 2, N - 1
                     B( I, J ) = B( I, J ) - DL( I-1 )*X( I-1, J ) - D( I )*X( I, J ) - DU( I )*X( I+1, J )
  110             CONTINUE
               }
  120       CONTINUE
         } else if ( LSAME( TRANS, 'T' ) ) {

            // Compute B := B - A**T*X

            DO 140 J = 1, NRHS
               if ( N.EQ.1 ) {
                  B( 1, J ) = B( 1, J ) - D( 1 )*X( 1, J )
               } else {
                  B( 1, J ) = B( 1, J ) - D( 1 )*X( 1, J ) - DL( 1 )*X( 2, J )                   B( N, J ) = B( N, J ) - DU( N-1 )*X( N-1, J ) - D( N )*X( N, J )
                  DO 130 I = 2, N - 1
                     B( I, J ) = B( I, J ) - DU( I-1 )*X( I-1, J ) - D( I )*X( I, J ) - DL( I )*X( I+1, J )
  130             CONTINUE
               }
  140       CONTINUE
         } else if ( LSAME( TRANS, 'C' ) ) {

            // Compute B := B - A**H*X

            DO 160 J = 1, NRHS
               if ( N.EQ.1 ) {
                  B( 1, J ) = B( 1, J ) - CONJG( D( 1 ) )*X( 1, J )
               } else {
                  B( 1, J ) = B( 1, J ) - CONJG( D( 1 ) )*X( 1, J ) - CONJG( DL( 1 ) )*X( 2, J )                   B( N, J ) = B( N, J ) - CONJG( DU( N-1 ) )* X( N-1, J ) - CONJG( D( N ) )*X( N, J )
                  DO 150 I = 2, N - 1
                     B( I, J ) = B( I, J ) - CONJG( DU( I-1 ) )* X( I-1, J ) - CONJG( D( I ) )* X( I, J ) - CONJG( DL( I ) )* X( I+1, J )
  150             CONTINUE
               }
  160       CONTINUE
         }
      }
      RETURN

      // End of CLAGTM

      }
