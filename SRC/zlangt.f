      double           FUNCTION ZLANGT( NORM, N, DL, D, DU );

*  -- LAPACK auxiliary routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             NORM;
      int                N;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         D( * ), DL( * ), DU( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                I;
      double             ANORM, SCALE, SUM, TEMP;
      // ..
      // .. External Functions ..
      bool               LSAME, DISNAN;
      // EXTERNAL LSAME, DISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZLASSQ
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, SQRT
      // ..
      // .. Executable Statements ..

      if ( N.LE.0 ) {
         ANORM = ZERO
      } else if ( LSAME( NORM, 'M' ) ) {

         // Find max(abs(A(i,j))).

         ANORM = ABS( D( N ) )
         DO 10 I = 1, N - 1
            IF( ANORM.LT.ABS( DL( I ) ) .OR. DISNAN( ABS( DL( I ) ) ) ) ANORM = ABS(DL(I))             IF( ANORM.LT.ABS( D( I ) ) .OR. DISNAN( ABS( D( I ) ) ) ) ANORM = ABS(D(I))             IF( ANORM.LT.ABS( DU( I ) ) .OR. DISNAN (ABS( DU( I ) ) ) ) ANORM = ABS(DU(I))
   10    CONTINUE
      } else if ( LSAME( NORM, 'O' ) .OR. NORM.EQ.'1' ) {

         // Find norm1(A).

         if ( N.EQ.1 ) {
            ANORM = ABS( D( 1 ) )
         } else {
            ANORM = ABS( D( 1 ) )+ABS( DL( 1 ) )
            TEMP = ABS( D( N ) )+ABS( DU( N-1 ) )
            IF( ANORM .LT. TEMP .OR. DISNAN( TEMP ) ) ANORM = TEMP
            DO 20 I = 2, N - 1
               TEMP = ABS( D( I ) )+ABS( DL( I ) )+ABS( DU( I-1 ) )
               IF( ANORM .LT. TEMP .OR. DISNAN( TEMP ) ) ANORM = TEMP
   20       CONTINUE
         }
      } else if ( LSAME( NORM, 'I' ) ) {

         // Find normI(A).

         if ( N.EQ.1 ) {
            ANORM = ABS( D( 1 ) )
         } else {
            ANORM = ABS( D( 1 ) )+ABS( DU( 1 ) )
            TEMP = ABS( D( N ) )+ABS( DL( N-1 ) )
            IF( ANORM .LT. TEMP .OR. DISNAN( TEMP ) ) ANORM = TEMP
            DO 30 I = 2, N - 1
               TEMP = ABS( D( I ) )+ABS( DU( I ) )+ABS( DL( I-1 ) )
               IF( ANORM .LT. TEMP .OR. DISNAN( TEMP ) ) ANORM = TEMP
   30       CONTINUE
         }
      } else if ( ( LSAME( NORM, 'F' ) ) .OR. ( LSAME( NORM, 'E' ) ) ) {

         // Find normF(A).

         SCALE = ZERO
         SUM = ONE
         zlassq(N, D, 1, SCALE, SUM );
         if ( N.GT.1 ) {
            zlassq(N-1, DL, 1, SCALE, SUM );
            zlassq(N-1, DU, 1, SCALE, SUM );
         }
         ANORM = SCALE*SQRT( SUM )
      }

      ZLANGT = ANORM
      RETURN

      // End of ZLANGT

      }
