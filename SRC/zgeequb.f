      SUBROUTINE ZGEEQUB( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, M, N;
      double             AMAX, COLCND, ROWCND;
      // ..
      // .. Array Arguments ..
      double             C( * ), R( * );
      COMPLEX*16         A( LDA, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                I, J;
      double             BIGNUM, RCMAX, RCMIN, SMLNUM, RADIX, LOGRDX;
      COMPLEX*16         ZDUM
      // ..
      // .. External Functions ..
      double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, LOG, DBLE, DIMAG
      // ..
      // .. Statement Functions ..
      double             CABS1;
      // ..
      // .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( DBLE( ZDUM ) ) + ABS( DIMAG( ZDUM ) )
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      if ( M.LT.0 ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -4
      }
      if ( INFO.NE.0 ) {
         xerbla('ZGEEQUB', -INFO );
         RETURN
      }

      // Quick return if possible.

      if ( M.EQ.0 .OR. N.EQ.0 ) {
         ROWCND = ONE
         COLCND = ONE
         AMAX = ZERO
         RETURN
      }

      // Get machine constants.  Assume SMLNUM is a power of the radix.

      SMLNUM = DLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM
      RADIX = DLAMCH( 'B' )
      LOGRDX = LOG( RADIX )

      // Compute row scale factors.

      DO 10 I = 1, M
         R( I ) = ZERO
   10 CONTINUE

      // Find the maximum element in each row.

      DO 30 J = 1, N
         DO 20 I = 1, M
            R( I ) = MAX( R( I ), CABS1( A( I, J ) ) )
   20    CONTINUE
   30 CONTINUE
      DO I = 1, M
         if ( R( I ).GT.ZERO ) {
            R( I ) = RADIX**INT( LOG(R( I ) ) / LOGRDX )
         }
      END DO

      // Find the maximum and minimum scale factors.

      RCMIN = BIGNUM
      RCMAX = ZERO
      DO 40 I = 1, M
         RCMAX = MAX( RCMAX, R( I ) )
         RCMIN = MIN( RCMIN, R( I ) )
   40 CONTINUE
      AMAX = RCMAX

      if ( RCMIN.EQ.ZERO ) {

         // Find the first zero scale factor and return an error code.

         DO 50 I = 1, M
            if ( R( I ).EQ.ZERO ) {
               INFO = I
               RETURN
            }
   50    CONTINUE
      } else {

         // Invert the scale factors.

         DO 60 I = 1, M
            R( I ) = ONE / MIN( MAX( R( I ), SMLNUM ), BIGNUM )
   60    CONTINUE

         // Compute ROWCND = min(R(I)) / max(R(I)).

         ROWCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
      }

      // Compute column scale factors.

      DO 70 J = 1, N
         C( J ) = ZERO
   70 CONTINUE

      // Find the maximum element in each column,
      // assuming the row scaling computed above.

      DO 90 J = 1, N
         DO 80 I = 1, M
            C( J ) = MAX( C( J ), CABS1( A( I, J ) )*R( I ) )
   80    CONTINUE
         if ( C( J ).GT.ZERO ) {
            C( J ) = RADIX**INT( LOG( C( J ) ) / LOGRDX )
         }
   90 CONTINUE

      // Find the maximum and minimum scale factors.

      RCMIN = BIGNUM
      RCMAX = ZERO
      DO 100 J = 1, N
         RCMIN = MIN( RCMIN, C( J ) )
         RCMAX = MAX( RCMAX, C( J ) )
  100 CONTINUE

      if ( RCMIN.EQ.ZERO ) {

         // Find the first zero scale factor and return an error code.

         DO 110 J = 1, N
            if ( C( J ).EQ.ZERO ) {
               INFO = M + J
               RETURN
            }
  110    CONTINUE
      } else {

         // Invert the scale factors.

         DO 120 J = 1, N
            C( J ) = ONE / MIN( MAX( C( J ), SMLNUM ), BIGNUM )
  120    CONTINUE

         // Compute COLCND = min(C(J)) / max(C(J)).

         COLCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
      }

      RETURN

      // End of ZGEEQUB

      }
