      SUBROUTINE CGEEQU( M, N, A, LDA, R, C, ROWCND, COLCND, AMAX, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, M, N;
      REAL               AMAX, COLCND, ROWCND
      // ..
      // .. Array Arguments ..
      REAL               C( * ), R( * )
      COMPLEX            A( LDA, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      // ..
      // .. Local Scalars ..
      int                I, J;
      REAL               BIGNUM, RCMAX, RCMIN, SMLNUM
      COMPLEX            ZDUM
      // ..
      // .. External Functions ..
      REAL               SLAMCH
      // EXTERNAL SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, MAX, MIN, REAL
      // ..
      // .. Statement Functions ..
      REAL               CABS1
      // ..
      // .. Statement Function definitions ..
      CABS1( ZDUM ) = ABS( REAL( ZDUM ) ) + ABS( AIMAG( ZDUM ) )
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
         CALL XERBLA( 'CGEEQU', -INFO )
         RETURN
      }

      // Quick return if possible

      if ( M.EQ.0 .OR. N.EQ.0 ) {
         ROWCND = ONE
         COLCND = ONE
         AMAX = ZERO
         RETURN
      }

      // Get machine constants.

      SMLNUM = SLAMCH( 'S' )
      BIGNUM = ONE / SMLNUM

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

         // Compute ROWCND = min(R(I)) / max(R(I))

         ROWCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
      }

      // Compute column scale factors

      DO 70 J = 1, N
         C( J ) = ZERO
   70 CONTINUE

      // Find the maximum element in each column,
      // assuming the row scaling computed above.

      DO 90 J = 1, N
         DO 80 I = 1, M
            C( J ) = MAX( C( J ), CABS1( A( I, J ) )*R( I ) )
   80    CONTINUE
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

         // Compute COLCND = min(C(J)) / max(C(J))

         COLCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
      }

      RETURN

      // End of CGEEQU

      }
