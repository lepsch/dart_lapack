      SUBROUTINE ZGBEQUB( M, N, KL, KU, AB, LDAB, R, C, ROWCND, COLCND, AMAX, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, KL, KU, LDAB, M, N;
      double             AMAX, COLCND, ROWCND;
      // ..
      // .. Array Arguments ..
      double             C( * ), R( * );
      COMPLEX*16         AB( LDAB, * )
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      // ..
      // .. Local Scalars ..
      int                I, J, KD;
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
      // INTRINSIC ABS, MAX, MIN, LOG, REAL, DIMAG
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
      } else if ( KL.LT.0 ) {
         INFO = -3
      } else if ( KU.LT.0 ) {
         INFO = -4
      } else if ( LDAB.LT.KL+KU+1 ) {
         INFO = -6
      }
      if ( INFO.NE.0 ) {
         xerbla('ZGBEQUB', -INFO );
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
      LOGRDX = LOG(RADIX)

      // Compute row scale factors.

      for (I = 1; I <= M; I++) { // 10
         R( I ) = ZERO
   10 CONTINUE

      // Find the maximum element in each row.

      KD = KU + 1
      for (J = 1; J <= N; J++) { // 30
         DO 20 I = MAX( J-KU, 1 ), MIN( J+KL, M )
            R( I ) = MAX( R( I ), CABS1( AB( KD+I-J, J ) ) )
   20    CONTINUE
   30 CONTINUE
      for (I = 1; I <= M; I++) {
         if ( R( I ).GT.ZERO ) {
            R( I ) = RADIX**INT( LOG( R( I ) ) / LOGRDX )
         }
      END DO

      // Find the maximum and minimum scale factors.

      RCMIN = BIGNUM
      RCMAX = ZERO
      for (I = 1; I <= M; I++) { // 40
         RCMAX = MAX( RCMAX, R( I ) )
         RCMIN = MIN( RCMIN, R( I ) )
   40 CONTINUE
      AMAX = RCMAX

      if ( RCMIN.EQ.ZERO ) {

         // Find the first zero scale factor and return an error code.

         for (I = 1; I <= M; I++) { // 50
            if ( R( I ).EQ.ZERO ) {
               INFO = I
               RETURN
            }
   50    CONTINUE
      } else {

         // Invert the scale factors.

         for (I = 1; I <= M; I++) { // 60
            R( I ) = ONE / MIN( MAX( R( I ), SMLNUM ), BIGNUM )
   60    CONTINUE

         // Compute ROWCND = min(R(I)) / max(R(I)).

         ROWCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
      }

      // Compute column scale factors.

      for (J = 1; J <= N; J++) { // 70
         C( J ) = ZERO
   70 CONTINUE

      // Find the maximum element in each column,
      // assuming the row scaling computed above.

      for (J = 1; J <= N; J++) { // 90
         DO 80 I = MAX( J-KU, 1 ), MIN( J+KL, M )
            C( J ) = MAX( C( J ), CABS1( AB( KD+I-J, J ) )*R( I ) )
   80    CONTINUE
         if ( C( J ).GT.ZERO ) {
            C( J ) = RADIX**INT( LOG( C( J ) ) / LOGRDX )
         }
   90 CONTINUE

      // Find the maximum and minimum scale factors.

      RCMIN = BIGNUM
      RCMAX = ZERO
      for (J = 1; J <= N; J++) { // 100
         RCMIN = MIN( RCMIN, C( J ) )
         RCMAX = MAX( RCMAX, C( J ) )
  100 CONTINUE

      if ( RCMIN.EQ.ZERO ) {

         // Find the first zero scale factor and return an error code.

         for (J = 1; J <= N; J++) { // 110
            if ( C( J ).EQ.ZERO ) {
               INFO = M + J
               RETURN
            }
  110    CONTINUE
      } else {

         // Invert the scale factors.

         for (J = 1; J <= N; J++) { // 120
            C( J ) = ONE / MIN( MAX( C( J ), SMLNUM ), BIGNUM )
  120    CONTINUE

         // Compute COLCND = min(C(J)) / max(C(J)).

         COLCND = MAX( RCMIN, SMLNUM ) / MIN( RCMAX, BIGNUM )
      }

      RETURN

      // End of ZGBEQUB

      }
