      SUBROUTINE DSYEQUB( UPLO, N, A, LDA, S, SCOND, AMAX, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, N;
      double             AMAX, SCOND;
      String             UPLO;
      // ..
      // .. Array Arguments ..
      double             A( LDA, * ), S( * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D0, ZERO = 0.0D0 ;
      int                MAX_ITER;
      const              MAX_ITER = 100 ;
      // ..
      // .. Local Scalars ..
      int                I, J, ITER;
      double             AVG, STD, TOL, C0, C1, C2, T, U, SI, D, BASE, SMIN, SMAX, SMLNUM, BIGNUM, SCALE, SUMSQ;
      bool               UP;
      // ..
      // .. External Functions ..
      double             DLAMCH;
      bool               LSAME;
      // EXTERNAL DLAMCH, LSAME
      // ..
      // .. External Subroutines ..
      // EXTERNAL DLASSQ, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, INT, LOG, MAX, MIN, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      if ( .NOT. ( LSAME( UPLO, 'U' ) .OR. LSAME( UPLO, 'L' ) ) ) {
         INFO = -1
      } else if ( N .LT. 0 ) {
         INFO = -2
      } else if ( LDA .LT. MAX( 1, N ) ) {
         INFO = -4
      }
      if ( INFO .NE. 0 ) {
         xerbla('DSYEQUB', -INFO );
         RETURN
      }

      UP = LSAME( UPLO, 'U' )
      AMAX = ZERO

      // Quick return if possible.

      if ( N .EQ. 0 ) {
         SCOND = ONE
         RETURN
      }

      for (I = 1; I <= N; I++) {
         S( I ) = ZERO
      END DO

      AMAX = ZERO
      if ( UP ) {
         for (J = 1; J <= N; J++) {
            DO I = 1, J-1
               S( I ) = MAX( S( I ), ABS( A( I, J ) ) )
               S( J ) = MAX( S( J ), ABS( A( I, J ) ) )
               AMAX = MAX( AMAX, ABS( A( I, J ) ) )
            END DO
            S( J ) = MAX( S( J ), ABS( A( J, J ) ) )
            AMAX = MAX( AMAX, ABS( A( J, J ) ) )
         END DO
      } else {
         for (J = 1; J <= N; J++) {
            S( J ) = MAX( S( J ), ABS( A( J, J ) ) )
            AMAX = MAX( AMAX, ABS( A( J, J ) ) )
            DO I = J+1, N
               S( I ) = MAX( S( I ), ABS( A( I, J ) ) )
               S( J ) = MAX( S( J ), ABS( A( I, J ) ) )
               AMAX = MAX( AMAX, ABS( A( I, J ) ) )
            END DO
         END DO
      }
      for (J = 1; J <= N; J++) {
         S( J ) = 1.0D0 / S( J )
      END DO

      TOL = ONE / SQRT( 2.0D0 * N )

      for (ITER = 1; ITER <= MAX_ITER; ITER++) {
         SCALE = 0.0D0
         SUMSQ = 0.0D0
         // beta = |A|s
         for (I = 1; I <= N; I++) {
            WORK( I ) = ZERO
         END DO
         if ( UP ) {
            for (J = 1; J <= N; J++) {
               DO I = 1, J-1
                  WORK( I ) = WORK( I ) + ABS( A( I, J ) ) * S( J )
                  WORK( J ) = WORK( J ) + ABS( A( I, J ) ) * S( I )
               END DO
               WORK( J ) = WORK( J ) + ABS( A( J, J ) ) * S( J )
            END DO
         } else {
            for (J = 1; J <= N; J++) {
               WORK( J ) = WORK( J ) + ABS( A( J, J ) ) * S( J )
               DO I = J+1, N
                  WORK( I ) = WORK( I ) + ABS( A( I, J ) ) * S( J )
                  WORK( J ) = WORK( J ) + ABS( A( I, J ) ) * S( I )
               END DO
            END DO
         }

         // avg = s^T beta / n
         AVG = 0.0D0
         for (I = 1; I <= N; I++) {
            AVG = AVG + S( I )*WORK( I )
         END DO
         AVG = AVG / N

         STD = 0.0D0
         DO I = N+1, 2*N
            WORK( I ) = S( I-N ) * WORK( I-N ) - AVG
         END DO
         dlassq(N, WORK( N+1 ), 1, SCALE, SUMSQ );
         STD = SCALE * SQRT( SUMSQ / N )

         IF ( STD .LT. TOL * AVG ) GOTO 999

         for (I = 1; I <= N; I++) {
            T = ABS( A( I, I ) )
            SI = S( I )
            C2 = ( N-1 ) * T
            C1 = ( N-2 ) * ( WORK( I ) - T*SI )
            C0 = -(T*SI)*SI + 2*WORK( I )*SI - N*AVG
            D = C1*C1 - 4*C0*C2

            if ( D .LE. 0 ) {
               INFO = -1
               RETURN
            }
            SI = -2*C0 / ( C1 + SQRT( D ) )

            D = SI - S( I )
            U = ZERO
            if ( UP ) {
               for (J = 1; J <= I; J++) {
                  T = ABS( A( J, I ) )
                  U = U + S( J )*T
                  WORK( J ) = WORK( J ) + D*T
               END DO
               DO J = I+1,N
                  T = ABS( A( I, J ) )
                  U = U + S( J )*T
                  WORK( J ) = WORK( J ) + D*T
               END DO
            } else {
               for (J = 1; J <= I; J++) {
                  T = ABS( A( I, J ) )
                  U = U + S( J )*T
                  WORK( J ) = WORK( J ) + D*T
               END DO
               DO J = I+1,N
                  T = ABS( A( J, I ) )
                  U = U + S( J )*T
                  WORK( J ) = WORK( J ) + D*T
               END DO
            }

            AVG = AVG + ( U + WORK( I ) ) * D / N
            S( I ) = SI
         END DO
      END DO

 999  CONTINUE

      SMLNUM = DLAMCH( 'SAFEMIN' )
      BIGNUM = ONE / SMLNUM
      SMIN = BIGNUM
      SMAX = ZERO
      T = ONE / SQRT( AVG )
      BASE = DLAMCH( 'B' )
      U = ONE / LOG( BASE )
      for (I = 1; I <= N; I++) {
         S( I ) = BASE ** INT( U * LOG( S( I ) * T ) )
         SMIN = MIN( SMIN, S( I ) )
         SMAX = MAX( SMAX, S( I ) )
      END DO
      SCOND = MAX( SMIN, SMLNUM ) / MIN( SMAX, BIGNUM )

      }
