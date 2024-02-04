      void csyequb(UPLO, N, A, LDA, S, SCOND, AMAX, WORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, N;
      REAL               AMAX, SCOND;
      String             UPLO;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * ), WORK( * );
      REAL               S( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      int                MAX_ITER;
      const              MAX_ITER = 100 ;
      // ..
      // .. Local Scalars ..
      int                I, J, ITER;
      REAL               AVG, STD, TOL, C0, C1, C2, T, U, SI, D, BASE, SMIN, SMAX, SMLNUM, BIGNUM, SCALE, SUMSQ;
      bool               UP;
      COMPLEX            ZDUM;
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH;
      //- bool               LSAME;
      // EXTERNAL LSAME, SLAMCH
      // ..
      // .. External Subroutines ..
      // EXTERNAL CLASSQ, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, AIMAG, INT, LOG, MAX, MIN, REAL, SQRT
      // ..
      // .. Statement Functions ..
      REAL               CABS1;
      // ..
      // .. Statement Function Definitions ..
      CABS1[ZDUM] = ( REAL( ZDUM ) ).abs() + ( AIMAG( ZDUM ) ).abs();
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0;
      if ( !( LSAME( UPLO, 'U' ) || LSAME( UPLO, 'L' ) ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -4;
      }
      if ( INFO != 0 ) {
         xerbla('CSYEQUB', -INFO );
         return;
      }

      UP = LSAME( UPLO, 'U' );
      AMAX = ZERO;

      // Quick return if possible.

      if ( N == 0 ) {
         SCOND = ONE;
         return;
      }

      for (I = 1; I <= N; I++) {
         S[I] = ZERO;
      }

      AMAX = ZERO;
      if ( UP ) {
         for (J = 1; J <= N; J++) {
            for (I = 1; I <= J-1; I++) {
               S[I] = max( S( I ), CABS1( A( I, J ) ) );
               S[J] = max( S( J ), CABS1( A( I, J ) ) );
               AMAX = max( AMAX, CABS1( A( I, J ) ) );
            }
            S[J] = max( S( J ), CABS1( A( J, J ) ) );
            AMAX = max( AMAX, CABS1( A( J, J ) ) );
         }
      } else {
         for (J = 1; J <= N; J++) {
            S[J] = max( S( J ), CABS1( A( J, J ) ) );
            AMAX = max( AMAX, CABS1( A( J, J ) ) );
            for (I = J+1; I <= N; I++) {
               S[I] = max( S( I ), CABS1( A( I, J ) ) );
               S[J] = max( S( J ), CABS1( A( I, J ) ) );
               AMAX = max( AMAX, CABS1( A( I, J ) ) );
            }
         }
      }
      for (J = 1; J <= N; J++) {
         S[J] = 1.0 / S( J );
      }

      TOL = ONE / sqrt( 2.0 * N );

      for (ITER = 1; ITER <= MAX_ITER; ITER++) {
         SCALE = 0.0;
         SUMSQ = 0.0;
         // beta = |A|s
         for (I = 1; I <= N; I++) {
             WORK[I] = ZERO;
         }
         if ( UP ) {
            for (J = 1; J <= N; J++) {
               for (I = 1; I <= J-1; I++) {
                  WORK[I] = WORK( I ) + CABS1( A( I, J ) ) * S( J );
                  WORK[J] = WORK( J ) + CABS1( A( I, J ) ) * S( I );
               }
               WORK[J] = WORK( J ) + CABS1( A( J, J ) ) * S( J );
            }
         } else {
            for (J = 1; J <= N; J++) {
               WORK[J] = WORK( J ) + CABS1( A( J, J ) ) * S( J );
               for (I = J+1; I <= N; I++) {
                  WORK[I] = WORK( I ) + CABS1( A( I, J ) ) * S( J );
                  WORK[J] = WORK( J ) + CABS1( A( I, J ) ) * S( I );
               }
            }
         }

         // avg = s^T beta / n
         AVG = 0.0;
         for (I = 1; I <= N; I++) {
            AVG = AVG + REAL( S( I )*WORK( I ) );
         }
         AVG = AVG / N;

         STD = 0.0;
         for (I = N+1; I <= 2*N; I++) {
            WORK[I] = S( I-N ) * WORK( I-N ) - AVG;
         }
         classq(N, WORK( N+1 ), 1, SCALE, SUMSQ );
         STD = SCALE * sqrt( SUMSQ / N );

         if (STD < TOL * AVG) GOTO 999;

         for (I = 1; I <= N; I++) {
            T = CABS1( A( I, I ) );
            SI = S( I );
            C2 = ( N-1 ) * T;
            C1 = REAL( N-2 ) * ( REAL( WORK( I ) ) - T*SI );
            C0 = -(T*SI)*SI + 2 * REAL( WORK( I ) ) * SI - N*AVG;
            D = C1*C1 - 4*C0*C2;

            if ( D <= 0 ) {
               INFO = -1;
               return;
            }
            SI = -2*C0 / ( C1 + sqrt( D ) );

            D = SI - S( I );
            U = ZERO;
            if ( UP ) {
               for (J = 1; J <= I; J++) {
                  T = CABS1( A( J, I ) );
                  U = U + S( J )*T;
                  WORK[J] = WORK( J ) + D*T;
               }
               for (J = I+1; J <= N; J++) {
                  T = CABS1( A( I, J ) );
                  U = U + S( J )*T;
                  WORK[J] = WORK( J ) + D*T;
               }
            } else {
               for (J = 1; J <= I; J++) {
                  T = CABS1( A( I, J ) );
                  U = U + S( J )*T;
                  WORK[J] = WORK( J ) + D*T;
               }
               for (J = I+1; J <= N; J++) {
                  T = CABS1( A( J, I ) );
                  U = U + S( J )*T;
                  WORK[J] = WORK( J ) + D*T;
               }
            }

            AVG = AVG + ( U + REAL( WORK( I ) ) ) * D / N;
            S[I] = SI;
         }
      }

      } // 999

      SMLNUM = SLAMCH( 'SAFEMIN' );
      BIGNUM = ONE / SMLNUM;
      SMIN = BIGNUM;
      SMAX = ZERO;
      T = ONE / sqrt( AVG );
      BASE = SLAMCH( 'B' );
      U = ONE / LOG( BASE );
      for (I = 1; I <= N; I++) {
         S[I] = BASE ** INT( U * LOG( S( I ) * T ) );
         SMIN = min( SMIN, S( I ) );
         SMAX = max( SMAX, S( I ) );
      }
      SCOND = max( SMIN, SMLNUM ) / min( SMAX, BIGNUM );

      }
