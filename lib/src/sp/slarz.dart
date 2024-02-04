      void slarz(SIDE, M, N, L, V, INCV, TAU, C, LDC, WORK ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             SIDE;
      int                INCV, L, LDC, M, N;
      REAL               TAU;
      // ..
      // .. Array Arguments ..
      REAL               C( LDC, * ), V( * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SAXPY, SCOPY, SGEMV, SGER
      // ..
      // .. External Functions ..
      //- bool               lsame;
      // EXTERNAL lsame
      // ..
      // .. Executable Statements ..

      if ( lsame( SIDE, 'L' ) ) {

         // Form  H * C

         if ( TAU != ZERO ) {

            // w( 1:n ) = C( 1, 1:n )

            scopy(N, C, LDC, WORK, 1 );

            // w( 1:n ) = w( 1:n ) + C( m-l+1:m, 1:n )**T * v( 1:l )

            sgemv('Transpose', L, N, ONE, C( M-L+1, 1 ), LDC, V, INCV, ONE, WORK, 1 );

            // C( 1, 1:n ) = C( 1, 1:n ) - tau * w( 1:n )

            saxpy(N, -TAU, WORK, 1, C, LDC );

            // C( m-l+1:m, 1:n ) = C( m-l+1:m, 1:n ) - ...
                                // tau * v( 1:l ) * w( 1:n )**T

            sger(L, N, -TAU, V, INCV, WORK, 1, C( M-L+1, 1 ), LDC );
         }

      } else {

         // Form  C * H

         if ( TAU != ZERO ) {

            // w( 1:m ) = C( 1:m, 1 )

            scopy(M, C, 1, WORK, 1 );

            // w( 1:m ) = w( 1:m ) + C( 1:m, n-l+1:n, 1:n ) * v( 1:l )

            sgemv('No transpose', M, L, ONE, C( 1, N-L+1 ), LDC, V, INCV, ONE, WORK, 1 );

            // C( 1:m, 1 ) = C( 1:m, 1 ) - tau * w( 1:m )

            saxpy(M, -TAU, WORK, 1, C, 1 );

            // C( 1:m, n-l+1:n ) = C( 1:m, n-l+1:n ) - ...
                                // tau * w( 1:m ) * v( 1:l )**T

            sger(M, L, -TAU, WORK, 1, V, INCV, C( 1, N-L+1 ), LDC );

         }

      }

      return;
      }
