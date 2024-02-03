      SUBROUTINE CLARZ( SIDE, M, N, L, V, INCV, TAU, C, LDC, WORK );

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             SIDE;
      int                INCV, L, LDC, M, N;
      COMPLEX            TAU;
      // ..
      // .. Array Arguments ..
      COMPLEX            C( LDC, * ), V( * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      COMPLEX            ONE, ZERO;
      const              ONE = ( 1.0, 0.0 ), ZERO = ( 0.0, 0.0 ) ;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CAXPY, CCOPY, CGEMV, CGERC, CGERU, CLACGV
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. Executable Statements ..

      if ( LSAME( SIDE, 'L' ) ) {

         // Form  H * C

         if ( TAU != ZERO ) {

            // w( 1:n ) = conjg( C( 1, 1:n ) )

            ccopy(N, C, LDC, WORK, 1 );
            clacgv(N, WORK, 1 );

            // w( 1:n ) = conjg( w( 1:n ) + C( m-l+1:m, 1:n )**H * v( 1:l ) )

            cgemv('Conjugate transpose', L, N, ONE, C( M-L+1, 1 ), LDC, V, INCV, ONE, WORK, 1 );
            clacgv(N, WORK, 1 );

            // C( 1, 1:n ) = C( 1, 1:n ) - tau * w( 1:n )

            caxpy(N, -TAU, WORK, 1, C, LDC );

            // C( m-l+1:m, 1:n ) = C( m-l+1:m, 1:n ) - ...
                                // tau * v( 1:l ) * w( 1:n )**H

            cgeru(L, N, -TAU, V, INCV, WORK, 1, C( M-L+1, 1 ), LDC );
         }

      } else {

         // Form  C * H

         if ( TAU != ZERO ) {

            // w( 1:m ) = C( 1:m, 1 )

            ccopy(M, C, 1, WORK, 1 );

            // w( 1:m ) = w( 1:m ) + C( 1:m, n-l+1:n, 1:n ) * v( 1:l )

            cgemv('No transpose', M, L, ONE, C( 1, N-L+1 ), LDC, V, INCV, ONE, WORK, 1 );

            // C( 1:m, 1 ) = C( 1:m, 1 ) - tau * w( 1:m )

            caxpy(M, -TAU, WORK, 1, C, 1 );

            // C( 1:m, n-l+1:n ) = C( 1:m, n-l+1:n ) - ...
                                // tau * w( 1:m ) * v( 1:l )**H

            cgerc(M, L, -TAU, WORK, 1, V, INCV, C( 1, N-L+1 ), LDC );

         }

      }

      return;

      // End of CLARZ

      }
