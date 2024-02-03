      SUBROUTINE SLARZ( SIDE, M, N, L, V, INCV, TAU, C, LDC, WORK )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             SIDE;
      int                INCV, L, LDC, M, N;
      REAL               TAU
      // ..
      // .. Array Arguments ..
      REAL               C( LDC, * ), V( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SAXPY, SCOPY, SGEMV, SGER
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. Executable Statements ..

      if ( LSAME( SIDE, 'L' ) ) {

         // Form  H * C

         if ( TAU.NE.ZERO ) {

            // w( 1:n ) = C( 1, 1:n )

            CALL SCOPY( N, C, LDC, WORK, 1 )

            // w( 1:n ) = w( 1:n ) + C( m-l+1:m, 1:n )**T * v( 1:l )

            CALL SGEMV( 'Transpose', L, N, ONE, C( M-L+1, 1 ), LDC, V, INCV, ONE, WORK, 1 )

            // C( 1, 1:n ) = C( 1, 1:n ) - tau * w( 1:n )

            CALL SAXPY( N, -TAU, WORK, 1, C, LDC )

            // C( m-l+1:m, 1:n ) = C( m-l+1:m, 1:n ) - ...
                               t // au * v( 1:l ) * w( 1:n )**T

            CALL SGER( L, N, -TAU, V, INCV, WORK, 1, C( M-L+1, 1 ), LDC )
         }

      } else {

         // Form  C * H

         if ( TAU.NE.ZERO ) {

            // w( 1:m ) = C( 1:m, 1 )

            CALL SCOPY( M, C, 1, WORK, 1 )

            // w( 1:m ) = w( 1:m ) + C( 1:m, n-l+1:n, 1:n ) * v( 1:l )

            CALL SGEMV( 'No transpose', M, L, ONE, C( 1, N-L+1 ), LDC, V, INCV, ONE, WORK, 1 )

            // C( 1:m, 1 ) = C( 1:m, 1 ) - tau * w( 1:m )

            CALL SAXPY( M, -TAU, WORK, 1, C, 1 )

            // C( 1:m, n-l+1:n ) = C( 1:m, n-l+1:n ) - ...
                               t // au * w( 1:m ) * v( 1:l )**T

            CALL SGER( M, L, -TAU, WORK, 1, V, INCV, C( 1, N-L+1 ), LDC )

         }

      }

      RETURN

      // End of SLARZ

      }
