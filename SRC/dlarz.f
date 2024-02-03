      SUBROUTINE DLARZ( SIDE, M, N, L, V, INCV, TAU, C, LDC, WORK )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             SIDE;
      int                INCV, L, LDC, M, N;
      double             TAU;
      // ..
      // .. Array Arguments ..
      double             C( LDC, * ), V( * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DAXPY, DCOPY, DGEMV, DGER
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. Executable Statements ..

      IF( LSAME( SIDE, 'L' ) ) THEN

         // Form  H * C

         IF( TAU.NE.ZERO ) THEN

            // w( 1:n ) = C( 1, 1:n )

            CALL DCOPY( N, C, LDC, WORK, 1 )

            // w( 1:n ) = w( 1:n ) + C( m-l+1:m, 1:n )**T * v( 1:l )

            CALL DGEMV( 'Transpose', L, N, ONE, C( M-L+1, 1 ), LDC, V, INCV, ONE, WORK, 1 )

            // C( 1, 1:n ) = C( 1, 1:n ) - tau * w( 1:n )

            CALL DAXPY( N, -TAU, WORK, 1, C, LDC )

            // C( m-l+1:m, 1:n ) = C( m-l+1:m, 1:n ) - ...
                               t // au * v( 1:l ) * w( 1:n )**T

            CALL DGER( L, N, -TAU, V, INCV, WORK, 1, C( M-L+1, 1 ), LDC )
         END IF

      ELSE

         // Form  C * H

         IF( TAU.NE.ZERO ) THEN

            // w( 1:m ) = C( 1:m, 1 )

            CALL DCOPY( M, C, 1, WORK, 1 )

            // w( 1:m ) = w( 1:m ) + C( 1:m, n-l+1:n, 1:n ) * v( 1:l )

            CALL DGEMV( 'No transpose', M, L, ONE, C( 1, N-L+1 ), LDC, V, INCV, ONE, WORK, 1 )

            // C( 1:m, 1 ) = C( 1:m, 1 ) - tau * w( 1:m )

            CALL DAXPY( M, -TAU, WORK, 1, C, 1 )

            // C( 1:m, n-l+1:n ) = C( 1:m, n-l+1:n ) - ...
                               t // au * w( 1:m ) * v( 1:l )**T

            CALL DGER( M, L, -TAU, WORK, 1, V, INCV, C( 1, N-L+1 ), LDC )

         END IF

      END IF

      RETURN

      // End of DLARZ

      }
