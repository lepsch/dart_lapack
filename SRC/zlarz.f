      SUBROUTINE ZLARZ( SIDE, M, N, L, V, INCV, TAU, C, LDC, WORK )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      String             SIDE;
      int                INCV, L, LDC, M, N;
      COMPLEX*16         TAU
      // ..
      // .. Array Arguments ..
      COMPLEX*16         C( LDC, * ), V( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         ONE, ZERO
      PARAMETER          ( ONE = ( 1.0D+0, 0.0D+0 ), ZERO = ( 0.0D+0, 0.0D+0 ) )
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZAXPY, ZCOPY, ZGEMV, ZGERC, ZGERU, ZLACGV
      // ..
      // .. External Functions ..
      bool               LSAME;
      // EXTERNAL LSAME
      // ..
      // .. Executable Statements ..

      IF( LSAME( SIDE, 'L' ) ) THEN

         // Form  H * C

         IF( TAU.NE.ZERO ) THEN

            // w( 1:n ) = conjg( C( 1, 1:n ) )

            CALL ZCOPY( N, C, LDC, WORK, 1 )
            CALL ZLACGV( N, WORK, 1 )

            // w( 1:n ) = conjg( w( 1:n ) + C( m-l+1:m, 1:n )**H * v( 1:l ) )

            CALL ZGEMV( 'Conjugate transpose', L, N, ONE, C( M-L+1, 1 ), LDC, V, INCV, ONE, WORK, 1 )
            CALL ZLACGV( N, WORK, 1 )

            // C( 1, 1:n ) = C( 1, 1:n ) - tau * w( 1:n )

            CALL ZAXPY( N, -TAU, WORK, 1, C, LDC )

            // C( m-l+1:m, 1:n ) = C( m-l+1:m, 1:n ) - ...
                               t // au * v( 1:l ) * w( 1:n )**H

            CALL ZGERU( L, N, -TAU, V, INCV, WORK, 1, C( M-L+1, 1 ), LDC )
         END IF

      ELSE

         // Form  C * H

         IF( TAU.NE.ZERO ) THEN

            // w( 1:m ) = C( 1:m, 1 )

            CALL ZCOPY( M, C, 1, WORK, 1 )

            // w( 1:m ) = w( 1:m ) + C( 1:m, n-l+1:n, 1:n ) * v( 1:l )

            CALL ZGEMV( 'No transpose', M, L, ONE, C( 1, N-L+1 ), LDC, V, INCV, ONE, WORK, 1 )

            // C( 1:m, 1 ) = C( 1:m, 1 ) - tau * w( 1:m )

            CALL ZAXPY( M, -TAU, WORK, 1, C, 1 )

            // C( 1:m, n-l+1:n ) = C( 1:m, n-l+1:n ) - ...
                               t // au * w( 1:m ) * v( 1:l )**H

            CALL ZGERC( M, L, -TAU, WORK, 1, V, INCV, C( 1, N-L+1 ), LDC )

         END IF

      END IF

      RETURN

      // End of ZLARZ

      }
