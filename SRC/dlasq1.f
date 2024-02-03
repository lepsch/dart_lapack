      SUBROUTINE DLASQ1( N, D, E, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, N;
      // ..
      // .. Array Arguments ..
      double             D( * ), E( * ), WORK( * );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ZERO;
      const              ZERO = 0.0D0 ;
      // ..
      // .. Local Scalars ..
      int                I, IINFO;
      double             EPS, SCALE, SAFMIN, SIGMN, SIGMX;
      // ..
      // .. External Subroutines ..
      // EXTERNAL DCOPY, DLAS2, DLASCL, DLASQ2, DLASRT, XERBLA
      // ..
      // .. External Functions ..
      double             DLAMCH;
      // EXTERNAL DLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SQRT
      // ..
      // .. Executable Statements ..

      INFO = 0
      if ( N.LT.0 ) {
         INFO = -1
         xerbla('DLASQ1', -INFO );
         RETURN
      } else if ( N.EQ.0 ) {
         RETURN
      } else if ( N.EQ.1 ) {
         D( 1 ) = ABS( D( 1 ) )
         RETURN
      } else if ( N.EQ.2 ) {
         dlas2(D( 1 ), E( 1 ), D( 2 ), SIGMN, SIGMX );
         D( 1 ) = SIGMX
         D( 2 ) = SIGMN
         RETURN
      }

      // Estimate the largest singular value.

      SIGMX = ZERO
      for (I = 1; I <= N - 1; I++) { // 10
         D( I ) = ABS( D( I ) )
         SIGMX = MAX( SIGMX, ABS( E( I ) ) )
      } // 10
      D( N ) = ABS( D( N ) )

      // Early return if SIGMX is zero (matrix is already diagonal).

      if ( SIGMX.EQ.ZERO ) {
         dlasrt('D', N, D, IINFO );
         RETURN
      }

      for (I = 1; I <= N; I++) { // 20
         SIGMX = MAX( SIGMX, D( I ) )
      } // 20

      // Copy D and E into WORK (in the Z format) and scale (squaring the
      // input data makes scaling by a power of the radix pointless).

      EPS = DLAMCH( 'Precision' )
      SAFMIN = DLAMCH( 'Safe minimum' )
      SCALE = SQRT( EPS / SAFMIN )
      dcopy(N, D, 1, WORK( 1 ), 2 );
      dcopy(N-1, E, 1, WORK( 2 ), 2 );
      dlascl('G', 0, 0, SIGMX, SCALE, 2*N-1, 1, WORK, 2*N-1, IINFO );

      // Compute the q's and e's.

      for (I = 1; I <= 2*N - 1; I++) { // 30
         WORK( I ) = WORK( I )**2
      } // 30
      WORK( 2*N ) = ZERO

      dlasq2(N, WORK, INFO );

      if ( INFO.EQ.0 ) {
         for (I = 1; I <= N; I++) { // 40
            D( I ) = SQRT( WORK( I ) )
         } // 40
         dlascl('G', 0, 0, SCALE, SIGMX, N, 1, D, N, IINFO );
      } else if ( INFO.EQ.2 ) {

      // Maximum number of iterations exceeded.  Move data from WORK
      // into D and E so the calling subroutine can try to finish

         for (I = 1; I <= N; I++) {
            D( I ) = SQRT( WORK( 2*I-1 ) )
            E( I ) = SQRT( WORK( 2*I ) )
         }
         dlascl('G', 0, 0, SCALE, SIGMX, N, 1, D, N, IINFO );
         dlascl('G', 0, 0, SCALE, SIGMX, N, 1, E, N, IINFO );
      }

      RETURN

      // End of DLASQ1

      }
