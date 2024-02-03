      void slasq1(N, D, E, WORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, N;
      // ..
      // .. Array Arguments ..
      REAL               D( * ), E( * ), WORK( * );
      // ..

// =====================================================================

      // .. Parameters ..
      REAL               ZERO;
      const              ZERO = 0.0 ;
      // ..
      // .. Local Scalars ..
      int                I, IINFO;
      REAL               EPS, SCALE, SAFMIN, SIGMN, SIGMX;
      // ..
      // .. External Subroutines ..
      // EXTERNAL SCOPY, SLAS2, SLASCL, SLASQ2, SLASRT, XERBLA
      // ..
      // .. External Functions ..
      REAL               SLAMCH;
      // EXTERNAL SLAMCH
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, SQRT
      // ..
      // .. Executable Statements ..

      INFO = 0;
      if ( N < 0 ) {
         INFO = -1;
         xerbla('SLASQ1', -INFO );
         return;
      } else if ( N == 0 ) {
         return;
      } else if ( N == 1 ) {
         D( 1 ) = ( D( 1 ) ).abs();
         return;
      } else if ( N == 2 ) {
         slas2(D( 1 ), E( 1 ), D( 2 ), SIGMN, SIGMX );
         D( 1 ) = SIGMX;
         D( 2 ) = SIGMN;
         return;
      }

      // Estimate the largest singular value.

      SIGMX = ZERO;
      for (I = 1; I <= N - 1; I++) { // 10
         D( I ) = ( D( I ) ).abs();
         SIGMX = max( SIGMX, ( E( I ) ) ).abs();
      } // 10
      D( N ) = ( D( N ) ).abs();

      // Early return if SIGMX is zero (matrix is already diagonal).

      if ( SIGMX == ZERO ) {
         slasrt('D', N, D, IINFO );
         return;
      }

      for (I = 1; I <= N; I++) { // 20
         SIGMX = max( SIGMX, D( I ) );
      } // 20

      // Copy D and E into WORK (in the Z format) and scale (squaring the
      // input data makes scaling by a power of the radix pointless).

      EPS = SLAMCH( 'Precision' );
      SAFMIN = SLAMCH( 'Safe minimum' );
      SCALE = sqrt( EPS / SAFMIN );
      scopy(N, D, 1, WORK( 1 ), 2 );
      scopy(N-1, E, 1, WORK( 2 ), 2 );
      slascl('G', 0, 0, SIGMX, SCALE, 2*N-1, 1, WORK, 2*N-1, IINFO );

      // Compute the q's and e's.

      for (I = 1; I <= 2*N - 1; I++) { // 30
         WORK( I ) = WORK( I )**2;
      } // 30
      WORK( 2*N ) = ZERO;

      slasq2(N, WORK, INFO );

      if ( INFO == 0 ) {
         for (I = 1; I <= N; I++) { // 40
            D( I ) = sqrt( WORK( I ) );
         } // 40
         slascl('G', 0, 0, SCALE, SIGMX, N, 1, D, N, IINFO );
      } else if ( INFO == 2 ) {

      // Maximum number of iterations exceeded.  Move data from WORK
      // into D and E so the calling subroutine can try to finish

         for (I = 1; I <= N; I++) {
            D( I ) = sqrt( WORK( 2*I-1 ) );
            E( I ) = sqrt( WORK( 2*I ) );
         }
         slascl('G', 0, 0, SCALE, SIGMX, N, 1, D, N, IINFO );
         slascl('G', 0, 0, SCALE, SIGMX, N, 1, E, N, IINFO );
      }

      return;
      }
