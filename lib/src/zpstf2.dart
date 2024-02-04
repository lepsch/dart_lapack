      void zpstf2(UPLO, N, A, LDA, PIV, RANK, TOL, WORK, INFO ) {

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double             TOL;
      int                INFO, LDA, N, RANK;
      String             UPLO;
      // ..
      // .. Array Arguments ..
      Complex         A( LDA, * );
      double             WORK( 2*N );
      int                PIV( N );
      // ..

// =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      Complex         CONE;
      const              CONE = ( 1.0, 0.0 ) ;
      // ..
      // .. Local Scalars ..
      Complex         ZTEMP;
      double             AJJ, DSTOP, DTEMP;
      int                I, ITEMP, J, PVT;
      bool               UPPER;
      // ..
      // .. External Functions ..
      //- double             DLAMCH;
      //- bool               LSAME, DISNAN;
      // EXTERNAL DLAMCH, LSAME, DISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZDSCAL, ZGEMV, ZLACGV, ZSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCONJG, MAX, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      INFO = 0;
      UPPER = LSAME( UPLO, 'U' );
      if ( !UPPER && !LSAME( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -4;
      }
      if ( INFO != 0 ) {
         xerbla('ZPSTF2', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Initialize PIV

      for (I = 1; I <= N; I++) { // 100
         PIV[I] = I;
      } // 100

      // Compute stopping value

      for (I = 1; I <= N; I++) { // 110
         WORK[I] = (A( I, I )).toDouble();
      } // 110
      PVT = MAXLOC( WORK( 1:N ), 1 );
      AJJ = (A( PVT, PVT )).toDouble();
      if ( AJJ <= ZERO || DISNAN( AJJ ) ) {
         RANK = 0;
         INFO = 1;
         GO TO 200;
      }

      // Compute stopping value if not supplied

      if ( TOL < ZERO ) {
         DSTOP = N * DLAMCH( 'Epsilon' ) * AJJ;
      } else {
         DSTOP = TOL;
      }

      // Set first half of WORK to zero, holds dot products

      for (I = 1; I <= N; I++) { // 120
         WORK[I] = 0;
      } // 120

      if ( UPPER ) {

         // Compute the Cholesky factorization P**T * A * P = U**H* U

         for (J = 1; J <= N; J++) { // 150

         // Find pivot, test for exit, else swap rows and columns
         // Update dot products, compute possible pivots which are
         // stored in the second half of WORK

            for (I = J; I <= N; I++) { // 130

               if ( J > 1 ) {
                  WORK[I] = WORK( I ) + DBLE( DCONJG( A( J-1, I ) )* A( J-1, I ) );
               }
               WORK[N+I] = (A( I, I )).toDouble() - WORK( I );

            } // 130

            if ( J > 1 ) {
               ITEMP = MAXLOC( WORK( (N+J):(2*N) ), 1 );
               PVT = ITEMP + J - 1;
               AJJ = WORK( N+PVT );
               if ( AJJ <= DSTOP || DISNAN( AJJ ) ) {
                  A[J, J] = AJJ;
                  GO TO 190;
               }
            }

            if ( J != PVT ) {

               // Pivot OK, so can now swap pivot rows and columns

               A[PVT, PVT] = A( J, J );
               zswap(J-1, A( 1, J ), 1, A( 1, PVT ), 1 );
               if (PVT < N) zswap( N-PVT, A( J, PVT+1 ), LDA, A( PVT, PVT+1 ), LDA );
               for (I = J + 1; I <= PVT - 1; I++) { // 140
                  ZTEMP = DCONJG( A( J, I ) );
                  A[J, I] = DCONJG( A( I, PVT ) );
                  A[I, PVT] = ZTEMP;
               } // 140
               A[J, PVT] = DCONJG( A( J, PVT ) );

               // Swap dot products and PIV

               DTEMP = WORK( J );
               WORK[J] = WORK( PVT );
               WORK[PVT] = DTEMP;
               ITEMP = PIV( PVT );
               PIV[PVT] = PIV( J );
               PIV[J] = ITEMP;
            }

            AJJ = sqrt( AJJ );
            A[J, J] = AJJ;

            // Compute elements J+1:N of row J

            if ( J < N ) {
               zlacgv(J-1, A( 1, J ), 1 );
               zgemv('Trans', J-1, N-J, -CONE, A( 1, J+1 ), LDA, A( 1, J ), 1, CONE, A( J, J+1 ), LDA );
               zlacgv(J-1, A( 1, J ), 1 );
               zdscal(N-J, ONE / AJJ, A( J, J+1 ), LDA );
            }

         } // 150

      } else {

         // Compute the Cholesky factorization P**T * A * P = L * L**H

         for (J = 1; J <= N; J++) { // 180

         // Find pivot, test for exit, else swap rows and columns
         // Update dot products, compute possible pivots which are
         // stored in the second half of WORK

            for (I = J; I <= N; I++) { // 160

               if ( J > 1 ) {
                  WORK[I] = WORK( I ) + DBLE( DCONJG( A( I, J-1 ) )* A( I, J-1 ) );
               }
               WORK[N+I] = (A( I, I )).toDouble() - WORK( I );

            } // 160

            if ( J > 1 ) {
               ITEMP = MAXLOC( WORK( (N+J):(2*N) ), 1 );
               PVT = ITEMP + J - 1;
               AJJ = WORK( N+PVT );
               if ( AJJ <= DSTOP || DISNAN( AJJ ) ) {
                  A[J, J] = AJJ;
                  GO TO 190;
               }
            }

            if ( J != PVT ) {

               // Pivot OK, so can now swap pivot rows and columns

               A[PVT, PVT] = A( J, J );
               zswap(J-1, A( J, 1 ), LDA, A( PVT, 1 ), LDA );
               if (PVT < N) zswap( N-PVT, A( PVT+1, J ), 1, A( PVT+1, PVT ), 1 );
               for (I = J + 1; I <= PVT - 1; I++) { // 170
                  ZTEMP = DCONJG( A( I, J ) );
                  A[I, J] = DCONJG( A( PVT, I ) );
                  A[PVT, I] = ZTEMP;
               } // 170
               A[PVT, J] = DCONJG( A( PVT, J ) );

               // Swap dot products and PIV

               DTEMP = WORK( J );
               WORK[J] = WORK( PVT );
               WORK[PVT] = DTEMP;
               ITEMP = PIV( PVT );
               PIV[PVT] = PIV( J );
               PIV[J] = ITEMP;
            }

            AJJ = sqrt( AJJ );
            A[J, J] = AJJ;

            // Compute elements J+1:N of column J

            if ( J < N ) {
               zlacgv(J-1, A( J, 1 ), LDA );
               zgemv('No Trans', N-J, J-1, -CONE, A( J+1, 1 ), LDA, A( J, 1 ), LDA, CONE, A( J+1, J ), 1 );
               zlacgv(J-1, A( J, 1 ), LDA );
               zdscal(N-J, ONE / AJJ, A( J+1, J ), 1 );
            }

         } // 180

      }

      // Ran to completion, A has full rank

      RANK = N;

      GO TO 200;
      } // 190

      // Rank is number of steps completed.  Set INFO = 1 to signal
      // that the factorization cannot be used to solve a system.

      RANK = J - 1;
      INFO = 1;

      } // 200
      return;
      }
