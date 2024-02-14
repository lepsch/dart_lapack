      void spstf2(final int UPLO, final int N, final Matrix<double> A_, final int LDA, final int PIV, final int RANK, final int TOL, final Array<double> _WORK_, final Box<int> INFO,) {
  final A = A_.dim();
  final _WORK = _WORK_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      double               TOL;
      int                INFO, LDA, N, RANK;
      String             UPLO;
      double               A( LDA, * ), WORK( 2*N );
      int                PIV( N );
      // ..

      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      double               AJJ, SSTOP, STEMP;
      int                I, ITEMP, J, PVT;
      bool               UPPER;
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH;
      //- bool               lsame, SISNAN;
      // EXTERNAL SLAMCH, lsame, SISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEMV, SSCAL, SSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, SQRT, MAXLOC

      // Test the input parameters

      INFO = 0;
      UPPER = lsame( UPLO, 'U' );
      if ( !UPPER && !lsame( UPLO, 'L' ) ) {
         INFO = -1;
      } else if ( N < 0 ) {
         INFO = -2;
      } else if ( LDA < max( 1, N ) ) {
         INFO = -4;
      }
      if ( INFO != 0 ) {
         xerbla('SPSTF2', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Initialize PIV

      for (I = 1; I <= N; I++) { // 100
         PIV[I] = I;
      } // 100

      // Compute stopping value

      PVT = 1;
      AJJ = A( PVT, PVT );
      for (I = 2; I <= N; I++) {
         if ( A( I, I ) > AJJ ) {
            PVT = I;
            AJJ = A( PVT, PVT );
         }
      }
      if ( AJJ <= ZERO || SISNAN( AJJ ) ) {
         RANK = 0;
         INFO = 1;
         GO TO 170;
      }

      // Compute stopping value if not supplied

      if ( TOL < ZERO ) {
         SSTOP = N * SLAMCH( 'Epsilon' ) * AJJ;
      } else {
         SSTOP = TOL;
      }

      // Set first half of WORK to zero, holds dot products

      for (I = 1; I <= N; I++) { // 110
         WORK[I] = 0;
      } // 110

      if ( UPPER ) {

         // Compute the Cholesky factorization P**T * A * P = U**T * U

         for (J = 1; J <= N; J++) { // 130

         // Find pivot, test for exit, else swap rows and columns
         // Update dot products, compute possible pivots which are
         // stored in the second half of WORK

            for (I = J; I <= N; I++) { // 120

               if ( J > 1 ) {
                  WORK[I] = WORK( I ) + A( J-1, I )**2;
               }
               WORK[N+I] = A( I, I ) - WORK( I );

            } // 120

            if ( J > 1 ) {
               ITEMP = MAXLOC( WORK( (N+J):(2*N) ), 1 );
               PVT = ITEMP + J - 1;
               AJJ = WORK( N+PVT );
               if ( AJJ <= SSTOP || SISNAN( AJJ ) ) {
                  A[J][J] = AJJ;
                  GO TO 160;
               }
            }

            if ( J != PVT ) {

               // Pivot OK, so can now swap pivot rows and columns

               A[PVT][PVT] = A( J, J );
               sswap(J-1, A( 1, J ), 1, A( 1, PVT ), 1 );
               if (PVT < N) sswap( N-PVT, A( J, PVT+1 ), LDA, A( PVT, PVT+1 ), LDA );
               sswap(PVT-J-1, A( J, J+1 ), LDA, A( J+1, PVT ), 1 );

               // Swap dot products and PIV

               STEMP = WORK( J );
               WORK[J] = WORK( PVT );
               WORK[PVT] = STEMP;
               ITEMP = PIV( PVT );
               PIV[PVT] = PIV( J );
               PIV[J] = ITEMP;
            }

            AJJ = sqrt( AJJ );
            A[J][J] = AJJ;

            // Compute elements J+1:N of row J

            if ( J < N ) {
               sgemv('Trans', J-1, N-J, -ONE, A( 1, J+1 ), LDA, A( 1, J ), 1, ONE, A( J, J+1 ), LDA );
               sscal(N-J, ONE / AJJ, A( J, J+1 ), LDA );
            }

         } // 130

      } else {

         // Compute the Cholesky factorization P**T * A * P = L * L**T

         for (J = 1; J <= N; J++) { // 150

         // Find pivot, test for exit, else swap rows and columns
         // Update dot products, compute possible pivots which are
         // stored in the second half of WORK

            for (I = J; I <= N; I++) { // 140

               if ( J > 1 ) {
                  WORK[I] = WORK( I ) + A( I, J-1 )**2;
               }
               WORK[N+I] = A( I, I ) - WORK( I );

            } // 140

            if ( J > 1 ) {
               ITEMP = MAXLOC( WORK( (N+J):(2*N) ), 1 );
               PVT = ITEMP + J - 1;
               AJJ = WORK( N+PVT );
               if ( AJJ <= SSTOP || SISNAN( AJJ ) ) {
                  A[J][J] = AJJ;
                  GO TO 160;
               }
            }

            if ( J != PVT ) {

               // Pivot OK, so can now swap pivot rows and columns

               A[PVT][PVT] = A( J, J );
               sswap(J-1, A( J, 1 ), LDA, A( PVT, 1 ), LDA );
               if (PVT < N) sswap( N-PVT, A( PVT+1, J ), 1, A( PVT+1, PVT ), 1 );
               sswap(PVT-J-1, A( J+1, J ), 1, A( PVT, J+1 ), LDA );

               // Swap dot products and PIV

               STEMP = WORK( J );
               WORK[J] = WORK( PVT );
               WORK[PVT] = STEMP;
               ITEMP = PIV( PVT );
               PIV[PVT] = PIV( J );
               PIV[J] = ITEMP;
            }

            AJJ = sqrt( AJJ );
            A[J][J] = AJJ;

            // Compute elements J+1:N of column J

            if ( J < N ) {
               sgemv('No Trans', N-J, J-1, -ONE, A( J+1, 1 ), LDA, A( J, 1 ), LDA, ONE, A( J+1, J ), 1 );
               sscal(N-J, ONE / AJJ, A( J+1, J ), 1 );
            }

         } // 150

      }

      // Ran to completion, A has full rank

      RANK = N;

      GO TO 170;
      } // 160

      // Rank is number of steps completed.  Set INFO = 1 to signal
      // that the factorization cannot be used to solve a system.

      RANK = J - 1;
      INFO = 1;

      } // 170
      }
