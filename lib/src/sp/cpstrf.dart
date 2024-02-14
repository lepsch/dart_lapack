      void cpstrf(final int UPLO, final int N, final Matrix<double> A_, final int LDA, final int PIV, final int RANK, final int TOL, final Array<double> _WORK_, final Box<int> INFO,) {
  final A = A_.dim();
  final _WORK = _WORK_.dim();

// -- LAPACK computational routine --
// -- LAPACK is a software package provided by Univ. of Tennessee,    --
// -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
      double               TOL;
      int                INFO, LDA, N, RANK;
      String             UPLO;
      Complex            A( LDA, * );
      double               WORK( 2*N );
      int                PIV( N );
      // ..

      double               ONE, ZERO;
      const              ONE = 1.0, ZERO = 0.0 ;
      Complex            CONE;
      const              CONE = ( 1.0, 0.0 ) ;
      Complex            CTEMP;
      double               AJJ, SSTOP, STEMP;
      int                I, ITEMP, J, JB, K, NB, PVT;
      bool               UPPER;
      // ..
      // .. External Functions ..
      //- REAL               SLAMCH;
      //- int                ILAENV;
      //- bool               lsame, SISNAN;
      // EXTERNAL SLAMCH, ILAENV, lsame, SISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMV, CHERK, CLACGV, CPSTF2, CSSCAL, CSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG, MAX, MIN, REAL, SQRT, MAXLOC

      // Test the input parameters.

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
         xerbla('CPSTRF', -INFO );
         return;
      }

      // Quick return if possible

      if (N == 0) return;

      // Get block size

      NB = ilaenv( 1, 'CPOTRF', UPLO, N, -1, -1, -1 );
      if ( NB <= 1 || NB >= N ) {

         // Use unblocked code

         cpstf2(UPLO, N, A( 1, 1 ), LDA, PIV, RANK, TOL, WORK, INFO );
         GO TO 230;

      } else {

      // Initialize PIV

         for (I = 1; I <= N; I++) { // 100
            PIV[I] = I;
         } // 100

      // Compute stopping value

         for (I = 1; I <= N; I++) { // 110
            WORK[I] = double( A( I, I ) );
         } // 110
         PVT = MAXLOC( WORK( 1:N ), 1 );
         AJJ = double( A( PVT, PVT ) );
         if ( AJJ <= ZERO || SISNAN( AJJ ) ) {
            RANK = 0;
            INFO = 1;
            GO TO 230;
         }

      // Compute stopping value if not supplied

         if ( TOL < ZERO ) {
            SSTOP = N * SLAMCH( 'Epsilon' ) * AJJ;
         } else {
            SSTOP = TOL;
         }


         if ( UPPER ) {

            // Compute the Cholesky factorization P**T * A * P = U**H * U

            for (K = 1; NB < 0 ? K >= N : K <= N; K += NB) { // 160

               // Account for last block not being NB wide

               JB = min( NB, N-K+1 );

               // Set relevant part of first half of WORK to zero,
               // holds dot products

               for (I = K; I <= N; I++) { // 120
                  WORK[I] = 0;
               } // 120

               for (J = K; J <= K + JB - 1; J++) { // 150

               // Find pivot, test for exit, else swap rows and columns
               // Update dot products, compute possible pivots which are
               // stored in the second half of WORK

                  for (I = J; I <= N; I++) { // 130

                     if ( J > K ) {
                        WORK[I] = WORK( I ) + double( CONJG( A( J-1, I ) )* A( J-1, I ) );
                     }
                     WORK[N+I] = double( A( I, I ) ) - WORK( I );

                  } // 130

                  if ( J > 1 ) {
                     ITEMP = MAXLOC( WORK( (N+J):(2*N) ), 1 );
                     PVT = ITEMP + J - 1;
                     AJJ = WORK( N+PVT );
                     if ( AJJ <= SSTOP || SISNAN( AJJ ) ) {
                        A[J][J] = AJJ;
                        GO TO 220;
                     }
                  }

                  if ( J != PVT ) {

                     // Pivot OK, so can now swap pivot rows and columns

                     A[PVT][PVT] = A( J, J );
                     cswap(J-1, A( 1, J ), 1, A( 1, PVT ), 1 );
                     if (PVT < N) cswap( N-PVT, A( J, PVT+1 ), LDA, A( PVT, PVT+1 ), LDA );
                     for (I = J + 1; I <= PVT - 1; I++) { // 140
                        CTEMP = CONJG( A( J, I ) );
                        A[J][I] = CONJG( A( I, PVT ) );
                        A[I][PVT] = CTEMP;
                     } // 140
                     A[J][PVT] = CONJG( A( J, PVT ) );

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

                  // Compute elements J+1:N of row J.

                  if ( J < N ) {
                     clacgv(J-1, A( 1, J ), 1 );
                     cgemv('Trans', J-K, N-J, -CONE, A( K, J+1 ), LDA, A( K, J ), 1, CONE, A( J, J+1 ), LDA );
                     clacgv(J-1, A( 1, J ), 1 );
                     csscal(N-J, ONE / AJJ, A( J, J+1 ), LDA );
                  }

               } // 150

               // Update trailing matrix, J already incremented

               if ( K+JB <= N ) {
                  cherk('Upper', 'Conj Trans', N-J+1, JB, -ONE, A( K, J ), LDA, ONE, A( J, J ), LDA );
               }

            } // 160

         } else {

         // Compute the Cholesky factorization P**T * A * P = L * L**H

            for (K = 1; NB < 0 ? K >= N : K <= N; K += NB) { // 210

               // Account for last block not being NB wide

               JB = min( NB, N-K+1 );

               // Set relevant part of first half of WORK to zero,
               // holds dot products

               for (I = K; I <= N; I++) { // 170
                  WORK[I] = 0;
               } // 170

               for (J = K; J <= K + JB - 1; J++) { // 200

               // Find pivot, test for exit, else swap rows and columns
               // Update dot products, compute possible pivots which are
               // stored in the second half of WORK

                  for (I = J; I <= N; I++) { // 180

                     if ( J > K ) {
                        WORK[I] = WORK( I ) + double( CONJG( A( I, J-1 ) )* A( I, J-1 ) );
                     }
                     WORK[N+I] = double( A( I, I ) ) - WORK( I );

                  } // 180

                  if ( J > 1 ) {
                     ITEMP = MAXLOC( WORK( (N+J):(2*N) ), 1 );
                     PVT = ITEMP + J - 1;
                     AJJ = WORK( N+PVT );
                     if ( AJJ <= SSTOP || SISNAN( AJJ ) ) {
                        A[J][J] = AJJ;
                        GO TO 220;
                     }
                  }

                  if ( J != PVT ) {

                     // Pivot OK, so can now swap pivot rows and columns

                     A[PVT][PVT] = A( J, J );
                     cswap(J-1, A( J, 1 ), LDA, A( PVT, 1 ), LDA );
                     if (PVT < N) cswap( N-PVT, A( PVT+1, J ), 1, A( PVT+1, PVT ), 1 );
                     for (I = J + 1; I <= PVT - 1; I++) { // 190
                        CTEMP = CONJG( A( I, J ) );
                        A[I][J] = CONJG( A( PVT, I ) );
                        A[PVT][I] = CTEMP;
                     } // 190
                     A[PVT][J] = CONJG( A( PVT, J ) );

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

                  // Compute elements J+1:N of column J.

                  if ( J < N ) {
                     clacgv(J-1, A( J, 1 ), LDA );
                     cgemv('No Trans', N-J, J-K, -CONE, A( J+1, K ), LDA, A( J, K ), LDA, CONE, A( J+1, J ), 1 );
                     clacgv(J-1, A( J, 1 ), LDA );
                     csscal(N-J, ONE / AJJ, A( J+1, J ), 1 );
                  }

               } // 200

               // Update trailing matrix, J already incremented

               if ( K+JB <= N ) {
                  cherk('Lower', 'No Trans', N-J+1, JB, -ONE, A( J, K ), LDA, ONE, A( J, J ), LDA );
               }

            } // 210

         }
      }

      // Ran to completion, A has full rank

      RANK = N;

      GO TO 230;
      } // 220

      // Rank is the number of steps completed.  Set INFO = 1 to signal
      // that the factorization cannot be used to solve a system.

      RANK = J - 1;
      INFO = 1;

      } // 230
      }
