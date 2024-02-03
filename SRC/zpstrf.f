      SUBROUTINE ZPSTRF( UPLO, N, A, LDA, PIV, RANK, TOL, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      double             TOL;
      int                INFO, LDA, N, RANK;
      String             UPLO;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * )
      double             WORK( 2*N );
      int                PIV( N );
      // ..

*  =====================================================================

      // .. Parameters ..
      double             ONE, ZERO;
      const              ONE = 1.0D+0, ZERO = 0.0D+0 ;
      COMPLEX*16         CONE
      const              CONE = ( 1.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      COMPLEX*16         ZTEMP
      double             AJJ, DSTOP, DTEMP;
      int                I, ITEMP, J, JB, K, NB, PVT;
      bool               UPPER;
      // ..
      // .. External Functions ..
      double             DLAMCH;
      int                ILAENV;
      bool               LSAME, DISNAN;
      // EXTERNAL DLAMCH, ILAENV, LSAME, DISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL ZDSCAL, ZGEMV, ZHERK, ZLACGV, ZPSTF2, ZSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC DBLE, DCONJG, MAX, MIN, SQRT, MAXLOC
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      if ( .NOT.UPPER && .NOT.LSAME( UPLO, 'L' ) ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -4
      }
      if ( INFO != 0 ) {
         xerbla('ZPSTRF', -INFO );
         RETURN
      }

      // Quick return if possible

      if (N == 0) RETURN;

      // Get block size

      NB = ILAENV( 1, 'ZPOTRF', UPLO, N, -1, -1, -1 )
      if ( NB.LE.1 || NB.GE.N ) {

         // Use unblocked code

         zpstf2(UPLO, N, A( 1, 1 ), LDA, PIV, RANK, TOL, WORK, INFO );
         GO TO 230

      } else {

      // Initialize PIV

         for (I = 1; I <= N; I++) { // 100
            PIV( I ) = I
         } // 100

      // Compute stopping value

         for (I = 1; I <= N; I++) { // 110
            WORK( I ) = DBLE( A( I, I ) )
         } // 110
         PVT = MAXLOC( WORK( 1:N ), 1 )
         AJJ = DBLE( A( PVT, PVT ) )
         if ( AJJ.LE.ZERO || DISNAN( AJJ ) ) {
            RANK = 0
            INFO = 1
            GO TO 230
         }

      // Compute stopping value if not supplied

         if ( TOL.LT.ZERO ) {
            DSTOP = N * DLAMCH( 'Epsilon' ) * AJJ
         } else {
            DSTOP = TOL
         }


         if ( UPPER ) {

            // Compute the Cholesky factorization P**T * A * P = U**H * U

            DO 160 K = 1, N, NB

               // Account for last block not being NB wide

               JB = MIN( NB, N-K+1 )

               // Set relevant part of first half of WORK to zero,
               // holds dot products

               for (I = K; I <= N; I++) { // 120
                  WORK( I ) = 0
               } // 120

               for (J = K; J <= K + JB - 1; J++) { // 150

               // Find pivot, test for exit, else swap rows and columns
               // Update dot products, compute possible pivots which are
               // stored in the second half of WORK

                  for (I = J; I <= N; I++) { // 130

                     if ( J.GT.K ) {
                        WORK( I ) = WORK( I ) + DBLE( DCONJG( A( J-1, I ) )* A( J-1, I ) )
                     }
                     WORK( N+I ) = DBLE( A( I, I ) ) - WORK( I )

                  } // 130

                  if ( J.GT.1 ) {
                     ITEMP = MAXLOC( WORK( (N+J):(2*N) ), 1 )
                     PVT = ITEMP + J - 1
                     AJJ = WORK( N+PVT )
                     if ( AJJ.LE.DSTOP || DISNAN( AJJ ) ) {
                        A( J, J ) = AJJ
                        GO TO 220
                     }
                  }

                  if ( J != PVT ) {

                     // Pivot OK, so can now swap pivot rows and columns

                     A( PVT, PVT ) = A( J, J )
                     zswap(J-1, A( 1, J ), 1, A( 1, PVT ), 1 );
                     if (PVT.LT.N) CALL ZSWAP( N-PVT, A( J, PVT+1 ), LDA, A( PVT, PVT+1 ), LDA );
                     for (I = J + 1; I <= PVT - 1; I++) { // 140
                        ZTEMP = DCONJG( A( J, I ) )
                        A( J, I ) = DCONJG( A( I, PVT ) )
                        A( I, PVT ) = ZTEMP
                     } // 140
                     A( J, PVT ) = DCONJG( A( J, PVT ) )

                     // Swap dot products and PIV

                     DTEMP = WORK( J )
                     WORK( J ) = WORK( PVT )
                     WORK( PVT ) = DTEMP
                     ITEMP = PIV( PVT )
                     PIV( PVT ) = PIV( J )
                     PIV( J ) = ITEMP
                  }

                  AJJ = SQRT( AJJ )
                  A( J, J ) = AJJ

                  // Compute elements J+1:N of row J.

                  if ( J.LT.N ) {
                     zlacgv(J-1, A( 1, J ), 1 );
                     zgemv('Trans', J-K, N-J, -CONE, A( K, J+1 ), LDA, A( K, J ), 1, CONE, A( J, J+1 ), LDA );
                     zlacgv(J-1, A( 1, J ), 1 );
                     zdscal(N-J, ONE / AJJ, A( J, J+1 ), LDA );
                  }

               } // 150

               // Update trailing matrix, J already incremented

               if ( K+JB.LE.N ) {
                  zherk('Upper', 'Conj Trans', N-J+1, JB, -ONE, A( K, J ), LDA, ONE, A( J, J ), LDA );
               }

            } // 160

         } else {

         // Compute the Cholesky factorization P**T * A * P = L * L**H

            DO 210 K = 1, N, NB

               // Account for last block not being NB wide

               JB = MIN( NB, N-K+1 )

               // Set relevant part of first half of WORK to zero,
               // holds dot products

               for (I = K; I <= N; I++) { // 170
                  WORK( I ) = 0
               } // 170

               for (J = K; J <= K + JB - 1; J++) { // 200

               // Find pivot, test for exit, else swap rows and columns
               // Update dot products, compute possible pivots which are
               // stored in the second half of WORK

                  for (I = J; I <= N; I++) { // 180

                     if ( J.GT.K ) {
                        WORK( I ) = WORK( I ) + DBLE( DCONJG( A( I, J-1 ) )* A( I, J-1 ) )
                     }
                     WORK( N+I ) = DBLE( A( I, I ) ) - WORK( I )

                  } // 180

                  if ( J.GT.1 ) {
                     ITEMP = MAXLOC( WORK( (N+J):(2*N) ), 1 )
                     PVT = ITEMP + J - 1
                     AJJ = WORK( N+PVT )
                     if ( AJJ.LE.DSTOP || DISNAN( AJJ ) ) {
                        A( J, J ) = AJJ
                        GO TO 220
                     }
                  }

                  if ( J != PVT ) {

                     // Pivot OK, so can now swap pivot rows and columns

                     A( PVT, PVT ) = A( J, J )
                     zswap(J-1, A( J, 1 ), LDA, A( PVT, 1 ), LDA );
                     if (PVT.LT.N) CALL ZSWAP( N-PVT, A( PVT+1, J ), 1, A( PVT+1, PVT ), 1 );
                     for (I = J + 1; I <= PVT - 1; I++) { // 190
                        ZTEMP = DCONJG( A( I, J ) )
                        A( I, J ) = DCONJG( A( PVT, I ) )
                        A( PVT, I ) = ZTEMP
                     } // 190
                     A( PVT, J ) = DCONJG( A( PVT, J ) )


                     // Swap dot products and PIV

                     DTEMP = WORK( J )
                     WORK( J ) = WORK( PVT )
                     WORK( PVT ) = DTEMP
                     ITEMP = PIV( PVT )
                     PIV( PVT ) = PIV( J )
                     PIV( J ) = ITEMP
                  }

                  AJJ = SQRT( AJJ )
                  A( J, J ) = AJJ

                  // Compute elements J+1:N of column J.

                  if ( J.LT.N ) {
                     zlacgv(J-1, A( J, 1 ), LDA );
                     zgemv('No Trans', N-J, J-K, -CONE, A( J+1, K ), LDA, A( J, K ), LDA, CONE, A( J+1, J ), 1 );
                     zlacgv(J-1, A( J, 1 ), LDA );
                     zdscal(N-J, ONE / AJJ, A( J+1, J ), 1 );
                  }

               } // 200

               // Update trailing matrix, J already incremented

               if ( K+JB.LE.N ) {
                  zherk('Lower', 'No Trans', N-J+1, JB, -ONE, A( J, K ), LDA, ONE, A( J, J ), LDA );
               }

            } // 210

         }
      }

      // Ran to completion, A has full rank

      RANK = N

      GO TO 230
      } // 220

      // Rank is the number of steps completed.  Set INFO = 1 to signal
      // that the factorization cannot be used to solve a system.

      RANK = J - 1
      INFO = 1

      } // 230
      RETURN

      // End of ZPSTRF

      }
