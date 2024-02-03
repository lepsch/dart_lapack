      SUBROUTINE SGEQPF( M, N, A, LDA, JPVT, TAU, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, M, N;
      // ..
      // .. Array Arguments ..
      int                JPVT( * );
      REAL               A( LDA, * ), TAU( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ZERO, ONE
      const              ZERO = 0.0E+0, ONE = 1.0E+0 ;
      // ..
      // .. Local Scalars ..
      int                I, ITEMP, J, MA, MN, PVT;
      REAL               AII, TEMP, TEMP2, TOL3Z
      // ..
      // .. External Subroutines ..
      // EXTERNAL SGEQR2, SLARF, SLARFG, SORM2R, SSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC ABS, MAX, MIN, SQRT
      // ..
      // .. External Functions ..
      int                ISAMAX;
      REAL               SLAMCH, SNRM2
      // EXTERNAL ISAMAX, SLAMCH, SNRM2
      // ..
      // .. Executable Statements ..

      // Test the input arguments

      INFO = 0
      if ( M.LT.0 ) {
         INFO = -1
      } else if ( N.LT.0 ) {
         INFO = -2
      } else if ( LDA.LT.MAX( 1, M ) ) {
         INFO = -4
      }
      if ( INFO.NE.0 ) {
         CALL XERBLA( 'SGEQPF', -INFO )
         RETURN
      }

      MN = MIN( M, N )
      TOL3Z = SQRT(SLAMCH('Epsilon'))

      // Move initial columns up front

      ITEMP = 1
      DO 10 I = 1, N
         if ( JPVT( I ).NE.0 ) {
            if ( I.NE.ITEMP ) {
               CALL SSWAP( M, A( 1, I ), 1, A( 1, ITEMP ), 1 )
               JPVT( I ) = JPVT( ITEMP )
               JPVT( ITEMP ) = I
            } else {
               JPVT( I ) = I
            }
            ITEMP = ITEMP + 1
         } else {
            JPVT( I ) = I
         }
   10 CONTINUE
      ITEMP = ITEMP - 1

      // Compute the QR factorization and update remaining columns

      if ( ITEMP.GT.0 ) {
         MA = MIN( ITEMP, M )
         CALL SGEQR2( M, MA, A, LDA, TAU, WORK, INFO )
         if ( MA.LT.N ) {
            CALL SORM2R( 'Left', 'Transpose', M, N-MA, MA, A, LDA, TAU, A( 1, MA+1 ), LDA, WORK, INFO )
         }
      }

      if ( ITEMP.LT.MN ) {

         // Initialize partial column norms. The first n elements of
         // work store the exact column norms.

         DO 20 I = ITEMP + 1, N
            WORK( I ) = SNRM2( M-ITEMP, A( ITEMP+1, I ), 1 )
            WORK( N+I ) = WORK( I )
   20    CONTINUE

         // Compute factorization

         DO 40 I = ITEMP + 1, MN

            // Determine ith pivot column and swap if necessary

            PVT = ( I-1 ) + ISAMAX( N-I+1, WORK( I ), 1 )

            if ( PVT.NE.I ) {
               CALL SSWAP( M, A( 1, PVT ), 1, A( 1, I ), 1 )
               ITEMP = JPVT( PVT )
               JPVT( PVT ) = JPVT( I )
               JPVT( I ) = ITEMP
               WORK( PVT ) = WORK( I )
               WORK( N+PVT ) = WORK( N+I )
            }

            // Generate elementary reflector H(i)

            if ( I.LT.M ) {
               CALL SLARFG( M-I+1, A( I, I ), A( I+1, I ), 1, TAU( I ) )
            } else {
               CALL SLARFG( 1, A( M, M ), A( M, M ), 1, TAU( M ) )
            }

            if ( I.LT.N ) {

               // Apply H(i) to A(i:m,i+1:n) from the left

               AII = A( I, I )
               A( I, I ) = ONE
               CALL SLARF( 'LEFT', M-I+1, N-I, A( I, I ), 1, TAU( I ), A( I, I+1 ), LDA, WORK( 2*N+1 ) )
               A( I, I ) = AII
            }

            // Update partial column norms

            DO 30 J = I + 1, N
               if ( WORK( J ).NE.ZERO ) {

                  // NOTE: The following 4 lines follow from the analysis in
                  // Lapack Working Note 176.

                  TEMP = ABS( A( I, J ) ) / WORK( J )
                  TEMP = MAX( ZERO, ( ONE+TEMP )*( ONE-TEMP ) )
                  TEMP2 = TEMP*( WORK( J ) / WORK( N+J ) )**2
                  if ( TEMP2 .LE. TOL3Z ) {
                     if ( M-I.GT.0 ) {
                        WORK( J ) = SNRM2( M-I, A( I+1, J ), 1 )
                        WORK( N+J ) = WORK( J )
                     } else {
                        WORK( J ) = ZERO
                        WORK( N+J ) = ZERO
                     }
                  } else {
                     WORK( J ) = WORK( J )*SQRT( TEMP )
                  }
               }
   30       CONTINUE

   40    CONTINUE
      }
      RETURN

      // End of SGEQPF

      }
