      SUBROUTINE CPSTF2( UPLO, N, A, LDA, PIV, RANK, TOL, WORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      REAL               TOL
      int                INFO, LDA, N, RANK;
      String             UPLO;
      // ..
      // .. Array Arguments ..
      COMPLEX            A( LDA, * )
      REAL               WORK( 2*N )
      int                PIV( N );
      // ..

*  =====================================================================

      // .. Parameters ..
      REAL               ONE, ZERO
      const              ONE = 1.0E+0, ZERO = 0.0E+0 ;
      COMPLEX            CONE
      const              CONE = ( 1.0E+0, 0.0E+0 ) ;
      // ..
      // .. Local Scalars ..
      COMPLEX            CTEMP
      REAL               AJJ, SSTOP, STEMP
      int                I, ITEMP, J, PVT;
      bool               UPPER;
      // ..
      // .. External Functions ..
      REAL               SLAMCH
      bool               LSAME, SISNAN;
      // EXTERNAL SLAMCH, LSAME, SISNAN
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEMV, CLACGV, CSSCAL, CSWAP, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC CONJG, MAX, REAL, SQRT
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      INFO = 0
      UPPER = LSAME( UPLO, 'U' )
      IF( .NOT.UPPER .AND. .NOT.LSAME( UPLO, 'L' ) ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, N ) ) THEN
         INFO = -4
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CPSTF2', -INFO )
         RETURN
      END IF

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      // Initialize PIV

      DO 100 I = 1, N
         PIV( I ) = I
  100 CONTINUE

      // Compute stopping value

      DO 110 I = 1, N
         WORK( I ) = REAL( A( I, I ) )
  110 CONTINUE
      PVT = MAXLOC( WORK( 1:N ), 1 )
      AJJ = REAL ( A( PVT, PVT ) )
      IF( AJJ.LE.ZERO.OR.SISNAN( AJJ ) ) THEN
         RANK = 0
         INFO = 1
         GO TO 200
      END IF

      // Compute stopping value if not supplied

      IF( TOL.LT.ZERO ) THEN
         SSTOP = N * SLAMCH( 'Epsilon' ) * AJJ
      } else {
         SSTOP = TOL
      END IF

      // Set first half of WORK to zero, holds dot products

      DO 120 I = 1, N
         WORK( I ) = 0
  120 CONTINUE

      IF( UPPER ) THEN

         // Compute the Cholesky factorization P**T * A * P = U**H * U

         DO 150 J = 1, N

         // Find pivot, test for exit, else swap rows and columns
         // Update dot products, compute possible pivots which are
         // stored in the second half of WORK

            DO 130 I = J, N

               IF( J.GT.1 ) THEN
                  WORK( I ) = WORK( I ) + REAL( CONJG( A( J-1, I ) )* A( J-1, I ) )
               END IF
               WORK( N+I ) = REAL( A( I, I ) ) - WORK( I )

  130       CONTINUE

            IF( J.GT.1 ) THEN
               ITEMP = MAXLOC( WORK( (N+J):(2*N) ), 1 )
               PVT = ITEMP + J - 1
               AJJ = WORK( N+PVT )
               IF( AJJ.LE.SSTOP.OR.SISNAN( AJJ ) ) THEN
                  A( J, J ) = AJJ
                  GO TO 190
               END IF
            END IF

            IF( J.NE.PVT ) THEN

               // Pivot OK, so can now swap pivot rows and columns

               A( PVT, PVT ) = A( J, J )
               CALL CSWAP( J-1, A( 1, J ), 1, A( 1, PVT ), 1 )
               IF( PVT.LT.N ) CALL CSWAP( N-PVT, A( J, PVT+1 ), LDA, A( PVT, PVT+1 ), LDA )
               DO 140 I = J + 1, PVT - 1
                  CTEMP = CONJG( A( J, I ) )
                  A( J, I ) = CONJG( A( I, PVT ) )
                  A( I, PVT ) = CTEMP
  140          CONTINUE
               A( J, PVT ) = CONJG( A( J, PVT ) )

               // Swap dot products and PIV

               STEMP = WORK( J )
               WORK( J ) = WORK( PVT )
               WORK( PVT ) = STEMP
               ITEMP = PIV( PVT )
               PIV( PVT ) = PIV( J )
               PIV( J ) = ITEMP
            END IF

            AJJ = SQRT( AJJ )
            A( J, J ) = AJJ

            // Compute elements J+1:N of row J

            IF( J.LT.N ) THEN
               CALL CLACGV( J-1, A( 1, J ), 1 )
               CALL CGEMV( 'Trans', J-1, N-J, -CONE, A( 1, J+1 ), LDA, A( 1, J ), 1, CONE, A( J, J+1 ), LDA )
               CALL CLACGV( J-1, A( 1, J ), 1 )
               CALL CSSCAL( N-J, ONE / AJJ, A( J, J+1 ), LDA )
            END IF

  150    CONTINUE

      } else {

         // Compute the Cholesky factorization P**T * A * P = L * L**H

         DO 180 J = 1, N

         // Find pivot, test for exit, else swap rows and columns
         // Update dot products, compute possible pivots which are
         // stored in the second half of WORK

            DO 160 I = J, N

               IF( J.GT.1 ) THEN
                  WORK( I ) = WORK( I ) + REAL( CONJG( A( I, J-1 ) )* A( I, J-1 ) )
               END IF
               WORK( N+I ) = REAL( A( I, I ) ) - WORK( I )

  160       CONTINUE

            IF( J.GT.1 ) THEN
               ITEMP = MAXLOC( WORK( (N+J):(2*N) ), 1 )
               PVT = ITEMP + J - 1
               AJJ = WORK( N+PVT )
               IF( AJJ.LE.SSTOP.OR.SISNAN( AJJ ) ) THEN
                  A( J, J ) = AJJ
                  GO TO 190
               END IF
            END IF

            IF( J.NE.PVT ) THEN

               // Pivot OK, so can now swap pivot rows and columns

               A( PVT, PVT ) = A( J, J )
               CALL CSWAP( J-1, A( J, 1 ), LDA, A( PVT, 1 ), LDA )
               IF( PVT.LT.N ) CALL CSWAP( N-PVT, A( PVT+1, J ), 1, A( PVT+1, PVT ), 1 )
               DO 170 I = J + 1, PVT - 1
                  CTEMP = CONJG( A( I, J ) )
                  A( I, J ) = CONJG( A( PVT, I ) )
                  A( PVT, I ) = CTEMP
  170          CONTINUE
               A( PVT, J ) = CONJG( A( PVT, J ) )

               // Swap dot products and PIV

               STEMP = WORK( J )
               WORK( J ) = WORK( PVT )
               WORK( PVT ) = STEMP
               ITEMP = PIV( PVT )
               PIV( PVT ) = PIV( J )
               PIV( J ) = ITEMP
            END IF

            AJJ = SQRT( AJJ )
            A( J, J ) = AJJ

            // Compute elements J+1:N of column J

            IF( J.LT.N ) THEN
               CALL CLACGV( J-1, A( J, 1 ), LDA )
               CALL CGEMV( 'No Trans', N-J, J-1, -CONE, A( J+1, 1 ), LDA, A( J, 1 ), LDA, CONE, A( J+1, J ), 1 )
               CALL CLACGV( J-1, A( J, 1 ), LDA )
               CALL CSSCAL( N-J, ONE / AJJ, A( J+1, J ), 1 )
            END IF

  180    CONTINUE

      END IF

      // Ran to completion, A has full rank

      RANK = N

      GO TO 200
  190 CONTINUE

      // Rank is number of steps completed.  Set INFO = 1 to signal
     t // hat the factorization cannot be used to solve a system.

      RANK = J - 1
      INFO = 1

  200 CONTINUE
      RETURN

      // End of CPSTF2

      }
