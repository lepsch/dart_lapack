      SUBROUTINE CGEBRD( M, N, A, LDA, D, E, TAUQ, TAUP, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LWORK, M, N;
      // ..
      // .. Array Arguments ..
      REAL               D( * ), E( * )
      COMPLEX            A( LDA, * ), TAUP( * ), TAUQ( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX            ONE
      PARAMETER          ( ONE = ( 1.0E+0, 0.0E+0 ) )
      // ..
      // .. Local Scalars ..
      bool               LQUERY;
      int                I, IINFO, J, LDWRKX, LDWRKY, LWKMIN, LWKOPT, MINMN, NB, NBMIN, NX, WS;
      // ..
      // .. External Subroutines ..
      // EXTERNAL CGEBD2, CGEMM, CLABRD, XERBLA
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN, REAL
      // ..
      // .. External Functions ..
      int                ILAENV;
      REAL               SROUNDUP_LWORK
      // EXTERNAL ILAENV, SROUNDUP_LWORK
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      INFO = 0
      MINMN = MIN( M, N )
      IF( MINMN.EQ.0 ) THEN
         LWKMIN = 1
         LWKOPT = 1
      ELSE
         LWKMIN = MAX( M, N )
         NB = MAX( 1, ILAENV( 1, 'CGEBRD', ' ', M, N, -1, -1 ) )
         LWKOPT = ( M+N )*NB
      END IF
      WORK( 1 ) = SROUNDUP_LWORK( LWKOPT )
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( N.LT.0 ) THEN
         INFO = -2
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -4
      ELSE IF( LWORK.LT.LWKMIN .AND. .NOT.LQUERY ) THEN
         INFO = -10
      END IF
      IF( INFO.LT.0 ) THEN
         CALL XERBLA( 'CGEBRD', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF

      // Quick return if possible

      IF( MINMN.EQ.0 ) THEN
         WORK( 1 ) = 1
         RETURN
      END IF

      WS = MAX( M, N )
      LDWRKX = M
      LDWRKY = N

      IF( NB.GT.1 .AND. NB.LT.MINMN ) THEN

         // Set the crossover point NX.

         NX = MAX( NB, ILAENV( 3, 'CGEBRD', ' ', M, N, -1, -1 ) )

         // Determine when to switch from blocked to unblocked code.

         IF( NX.LT.MINMN ) THEN
            WS = LWKOPT
            IF( LWORK.LT.WS ) THEN

               // Not enough work space for the optimal NB, consider using
               // a smaller block size.

               NBMIN = ILAENV( 2, 'CGEBRD', ' ', M, N, -1, -1 )
               IF( LWORK.GE.( M+N )*NBMIN ) THEN
                  NB = LWORK / ( M+N )
               ELSE
                  NB = 1
                  NX = MINMN
               END IF
            END IF
         END IF
      ELSE
         NX = MINMN
      END IF

      DO 30 I = 1, MINMN - NX, NB

         // Reduce rows and columns i:i+ib-1 to bidiagonal form and return
        t // he matrices X and Y which are needed to update the unreduced
         // part of the matrix

         CALL CLABRD( M-I+1, N-I+1, NB, A( I, I ), LDA, D( I ), E( I ), TAUQ( I ), TAUP( I ), WORK, LDWRKX, WORK( LDWRKX*NB+1 ), LDWRKY )

         // Update the trailing submatrix A(i+ib:m,i+ib:n), using
         // an update of the form  A := A - V*Y**H - X*U**H

         CALL CGEMM( 'No transpose', 'Conjugate transpose', M-I-NB+1, N-I-NB+1, NB, -ONE, A( I+NB, I ), LDA, WORK( LDWRKX*NB+NB+1 ), LDWRKY, ONE, A( I+NB, I+NB ), LDA )
         CALL CGEMM( 'No transpose', 'No transpose', M-I-NB+1, N-I-NB+1, NB, -ONE, WORK( NB+1 ), LDWRKX, A( I, I+NB ), LDA, ONE, A( I+NB, I+NB ), LDA )

         // Copy diagonal and off-diagonal elements of B back into A

         IF( M.GE.N ) THEN
            DO 10 J = I, I + NB - 1
               A( J, J ) = D( J )
               A( J, J+1 ) = E( J )
   10       CONTINUE
         ELSE
            DO 20 J = I, I + NB - 1
               A( J, J ) = D( J )
               A( J+1, J ) = E( J )
   20       CONTINUE
         END IF
   30 CONTINUE

      // Use unblocked code to reduce the remainder of the matrix

      CALL CGEBD2( M-I+1, N-I+1, A( I, I ), LDA, D( I ), E( I ), TAUQ( I ), TAUP( I ), WORK, IINFO )
      WORK( 1 ) = SROUNDUP_LWORK( WS )
      RETURN

      // End of CGEBRD

      }
