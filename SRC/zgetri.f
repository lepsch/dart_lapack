      SUBROUTINE ZGETRI( N, A, LDA, IPIV, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LWORK, N;
      // ..
      // .. Array Arguments ..
      int                IPIV( * );
      COMPLEX*16         A( LDA, * ), WORK( * )
      // ..

*  =====================================================================

      // .. Parameters ..
      COMPLEX*16         ZERO, ONE
      const              ZERO = ( 0.0D+0, 0.0D+0 ), ONE = ( 1.0D+0, 0.0D+0 ) ;
      // ..
      // .. Local Scalars ..
      bool               LQUERY;
      int                I, IWS, J, JB, JJ, JP, LDWORK, LWKOPT, NB, NBMIN, NN;
      // ..
      // .. External Functions ..
      int                ILAENV;
      // EXTERNAL ILAENV
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZGEMM, ZGEMV, ZSWAP, ZTRSM, ZTRTRI
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input parameters.

      INFO = 0
      NB = ILAENV( 1, 'ZGETRI', ' ', N, -1, -1, -1 )
      LWKOPT = MAX( 1, N*NB )
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      if ( N.LT.0 ) {
         INFO = -1
      } else if ( LDA.LT.MAX( 1, N ) ) {
         INFO = -3
      } else if ( LWORK.LT.MAX( 1, N ) .AND. .NOT.LQUERY ) {
         INFO = -6
      }
      if ( INFO.NE.0 ) {
         xerbla('ZGETRI', -INFO );
         RETURN
      } else if ( LQUERY ) {
         RETURN
      }

      // Quick return if possible

      IF( N.EQ.0 ) RETURN

      // Form inv(U).  If INFO > 0 from ZTRTRI, then U is singular,
      // and the inverse is not computed.

      ztrtri('Upper', 'Non-unit', N, A, LDA, INFO );
      IF( INFO.GT.0 ) RETURN

      NBMIN = 2
      LDWORK = N
      if ( NB.GT.1 .AND. NB.LT.N ) {
         IWS = MAX( LDWORK*NB, 1 )
         if ( LWORK.LT.IWS ) {
            NB = LWORK / LDWORK
            NBMIN = MAX( 2, ILAENV( 2, 'ZGETRI', ' ', N, -1, -1, -1 ) )
         }
      } else {
         IWS = N
      }

      // Solve the equation inv(A)*L = inv(U) for inv(A).

      if ( NB.LT.NBMIN .OR. NB.GE.N ) {

         // Use unblocked code.

         DO 20 J = N, 1, -1

            // Copy current column of L to WORK and replace with zeros.

            DO 10 I = J + 1, N
               WORK( I ) = A( I, J )
               A( I, J ) = ZERO
   10       CONTINUE

            // Compute current column of inv(A).

            IF( J.LT.N ) CALL ZGEMV( 'No transpose', N, N-J, -ONE, A( 1, J+1 ), LDA, WORK( J+1 ), 1, ONE, A( 1, J ), 1 )
   20    CONTINUE
      } else {

         // Use blocked code.

         NN = ( ( N-1 ) / NB )*NB + 1
         DO 50 J = NN, 1, -NB
            JB = MIN( NB, N-J+1 )

            // Copy current block column of L to WORK and replace with
            // zeros.

            DO 40 JJ = J, J + JB - 1
               DO 30 I = JJ + 1, N
                  WORK( I+( JJ-J )*LDWORK ) = A( I, JJ )
                  A( I, JJ ) = ZERO
   30          CONTINUE
   40       CONTINUE

            // Compute current block column of inv(A).

            IF( J+JB.LE.N ) CALL ZGEMM( 'No transpose', 'No transpose', N, JB, N-J-JB+1, -ONE, A( 1, J+JB ), LDA, WORK( J+JB ), LDWORK, ONE, A( 1, J ), LDA )
            ztrsm('Right', 'Lower', 'No transpose', 'Unit', N, JB, ONE, WORK( J ), LDWORK, A( 1, J ), LDA );
   50    CONTINUE
      }

      // Apply column interchanges.

      DO 60 J = N - 1, 1, -1
         JP = IPIV( J )
         IF( JP.NE.J ) CALL ZSWAP( N, A( 1, J ), 1, A( 1, JP ), 1 )
   60 CONTINUE

      WORK( 1 ) = IWS
      RETURN

      // End of ZGETRI

      }
