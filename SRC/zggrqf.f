      SUBROUTINE ZGGRQF( M, P, N, A, LDA, TAUA, B, LDB, TAUB, WORK, LWORK, INFO )

*  -- LAPACK computational routine --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--

      // .. Scalar Arguments ..
      int                INFO, LDA, LDB, LWORK, M, N, P;
      // ..
      // .. Array Arguments ..
      COMPLEX*16         A( LDA, * ), B( LDB, * ), TAUA( * ), TAUB( * ), WORK( * )
      // ..

*  =====================================================================

      // .. Local Scalars ..
      bool               LQUERY;
      int                LOPT, LWKOPT, NB, NB1, NB2, NB3;
      // ..
      // .. External Subroutines ..
      // EXTERNAL XERBLA, ZGEQRF, ZGERQF, ZUNMRQ
      // ..
      // .. External Functions ..
      int                ILAENV;
      // EXTERNAL ILAENV
      // ..
      // .. Intrinsic Functions ..
      // INTRINSIC INT, MAX, MIN
      // ..
      // .. Executable Statements ..

      // Test the input parameters

      INFO = 0
      NB1 = ILAENV( 1, 'ZGERQF', ' ', M, N, -1, -1 )
      NB2 = ILAENV( 1, 'ZGEQRF', ' ', P, N, -1, -1 )
      NB3 = ILAENV( 1, 'ZUNMRQ', ' ', M, N, P, -1 )
      NB = MAX( NB1, NB2, NB3 )
      LWKOPT = MAX( 1, MAX( N, M, P )*NB )
      WORK( 1 ) = LWKOPT
      LQUERY = ( LWORK.EQ.-1 )
      IF( M.LT.0 ) THEN
         INFO = -1
      ELSE IF( P.LT.0 ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -3
      ELSE IF( LDA.LT.MAX( 1, M ) ) THEN
         INFO = -5
      ELSE IF( LDB.LT.MAX( 1, P ) ) THEN
         INFO = -8
      ELSE IF( LWORK.LT.MAX( 1, M, P, N ) .AND. .NOT.LQUERY ) THEN
         INFO = -11
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'ZGGRQF', -INFO )
         RETURN
      ELSE IF( LQUERY ) THEN
         RETURN
      END IF

      // RQ factorization of M-by-N matrix A: A = R*Q

      CALL ZGERQF( M, N, A, LDA, TAUA, WORK, LWORK, INFO )
      LOPT = INT( WORK( 1 ) )

      // Update B := B*Q**H

      CALL ZUNMRQ( 'Right', 'Conjugate Transpose', P, N, MIN( M, N ), A( MAX( 1, M-N+1 ), 1 ), LDA, TAUA, B, LDB, WORK, LWORK, INFO )
      LOPT = MAX( LOPT, INT( WORK( 1 ) ) )

      // QR factorization of P-by-N matrix B: B = Z*T

      CALL ZGEQRF( P, N, B, LDB, TAUB, WORK, LWORK, INFO )
      WORK( 1 ) = MAX( LOPT, INT( WORK( 1 ) ) )

      RETURN

      // End of ZGGRQF

      }
